"""Read original P301 CRF PDF bytes. No adjudication or CSV findings are inputs.
Python 3.10+, PyMuPDF. Run from any directory; see CRF_EVIDENCE.md."""
from pathlib import Path
import argparse, collections, concurrent.futures, csv, gzip, hashlib, json, os, re, shutil, tempfile, zlib
from lib import extract, fieldscan, front, illscan, scan

ROOT=Path(__file__).resolve().parents[1]
def csvout(path,cols,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(cols);w.writerows(rows)
def work(job):
    f,name,size,crc,out=job; out=Path(out);f=Path(f);doc=Path(name).stem
    with f.open('rb') as h:
        c=0
        for b in iter(lambda:h.read(1024*1024),b''):c=zlib.crc32(b,c)
    if f.stat().st_size!=size or '%08x'%c!=crc:raise ValueError(name+': archive CRC32 or size mismatch')
    for mod in (extract,scan,fieldscan,front):mod.w=out
    fieldscan.out=out/'hist';front.out=out/'front';alias=out/'aliases'/name
    try:
        try:os.link(f,alias)
        except OSError:shutil.copyfile(f,alias)
        extract.extract(alias)
    finally:alias.unlink(missing_ok=True)
    gz=out/'text'/(doc+'.gz')
    # Front uses the dark filled mark, not the printed options. Unreadable stays blank.
    front.parse(gz,f);fieldscan.go(gz)
    raw=gzip.decompress(gz.read_bytes()).decode().split('\f')[:-1]
    exposures=[]; dosechanges=[]
    for i,s in enumerate(raw,1):
        fd=re.search(r'^Folder: (.+)$',s,re.M); fm=re.search(r'^Form: (.+)$',s,re.M)
        if not(fd and fm) or '\nAudit\n' not in s:continue
        bt=re.search(r'FDA-CBER-[0-9-]+',s);bt=bt.group() if bt else ''
        if fd[1].startswith('Visit 2 Day 29') and fm[1]=='Exposure' and re.search(r'Was study treatment given\?\s*\nAudit',s):
            vals=re.findall(r"User\s+entered\s+(?:'([^']{1,120})'|(empty))",s,re.I)
            exposures += [[i,bt,' '.join(v.split()) if v else '<EMPTY>'] for v,_ in vals]
        if fm[1]=='Unblinding' and re.search(r'Actual Dose 2\s*\nAudit',s):
            # Store the source locator and flags, not staff names. Date/value histories are in hist/.
            dosechanges.append([i,bt,bool(re.search('in order to populate', ' '.join(s.split()),re.I))])
    visit_hits=[[pg,re.search(r'FDA-CBER-[0-9-]+',raw[pg-1]).group() if re.search(r'FDA-CBER-[0-9-]+',raw[pg-1]) else '',snip] for _,pg,snip in illscan.look(gz)]
    (out/'source'/(doc+'.json')).write_text(json.dumps(dict(day29_given=exposures,dose2_audit=dosechanges,illness_candidates=visit_hits)),encoding='utf8')
    if not (out/'KEEP_TEXT').exists():gz.unlink()
    return name
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True,type=Path,help='Tree of locally extracted original PDFs (nested directories supported)')
    p.add_argument('--out',required=True,type=Path,help='NEW, separate scan directory');p.add_argument('--jobs',default=2,type=int);p.add_argument('--allow-partial',action='store_true',help='Allow a subset; exclusions are listed, not counted as negatives');p.add_argument('--keep-text',action='store_true',help='Save large, local-only compressed page texts for manual inspection')
    a=p.parse_args();out=a.out.resolve();root=a.root.resolve()
    if not root.is_dir() or a.jobs<1:p.error('existing source directory and positive --jobs required')
    if out.exists() and any(out.iterdir()):p.error('--out must be new or empty (no mixed versions/cached runs)')
    out.mkdir(parents=True,exist_ok=True)
    for k in ('keys','text','hist','front','source','aliases'):(out/k).mkdir()
    if a.keep_text:(out/'KEEP_TEXT').touch()
    manifest=ROOT/'evidence/corpus_manifest.csv'
    with manifest.open(encoding='utf-8-sig') as h:ms={r['name'].lower():r for r in csv.DictReader(h) if r['class']=='full'}
    pat=re.compile(r'((?:125752|27073)_s\d+_m5_crf_mrna-1273-p301-us\d{7}\.pdf)$',re.I);found=collections.defaultdict(list)
    for f in root.rglob('*'):
        match=pat.search(f.name)
        if match and match[1].lower() in ms and f.is_file() and out not in f.resolve().parents:found[match[1].lower()].append(f)
    chosen={};bad=[]
    for k,fs in found.items():
        good=[f for f in sorted(fs) if f.stat().st_size==int(ms[k]['us']) and ('%08x'%zlib.crc32(f.read_bytes())==ms[k]['crc32_hex'])]
        if not good:bad.append(k)
        else:chosen[k]=good[0]
    csvout(out/'coverage.csv',['file','status'],[(r['name'],'selected' if k in chosen else 'invalid bytes' if k in bad else 'missing locally') for k,r in ms.items()])
    if bad or not chosen or (not a.allow_partial and len(chosen)!=len(ms)):p.error('Incomplete or incorrect source set. See coverage.csv; partial mode never accepts invalid bytes.')
    jobs=[(str(f),ms[k]['name'],int(ms[k]['us']),ms[k]['crc32_hex'],str(out)) for k,f in sorted(chosen.items())]
    with concurrent.futures.ProcessPoolExecutor(a.jobs) as ex:
        for i,_ in enumerate(ex.map(work,jobs),1):
            if i%25==0:print(f'{i}/{len(jobs)}',flush=True)
    (out/'FINISHED.json').write_text(json.dumps({'read':len(jobs),'expected':len(ms),'partial':len(jobs)!=len(ms),'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest()}),encoding='utf8')
    print(f'{len(jobs)}/{len(ms)} substantive exports read. Next: reconcile_fields.py.')
if __name__=='__main__':main()

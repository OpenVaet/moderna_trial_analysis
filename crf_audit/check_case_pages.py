"""Verify page identity for the 19 reviewed examples and the 23 targeted EUA serious forms.
This is a targeted source check, not automated adjudication or a search for additional cases."""
import argparse,collections,re,zlib
from pathlib import Path
import fitz
from reconcile_fields import ROOT,read,save,date,normal
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--allow-partial',action='store_true');a=p.parse_args()
 e=ROOT/'evidence';manifest=read(e/'corpus_manifest.csv').set_index('name');cur=read(e/'reviewed_cases.csv');sae=read(e/'previous_23_cases_CRF_crosswalk.csv')
 want=set(cur.filename)|set(sae.EUA_file);local={};pat=re.compile(r'((?:27073|125752)_s\d+_m5_crf_mrna-1273-p301-us\d{7}\.pdf)$',re.I);bylower={n.lower():n for n in want}
 for f in a.root.rglob('*'):
  m=pat.search(f.name)
  if not f.is_file() or not m or m[1].lower() not in bylower:continue
  n=bylower[m[1].lower()];r=manifest.loc[n]
  if f.stat().st_size==int(r.us) and '%08x'%zlib.crc32(f.read_bytes())==r.crc32_hex:local[n]=f
 missing=sorted(want-local.keys())
 if (missing and not a.allow_partial) or not local:raise SystemExit('Missing or invalid documents: '+', '.join(missing))
 if a.out.exists() and any(a.out.iterdir()):p.error('Use an empty output directory')
 a.out.mkdir(exist_ok=True,parents=True);notes=[];fatal=[]
 for n,f in sorted(local.items()):
  with fitz.open(f) as d:
   for r in cur[cur.filename==n].itertuples():
    pages=[int(s.strip()) for s in r.PDF_pages_1based.split(',')];bates=r.FDA_Bates.split(',');assert len(pages)==len(bates)
    for pg,bt in zip(pages,bates):
     t=d[pg-1].get_text();ok=(r.subject in t and bt.strip() in t);notes.append([r.finding,r.subject,n,pg,bt.strip(),ok])
   for r in sae[sae.EUA_file==n].itertuples():
    pg=int(r.EUA_form_page);sp=int(r.EUA_CRF_seriousness_page);t=d[pg-1].get_text();u=d[sp-1].get_text()
    term=re.search(r'\nAdverse event\s*\n(.+?)\nWas this a medically-attended AE\?',t,re.S);st=re.search(r'Start date\s*\(dd MMM yyyy\)\s*\n(\d{1,2} \w{3} 20\d\d)',t)
    val=re.search(r"^User\s+entered\s+(?:'([^']+)'|(empty))",u,re.M);form='Form: Adverse Events ('+str(int(r.AESPID.split('-')[1]))+')'
    ok=bool(term and st and val and val[1]=='Yes (Y)' and form in u and form in t and r.EUA_Bates in t and r.subject in t and r.subject in u and normal(term[1])==normal(r.eua_term) and date(st[1])==r.eua_onset[:10] and 'Is the adverse event serious?\nAudit' in u)
    fatal.append([r.subject,r.AESPID,n,pg,r.EUA_Bates,sp,r.eua_term,r.eua_onset,val[1] if val else '',ok])
 save(a.out,'reviewed_source_page_checks','finding,subject,filename,pdf_page,FDA_Bates,verified',notes)
 save(a.out,'serious_source_page_checks','subject,AESPID,filename,form_page,FDA_Bates,seriousness_page,EUA_term,EUA_onset,source_seriousness,verified',fatal)
 save(a.out,'missing_documents','filename',[[n] for n in missing])
 if not all(x[-1] for x in notes+fatal):raise SystemExit('At least one source check failed. Inspect output; do not treat a missing/mismatching source as corroboration.')
 print(f'{len(notes)} page references for reviewed cases and {len(fatal)} serious EUA forms verified; {len(missing)} documents missing.')
if __name__=='__main__':main()

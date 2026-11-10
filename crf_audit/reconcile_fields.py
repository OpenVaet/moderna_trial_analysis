"""Join newly extracted CRF observations to the named release comparators; never read archived findings to generate observations."""
from pathlib import Path
import argparse, collections, csv, datetime as dt, gzip, hashlib, json, re
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def sha(f):
    h=hashlib.sha256()
    with f.open('rb') as k:
        for block in iter(lambda:k.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def date(s):
    for fmt in ('%d %b %Y','%Y-%m-%d'):
        try:return dt.datetime.strptime(s,fmt).date().isoformat()
        except ValueError:pass
    return ''
def normal(t):return ''.join(re.findall(r'[A-Z0-9]',str(t).upper()))
def yes(t):return str(t).upper() in ('Y','TRUE','1')
def save(base,f,cols,rows):
    with (base/(f+'.csv')).open('w',encoding='utf8',newline='') as h:
        w=csv.writer(h);w.writerow(cols.split(','));w.writerows(rows)
def read(f):
    return pd.read_sas(f,encoding='latin1').fillna('') if f.suffix.lower()=='.xpt' else pd.read_csv(f,dtype=str,keep_default_na=False)
def loadjson(path):return json.loads(path.read_text(encoding='utf8'))
def exp_arm(val):return {'PLACEBO':'Placebo','MRNA-1273':'mRNA-1273'}.get(val.upper(),'')
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('scan','out','ex','ae-s5','ae-s3','ae-s27'):p.add_argument('--'+k,type=Path,required=True)
    compact=ROOT/'evidence/crf_source/comparator_rows'; compact_hash=loadjson(compact/'SHA256.json')
    for k in ('adsl','pcr_eua','pcr_prod'):p.add_argument('--'+k.replace('_','-'),type=Path,default=compact/(k+'.csv'))
    p.add_argument('--equivalent-inputs',action='store_true',help='Explicit opt-in to re-exported CSV/alternate equivalent bytes; release identity must be established by the user')
    a=p.parse_args();out=a.out
    if out.exists() and any(out.iterdir()):p.error('Use a fresh output directory')
    cov=loadjson(a.scan/'FINISHED.json');ms=ROOT/'evidence/corpus_manifest.csv'
    if sha(ms)!=cov['manifest_sha256']:p.error('Mismatch between scan and corpus manifest')
    manifest=read(ms);manifest['doc']=manifest.name.str[:-4];manifest=manifest.set_index('doc')
    cmp=read(ROOT/'evidence/crf_source/comparators.csv').set_index('file').checksum
    raw={'adsl':a.adsl,'ex':a.ex,'ae_s5':a.ae_s5,'ae_s3':a.ae_s3,'ae_s27':a.ae_s27,'pcr_eua':a.pcr_eua,'pcr_prod':a.pcr_prod}
    expected={'adsl':'BLA_ADSL_selected.csv','pcr_eua':'EUA/pcr_record_audit.csv','pcr_prod':'Prod/pcr_record_audit.csv'}
    for nm,snap,dom in [('ex','S3','ex'),('ae_s5','S5','ae'),('ae_s3','S3','ae'),('ae_s27','S27','ae')]:
        expected[nm]=next(n for n in cmp.index if '_'+snap+'_' in n and n.endswith('_'+dom+'.xpt'))
    inputs=[]
    for k,f in raw.items():
        checksum=sha(f);matched=checksum in (cmp.loc[expected[k]],compact_hash.get(k,''))
        if not matched and not a.equivalent_inputs:p.error('Input checksum mismatch: '+str(f)+'; see --equivalent-inputs')
        inputs.append([k,f.name,checksum,matched])
    out.mkdir(parents=True,exist_ok=True);save(out,'input_checksums','role,filename,sha256,reference_bytes_match',inputs)
    ad=read(a.adsl).set_index('SUBJID',verify_integrity=True)
    ex=read(a.ex);ex=ex[(ex.EPOCH=='BLINDED TREATMENT') & pd.to_numeric(ex.EXDOSE,errors='coerce').gt(0) & ex.EXTRT.map(exp_arm).ne('')]
    ex['sid']=ex.USUBJID.str.replace('-','',regex=False).str[-9:];xp={}
    for key,g in ex.groupby(['sid','VISIT']):
        v=set(g.EXTRT.map(exp_arm));xp[key]=v.pop() if len(v)==1 else 'Ambiguous'
    ae={}
    for sn,f in [('S5',a.ae_s5),('S3',a.ae_s3),('S27',a.ae_s27)]:
        b=read(f);b['sid']=b.USUBJID.str.replace('-','',regex=False).str[-9:];b['n']=b.AESPID.str.extract(r'^AE-(\d+)$')[0];b=b[b.n.notna()];b['n']=b.n.astype(int);ae[sn]=b.set_index(['sid','n'],verify_integrity=True)
    lab={}
    for sn,f in [('S5',a.pcr_eua),('S3',a.pcr_prod)]:
        b=read(f);usable=b.usable_dated.map(yes)&b.binary.map(yes)
        for x in b[usable].itertuples():
            key=(sn,x.subject,x.date);v=lab.setdefault(key,set())
            if yes(x.positive):v.add('POSITIVE')
            if yes(x.negative):v.add('NEGATIVE')
    arms=[];flags=[];events=[];pcr=[];contacts=[];leads=[];details=[];rejects=[];d2audit=[]
    docs=sorted((a.scan/'front').glob('*.json'));assert len(docs)==cov['read']
    for path in docs:
        doc=path.stem;m=manifest.loc[doc];sid=m.subject;seq=m.submission;sn='S5' if seq=='S6' else 'S3';cut='2020-11-25' if sn=='S5' else '2021-03-26'
        d=loadjson(path);src=loadjson(a.scan/'source'/(doc+'.json'))
        with gzip.open(a.scan/'hist'/(doc+'.gz'),'rt',encoding='utf8') as f:hist={(folder,form,k):v for folder,form,k,v in json.load(f)}
        for pg,bt,snip in src['illness_candidates']:leads.append([doc,pg]);details.append([doc,sid,seq,pg,bt,snip])
        for pg,bt,workaround in src['dose2_audit']:d2audit.append([doc,sid,seq,pg,bt,workaround])
        for pg,folder,bt,assigned,d1,d2 in d['arm']:
            if sid not in ad.index:raise ValueError('ADSL missing participant: '+sid)
            s=ad.loc[sid];e1=xp.get((sid,'Visit 1 Day 1'),'Not Administered');e2=xp.get((sid,'Visit 2 Day 29'),'Not Administered')
            ar=[doc,sid,seq,pg,bt,assigned,d1,d2,s.TRT01P,s.DOSE1FL,s.DOSE2FL,e1,e2];arms.append(ar)
            if d2 in ('Placebo','mRNA-1273') and s.DOSE2FL=='N' and e2=='Not Administered':
                if len(set(ex.loc[ex.sid.eq(sid),'VISIT']))!=1 or len(ex.loc[ex.sid.eq(sid)])!=1:raise ValueError('Unexpected blinded administrations for '+sid)
                au=src['dose2_audit'];day=src['day29_given'];flags.append(ar+[au[0][0] if au else '',day[0][2] if day else 'No extracted Yes/No entry',day[0][0] if day else '',day[0][1] if day else ''])
        for pg,bt,form,term,start,_ in d['ae']:
            start=date(start);n=int(re.search(r'\((\d+)\)',form)[1])
            if not start or start>cut or not term:rejects.append([doc,pg,'AE','blank term, incomplete date or after cutoff']);continue
            for table in (['S5'] if sn=='S5' else ['S3','S27']):
                z=ae[table];found=(sid,n) in z.index;x=z.loc[(sid,n)] if found else {};t=x.get('AETERM','');st=x.get('AESTDTC','')
                events.append([doc,sid,seq,pg,bt,n,term,start,x.get('AESEQ',''),t,st,x.get('AESER',''),x.get('AEOUT',''),table,found,term==t,start==st[:10],True,normal(term)==normal(t)])
        for (fd,fo,k),v in hist.items():
            if not v:continue
            if fo=='COVID-19 Contact' and k=='contactdate':
                # Timeline is printed newest to oldest. A correction can make the earliest contact differ from the current field.
                valid=[(date(val),ts[:11].strip(),pg) for pg,val,ts in v if date(val) and ts]
                if not valid:continue
                old,_,_=valid[-1];oldraw=next(z for z in reversed(v) if date(z[1]) and z[2]);oldp=oldraw[0]
                first=date(oldraw[2].rsplit(' ',1)[0]);current=v[0][1]
                if not first:continue
                lag=(dt.date.fromisoformat(first)-dt.date.fromisoformat(old)).days
                contacts.append([doc,sid,seq,fd,oldp,old,first,lag,current,date(current)==old and lag>=0])
            if fo=='COVID Diagnostic Test' and k=='positive':
                tv=hist.get((fd,fo,'testdate'),[]);dd=date(tv[0][1]) if tv else '';val=v[0][1];bin={'Yes (Y)':'POSITIVE','No (N)':'NEGATIVE'}.get(val)
                if not dd or dd>cut or not bin:rejects.append([doc,v[0][0],'PCR','nonbinary, incomplete date or after cutoff']);continue
                res=lab.get((sn,sid,dd),set());status=next(iter(res)) if len(res)==1 else 'MIXED' if res else 'NO SAME-DATE BINARY'
                pcr.append([doc,sid,seq,fd,v[0][0],dd,bin,status,bin==status])
    ac='doc,sid,seq,pg,bates,assigned,dose1,dose2,rand,DOSE1FL,DOSE2FL,EX1,EX2'
    for name,cols,rows in [('arm_checkbox',ac,arms),('dose2_flags',ac+',audit_unblind_dose_pg,exposure_crf_response,exposure_crf_audit_page,exposure_bates',flags),('AE_term_onset_checks','doc,sid,seq,pg,bates,form_no,crf_term,crf_start,AESEQ,AETERM,AESTDTC,AESER,AEOUT,snapshot,present,term_equal,start_equal,onset_through_cutoff,term_equal_normalized',events),('pcr_comparison','doc,sid,seq,folder,positive_field_pg,test_date,CRF_result,same_date_MB,agree',pcr),('contact_lags','doc,sid,seq,folder,pg,old_contact,first_recorded_entry_gmt_date,entry_lag,current,eligible',contacts),('illness_visit_search_hits','doc,pg',leads),('illness_candidates_with_context','doc,sid,seq,pg,bates,local_text_excerpt',details),('source_exclusions','doc,pg,domain,reason',rejects),('dose2_field_history_locators','doc,sid,seq,pg,bates,explicit_page_population_phrase',d2audit)]:save(out,name,cols,rows)
    (out/'coverage.json').write_text(json.dumps(cov,indent=2),encoding='utf8')
    print(f"{len(docs)} exports; {len(arms)} assignment forms; {len(flags)} second-dose discrepancies; {len(events)} AE comparisons; {len(pcr)} PCR comparisons; {len(leads)} illness candidate pages. Coverage: "+('PARTIAL' if cov['partial'] else 'complete manifest'))
if __name__=='__main__': main()

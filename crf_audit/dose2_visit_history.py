"""Reopen original PDFs for visits/dates and field-entry history in the readable one-dose cohort. Run after reconcile_fields.py."""
import argparse, csv, datetime as dt, hashlib, re, zlib
from pathlib import Path
import fitz
import pandas as pd
from reconcile_fields import ROOT,read,save,loadjson,sha
STAMP=r'\d{1,2} [A-Za-z]{3} 20\d\d \d\d:\d\d:\d\d'
def get(p,t):
 m=re.search(p,t,re.M);return m[1].strip() if m else ''
def day(t,f):
 v=get(re.escape(f)+r'\s*\n(\d{1,2} [A-Za-z]{3} 20\d\d)',t);return dt.datetime.strptime(v,'%d %b %Y').date().isoformat() if v else ''
def box(p,title,opts):
 loc=p.search_for(title)
 if len(loc)!=1:return ''
 y=loc[0].y0;r=[j['rect'] for j in p.get_drawings() if j['type'] in ('f','fs') and j['fill'] and max(j['fill'])<.65 and j['rect'].x0>490 and 7<=j['rect'].width<=17 and 7<=j['rect'].height<=17]
 hit=[k for k in opts for a in p.search_for(k) if a.x0>340 and y-2<=a.y0<=y+16*(len(opts)-1)+5 and any(abs((a.y0+a.y1)/2+3.5-(b.y0+b.y1)/2)<3 for b in r)]
 return hit[0] if len(hit)==1 else ''
def hist(t):
 vals=[]
 for m in re.finditer(r"^(Amendment Manager:\s*)?User\s+entered\s+(?:'([^']{1,170})'|(empty)[.;])(.{0,420}?)("+STAMP+')',t,re.M|re.S):
  if m[1] or re.search(r'\bSystem\s*$',m[4]) or m[3]:continue
  if re.search(r'\n(?:Query |User (?:entered|opened|closed|signature)|DataPoint |Comment added|Data hard locked)',m[4]):raise ValueError('Uncertain timestamp boundary')
  vals.append((dt.datetime.strptime(m[5],'%d %b %Y %H:%M:%S').isoformat(' '),' '.join(m[2].split())))
 return vals
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('root','joins','ex','out'):p.add_argument('--'+k,required=True,type=Path)
 a=p.parse_args();ids=read(a.joins/'arm_checkbox.csv');ids=ids[(ids.DOSE2FL=='N')&(ids.dose2!='') & (ids.EX2=='Not Administered')]
 if ids.sid.duplicated().any():raise ValueError('Duplicate participants; adjudicate version differences before counting')
 files=read(ROOT/'evidence/corpus_manifest.csv').set_index('name')
 ebase=ROOT/'evidence/crf_source';ad=read(ebase/'comparator_rows/adsl.csv').set_index('SUBJID')
 xm=read(ebase/'comparators.csv');allowed=xm[xm.file.str.contains('_S3_')&xm.file.str.endswith('_ex.xpt')].checksum.item()
 if sha(a.ex)!=allowed:raise ValueError('Unexpected S3 EX bytes')
 ex=read(a.ex);ex=ex[ex.EPOCH.str.contains('OPEN LABEL') & ex.EXTRT.eq('mRNA-1273') & pd.to_numeric(ex.EXDOSE,errors='coerce').gt(0)]
 paths={}
 for f in a.root.rglob('*'):
  if f.is_file() and f.suffix.lower()=='.pdf':
   m=re.search(r'((?:125752|27073)_s\d+_m5_crf_mrna-1273-p301-us\d{7}\.pdf)$',f.name,re.I)
   if m:paths.setdefault(m[1].lower(),[]).append(f)
 rows=[]
 for x in ids.itertuples():
  name=x.doc+'.pdf';src=files.loc[name];orig=next((f for f in paths[name.lower()] if f.stat().st_size==int(src.us) and '%08x'%zlib.crc32(f.read_bytes())==src.crc32_hex),None)
  if not orig:raise ValueError(name)
  with fitz.open(orig) as d:
   pg=int(x.pg);tf=d[pg-1].get_text();folder=get(r'^Folder: (.+)',tf);fb=folder.replace('(','\\(').replace(')','\\)').encode()
   field=[];assign=[];vp='';when='';loc='';workaround=False
   for i,page in enumerate(d):
    raw=page.read_contents()
    if fb not in raw or (b'(Form: Unblinding)' not in raw and b'(Form: Visit Date)' not in raw):continue
    t=page.get_text()
    if get(r'^Folder: (.+)',t)!=folder:continue
    if 'Form: Visit Date\n' in t and '\nAudit\n' not in t:
     vp=i+1;when=day(t,'Visit date (dd MMM yyyy)');loc=box(page,"Was visit performed at the participant's home or at the clinic?",['Home','Clinic'])
    if '\nActual Dose 2\nAudit\n' in t:field.extend((ts,v,i+1,get(r'(FDA-CBER-[\d-]+)',t)) for ts,v in hist(t));workaround|='in order to populate' in ' '.join(t.lower().split())
    if '\nParticipant randomization assignment\nAudit\n' in t:assign+=hist(t)
   target=sorted(i for i in field if i[1].split(' (')[0]==x.dose2); aa=sorted(k for k in assign if k[1].split(' (')[0]==x.rand)
   if not target:raise ValueError(name+': cannot establish non-system entry for displayed value')
   fi=sorted(field)[0];end=target[0];ass=aa[0][0] if aa else '';ub=day(tf,'Date of unblinding (dd MMM yyyy)')
   adm=ex[ex.USUBJID.str.replace('-','',regex=False).str.endswith(x.sid)].sort_values('EXSTDTC');a1=adm.EXSTDTC.iloc[0][:10] if len(adm) else '';dsub=ad.loc[x.sid]
   sg=get(r'Data signed:\s*\n\s*('+STAMP+')',tf);sg=dt.datetime.strptime(sg,'%d %b %Y %H:%M:%S').isoformat(' ') if sg else ''
   rows.append([x.sid,x.rand,x.dose2!='Not Administered',x.dose2,folder,loc,when,ub,a1,len(adm),end[0],fi[0],fi[1],ass,bool(ass and ass==end[0]),sg,dsub.PPROTFL,dsub.SAFFL,dsub.FASFL,dsub.PPEREAS,workaround,name,pg,vp,end[2],end[3]])
 cols='subject,arm,discrepant_final_D2,final_D2,folder,decision_visit_location_as_checked,decision_visit_date,CRF_unblinding_date,first_open_label_mRNA_date,recorded_OL_doses,first_non_system_entry_of_final_D2_GMT,first_nonblank_field_GMT,first_nonblank_field_value,first_assigned_arm_GMT,assigned_arm_same_second,final_form_signature_GMT,PPROTFL,SAFFL,FASFL,PPEREAS,US331_explicit_workaround,filename,CRF_page,visit_page,audit_field_page,audit_field_Bates'
 if a.out.exists() and any(a.out.iterdir()):p.error('Output directory must be new or empty.')
 a.out.mkdir(parents=True,exist_ok=True);save(a.out,'dose2_visit_history_rebuilt',cols,rows)
 new=read(a.out/'dose2_visit_history_rebuilt.csv').set_index('subject')
 ref=read(ebase/'Moderna_171_dose2_visit_comparison.csv').set_index('subject')
 scanned=set(read(a.joins.parent/'scan/coverage.csv').query("status=='selected'").file) if (a.joins.parent/'scan/coverage.csv').exists() else set(read(a.joins/'arm_checkbox.csv').doc+'.pdf')
 expected=set(ref[ref.filename.isin(scanned)].index)
 if expected!=set(new.index):raise ValueError('Subject sets differ')
 bad=[]
 for sid,row in new.iterrows():
  for col in new.columns:
   a1=row[col];b=ref.loc[sid,col]
   if col in ['visit_page','audit_field_page','CRF_page']:a1=a1.split('.')[0];b=b.split('.')[0]
   if a1!=b:bad.append([sid,col,a1,b])
 save(a.out,'workflow_differences','subject,field,new,reference',bad)
 if bad:raise SystemExit('Dose workflow differs: see workflow_differences.csv')
 print(f'{len(rows)} participants verified against the published workflow fields.')
if __name__=='__main__':main()

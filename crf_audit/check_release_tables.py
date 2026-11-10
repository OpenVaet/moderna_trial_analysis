"""Check published CRF rows against original S5/S3/S27 AE and S3 EX, plus independent PCR comparator rows.
This verifies the study-table side of the joins; it is not another PDF extraction."""
import argparse, collections
from pathlib import Path
import pandas as pd
from reconcile_fields import ROOT,read,sha,yes,normal,exp_arm
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('ex','ae-s5','ae-s3','ae-s27'):p.add_argument('--'+k,type=Path,required=True)
 a=p.parse_args();e=ROOT/'evidence';fp=read(e/'crf_source/comparators.csv').set_index('file')
 def checked(f,sn,domain):
  match=next(n for n in fp.index if '_'+sn+'_' in n and n.endswith('_'+domain+'.xpt'))
  if sha(f)!=fp.loc[match,'checksum']:raise ValueError('Unexpected input '+str(f))
  g=read(f);g['sid']=g.USUBJID.str.replace('-','',regex=False).str[-9:];return g
 ae={sn:checked(f,sn,'ae').set_index(['sid','AESPID'],verify_integrity=True) for sn,f in [('S5',a.ae_s5),('S3',a.ae_s3),('S27',a.ae_s27)]}
 obs=read(e/'AE_term_onset_checks.csv')
 for r in obs.itertuples():
  x=ae[r.snapshot].loc[(r.sid,'AE-%03d'%int(r.form_no))]
  assert normal(r.crf_term)==normal(x.AETERM) and r.crf_start==x.AESTDTC[:10] and r.AESER==x.AESER and str(float(r.AESEQ))==str(float(x.AESEQ)),(r.doc,r.pg)
 assert obs.groupby('snapshot').size().to_dict()=={'S5':1118,'S3':4362,'S27':4362}
 # Reconstruct the 566 keys from whole tables, not the previous selected 287 people.
 old=ae['S5'].reset_index();new=ae['S3'].reset_index();am=ae['S27'].reset_index();assert set(zip(new.sid,new.AESPID))==set(zip(am.sid,am.AESPID))
 new['n']=new.AESPID.str[3:].astype(int);keyset=set(zip(new.sid,new.AESPID))
 bounds=new.groupby('sid').n.agg(['min','max']).to_dict('index')
 def gap(r):
  b=bounds.get(r.sid);return b and b['min']<int(r.AESPID[3:])<b['max'] and (r.sid,r.AESPID) not in keyset
 gaps={r.sid+':'+r.AESPID for r in old.itertuples() if gap(r)};reg=read(e/'form_dispositions_566.csv');assert len(gaps)==566 and set(reg.key)==gaps
 assert set(read(e/'serious_12.csv').key)==set(reg.loc[reg.EUA_serious=='Y','key'])
 x=read(e/'arm_checkbox.csv');roster=read(e/'crf_source/comparator_rows/adsl.csv').set_index('SUBJID',verify_integrity=True)
 ex=checked(a.ex,'S3','ex');ex=ex[(ex.EPOCH=='BLINDED TREATMENT')&pd.to_numeric(ex.EXDOSE,errors='coerce').gt(0)&ex.EXTRT.map(exp_arm).ne('')]
 ix={(s,v):set(g.EXTRT.map(exp_arm)) for (s,v),g in ex.groupby(['sid','VISIT'])}
 for r in x.itertuples():
  sub=roster.loc[r.sid];assert (r.rand,r.DOSE1FL,r.DOSE2FL)==(sub.TRT01P,sub.DOSE1FL,sub.DOSE2FL)
  for v,expected in [('Visit 1 Day 1',r.EX1),('Visit 2 Day 29',r.EX2)]:
   assert ix.get((r.sid,v),{'Not Administered'})=={expected}
 def pcrlookup(name):
  q=read(e/'crf_source/comparator_rows'/name); g=q[q.usable_dated.map(yes)&q.binary.map(yes)]
  res=collections.defaultdict(set)
  for r in g.itertuples():
   if yes(r.positive):res[r.subject,r.date].add('POSITIVE')
   if yes(r.negative):res[r.subject,r.date].add('NEGATIVE')
  return q,res
 prod,P=pcrlookup('pcr_prod.csv');_,E=pcrlookup('pcr_eua.csv');q=read(e/'pcr_comparison.csv')
 for r in q.itertuples():
  z=(E if r.seq=='S6' else P).get((r.sid,r.test_date),set())
  assert r.same_date_MB==(next(iter(z)) if len(z)==1 else 'MIXED' if z else 'NO SAME-DATE BINARY')
 # This 55-row exhibit is a focused selection, not all tests in these participants.
 for r in read(e/'case_PCR_records.csv').itertuples():
  y=prod[(prod.subject==r.subject)&(prod.MBSEQ==r.MBSEQ)&(prod.date==r.date)]
  assert any(t.raw_result==r.raw_result and t.route==r.route and t.visits==r.visits and yes(t.usable_dated)==yes(r.usable_dated) for t in y.itertuples()),(r.subject,r.MBSEQ)
 print('Original XPTs verified: 1,118 + 4,362 + 4,362 AE transcription comparisons, 566 original form keys, and exposure/ADSL values on all 1,071 unblinding form observations. PCR comparator: 640 dated comparisons, 55 focused rows verified.')
if __name__=='__main__':main()

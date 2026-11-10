"""Rebuild the focused deviation exhibits for the published, manually selected illness cases.
Run prepare_source_columns.R for the linked SUPPDV decision fields first."""
import argparse
from pathlib import Path
from reconcile_fields import ROOT,read,sha
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ['dv','addv','suppdv-csv','out']:p.add_argument('--'+n,required=True,type=Path)
 a=p.parse_args();expected=read(ROOT/'evidence/crf_source/comparators.csv')
 for n,f in [('S_dv.xpt',a.dv),('A_D_addv.xpt',a.addv)]:
  if sha(f)!=expected[expected.file.str.endswith(n)].checksum.item():raise ValueError('Unexpected '+n+' snapshot')
 dv=read(a.dv);ad=read(a.addv);sup=read(a.suppdv_csv)
 if a.out.exists() and any(a.out.iterdir()):p.error('Use an empty output directory.')
 a.out.mkdir(exist_ok=True,parents=True)
 sub=dv[dv.SUBJID.isin(['US3932018','US3932147','US3632127'])][['SUBJID','DVSEQ','DVSTDTC','DVTERM','DVSPID']].sort_values(['SUBJID','DVSEQ'])
 c=ad[ad.SUBJID=='US3932018'][['SUBJID','DVSEQ','DVSTDTC','DVTERM','DVTERM1','DVSIG']].copy()
 for k in ['DECSRTNL','DECSRTN1','DECSRTN2','DECSRTN3']:
  ss=sup[sup.QNAM==k];kv={float(r.IDVARVAL):r.QVAL for r in ss.itertuples()}
  if len(kv)!=len(ss):raise ValueError('Repeated qualifier '+k)
  c[k]=c.DVSEQ.map(kv).fillna('')
 for n,x in [('illness_DV',sub),('dv_links',c)]:
  f=a.out/(n+'.csv');x.to_csv(f,index=False)
  cols=list(x.columns);k=['SUBJID','DVSEQ'];ref=read(ROOT/'evidence'/f.name).sort_values(k).reset_index(drop=True)
  if not read(f).sort_values(k).reset_index(drop=True).equals(ref):raise SystemExit('Source comparison differs for '+n+'; output retained for review')
  print(n,len(x),'rows verified')
if __name__=='__main__':main()

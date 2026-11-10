"""Make compact, independent study-table inputs; never read reported CRF results.
Use the complete selected ADSL CSV (or export the required columns from licensing ADSL)
and both pcr_record_audit.csv files produced by analyse_pcr_testing_by_lab.R."""
import argparse,csv,hashlib
from pathlib import Path
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('adsl','pcr-eua','pcr-prod','out'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();m=pd.read_csv(BASE/'evidence/corpus_manifest.csv',dtype=str);ids=set(m.loc[m['class']=='full','subject']); a.out.mkdir(exist_ok=False,parents=True)
 for name,f,cols,idcol in [('adsl',a.adsl,['SUBJID','TRT01P','DOSE1FL','DOSE2FL','PPROTFL','SAFFL','FASFL','PPEREAS','UBLDDT','PDVDT'],'SUBJID'),('pcr_eua',a.pcr_eua,['subject','date','MBSEQ','visits','raw_result','route','usable_dated','binary','positive','negative','illness_labelled','baseline'],'subject'),('pcr_prod',a.pcr_prod,['subject','date','MBSEQ','visits','raw_result','route','usable_dated','binary','positive','negative','illness_labelled','baseline'],'subject')]:
  x=pd.read_csv(f,dtype=str,keep_default_na=False,usecols=cols);x=x[x[idcol].isin(ids)].copy()
  if name=='adsl':assert x.SUBJID.is_unique and set(x.SUBJID)==ids
  x=x.sort_values([idcol]+(['date','MBSEQ'] if name!='adsl' else []),kind='stable')[cols]
  x.to_csv(a.out/(name+'.csv'),index=False,lineterminator='\n')
if __name__=='__main__':main()

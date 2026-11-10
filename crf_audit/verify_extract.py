"""Compare a new full or partial run to the cited ledgers. Exit nonzero on differences.
Select by the scanned document manifest, never by the set of successful finding rows."""
import argparse,csv,collections,decimal,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
NUM={'pg','AESEQ','form_no','positive_field_pg','audit_unblind_dose_pg','exposure_crf_audit_page','entry_lag'}
NAMES=('arm_checkbox','dose2_flags','contact_lags','AE_term_onset_checks','pcr_comparison','illness_visit_search_hits')
def read(path):
 with path.open(encoding='utf-8-sig',newline='') as h:
  rd=csv.DictReader(h);return rd.fieldnames,list(rd)
def canonical(k,v):
 if k in NUM and v:return str(decimal.Decimal(v).normalize())
 if v.upper() in ('TRUE','FALSE'):return v.upper()
 return v
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--scan',required=True,type=Path);p.add_argument('--joins',required=True,type=Path);a=p.parse_args()
 info=json.loads((a.scan/'FINISHED.json').read_text());assert json.loads((a.joins/'coverage.json').read_text())==info
 docs={Path(row['file']).stem for row in read(a.scan/'coverage.csv')[1] if row['status']=='selected'};assert len(docs)==info['read']
 result=[];fails=[]
 for n in NAMES:
  f=n+'.csv';ca,aa=read(ROOT/'evidence'/f);cb,bb=read(a.joins/f)
  if ca!=cb:raise ValueError(f+' schema mismatch')
  A=collections.Counter(tuple(canonical(k,r[k]) for k in ca) for r in aa if r['doc'] in docs)
  B=collections.Counter(tuple(canonical(k,r[k]) for k in ca) for r in bb)
  missing,extra=A-B,B-A;result.append((n,sum(A.values()),sum(B.values()),sum(missing.values()),sum(extra.values())))
  for kind,dd in [('reference_only',missing),('new_only',extra)]:
   for row,count in dd.items():fails.append([n,kind,count,json.dumps(dict(zip(ca,row)))])
 for fn,cols,data in [('comparison_results.csv',['ledger','reference_in_scanned_documents','new_rows','missing','extra'],result),('comparison_differences.csv',['ledger','difference','multiplicity','row'],fails)]:
  with (a.joins/fn).open('w',newline='',encoding='utf8') as h:w=csv.writer(h);w.writerow(cols);w.writerows(data)
 print(f"{len(docs)}/{info['expected']} PDFs; "+'; '.join(f'{n}: {num} rows, -{miss}/+{extra}' for n,_,num,miss,extra in result))
 if fails:raise SystemExit('Differences found; review comparison_differences.csv before using new counts.')
if __name__=='__main__': main()

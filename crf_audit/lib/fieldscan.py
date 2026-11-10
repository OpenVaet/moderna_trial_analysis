from pathlib import Path
import gzip,json,re,sys,time,concurrent.futures as cf,collections,os
from .scan import wanted
w=Path('.');out=w/'hist'
def go(f):
 d=collections.defaultdict(list);pages=gzip.decompress(f.read_bytes()).decode().split('\f')
 for i,s in enumerate(pages,1):
  if '\nAudit' not in s or 'User entered' not in s:continue
  a=re.search(r'([^\n]+)\nAudit\s*\nUser\s*\nTime \(GMT\)',s)
  if not a:continue
  typ=wanted.get(a[1].strip())
  if not typ:continue
  fd=re.search(r'^Folder: (.*)$',s,re.M);fm=re.search(r'^Form: (.*)$',s,re.M)
  if not fd or not fm:continue
  fo=fm[1].strip();folder=fd[1].strip()
  if not (re.fullmatch(r'Adverse Events \(\d+\)',fo) or fo in ('Unblinding','COVID Diagnostic Test','COVID-19 Contact','End of Study')):continue
  body=s[a.end():];ms=list(re.finditer(r"User\s+entered\s+(?:'([^']{1,300})'|(empty)\.)",body,re.I))
  values=d[(folder,fo,typ)]
  for m in ms:
   val=' '.join((m[1] or '<EMPTY>').split())
   # The displayed GMT audit-action time must follow this entry before the next action.
   tail=re.split(r"\n(?:User (?:entered|opened|closed|coded|signature)|Query '|DataPoint |Data hard locked|Data entry locked|External Audit Record)",body[m.end():],maxsplit=1)[0]
   date=re.search(r'\b(\d{1,2} [A-Za-z]{3} 20\d\d \d\d:\d\d:\d\d)\b',tail)
   values.append([i,val,date[1] if date else ''])
 for k in list(d):
  if not d[k]:del d[k]
 dest=out/(f.stem+'.gz');tmp=dest.with_suffix('.part')
 with gzip.open(tmp,'wt',compresslevel=1) as fp:json.dump([[*k,v] for k,v in d.items()],fp,separators=(',',':'))
 tmp.replace(dest)

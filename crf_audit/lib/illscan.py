from pathlib import Path
import re,gzip,concurrent.futures as cf,csv
rx=re.compile(r"(?:not (?:brought|called|seen|brought back)|(?:did not|failed to|didn.t|unable to) (?:bring|conduct|schedule|come in|complete)|inadvertently (?:missed|not)|(?:missed|no) (?:opportunity for|opportunity to)).{0,70}?(?:illness visit|illness assessment)|(?:illness visit|illness assessment).{0,32}?(?:was not (?:done|completed|performed|conducted)|not (?:done|completed|performed)|inadvertently missed)",re.I)
def look(p):
 out=[]
 for i,s in enumerate(gzip.decompress(p.read_bytes()).decode().split('\f'),1):
  if not re.search('illness (?:visit|assessment)',s,re.I):continue
  t=' '.join(s.split())
  for m in list(rx.finditer(t))[:1]:out.append((p.stem,i,t[max(0,m.start()-100):m.end()+80]))
 return out

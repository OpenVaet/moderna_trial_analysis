from pathlib import Path
import re,gzip,json,collections,sys,time,concurrent.futures as cf,fitz
w=Path('.');out=w/'front'
def day(t,name):
 mm=re.search(name+r'\s*\(dd MMM yyyy\)\s*\n(\d{1,2} [A-Za-z]{3} 20\d\d)',t,re.I);return mm[1] if mm else ''
def armread(pg,t,ttl):
 rect=pg.search_for(ttl); 
 if len(rect)!=1:return ''
 r=rect[0];prods=[]
 for key in ['mRNA-1273','Placebo','Not Administered']:
  for x in pg.search_for(key):
   if x.x0>340 and r.y0-2<=x.y0<=r.y0+(21 if ttl=='Participant randomization assignment' else 37):
    prods.append((key,x))
 fills=[z['rect'] for z in pg.get_drawings() if z['type'] in ['f','fs'] and z['fill'] and max(z['fill'])<.65 and z['rect'].x0>490 and 7<=z['rect'].width<=17 and 7<=z['rect'].height<=17]
 checked=[k for k,x in prods if any(abs((x.y0+x.y1)/2+3.5-(b.y0+b.y1)/2)<3 for b in fills)]
 return checked[0] if len(checked)==1 else ''
def parse(f, source=None):
 ss=gzip.decompress(f.read_bytes()).decode().split('\f');d=None;ae=[];arm=[];other=[]
 for i,s in enumerate(ss,1):
  if '\nAudit\n' in s or not s:continue
  mini=' '.join(s[:460].split())
  if not any(k in mini for k in ['Form: Adverse Events (','Form: Unblinding','Form: COVID Diagnostic Test','Form: COVID-19 Contact']):continue
  if not ('Was this a medically-attended AE?' in s or 'Is the adverse event serious?' in s or 'Actual Dose 2' in s or 'Did Subject Test Positive' in s or 'Has the subject reported' in s):continue
  if d is None:d=fitz.open(source or w/'pdf'/(f.stem+'.pdf'))
  pg=d[i-1];t=pg.get_text();fo=re.search(r'^Form: (.*)$',t,re.M);fd=re.search(r'^Folder: (.*)$',t,re.M)
  if not fo or not fd:continue
  b=re.search(r'FDA-CBER-[\d-]+',t);b=b[0] if b else ''
  if fo[1]=='Unblinding' and 'Actual Dose 2' in t:
   arm.append([i,fd[1],b,*[armread(pg,t,k) for k in ['Participant randomization assignment','Actual Dose 1','Actual Dose 2']]])
  elif fo[1].startswith('Adverse Events (') and 'Start date' in t:
   m=re.search(r'\nAdverse event\s*\n(.{0,600}?)\nWas this a medically-attended AE\?',t,re.S)
   if m:ae.append([i,b,fo[1],' '.join(m[1].split()),day(t,'Start date'),day(t,'If not Ongoing, end date')])
  elif fo[1] in ['COVID Diagnostic Test','COVID-19 Contact']:other.append([i,fd[1],fo[1],b])
 if d:d.close()
 dest=out/(f.stem+'.json');tmp=dest.with_suffix('.tmp');tmp.write_text(json.dumps(dict(ae=ae,arm=arm,other=other)));tmp.replace(dest)

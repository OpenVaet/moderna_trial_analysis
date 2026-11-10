from pathlib import Path
import re,fitz,zlib,gzip,time,concurrent.futures as cf,sys,subprocess,os,json
w=Path('.')
lit=re.compile(rb'\(((?:\\.|[^()\\])*)\)')
esc=re.compile(rb'\\([0-7]{1,3}|[()\\nrtbf])')
m={b'n':b'\n',b'r':b'\r',b't':b'\t',b'b':b'\b',b'f':b'\f'}
def trans(e):
 a=e[1];return bytes([int(a,8)%256]) if a[:1].isdigit() else m.get(a,a)
def extract(f):
 out=w/'text'/(f.stem+'.gz');d=fitz.open(f);tmp=out.with_suffix('.tmp')
 fallback=blank=img=0;chars=0;size=len(d)
 with gzip.open(tmp,'wt',encoding='utf8',compresslevel=2) as g:
  for i,pg in enumerate(d):
   raw=pg.read_contents()
   fields=lit.findall(raw)
   # fast path only on plain Rave printable content with native subject/form labels and Bates.
   if not(any(v.startswith(b'Form:') for v in fields) and any(v.startswith(b'FDA-CBER-') for v in fields) and b'\x00' not in b''.join(fields) and not re.search(rb'<[a-fA-F\d]{6,}>\s*Tj',raw)):
    s=pg.get_text();fallback+=1
   else:s=b'\n'.join(esc.sub(trans,v) if b'\\' in v else v for v in fields).decode('cp1252','replace')
   if len(re.sub('FDA-CBER-[0-9-]+|\\s+','',s))<30:
    blank+=1
    if pg.get_images():img+=1
   chars+=len(s);g.write(s.replace('\f',' ')+'\f')
 tmp.replace(out);(w/'keys'/(f.stem+'.json')).write_text(json.dumps(dict(name=f.name,pages=size,fallback=fallback,short_pages=blank,short_with_image=img,chars=chars)))

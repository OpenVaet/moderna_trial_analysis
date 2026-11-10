from pathlib import Path
import gzip,re,json,os,time,concurrent.futures as cf,collections
w=Path('.')
patterns={
'ae_history':r'(?:mov\w*|transfer\w*|record\w*|report\w*|captur\w*|enter\w*|add\w*|remov\w*|delet\w*|inactivat\w*).{0,50}(?:medical history|med(?:ical)?\s+hx)|(?:medical history|med(?:ical)?\s+hx).{0,65}(?:not an? (?:AE|adverse event)|instead of|remove|inactivate|delet|should not be)',
'ae_remove':r'(?:please|should|need(?:s)? to|instructed to|requested to|ask(?:ed)? to).{0,32}(?:inactivat\w*|remov\w*|delet\w*).{0,60}(?:\bAE\b|adverse event|\bSAE\b)|(?:\bAE\b|adverse event|\bSAE\b).{0,55}(?:should be (?:inactivated|removed|deleted)|does not (?:need|require) (?:to be )?report)',
'late_covid_visit':r'(?:retrospect\w*|retroactiv\w*|backdat\w*).{0,75}(?:illness|covid|symptom|test|swab|visit)|(?:illness|covid|symptom|swab|visit).{0,75}(?:retrospect\w*|retroactiv\w*|backdat\w*)',
'covid_redirect':r'(?:inactivat\w*|remov\w*|delet\w*|instead of|do not (?:report|enter)).{0,75}(?:covid-?19|illness visit|symptom log)|(?:covid-?19|illness visit|symptom log).{0,75}(?:inactivat\w*|remov\w*|delet\w*|instead of|do not (?:report|enter))',
'missed_swab':r'(?:swab|PCR|nasopharyngeal).{0,50}(?:not (?:obtained|taken|collected|performed|done)|missed|not available|not ordered|not tested|unable to collect|refused|declined)|(?:did not|unable to|declined|refused|forgot to|not able to|failed to|could not).{0,36}(?:swab|PCR)',
'deviation_pressure':r'(?:avoid|prevent|not (?:to|have to) (?:record|report)).{0,45}(?:protocol deviation|\bdeviation)|(?:otherwise|if not|or (?:else|a)|would result in).{0,25}(?:protocol deviation|deviation)',
'treatment_error':r'(?:incorrect|wrong|inadvertent\w*|in error|accidentally).{0,42}(?:product|placebo|vaccine|randomiz\w*|assignment|treatment arm)|(?:product|placebo|vaccine|randomiz\w*|assignment|treatment arm).{0,45}(?:incorrect|inadvertent|in error|wrong arm|accidentally)',
'bleeding':r'(?:(?:heavy|irregular|abnormal|vaginal|uterine|menstrual|intermenstrual|breakthrough|postmenopausal) (?:menstrual |vaginal |uterine )?(?:bleeding|hemorrhage)|menorrhagia|metrorrhagia|amenorrh\w*|abnormal menstru\w*)',
'safety_db':r'(?:safety (?:database|team|report|form)|CIOMS|pharmacovigilance).{0,60}(?:missing|not (?:reported|recorded|entered)|not in|inconsisten\w*|not match)|(?:missing|not (?:reported|recorded|entered)|does not match).{0,60}(?:safety database|safety report|CIOMS|SAE form)',
'diary_ae':r'(?:diary|solicited).{0,70}(?:should (?:be|not be) (?:recorded|reported|entered)|inactivat\w*|must be reported|please (?:add|record|enter)|not reported as an AE)|(?:should (?:be|not be) (?:recorded|reported|entered)|inactivat\w*|must be reported|please (?:add|record|enter)).{0,70}(?:diary|solicited)',
}
pat={k:re.compile(v,re.I) for k,v in patterns.items()}
wanted={
'Is the adverse event serious?':'sae','Is this event a confirmed diagnosis of Symptomatic Covid-19?':'covidAE','Outcome':'outcome','Relationship to investigational product':'rel','Adverse event':'aeterm',
'Date of unblinding (dd MMM yyyy)':'ubdate','Participant randomization assignment':'assigned','Actual Dose 1':'dose1','Actual Dose 2':'dose2','Was the participant unblinded?':'unblind','Did Subject Test Positive For SARS-CoV-2 by RT-PCR?':'positive','Was the Subject Tested For SARS-CoV-2 by RT-PCR?':'tested','Date of Test':'testdate','Date of Visit':'visitdate','Date of Contact':'contactdate','Has the subject reported symptoms of SARS-COV-2?':'reported','Was this diagnostic test performed at a lab other than the Study Central Lab?':'local','Date of death (dd MMM yyyy)':'deathdate'}
def scan(f):
 pg=gzip.decompress(f.read_bytes()).decode().split('\f')[:-1];fields=[];aes=[];hits=[];cnt=collections.Counter();forms=set()
 for page,s in enumerate(pg,1):
  fld=re.search(r'^Folder: (.*)$',s,re.M);fm=re.search(r'^Form: (.*)$',s,re.M)
  if not (fld and fm):continue
  fd=fld[1].strip();fo=fm[1].strip();forms.add((fd,fo))
  # Find only explicit audit-field headers; values marked by the EDC, latest first.
  a=re.search(r'([^\n]+)\nAudit\s*\nUser\s*\nTime \(GMT\)',s)
  if a:
   q=a[1].strip();typ=wanted.get(q,'')
   if typ and (re.match(r'^Adverse Events \(\d+\)$',fo) or fo in ('Unblinding','COVID Diagnostic Test','COVID-19 Contact','End of Study')):
    vals=re.findall(r"User entered (?:'([^']{1,220})'|(empty)\.)",s[a.end():],re.I)
    vals=[re.sub(r'\s+',' ',str(v or '<EMPTY>')).strip() for v,b in vals]
    if vals:fields.append([page,fd,fo,typ,vals])
  if fd=='Adverse Events' and re.fullmatch(r'Adverse Events \(\d+\)',fo) and not a and re.search(r'\nAdverse event\s*\n',s) and 'Start date (dd MMM yyyy)' in s and 'Was this a medically-attended AE?' in s:
   t=re.search(r'\nAdverse event\s*\n(.{1,175}?)\nWas this a medically-attended AE\?',s,re.S);dt=re.search(r'Start date \(dd MMM yyyy\)\s*\n(\d{1,2} [A-Z]{3} \d{4})',s,re.I);end=re.search(r'If not Ongoing, end date \(dd MMM yyyy\)\s*\n(\d{1,2} [A-Z]{3} \d{4})',s,re.I)
   if t:aes.append([page,fo,re.sub(r'\s+',' ',t[1]).strip(),dt[1] if dt else '',end[1] if end else ''])
  # search audit narratives incl continuation pages, not occurrence of a question/template.
  if not any(x in s for x in ('query','Query','reason for change','Reason for change','User entered')):continue
  if not any(x in s.lower() for x in ('medical history','med hx','inactivat','remov','delet','retro','backdat','swab','pcr','deviation','incorrect','wrong','error','bleed','menorr','metrorr','amenorr','safety database','safety report','cioms','diary','solicited','inadvertent','instead of','not reported','pharmacovigilance','menstr','hemorrhag')):continue
  ss=' '.join(s.split())
  # Scientific/rule flags stored with source page; do not call flags deviations.
  for k,r in pat.items():
   h=r.search(ss)
   if h:
    cnt[k]+=1
    hits.append([page,k,fd,fo,ss[max(0,h.start()-70):h.end()+125]])
 result={'fields':fields,'aes':aes,'hits':hits,'forms':sorted(forms),'n':len(pg),'counts':cnt}
 path=w/'scans'/(f.stem+'.json.gz');tmp=path.with_suffix('.part')
 with gzip.open(tmp,'wt',encoding='utf8',compresslevel=2) as fout:json.dump(result,fout,separators=(',',':'),ensure_ascii=False)
 tmp.replace(path)

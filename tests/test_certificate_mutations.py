"""Reject corrupted certificates for the intended reason, without a success claim."""
import json,subprocess,sys,tempfile,copy
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
CERT=ROOT/'certificates'/'sharp_recovery.json'
VERIFY=ROOT/'proofs'/'verify_sharp_recovery.py'
src=json.loads(CERT.read_text())
def run(data,code=None):
 with tempfile.TemporaryDirectory() as t:
  p=Path(t)/'certificate.json';p.write_text(json.dumps(data))
  r=subprocess.run([sys.executable,str(VERIFY),str(p)],capture_output=True,text=True)
  if code is None:assert r.returncode==0,r.stderr
  else:
   assert r.returncode!=0 and code in r.stderr,(code,r.stdout,r.stderr)
   assert 'CERTIFICATE VALID' not in r.stdout,'FALSE_SUCCESS'
run(src);print('PASS unmodified certificate')
x=copy.deepcopy(src);x['t']='1001/1000';run(x,'CLAIMED_CONSTANT');print('PASS wrong constant rejected')
x=copy.deepcopy(src);x['reduced_G'][0][0]='-1';run(x,'POSITIVE_DEFINITE');print('PASS broken positivity rejected')
x=copy.deepcopy(src)
from fractions import Fraction as F
x['reduced_G'][0][0]=str(F(x['reduced_G'][0][0])+1);run(x,'SOS_IDENTITY');print('PASS changed polynomial rejected')
x=copy.deepcopy(src)
# Add the identical positive identity summand to Q and S. Difference and PSD survive.
# Their trace normalization does not. This is an isolated normalization test.
for key in ['G','J']:
 B=x['basis_'+key];R=x['reduced_'+key];n=len(R)
 for i,row in enumerate(B):row.append('1' if i==0 else '0')
 for row in R:row.append('0')
 R.append(['0']*n+['1'])
run(x,'PARTIAL_TRACE');print('PASS isolated normalization corruption rejected')
r=subprocess.run([sys.executable,'-O',str(VERIFY)],capture_output=True,text=True)
assert r.returncode!=0 and 'requires assertions' in r.stderr
print('PASS -O refused')

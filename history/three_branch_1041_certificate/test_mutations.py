from pathlib import Path
from fractions import Fraction as F
from tempfile import TemporaryDirectory
import copy,json
from verify import verify
base=json.loads(Path(__file__).with_name('certificate.json').read_text())
def reject(z,expected):
 with TemporaryDirectory() as d:
  p=Path(d)/'mutant.json';p.write_text(json.dumps(z))
  try:verify(p)
  except AssertionError as e:assert str(e)==expected,(str(e),expected)
  else:raise AssertionError('MUTATION_ACCEPTED')
verify(Path(__file__).with_name('certificate.json'))
z=copy.deepcopy(base);z['t']='1';reject(z,'CLAIMED_CONSTANT')
z=copy.deepcopy(base);z['G'][0][0]='-1';reject(z,'POSITIVE_DEFINITE')
z=copy.deepcopy(base);z['J'][0][0]=str(F(z['J'][0][0])+1);reject(z,'SOS_IDENTITY')
z=copy.deepcopy(base)
for name in ['G','J']:z[name][0][0]=str(F(z[name][0][0])+1)
reject(z,'PARTIAL_TRACE')
print('PASS: original accepted; four mutations rejected for the intended reasons.')

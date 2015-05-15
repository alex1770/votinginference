import sys,shlex

from parsedatafiles import d,cc,cs

years=[2010,2015]

for c in sorted(d.keys()):
  assert list(d[c])==years
  tt=[sum([x[2] for x in d[c][y]]) for y in years]
  print "%8d %8d"%tuple(tt),
  print "%8d %s %7.3f"%(cs[c],cc[c],tt[1]/float(tt[0]))

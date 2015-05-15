import sys,shlex

from parsedatafiles import d

years=[2010,2015]

swing={}
r={}
for c in sorted(d.keys()):
  assert list(d[c])==years
  w=[]
  for y in years:
    d[c][y].sort(key=lambda x:-x[2])
    p=d[c][y][0][0]# winning party
    w.append(p)
  k=tuple(w)
  swing[k]=swing.get(k,0)+1

for k in sorted(list(swing),key=lambda x:-swing[x]):
  print k,swing[k]

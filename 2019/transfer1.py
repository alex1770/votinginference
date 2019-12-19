#!/usr/bin/env python3

# Simple voting transfer to the major two parties
# E.g., python3 transfer1.py 0.7 0.3 0.3 0.7

from voting2019dict import vote2019

import sys
x0=float(sys.argv[1])# Transfer from Lib Dems to Labour
x1=float(sys.argv[2])# Transfer from Lib Dems to Conservative
x2=float(sys.argv[3])# Transfer from TBP to Labour
x3=float(sys.argv[4])# Transfer from TBP to Conservative

d={}
cw={}# map from constituency to (winning party, winning vote)
for const in vote2019:
  v=vote2019[const].copy()
  for (donor,x,y) in [("Liberal Democrat",x0,x1), ("The Brexit Party",x2,x3)]:
    if donor in v:
      t=v[donor][1]
      v[donor]=(v[donor][0],t*(1-x-y))
      for (donee,z) in [("Labour",x), ("Conservative",y)]:
        y=0
        v[donee]=(v[donee][0], v[donee][1]+z*t)
  best=(None,-1)
  for party in v:
    votes=v[party][1]
    d[party]=d.get(party,0)+votes
    if votes>best[1]: best=(party,votes)
  cw[const]=best

l=list(d)
l.sort(key=lambda x:-d[x])
for x in l:
  if d[x]<100000: break
  print("%9d"%d[x],x)
print()

e={}
for const in cw:
  party=cw[const][0]
  e[party]=e.get(party,0)+1
l=list(e)
l.sort(key=lambda x:-e[x])
for x in l: print("%4d"%e[x],x)
print()

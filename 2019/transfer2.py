#!/usr/bin/env python3

# Constituency-dependent voting transfer to one of the top two parties
# E.g., python3 transfer2.py 0.7 0.7 0.3 0.7 0.3 0 0.7 0.7
#       python3 transfer2.py 1 1 0 1 0.5 0 0.5 1

from voting2019dict import vote2019

import sys

xfer={}
xfer[("Liberal Democrat","Labour")]           = float(sys.argv[1]) # transfer proportion of Lib Dem to Labour
xfer[("Labour","Liberal Democrat")]           = float(sys.argv[2]) # transfer proportion of Labour to Lib Dem
xfer[("Liberal Democrat","Conservative")]     = float(sys.argv[3]) # etc
xfer[("Conservative","Liberal Democrat")]     = float(sys.argv[4])
xfer[("The Brexit Party","Labour")]           = float(sys.argv[5])
xfer[("Labour","The Brexit Party")]           = float(sys.argv[6])
xfer[("The Brexit Party","Conservative")]     = float(sys.argv[7])
xfer[("Conservative","The Brexit Party")]     = float(sys.argv[8])

wins={}
totvotes={}
for const in vote2019:
  v=vote2019[const]
  vl=[(v[party][1],party) for party in v]
  vl.sort(reverse=True)
  # Try to transfer votes from rank>=2 parties to rank=0,1 parties
  votes=[v for (v,p) in vl]
  for i in range(2):
    for j in range(2,len(vl)):
      t=(vl[j][1],vl[i][1])
      if t in xfer: x=xfer[t]*vl[j][0];votes[i]+=x;votes[j]-=x
  for j in range(2,len(vl)):
    if votes[j]<-1e-8: print("Error: transfered more than 100% of votes from",vl[j][1]);sys.exit(1)
  i=(votes[0]<votes[1])# index of winning party
  party=vl[i][1]
  wins[party]=wins.get(party,0)+1
  for i in range(len(vl)): totvotes[vl[i][1]]=totvotes.get(vl[i][1],0)+votes[i]

l=list(totvotes)
l.sort(key=lambda x:-totvotes[x])
for x in l:
  if totvotes[x]<100000: break
  print("%9d"%totvotes[x],x)
print()

l=list(wins)
l.sort(key=lambda x:-wins[x])
for x in l: print("%4d"%wins[x],x)

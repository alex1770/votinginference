import sys,shlex

cc={}# map from constituency to country
cs={}# map from constituency to electorate size
f=open('constinfo','r')
for x in f:
  y=shlex.split(x)
  cc[y[0]]=y[2][0]
  cs[y[0]]=int(y[1])
f.close()

d={}# d[constituency][year]=[(party,candidate,vote),...]
f=open('votingdata','r')
for x in f:
  y=shlex.split(x)
  y[1]=int(y[1]);y[4]=int(y[4])
  if y[0] not in d: d[y[0]]={}
  if y[1] not in d[y[0]]: d[y[0]][y[1]]=[]
  d[y[0]][y[1]].append(tuple(y[2:]))
f.close()

# -*- coding: utf-8 -*-

import sys,shlex
from numpy import array,zeros,dot
from numpy.linalg import lstsq,det,norm
from scipy.optimize import fmin_slsqp
opmode=1;o1pow=1

countries=['E']
#parties=['Conservative','Other']
#parties=['Labour','Conservative','Other']
#parties=['Labour','Conservative','Other','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','Other','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Other','NoVote']
parties=['Labour','Conservative','Liberal Democrat','UKIP','Green','Other','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Scottish Green','SNP','Other','NoVote'];countries=['S']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Scottish Green','SNP','Other'];countries=['S']
#parties=['Labour','Conservative','Liberal Democrat','SNP','Other','NoVote'];countries=['S']
#parties=['Labour','Conservative','Liberal Democrat','Green','Plaid Cymru','UKIP','Other','NoVote'];countries=['W']
# There is a hard-coded assumption that "Other" is on the list of parties
# (and that there are no real parties called "Other" or "NoVote").

shortname={'Conservative':'CON', 'Labour':'LAB', 'Liberal Democrat':'LD', 'UKIP':'UKIP', 'Green':'GRN', 'SNP':'SNP', 'Plaid Cymru':'PC', 'Independent':'IND', 'BNP':'BNP', 'TUSC':'TUSC', 'Sinn Féin':'SF', 'Scottish Green':'SCG', 'Other':'OTH', 'NoVote':'NOV'}

n=len(parties)
N=n*n# Number of variables (ignoring the fact that there are 2n-1 constraints)
colconstr=True

years=[2010,2015]

seed=1;prlevel=1

print "Using countries:",', '.join(countries)
print "Using parties:",', '.join(parties)
print "seed:",seed
print "opmode, o1pow:",opmode,o1pow
print "colconstraints:",colconstr

# Get voting data, country codes, constituency sizes
from parsedatafiles import d,cc,cs

# Simplify and reduce data to allowable parties and countries, amalgamating 'Others'
sid={}# Simplified voting data by constituency and year
ptn={}# party -> number
for i in range(n): ptn[parties[i]]=i
for c in sorted(list(d)):
  if cc[c] not in countries: continue
  sp=0
  for y in years:
    for (p,cand,v) in d[c][y]:
      if p=='Speaker': sp=1
  if sp: continue# Exclude the speaker's constituency, because it is very different
  sid[c]={}
  for y in years:
    sid[c][y]=zeros(n)
    t=sum([x[2] for x in d[c][y]])# turnout
    for (p,cand,v) in d[c][y]+[('NoVote','Russell',cs[c]-t)]:# p=party, v=#votes
      sid[c][y][ptn[p if p in parties else 'Other']]+=v
print "Size of reducted constituency list:",len(sid)

tot=[]# Total votes for each party, by year
for y in years:
  tot.append([0]*n)
  for c in sid:
    for i in range(n): tot[-1][i]+=sid[c][y][i]

def prprob(prob):
  prec=3
  print "      ",
  for i in range(n): print " "*(prec-1)+"%4s"%shortname[parties[i]],
  print
  for i in range(n):
    print " %4s:"%shortname[parties[i]],
    for j in range(n):
      print "%*.*f"%(prec+3,prec,prob[i,j]),
    print

def prtot(tot):
  print "     ",
  for y in years: print "%8d"%y,
  print
  for i in range(n):
    print " %4s"%shortname[parties[i]],
    for y in range(2): print "%8d"%tot[y][i],
    print

def prflow(prob,tot):
  un=1000
  print "Flows in units of",un
  tt=[0]*n
  print "      ",
  for i in range(n): print "    %4s"%shortname[parties[i]],
  print
  for i in range(n):
    print " %4s:"%shortname[parties[i]],
    for j in range(n):
      flow=tot[0][i]*prob[i,j]
      print "%8.0f"%(flow/un),
      tt[j]+=flow
    print
  print " -----",
  for i in range(n): print "--------",
  print
  print "  Tot:",
  for i in range(n): print "%8.0f"%(tt[i]/un),
  print
  print " Real:",
  for i in range(n): print "%8.0f"%(tot[1][i]/un),
  print

def decodeparams(params):
  prob=zeros((n,n))
  for i in range(n):
    for j in range(n):
      prob[i,j]=params[i*n+j]
  return prob

def encodeparams(prob):
  params=[0]*N
  for i in range(n):
    for j in range(n):
      params[i*n+j]=prob[i,j]
  return params

bestv=1e30;its=0
def negprob(params):
  global bestv,its,resid
  prob=decodeparams(params)
  negp=0.;resid=[]
  if prlevel>=3: print params;prprob(prob)
  for c in sorted(list(sid)):
    #print c
    mu=zeros(n);cov=zeros((n,n))
    s0=sid[c][years[0]]
    s1=sid[c][years[1]]
    for p0,v0 in enumerate(s0):
      for p1,v1 in enumerate (s1):
        mu[p1]+=v0*prob[p0,p1]
        if opmode==2:
          for p2,v2 in enumerate (s1):
            cov[p1,p2]+=v0*prob[p0,p1]*((p1==p2)-prob[p0,p2])
    z=zeros(n)
    for p,v in enumerate(s1): z[p]=v-mu[p]
    if opmode==0: negp+=dot(z,z)# simple least-squares
    elif opmode==1:
      d0=dot(z,z)/dot(s0,s0)
      z/=dot(s0,s0)
      qrob=zeros((n,n))
      for p0 in range(n):
        for p1 in range(n):
          qrob[p0,p1]=prob[p0,p1]+s0[p0]*z[p1]
          if qrob[p0,p1]<0: d0+=qrob[p0,p1]**2
          elif qrob[p0,p1]>1: d0+=(qrob[p0,p1]-1)**2
      #print c;prprob(qrob);print
      #di=prob-qrob;print d0,norm(di)**2
      d0**=o1pow/2.
    else:
      z0=z;cov0=cov
      z=z[0:n-1];cov=cov[0:n-1,0:n-1]# There's a sum=0 constraint, so use subspace
      try:
        d0=dot(z,lstsq(cov,z)[0])/2+log(det(cov))/2
      except:
        print c
        print z0
        print cov0
        print z
        print cov
        assert 0
    resid.append((d0,c))
    negp+=d0
  its+=1
  if (negp<bestv)+prlevel>=2:
    print "%d: negprob/const ="%its,negp/len(sid)
    prtot(tot)
    prprob(prob)
    prflow(prob,tot)
    print
    sys.stdout.flush()
  if negp<bestv: bestv=negp
  return negp

prob=array([[.8*(j==i)+0.2/(n-1)*(j!=i) for j in range(n)] for i in range(n)])
params=encodeparams(prob)

if 0:# Some starting point
  prob=array([
    [0.945,  0.002,  0.0005,  0.019,  0.0005,  0.0005,  0.034],
    [0.001,  0.942,  0.0005,  0.045,  0.0005,  0.0005,  0.012],
    [0.124,  0.144,  0.369,  0.104,  0.119,  0.032,  0.108],
    [0.001,  0.002,  0.0005,  0.994,  0.001,  0.001,  0.002],
    [0.059,  0.002,  0.001,  0.001,  0.933,  0.002,  0.001],
    [0.075,  0.001,  0.0005,  0.560,  0.0005,  0.022,  0.341],
    [0.040,  0.015,  0.0005,  0.076,  0.005,  0.001,  0.863]
  ])
  params=encodeparams(prob)

def fmin(np,ar): return fmin_slsqp(np,ar,epsilon=1e-4,f_eqcons=lambda params:[sum([params[i*n+j] for j in range(n)])-1 for i in range(n)]+[sum([params[i*n+j]*tot[0][i] for i in range(n)])-tot[1][j] for j in range(n-1)]*colconstr,bounds=[(1e-4,1-1e-4)]*N)

best=fmin(negprob,params)
print best
resid.sort()
for (v,c) in resid: print "%12g"%v,c
print "Done"

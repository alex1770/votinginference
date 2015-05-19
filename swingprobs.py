# -*- coding: utf-8 -*-

import sys,shlex
from numpy import array,zeros,ones,dot,outer
from numpy.linalg import lstsq,det,norm,eigvalsh,eigvals
from scipy.optimize import fmin_slsqp
from math import log,exp,pi,sqrt
from random import normalvariate
opmode=3;o1pow=2;lam=2.0
colconstr=True
rap=True

countries=['E']
#parties=['Conservative','NoVote']
#parties=['Labour','Conservative','NoVote']
#parties=['Labour','Conservative','Other','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','Other','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Other','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Green','Other','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Green','NoVote']
parties=['Conservative','Labour','Liberal Democrat','UKIP','Green','NoVote']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Scottish Green','SNP','Other','NoVote'];countries=['S']
#parties=['Labour','Conservative','Liberal Democrat','UKIP','Scottish Green','SNP','NoVote'];countries=['S']
#parties=['Labour','SNP','Liberal Democrat','Conservative','Other','NoVote'];countries=['S']
#parties=['Labour','Conservative','Liberal Democrat','Green','Plaid Cymru','UKIP','Other','NoVote'];countries=['W']
#parties=['SNP','Labour','Liberal Democrat','Conservative','NoVote'];countries=['S']
# There is a hard-coded assumption that "NoVote" is on the list of parties
# (and that there are no real parties called "Other" or "NoVote").
assert 'NoVote' in parties

shortname={'Conservative':'CON', 'Labour':'LAB', 'Liberal Democrat':'LD', 'UKIP':'UKIP', 'Green':'GRN', 'SNP':'SNP', 'Plaid Cymru':'PC', 'Independent':'IND', 'BNP':'BNP', 'TUSC':'TUSC', 'Sinn Féin':'SF', 'Scottish Green':'SCG', 'Other':'OTH', 'NoVote':'NOV'}

n=len(parties)
N=n*n# Number of variables (ignoring the fact that there are 2n-1 constraints)
one=ones(n)

years=[2010,2015]

seed=1;prlevel=1

print "Using countries:",', '.join(countries)
print "Using parties:",', '.join(parties)
print "Seed:",seed
print "opmode, o1pow, lam:",opmode,o1pow,lam
print "colconstraints:",colconstr
print "Require all (destination) parties present:",rap

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
      if p not in parties:
        if 'Other' in parties: p='Other'
        else: p='NoVote'
      sid[c][y][ptn[p]]+=v
  if rap and sum([x>0 for x in sid[c][years[1]]])<n: del sid[c]
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
  for y in years: print "                   %8d     "%y,
  print
  print "     ",
  for y in years: print "    Number   %exclNOV   %inclNOV",
  print
  tt={}
  for y in range(2):
    s=float(sum(tot[y]));tt[y]=[s-tot[y][ptn['NoVote']],s]
  for i in range(n):
    print " %4s"%shortname[parties[i]],
    for y in range(2):
      print "  %8d"%tot[y][i],
      if parties[i]=='NoVote': print "       ---",
      else: print "     %5.1f"%(tot[y][i]/tt[y][0]*100),
      print "     %5.1f"%(tot[y][i]/tt[y][1]*100),
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

def prod(l): return reduce(lambda x,y:x*y,l)

def op3average(ns,prob):
  trob=zeros((n,n))
  for i in xrange(ns):
    D1=array([exp(normalvariate(0,1/sqrt(lam))) for j in range(n)])
    D0=1/dot(prob,D1)
    trob+=outer(D0,one)*prob*outer(one,D1)
  trob/=ns
  return trob

bestv=1e30;its=0
def negprob(params):
  global bestv,its,resid
  prob=decodeparams(params)
  if (prob<0).any(): print "BAD";return 1e30# Workaround bug in slsqp
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
        if opmode==0:
          for p2,v2 in enumerate (s1):
            cov[p1,p2]+=v0*prob[p0,p1]*((p1==p2)-prob[p0,p2])
    z=zeros(n)
    for p,v in enumerate(s1): z[p]=v-mu[p]
    if opmode==0: 
      z0=z;cov0=cov.copy()
      z=z[0:n-1];cov=cov[0:n-1,0:n-1]# There's a sum=0 constraint, so use subspace
      try:
        d0=dot(z,lstsq(cov,z)[0])/2+log(det(cov))/2+(n-1)/2.*log(2*pi)
      except:
        print c,"\n"
        prprob(prob)
        print z0,"\n"
        print cov0,"\n"
        print eigvalsh(cov0),"\n"
        print z,"\n"
        print cov,"\n"
        print eigvalsh(cov),"\n"
        assert 0
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
    elif opmode==2 or opmode==3:
      deb=0
      qrob=prob.copy()
      if deb>=1: prprob(qrob);print s0;print s1;print
      D0=one.copy();D1=one.copy()
      # (Add something small to s1, or otherwise tweak, if zero components are allowed,
      # otherwise the SK-iteration below may not terminate.)
      while 1:
        m1=s1/dot(s0,qrob)
        qrob=qrob*outer(one,m1);D1=D1*m1
        if deb>=2: print "m1",m1;prprob(qrob);print
        m0=one/dot(qrob,one)
        qrob=qrob*outer(m0,one);D0=D0*m0
        if deb>=2: print "m0",m0;prprob(qrob);print
        err=dot(m0-one,m0-one)+dot(m1-one,m1-one)
        if err<1e-12: break
      if deb>=1:
        prprob(qrob)
        prprob(outer(D0,one)*prob*outer(one,D1))
        print D0;print D1
        print
      assert min(D1)>1e-6
      # assert D0 == 1/dot(prob,D1)
      if opmode==2:
        rtrnorm=sum(D1)/n;D0*=rtrnorm;D1/=rtrnorm
        #lftnorm=sum(D0)/n;D0/=lftnorm;D1*=lftnorm
      elif opmode==3:
        rdetnorm=prod(D1)**(1./n);D0*=rdetnorm;D1/=rdetnorm
      # Model is: Constituency swing matrix, A = diag((Pd)^{-1}).P.diag(d)
      # d=D1, P=prob=overall swing matrix
      # s0^t.A = vector of predicted 2015-votes for each party
      # jac = Jacobian d(predicted 2015-vote for party r)/d(D1[s])
      jac=zeros((n,n))
      for r in range(n):
        for s in range(n):
          for i in range(n):
            jac[r,s]+=(r==s)*s0[i]*D0[i]*prob[i,r]-s0[i]*D0[i]**2*prob[i,r]*prob[i,s]*D1[r]
          if 0:
            eps=1e-5
            epsv=array([(i==s)*eps for i in range(n)])
            D1a=D1+epsv;D0a=1/dot(prob,D1a);s1a=dot(s0,outer(D0a,one)*prob*outer(one,D1a))
            D1b=D1-epsv;D0b=1/dot(prob,D1b);s1b=dot(s0,outer(D0b,one)*prob*outer(one,D1b))
            der=(s1a-s1b)/(2*eps)
            print jac[r,s],der[r]
      if opmode==2:
        #print det(jac),det(jac[0:n-1,0:n-1]),det(jac[1:n,1:n])
        dj=log(det(jac[0:n-1,0:n-1]))
        #D0=D0/sum(D0)*n;D1=D1/sum(D1)*n
        #D0=D0/prod(D0)**(1./n);D1=D1/prod(D1)**(1./n)
        #print sum(D0),sum(D1)
        #d0=dot(D0-one,D0-one)**(o1pow/2.)
        #d0=dot(D1-one,D1-one)**(o1pow/2.)
        d0=sum([abs(x-1) for x in D1])+dj
      elif opmode==3:
        dj=log(det(jac[0:n-1,0:n-1])/D1[n-1])# prod(D1[:n-1])=1/D1[n-1]
        d0=lam*sum([log(x)**2/2. for x in D1])+(n-1)*log(2*pi/lam)/2+log(n)/2+dj
    else: assert 0
    resid.append((d0,c))
    negp+=d0
  its+=1
  if (negp<bestv)+prlevel>=2:
    print "\n%d: negprob/const ="%its,negp/len(sid)
    print;prtot(tot)
    print;prprob(prob)
    print;prflow(prob,tot)
    if opmode==3:
      ns=10000;trob=op3average(ns,prob)
      print;print ns,"sample"+"s"*(ns!=1),"averaged swings";prprob(trob)
      print;prflow(trob,tot)
    print
    sys.stdout.flush()
  if negp<bestv: bestv=negp
  return negp/(100000. if opmode==0 else 1.)# Workaround buggy scale requirement in slsqp

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

if 1 and 'E' in countries:
  prob=array([
    [0.956,  0.006,  0.001,  0.032,  0.001,  0.004],
    [0.020,  0.955,  0.002,  0.016,  0.001,  0.006],
    [0.107,  0.123,  0.344,  0.292,  0.116,  0.018],
    [0.001,  0.001,  0.001,  0.995,  0.001,  0.001],
    [0.011,  0.001,  0.001,  0.001,  0.982,  0.004],
    [0.035,  0.015,  0.001,  0.035,  0.010,  0.904]
  ])
  params=encodeparams(prob)

if 1 and 'S' in countries:
  prob=array([
    [0.968,  0.001,  0.001,  0.029,  0.001],
    [0.321,  0.660,  0.017,  0.001,  0.001],
    [0.483,  0.040,  0.384,  0.092,  0.001],
    [0.021,  0.009,  0.057,  0.912,  0.001],
    [0.269,  0.001,  0.001,  0.001,  0.728]
  ])
  params=encodeparams(prob)

def fmin(np,ar): return fmin_slsqp(np,ar,iter=1000000,epsilon=1e-4,f_eqcons=lambda params:[sum([params[i*n+j] for j in range(n)])-1 for i in range(n)]+[sum([params[i*n+j]*tot[0][i] for i in range(n)])-tot[1][j] for j in range(n-1)]*colconstr,bounds=[(1e-4,1-1e-4)]*N)

# fprime_eqcons=lambda x:array([[int(j>=i*n and j<(i+1)*n) for j in range(N)] for i in range(n)]),

best=fmin(negprob,params)
resid.sort()
for (v,c) in resid: print "%12g"%v,c
prob=decodeparams(best)
print;prtot(tot)
print;prprob(prob)
print;prflow(prob,tot)
if opmode==3:
  ns=1000000;trob=op3average(ns,prob)
  print;print ns,"sample"+"s"*(ns!=1),"averaged swings";prprob(trob)
  print;prflow(trob,tot)
print "Done"

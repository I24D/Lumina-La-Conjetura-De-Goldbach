#!/usr/bin/env python3
"""Independent finite test: localized conductor identity in the fixed-switch support.
Proxy parameters are deliberately independent of asymptotic X; no distribution theorem tested.
"""
import math
from collections import defaultdict

def primes_upto(n):
    flags=bytearray(b"\x01")*(n+1); flags[:2]=b"\x00\x00"
    for p in range(2,math.isqrt(n)+1):
        if flags[p]:
            flags[p*p:n+1:p]=b"\x00"*(((n-p*p)//p)+1)
    return [p for p in range(2,n+1) if flags[p]]

def divisors(fs):
    out=[1]
    for p in fs: out += [d*p for d in out]
    return out

def vplus(fs,L):
    t=1
    for j,p in enumerate(sorted(fs,reverse=True),1):
        t*=p
        if j%2 and p*p*t>=L: return False
    return t<L

def support(N,A,z,y,D,primes):
    ps=[p for p in primes if 2<p<z and math.gcd(p,A*N)==1]
    out=[]
    for s in primes:
        if s>=y: break
        if s<z or math.gcd(s,A*N)!=1: continue
        for mask in range(1<<len(ps)):
            fs=[p for i,p in enumerate(ps) if mask>>i&1]
            m=s*math.prod(fs)
            if vplus(fs,D/s):
                out.append((m,(-1)**len(fs),[s]+fs))
    return out

def run(N,A,z,y,D,lam,primes):
    rows=support(N,A,z,y,D,primes)
    assert rows and math.gcd(N,A)==1
    afs=[p for p in primes if A%p==0]
    phiA=A
    for p in afs: phiA=phiA//p*(p-1)
    seq=[(n,lam[n]) for n in range(N,0,-A) if lam[n]]
    psiA=math.fsum(w for n,w in seq)
    psi=math.fsum(lam)
    EA=psiA-psi/phiA
    K=defaultdict(float); fact={1:[]}
    for m,b,fs in rows:
        ph=math.prod(p-1 for p in fs)
        for h in divisors(fs):
            K[h]+=b/ph
            fact[h]=[p for p in fs if h%p==0]
    # Primitive character kernel on squarefree h, computed by local
    # orthogonality, with zero for non-units. S[1] is the exterior progression sum.
    S={}
    for h,fs in fact.items():
        terms=[]
        for n,w in seq:
            if math.gcd(n,h)!=1: continue
            k=1
            for p in fs:
                k*=p-2 if (n-N)%p==0 else -1
            terms.append(w*k)
        S[h]=math.fsum(terms)
    R=[]; induced=[]; maxerr=0.0
    for m,b,fs in rows:
        ph=math.prod(p-1 for p in fs)
        psim=math.fsum(w for n,w in seq if (N-n)%m==0)
        R.append(b*(psim-psi/(phiA*ph)))
        Im=0.0
        for p in fs:
            modulus=A*(m//p)
            power=p
            while power<=N:
                if (power-N)%modulus==0:
                    Im+=math.log(p)/(p-1)
                if power>N//p: break
                power*=p
        induced.append(b*Im)
        reconstructed=math.fsum(S[h] for h in divisors(fs))/ph-Im
        maxerr=max(maxerr,abs(psim-reconstructed))
    R=math.fsum(R); I=math.fsum(induced)
    rebuilt=K[1]*EA+math.fsum(K[h]*S[h] for h in K if h>1)-I
    err=abs(R-rebuilt)
    assert err<=1e-7*max(1,abs(R)) and maxerr<=1e-7
    assert any(len(fs)>1 for _,_,fs in rows), "test must contain composites"
    B=math.fsum(abs(S[h])**2/math.prod(p-1 for p in fact[h]) for h in K if h>1)
    G=math.fsum(math.prod(p-1 for p in fact[h])*abs(K[h])**2 for h in K)
    C=math.fsum(K[h]*S[h] for h in K if h>1)
    assert C*C<=G*B+1e-6*max(1,G*B)
    print(f"N={N} A={A} z={z} y={y} D={D} rows={len(rows)} H={len(K)} "
          f"R={R:.9f} rebuilt={rebuilt:.9f} err={err:.3e} max_module_err={maxerr:.3e} "
          f"EA={EA:.6f} I={I:.6f} B_local={B:.6e} G_local={G:.6e}")
    return err,maxerr

if __name__=="__main__":
    N=200002
    primes=primes_upto(N)
    lam=[0.0]*(N+1)
    for p in primes:
        x=p
        while x<=N:
            lam[x]=math.log(p)
            if x>N//p: break
            x*=p
    checks=[run(N,A,11,37,6000,lam,primes) for A in (3,9,27)]
    print("PASS 3 proxy selectors; exact local regrouping and Cauchy checked; no asymptotic claim")

#!/usr/bin/env python3
"""Finite support diagnostic for the corrected second Cauchy factor.
N=1000000 (even), A=3, D=floor(sqrt(N/A), no log loss; small proxy only.
Direct primitive-character sums are compared with M_r + Delta_r.
Standard library only; no Goldbach proof.
"""
import itertools
import math

N = 1000000
A = 3
KAPPAS = (1.0, 1.5, 2.0)

def primes_upto(n):
    a=bytearray(b"\x01")*(n+1)
    a[:2]=b"\x00\x00"
    for p in range(2,math.isqrt(n)+1):
        if a[p]:
            j=p*p
            a[j:n+1:p]=b"\x00"*((n-j)//p+1)
    return [p for p in range(2,n+1) if a[p]]

def factor_powers(n):
    out=[]
    p=2
    while p*p<=n:
        if n%p==0:
            a=0
            while n%p==0:
                n//=p; a+=1
            out.append((p,a))
        p+=1 if p==2 else 2
    if n>1: out.append((n,1))
    return out

def divisors_all(factors):
    out=[1]
    for p,a in factors:
        old=list(out); power=1
        for _ in range(a):
            power*=p
            out += [x*power for x in old]
    return out

def phi(n):
    out=n
    for p,_ in factor_powers(n):
        out-=out//p
    return out

def mobius(n):
    f=factor_powers(n)
    if any(a>1 for _,a in f): return 0
    return -1 if len(f)%2 else 1

def primitive_root_odd_prime_power(p,a):
    m=p**a; ph=p**(a-1)*(p-1)
    qfac=[q for q,_ in factor_powers(ph)]
    for g in range(2,m):
        if math.gcd(g,m)==1 and all(pow(g,ph//q,m)!=1 for q in qfac):
            return g
    raise RuntimeError(f"primitive root not found for {m}")

def local_primitive_characters(p,a):
    """List (label, character) for primitive chars modulo p^a."""
    m=p**a
    if p==2:
        if a==1:
            return []  # U(2) has no character of conductor 2.
        if a==2:
            return [("odd", lambda n: (1+0j) if n%4==1 else (-1+0j))]
        order=2**(a-2)
        logs={}; x=1
        for t in range(order):
            logs[x]=t
            x=x*5%m
        out=[]
        for sign in (0,1):
            for k in range(1,order,2):
                def val(n,m=m,logs=logs,order=order,sign=sign,k=k):
                    n%=m
                    u=0 if n%4==1 else 1
                    x=n if u==0 else (-n)%m
                    t=logs[x]
                    angle=2*math.pi*k*t/order
                    z=complex(math.cos(angle),math.sin(angle))
                    return (-1 if sign and u else 1)*z
                out.append(((sign,k),val))
        return out
    ph=p**(a-1)*(p-1)
    g=primitive_root_odd_prime_power(p,a)
    logs={}; x=1
    for t in range(ph):
        logs[x]=t
        x=x*g%m
    out=[]
    for k in range(1,ph):
        if k%p==0: continue  # factors through p^(a-1)
        def val(n,m=m,logs=logs,ph=ph,k=k):
            t=logs[n%m]
            angle=2*math.pi*k*t/ph
            return complex(math.cos(angle),math.sin(angle))
        out.append((k,val))
    return out

def primitive_characters(modulus):
    if modulus<=1:
        return
    factors=factor_powers(modulus)
    local=[local_primitive_characters(p,a) for p,a in factors]
    if any(not x for x in local):
        return
    for choice in itertools.product(*local):
        labels=tuple(item[0] for item in choice)
        funcs=tuple(item[1] for item in choice)
        def val(n,funcs=funcs):
            z=1+0j
            for f in funcs: z*=f(n)
            return z
        yield labels,val

def iter_prime_powers(p):
    x=p
    while x<=N:
        yield x
        if x>N//p: break
        x*=p

def build_lambda(N):
    out=[]
    for p in primes_upto(N):
        x=p
        while x<=N:
            out.append((x,math.log(p)))
            if x>N//p: break
            x*=p
    return out

def vplus(fs_desc,level):
    product=1
    for j,p in enumerate(fs_desc,1):
        product*=p
        if j%2==1 and p*p*product>=level:
            return False
    return True

def squarefree_divisors(fs):
    out=[1]
    for p in fs:
        out += [x*p for x in out]
    return out

def delta_r(r):
    factors=factor_powers(r)
    delta=0j
    for p,a in factors:
        q=r//p
        if a==1 and q>1:
            subtotal=0j
            for _,chi in primitive_characters(q):
                pp=0j
                x=p
                j=1
                while x<=N:
                    pp += chi(p)**j*math.log(p)
                    if x>N//p: break
                    x*=p
                    j+=1
                subtotal += chi(N).conjugate()*pp
            delta+=subtotal
    if len(factors)==1 and factors[0][1]==1:
        p=factors[0][0]
        delta += sum(math.log(p) for _ in iter_prime_powers(p))
    return delta

def main():
    X=N//A
    D=math.isqrt(X)
    z=int(X**0.1)
    y=int(X**(1/3))
    primes=primes_upto(y-1)
    small=[p for p in primes if p<z and math.gcd(p,A*N)==1]
    rows=[]
    H=set()
    for s in primes:
        if not (z<=s<y) or math.gcd(s,A*N)!=1: continue
        L=D//s
        for mask in range(1<<len(small)):
            fs=[small[i] for i in range(len(small)) if (mask>>i)&1]
            d=math.prod(fs)
            if d<L and vplus(sorted(fs,reverse=True),L):
                rows.append((s,d))
                H.update(squarefree_divisors([s]+fs))
    adivs=[a for a in divisors_all(factor_powers(A)) if math.gcd(a,1)==1]
    conductors=sorted({a*h for a in adivs for h in H
                       if math.gcd(a,h)==1 and a*h>1})
    lam=build_lambda(N)
    psi_total=sum(w for _,w in lam)
    # Needed divisors for all conductors.
    divmap={r:divisors_all(factor_powers(r)) for r in conductors}
    needed=sorted({d for ds in divmap.values() for d in ds})
    progress={}
    for d in needed:
        target=N%d
        progress[d]=sum(w for n,w in lam if n%d==target)
    # Direct residue sums, then direct primitive-character S_r.
    direct={}
    M={}
    Delta={}
    errors={}
    for r in conductors:
        counts=[0.0]*r
        for n,w in lam:
            counts[n%r]+=w
        Sr=0j
        for _,chi in primitive_characters(r):
            psichi=sum(counts[a]*chi(a) for a in range(r)
                       if math.gcd(a,r)==1)
            Sr += chi(N).conjugate()*psichi
        direct[r]=Sr
        mr=0.0
        for d in divmap[r]:
            E=progress[d]-psi_total/phi(d)
            mr += mobius(r//d)*phi(d)*E
        M[r]=mr
        Delta[r]=delta_r(r)
        errors[r]=abs(Sr-(mr+Delta[r]))

    print(f"PARAMS N={N} A={A} X={X} D={D} z={z} y={y} small={small}")
    print(f"SUPPORT rows={len(rows)} H={len(H)} conductors={len(conductors)}")
    print(f"PRIMITIVE_SUMS nonzero={sum(abs(v)>1e-8 for v in direct.values())}")
    max_abs=max(errors.values())
    max_rel=max(errors[r]/max(1.0,abs(direct[r])) for r in conductors)
    print(f"IDENTITY max_abs_error={max_abs:.6e} max_rel_error={max_rel:.6e}")
    for kappa in KAPPAS:
        R0=math.log(N)**kappa
        high=[r for r in conductors if r>R0]
        BM=sum(abs(M[r])**2/phi(r) for r in high)
        BD=sum(abs(Delta[r])**2/phi(r) for r in high)
        BC=sum(2*(M[r]*Delta[r].conjugate()).real/phi(r) for r in high)
        Bcorr=BM+BD+BC
        Bdirect=sum(abs(direct[r])**2/phi(r) for r in high)
        ratio=BD/Bcorr if Bcorr else float("nan")
        print(
            f"kappa={kappa:.1f} R0={R0:.3f} high={len(high)} "
            f"BM={BM:.9e} Bcross={BC:.9e} Bdelta={BD:.9e} "
            f"Bcorr={Bcorr:.9e} Bdirect={Bdirect:.9e} "
            f"cross_over_BM={BC/BM:.6e} delta_over_BM={BD/BM:.6e} "
            f"net_rel_change={(Bcorr-BM)/BM:.6e} "
            f"delta_energy_share={ratio:.6e} direct_error={abs(Bcorr-Bdirect):.3e}"
        )
    # List high-conductor support for reproducibility.
    print("HIGH_CONDUCTORS_kappa1",
          ",".join(str(r) for r in conductors if r>math.log(N)))

if __name__=="__main__":
    main()
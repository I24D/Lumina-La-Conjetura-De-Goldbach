#!/usr/bin/env python3
"""Independent check of the imprimitive prime-power correction in doc 96.
Direct prime-power sums and explicit primitive characters for odd conductors.
N=999983. Standard library only. Floating-point diagnostics, not a proof.
"""
import itertools
import math

N = 999983
TEST_R = (3, 5, 7, 11, 15, 21, 35, 55, 105, 385, 9, 25, 45, 75, 225, 315)

def primes_upto(n):
    a = bytearray(b"\x01") * (n + 1)
    a[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if a[p]:
            start = p*p
            a[start:n+1:p] = b"\x00" * ((n-start)//p+1)
    return [p for p in range(2,n+1) if a[p]]

def factor_powers(n):
    out=[]
    p=2
    while p*p<=n:
        if n%p==0:
            a=0
            while n%p==0:
                n//=p
                a+=1
            out.append((p,a))
        p += 1 if p==2 else 2
    if n>1: out.append((n,1))
    return out

def divisors_all(factors):
    out=[1]
    for p,a in factors:
        old=list(out)
        power=1
        for _ in range(a):
            power*=p
            out += [x*power for x in old]
    return out

def phi(n):
    x=n
    out=n
    for p,_ in factor_powers(n):
        out -= out//p
    return out

def mobius(n):
    factors=factor_powers(n)
    if any(a>1 for _,a in factors):
        return 0
    return -1 if len(factors)%2 else 1

def primitive_root_odd_prime_power(p,a):
    m=p**a
    ph=p**(a-1)*(p-1)
    fac=[q for q,_ in factor_powers(ph)]
    for g in range(2,m):
        if math.gcd(g,m)==1 and all(pow(g,ph//q,m)!=1 for q in fac):
            return g
    raise RuntimeError(f"primitive root not found for {m}")

PRIMES = primes_upto(N)
LAM = []
for p in PRIMES:
    x=p
    while x<=N:
        LAM.append((x, math.log(p)))
        if x > N//p: break
        x*=p
PSI_TOTAL = sum(w for _,w in LAM)

def local_log_table(p,a):
    m=p**a
    ph=p**(a-1)*(p-1)
    g=primitive_root_odd_prime_power(p,a)
    out={}
    x=1
    for t in range(ph):
        out[x]=t
        x=x*g%m
    return m,ph,out

def primitive_characters(factors):
    """Yield explicit primitive characters modulo product p^a (all p odd)."""
    local=[(p,a,*local_log_table(p,a)) for p,a in factors]
    choices=[tuple(k for k in range(1,ph) if k%p!=0)
             for p,a,m,ph,logs in local]
    for ks in itertools.product(*choices):
        def value(n, local=local, ks=ks):
            z=1+0j
            for (p,a,m,ph,logs),k in zip(local,ks):
                t=logs[n%m]
                angle=2*math.pi*k*t/ph
                z *= complex(math.cos(angle), math.sin(angle))
            return z
        yield ks,value

def theta_prime(p):
    total=0.0
    x=p
    while x<=N:
        total += math.log(p)
        if x>N//p: break
        x*=p
    return total

def run_r(r):
    factors=factor_powers(r)
    assert all(p%2==1 for p,_ in factors), "tests use odd conductors"
    divs=divisors_all(factors)
    progression={d:0.0 for d in divs}
    residue_counts=[0.0]*r
    for n,w in LAM:
        residue_counts[n%r]+=w
        for d in divs:
            if n%d==N%d:
                progression[d]+=w

    # S_r from explicit primitive characters modulo r.
    S=0j
    for _,chi in primitive_characters(factors):
        S += chi(N).conjugate() * sum(
            residue_counts[a]*chi(a)
            for a in range(r) if math.gcd(a,r)==1
        )

    # M_r = sum_{d|r} mu(r/d) phi(d) E(N;d,N).
    M=0.0
    for d in divs:
        E=progression[d]-PSI_TOTAL/phi(d)
        M += mobius(r//d)*phi(d)*E

    # Delta: for p||r, primitive characters of conductor r/p
    # lose the p^j terms when induced to moduli divisible by p.
    Delta=0j
    terms=[]
    for p,a in factors:
        q=r//p
        if a==1 and q>1:
            subtotal=0j
            for _,chi_q in primitive_characters(factor_powers(q)):
                prime_power_sum=0j
                x=p
                j=1
                while x<=N:
                    prime_power_sum += (chi_q(p)**j)*math.log(p)
                    if x>N//p: break
                    x*=p
                    j+=1
                subtotal += chi_q(N).conjugate()*prime_power_sum
            Delta += subtotal
            terms.append((p,q,subtotal))
    # The principal-character correction survives only when r is prime.
    if len(factors)==1 and factors[0][1]==1:
        Delta += theta_prime(factors[0][0])

    err=abs(S-(M+Delta))
    scale=max(1.0,abs(S))
    return {
        "r":r,"S":S,"M":M,"Delta":Delta,"err":err,
        "rel":err/scale,"delta_over_S":abs(Delta)/max(abs(S),1e-300),
        "terms":terms
    }

def main():
    rows=[run_r(r) for r in TEST_R]
    for x in rows:
        print(
            f"r={x['r']:4d} |S|={abs(x['S']):.9e} "
            f"|Delta|/|S|={x['delta_over_S']:.6%} "
            f"abs_error={x['err']:.3e} rel_error={x['rel']:.3e}"
        )
    print("MAX_ABSOLUTE_ERROR",max(x["err"] for x in rows))
    print("MAX_RELATIVE_ERROR",max(x["rel"] for x in rows))
    print("DELTA_RATIO_RANGE",
          min(x["delta_over_S"] for x in rows),
          max(x["delta_over_S"] for x in rows))

if __name__=="__main__":
    main()

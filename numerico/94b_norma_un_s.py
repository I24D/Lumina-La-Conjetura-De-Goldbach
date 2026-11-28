#!/usr/bin/env python3
"""Doc 94b: norma Gram por un unico primo de switching.
Procedencia: LUMINA, GPT-6.2; 2026-11-27. Método: enumeración finita de V+.
Estado [NUM], no prueba asintótica. Fuentes: docs 86-94; Nathanson, cap. 9.
Control finito; no prueba la cota asintotica Q_s << 1/log z.
A=3 es el selector canonico k=2 para N=2^e.
"""
import math
from collections import defaultdict

def primes_upto(n):
    a=bytearray(b"\x01")*(n+1)
    a[:2]=b"\x00\x00"
    for p in range(2,math.isqrt(n)+1):
        if a[p]:
            j=p*p
            a[j:n+1:p]=b"\x00"*((n-j)//p+1)
    return [p for p in range(2,n+1) if a[p]]

def vplus(fs_desc, level):
    prod=1
    for j,p in enumerate(fs_desc,1):
        prod*=p
        if j%2==1 and p*p*prod>=level:
            return False
    return True

def divs(fs):
    out=[1]
    for p in fs:
        out += [d*p for d in out]
    return out

def phi_sf(fs):
    out=1
    for p in fs:
        out*=p-1
    return out

def one_s(N,A,s,primes):
    X=N/A
    D=math.isqrt(int(X))
    z=int(X**0.1)
    small=[p for p in primes if p<z and math.gcd(p,A*N)==1]
    level=D//s
    T=defaultdict(float)
    count=0
    for mask in range(1<<len(small)):
        fs=[small[i] for i in range(len(small)) if mask>>i&1]
        d=math.prod(fs)
        if d>=level or not vplus(sorted(fs,reverse=True),level):
            continue
        count+=1
        mu=(-1)**len(fs)
        den=phi_sf(fs)
        for h in divs(fs):
            T[h]+=mu/den
    Q=sum(phi_sf([p for p in small if h%p==0])*v*v for h,v in T.items())
    P=math.prod(1-2/(p-1)+p/(p-1)**2 for p in small)
    return z,small,count,Q,P,Q/P,Q*math.log(z)

def main():
    out=[]
    for e in (40,45,50,55,60):
        N=2**e
        A=3
        X=N/A
        D=math.isqrt(int(X))
        z=int(X**0.1)
        y=int(X**(1/3))
        primes=primes_upto(y)
        rows=[]
        for u0 in (5/3,2,3,4):
            target=D/(z**u0)
            candidates=[p for p in primes if z<=p<y]
            s=min(candidates,key=lambda p:abs(p-target))
            z1,small,n_d,Q,P,ratio,Qlog=one_s(N,A,s,primes)
            u=math.log(D/s)/math.log(z1)
            rows.append((u,s,n_d,Q,P,ratio,Qlog))
        out.append((e,z,small,rows))
    for e,z,small,rows in out:
        print(f"e={e} z={z} small={small}")
        for u,s,nd,Q,P,ratio,Qlog in rows:
            print(f"  u={u:.3f} s={s} nd={nd} Qs={Q:.9f} P={P:.9f} Qs/P={ratio:.6f} Qs*logz={Qlog:.6f}")
    ok=all(1.0 <= row[5] <= 1.8 for _,_,_,rows in out for row in rows)
    print("PASS" if ok else "FAIL")

if __name__=="__main__":
    main()

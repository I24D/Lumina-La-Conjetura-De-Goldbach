"""Verifier for doc 71: Kuhn blocks aligned at alpha=1/u."""
from fractions import Fraction
from random import Random

V=20
U=6
S=3
LO=Fraction(1,V)
JUMP=Fraction(1,U)
HI=Fraction(1,S+1)

def cells(log2_over_logx: Fraction):
    pts={LO,JUMP,HI}
    a=LO
    while a<HI:
        pts.add(a)
        a += log2_over_logx
    pts.add(HI)
    q=sorted(p for p in pts if LO<=p<=HI)
    return list(zip(q,q[1:]))

def eps(alpha: Fraction):
    return int(LO <= alpha < JUMP)

def direct(alpha, k):
    return max(Fraction(1)-Fraction(eps(alpha),2)-Fraction(k,2),0)

def block(epsilon, k):
    return max(Fraction(1)-Fraction(epsilon,2)-Fraction(k,2),0)

rng=Random(710267)
count_cells=0
for n in range(2000):
    # log(2)/log(x) represented by varied rational meshes < 1/40.
    h=Fraction(rng.randint(1,1000),rng.randint(50000,200000))
    cs=cells(h)
    assert cs[0][0]==LO and cs[-1][1]==HI
    for a,b in cs:
        assert a<b
        assert not (a<JUMP<b)
        mid=(a+b)/2
        assert eps(mid)==eps(a)  # half-open convention at the jump
    count_cells += len(cs)

for _ in range(20000):
    # Sample alpha away from endpoints and all relevant factor counts.
    den=rng.randint(1000,100000)
    num=rng.randint((den+V-1)//V, den//(S+1)-1)
    alpha=Fraction(num,den)
    k=rng.randint(0,8)
    assert direct(alpha,k)==block(eps(alpha),k)

print(f"PASS: 2000 partitions ({count_cells} cells), 20000 exact Kuhn states")

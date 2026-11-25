"""Verifica las identidades y contraejemplos del doc 36; no prueba Goldbach."""


def mu_squarefree(primes):
    return -1 if len(primes) % 2 else 1


def coprime_factor_coeffs(primes):
    omega = len(primes)
    mu = mu_squarefree(primes)
    den = 1 << omega
    # Cada subconjunto de los primos define una factorizacion coprima ordenada.
    return [mu / den for _ in range(1 << omega)]


checked = 0
for omega in range(0, 13):
    primes = list(range(omega))  # solo importa el numero de primos distintos
    coeffs = coprime_factor_coeffs(primes)
    assert sum(coeffs) == mu_squarefree(primes)
    assert sum(abs(c) for c in coeffs) == 1
    checked += 1

for m in range(2, 101):
    # a_d=mu(d): cada termino mu(d)*a_d vale 1.
    terms = [1] * m
    E = sum(terms)
    diagonal = sum(t * t for t in terms)
    cross = sum(terms[i] * terms[j] for i in range(m) for j in range(m) if i != j)
    assert E * E == diagonal + cross
    assert diagonal == m and cross == m * (m - 1)
    checked += 1

# Dos terminos de signos opuestos: X<0.
terms = [1, -1]
E = sum(terms)
diagonal = sum(t * t for t in terms)
cross = sum(terms[i] * terms[j] for i in range(2) for j in range(2) if i != j)
assert E == 0 and diagonal == 2 and cross == -2
assert E * E == diagonal + cross
checked += 1

print(f"PASS: {checked} familias; colapso por fibras y ambos signos de X exactos")

# 30 - Serie de Dirichlet del cuello Mertens del termino principal B3

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 07:36 ET, bucle 30 min.
- **Metodo:** factorizacion euleriana exacta del sumando pendiente en (28.5), seguida de auditoria de uniformidad.
- **Estado:** **[T] elemental** para la identidad de series y su dato local en `s=0`; **[P]** para la estimacion uniforme de colas. No es una cota de `B3` y no demuestra Goldbach.

## 1. La serie exacta

Sea `q` par y, para `Re(s)>0`,

```text
A_q(s) = sum_{(e,q)=1} mu(e)/(phi(e) e^s).
```

La serie converge absolutamente y su producto de Euler es

```text
A_q(s) = prod_{p no divide q} (1 - p^(-s)/(p-1)).             (30.1)
```

Separando el polo de zeta,

```text
A_q(s) = zeta(1+s)^(-1) H_q(s),                               (30.2)

H_q(s) = prod_{p no divide q}
           (1 - p^(-s)/(p-1))/(1-p^(-1-s))
         * prod_{p divide q} (1-p^(-1-s))^(-1).               (30.3)
```

Para `q` par, el factor problematico `p=2` esta en el segundo producto. Los factores del primer producto son `1+O(p^(-2))` cerca de `s=0`, de modo que `H_q` es holomorfa en un entorno de cero para cada `q` fijo. Por tanto

```text
A_q(s) = H_q(0) s + O_q(s^2),                                 (30.4)

H_q(0) = 2 C_2 prod_{p|q, p>2} (p-1)/(p-2),                  (30.5)

C_2 = prod_{p>2} (1 - 1/(p-1)^2).
```

Estas identidades son algebraicas a nivel de productos absolutamente convergentes para `Re(s)>0`; (30.4)--(30.5) siguen por continuacion local de los productos de (30.3) y `zeta(1+s)^(-1)=s+O(s^2)`.

## 2. Dependencia en los primos excluidos [T]

El factor local no crece exponencialmente con `omega(q)`. En efecto,

```text
prod_{p|q,p>2} (p-1)/(p-2)
 = prod_{p|q,p>2} p/(p-1)
   * prod_{p|q,p>2} (1-1/(p-1)^2)^(-1)
 << q/phi(q) << log log(3q).                                  (30.6)
```

El segundo producto esta acotado por `C_2^(-1)` y la ultima cota es la cota maximal clasica para `q/phi(q)`. Asi, el valor local en cero solo paga un factor `O(log log q)`. La mera exclusion de los primos de `rN` no crea una perdida de tipo `exp(omega(rN))`.

## 3. Lectura para la suma pendiente (28.5)

Formalmente, derivar (30.2) muestra que la suma amortiguada de
`mu(e) log(e)/phi(e)` tiene limite Abeliano `-H_q(0)`, mientras la de
`mu(e)/phi(e)` tiene limite Abeliano cero. Por ello, para

```text
log(re) = log r + log e,
```

la constante `-H_q(0)` aparece en cada suma inicial pero se cancela al tomar una diferencia diadica `E<e<=y`. No hay un termino principal constante que impida por si mismo cerrar (28.5).

Esto no autoriza a sustituir el limite Abeliano por una cota uniforme de sumas parciales. El lema que realmente falta puede formularse como control uniforme de las colas

```text
S_q(x) = sum_{e>x,(e,q)=1} mu(e)/phi(e),
T_q(x) = sum_{e>x,(e,q)=1} mu(e) log(e)/phi(e),               (30.7)
```

entendidas tras justificar la convergencia ordinaria, para `q=rN`,
`r<=N^delta` y `x>=N^(1/2-delta)/log(N)^B`. Bastaria una cota del tipo

```text
|T_q(x)| + log(2N)|S_q(x)|
  <<_A (log log(3q))^O(1) / log(x)^A                         (30.8)
```

para todo `A`, con constantes uniformes en ese rango. Entonces sumacion parcial junto con (28.4) cerraria la parte `r<=N^delta` del termino principal. (30.8) no se demuestra aqui.

## 4. Avance y limite honesto

- **[T]** Se identifico exactamente la serie de Dirichlet del cuello de (28.5), su cero simple en `s=0` y la dependencia local `2 C_2 prod_{p|q,p>2}(p-1)/(p-2)`.
- **[T]** Esa dependencia es a lo sumo `O(log log q)`; no es el obstaculo principal.
- **[P]** Falta trasladar una region libre de ceros para zeta a (30.8) con uniformidad cuantitativa en los factores locales de `q`. No se invoca esa conclusion sin prueba.
- Aunque se cerrara este termino principal, seguiria abierto el error de dispersion de `B3`, que contiene la barrera de paridad.

## 5. Siguiente paso

Probar o refutar (30.8) mediante un lema de Perron/region libre de ceros con las derivadas de los factores locales controladas uniformemente en `q`; auditar en particular el costo de `sum_{p|q} log(p)/(p-2)`, no solo el valor (30.6).

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto en esta ejecucion. El resultado queda respaldado localmente; no se afirma persistencia remota.
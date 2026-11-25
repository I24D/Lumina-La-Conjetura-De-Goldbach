# 53 - Instanciacion exacta de la sucesion switched para Goldbach: A5 se normaliza, A4 es un BV de casi-primos

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 18:44 ET, bucle 30 min.
- **Metodo:** especializacion literal de Matomaki--Zuniga Alterman, ecuaciones (3.6), (3.11), (3.12), Definicion 2.4 y Assumption 3.1, a la reflexion aditiva de Goldbach.
- **Estado:** **[T-lit]** para las definiciones de la fuente; **[T]** para las identidades finitas 53.1--53.6; **[P-rev]** para delimitar exactamente lo que Bombieri--Vinogradov clasico no entrega. No se prueba A4, una cota de capa ni Goldbach.

Fuente primaria: K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*, Math. Proc. Camb. Phil. Soc. 179 (2025), 351--372, <https://doi.org/10.1017/S0305004125000258> (arXiv:2405.19063v3).

## 1. La pareja reflejada y la sucesion original

Fijemos un par `N>=4` y tomemos `x=N` y

```text
A_N = {(ell,n) in [1,N]^2 : ell+n=N}.                       (53.1)
```

La ecuacion (3.6) de la fuente se especializa sin aproximacion a

```text
a_n = sum_{p<=N, (p,n) in A_N} 1 = 1_{N-n es primo}.        (53.2)
```

Para `N` par, `sum_n a_n = pi(N-1)=pi(N)`. Por tanto podemos elegir
`X_1=pi(N)`. Si `(d,N)=1`, la progresion `d|n` equivale a
`p=N-n == N (mod d)`, y el termino principal natural es
`pi(N)/phi(d)`. Esto corresponde a

```text
g_N(d)/d = 1/phi(d)  si (d,N)=1,
g_N(p)=0             si p|N.                                (53.3)
```

Los primos que dividen `N` producen solo las excepciones locales habituales.

## 2. Formula exacta de b_{P,ell}

Escribamos

```text
V_P(m)=|mu(m)| W_P^+(m) + C_w/log N,
E_P(t)=V_P(m) si t=p_1 m cumple:
       P<p_1<=2P, p|m => p>p_1, S<=Omega(m)<=R-1;
       0 en otro caso.                                      (53.4)
```

La unicidad de `ell=N-t` en (53.1) convierte (3.11) exactamente en

```text
b_{P,ell}=E_P(N-ell),  1<=ell<=N-1.                         (53.5)
```

En particular, para todo `d>=1`,

```text
sum_{d|ell} b_{P,ell}
 = sum_{t<=N-1, t==N (mod d)} E_P(t).                        (53.6)
```

Esta es la forma correcta del problema switched para Goldbach: distribuir
productos casi-primos ponderados `t=p_1m` en la clase `N mod d`. El script
`numerico/verifica_bPell_goldbach.py` comprueba (53.5)--(53.6) por enumeracion
exacta en 252 configuraciones finitas. Es un test algebraico, no evidencia de
distribucion asintotica.

## 3. A5 no es el cuello para la reflexion Goldbach

En A5 aparece el factor

```text
X_1 / sum_{p<=N} 1 = pi(N)/pi(N)=1.                         (53.7)
```

Luego el `X_2(alpha)` requerido es simplemente la masa total de `E_P`, salvo
el posible punto frontera `t=N`, que corresponderia a `ell=0` y aporta a lo
sumo un termino divisor-acotado. Siempre que la masa principal tienda a
infinito, esa frontera es `o(X_2)`. Asi, **la comparabilidad A5 se normaliza
tautologicamente** en esta familia; no requiere cancelacion de primos.

Esto no verifica por si solo todas las hipotesis de Assumption 3.1: aun hay
que fijar uniformemente los terminos de error y las excepciones `p|N`.

## 4. A4: el insumo analitico que queda

Usando (53.6), A4 pide, uniformemente en el bloque `P=N^alpha`, una cota del
tipo

```text
sum_{d<=N^theta} mu(d)^2
 | sum_{t<=N-1, t==N (mod d)} E_P(t)
   - g_N(d) X_2(alpha)/d |
 << X_2(alpha)/(log N)^100.                                 (53.8)
```

**[P-rev]** Bombieri--Vinogradov clasico para la sucesion de primos demuestra
la version de (53.8) para `a_n`, es decir A1 hasta nivel `<1/2`; no demuestra
automaticamente (53.8), porque `E_P` es una convolucion restringida de un
primo pequeno con varios factores primos ordenados y un peso `W_P^+`.

La aparente condicion supra-raiz `q r d > N^{1/2}` de los docs 39--43 ya no
aparece como un unico modulo despues del switching: el modulo externo en
(53.8) es solo `d<=N^theta`. El precio es que la complejidad `q r` queda
dentro de los coeficientes `E_P`. Por tanto, el cuello correcto no es
"BV con modulo qrd", sino un **teorema BV para esta convolucion de casi-primos,
uniforme en alpha y en la condicion de menor factor**. Invocar solo BV para
primos confunde ambas afirmaciones.

Para aplicar el Teorema 3.2 con sus parametros haria falta `theta_2(alpha)=
0.267` en todo `alpha in [1/20,1/4]`, junto con el control uniforme de los
factores locales. Este documento no prueba esa estimacion.

## 5. Veredicto y siguiente paso

- **[T]** `b_{P,ell}` es exactamente la reflexion `E_P(N-ell)` y sus sumas
  divisibles son clases `t==N mod d`.
- **[T]** Con `X_1=pi(N)`, el cociente de A5 es exactamente uno; la unica
  discrepancia estructural es el extremo `t=N`.
- **[P-rev]** A4 no es una consecuencia literal de BV para primos. Es un BV
  uniforme para la convolucion restringida `E_P`; el switching elimina el
  modulo compuesto supra-raiz solo trasladando su dificultad a coeficientes.
- **Siguiente:** descomponer `E_P` en convoluciones de numero fijo de primos
  mediante Buchstab y localizar una version publicada del principio de
  induccion de Bombieri para verificar (o limitar) el nivel `0.267`.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto. No se afirma persistencia remota.

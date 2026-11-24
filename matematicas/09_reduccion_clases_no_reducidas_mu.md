# 09 — Reduccion exacta de clases no reducidas para sumas de Mobius

## Procedencia y estado

- Agente: LUMINA/OpenClaw.
- Fecha: 2026-11-24, America/New_York.
- Metodo: revision adversarial de 06 y 08; derivacion algebraica exacta.
- Estado: **[P]** para la identidad de reduccion. La estimacion agregada final
  sigue **ABIERTA**; no se eleva el Lema 6.2 a demostrado.
- Alcance: aislar con precision el hueco detectado por Codex en 08, seccion 2.3.

## 1. El problema

El paso cuestionado en 06 necesita controlar, para clases posiblemente no
reducidas,

\[
M(y;q,r)=\sum_{n\le y,\ n\equiv r\pmod q}\mu(n).
\]

Los enunciados de Bombieri--Vinogradov localizados para `mu` controlan clases
`(a,q)=1`. No es valido sustituir sin prueba el maximo reducido por un maximo
sobre todo `r`.

## 2. Identidad exacta de descenso de conductor [P]

Ponga `g=(r,q)`, `q=gh`, `r=gb`, de modo que `(b,h)=1`.

**Lema 9.1.** Si `g` no es libre de cuadrados, `M(y;q,r)=0`. Si `g` es libre
de cuadrados, entonces

\[
M(y;q,r)=\mu(g)
\sum_{\substack{k\ge1\\p\mid k\Rightarrow p\mid g\\(k,h)=1\\gk\le y}}
\sum_{\substack{u\le y/(gk)\\u\equiv b\bar k\pmod h}}\mu(u). \tag{9.1}
\]

Aqui `bar k` es el inverso de `k` modulo `h`; la clase interior es reducida.

**Prueba.** Todo `n` de la clase es `n=gm`, con `m=b (mod h)`. Si `g` tiene
un factor cuadrado, `mu(n)=0`. Para `g` libre de cuadrados,

\[
\mu(gm)=\mu(g)\mu(m)1_{(m,g)=1}.
\]

La identidad de series de Dirichlet finitas

\[
\mu(m)1_{(m,g)=1}
=\sum_{\substack{k\mid m\\p\mid k\Rightarrow p\mid g}}\mu(m/k) \tag{9.2}
\]

se verifica prima por prima: para `p|g`, `(1-X)(1+X+X^2+...)=1`, y para
`p` fuera de `g` queda el factor `1-X`. Escribiendo `m=ku`, la congruencia
`ku=b (mod h)` no tiene soluciones si `(k,h)>1`, porque `(b,h)=1`; si
`(k,h)=1` equivale a `u=b bar(k) (mod h)`. Esto prueba (9.1). □

La condicion `(m,g)=1` no puede borrarse. Por ejemplo, para `g=m=2`,
`mu(gm)=mu(4)=0`, mientras `mu(g)mu(m)=1`.

## 3. Lo que falta para cerrar el Lema 6.2

La identidad (9.1) convierte cada clase no reducida en clases reducidas, pero
produce dos parametros nuevos: `g` y un numero `g`-liso `k`. El paso restante
no es algebraico sino una estimacion agregada uniforme.

Una entrada suficiente, escrita sin esconder cuantificadores, es la siguiente.
Para cada `A` debe existir `B` tal que, con
`Q<=sqrt(x)/(log x)^B`, se tenga

\[
\sum_{g\le Q}^{\squarefree}
\sum_{h\le Q/g}
\sum_{\substack{k:\ p|k\Rightarrow p|g\\(k,h)=1}}
\max_{z\le x/(gk)}\max_{(a,h)=1}
\left|\sum_{u\le z,\ u\equiv a\pmod h}\mu(u)\right|
\ll_A \frac{x}{(\log x)^A}. \tag{WNR}
\]

El rango de `k` es finito implicitamente (`gk<=x`). **[C: WNR]** La entrada
(WNR), aplicada con `x` del orden de `K`, repara exactamente la linea 168 de
06 y valida sus terminos Tipo I, despues de la sumacion parcial ya escrita.

## 4. Por que no se declara (WNR) demostrado

Para `k` pequeno, (9.1) sugiere aplicar BV reducido a longitud `x/(gk)` y
modulo `h`. Pero el nivel permitido cambia con `gk`; no basta aplicar una sola
vez BV a longitud `x`. Para `k` grande hay que demostrar una cola uniforme de
numeros `g`-lisos y sumar tambien los casos donde una progresion contiene cero
o un termino.

Los pesos armonicos son prometedores: localmente,

\[
\sum_{k:\ p|k\Rightarrow p|g}\frac1k
=\prod_{p|g}(1-1/p)^{-1},
\]

y al combinar la eleccion de `g|q` aparecen productos de orden polilogaritmico.
Sin embargo, convertir esta observacion en (WNR) exige una particion diadica en
`gk`, verificar el nivel `h<=sqrt(x/(gk))/log^B`, y controlar la cola sin perder
un factor `x^{o(1)}`. Ese libro de cuentas no se completa aqui.

## 5. Dictamen

1. La objecion de Codex a 06 es correcta: la frase "se absorbe con un divisor
   de q" no era una prueba.
2. La reduccion a clases reducidas si existe y queda demostrada en (9.1).
3. El hueco se reduce ahora a la desigualdad explicita (WNR), no a una vaga
   version de BV para clases arbitrarias.
4. Aun cerrando (WNR), el termino Tipo II `B_3` de 06 permanece abierto y
   contiene la barrera principal. Por tanto este avance no demuestra Goldbach.

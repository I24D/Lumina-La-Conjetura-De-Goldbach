# 85 - Primera suma de switching de Chen con factor exterior Q
## Procedencia y estado
- **Agente:** LUMINA/Muse. **Fecha:** 2026-11-26 13:55 ET.
- **Metodo:** escritura explicita de la suma de switching para la subincidencia `u=1` con selector `Q_N` (docs 82-84) desde la presentacion primaria de Chen, y auditoria de la estimacion bilineal que realmente necesita.
- **Fuente primaria:** Chen Jing-run, "On the representation of a larger even integer as the sum of a prime and the product of at most two primes", *Sci. Sinica* 16 (1973), 157-176; presentacion moderna en Nathanson, *Additive Number Theory: The Classical Bases* (1996), cap. 10 (teorema de Chen y switching); criba lineal en Halberstam-Richert, *Sieve Methods* (1974).
- **Estado:** **[T]** identidades finitas y lema de pocos divisores (sec. 3); **[P-rev]** transferencia de pesos y constantes de Chen a `a_Q`; **[N]** el Lema SW-Q (sec. 5) no se tiene. No se demuestra Goldbach; no hay victoria parcial.
## 1. Marco
`N` par grande, `k>=2` fijo, `Q=Q_N` selector de (83.1) (`Q<<_k (log N)^(k-1)`, `(Q,N)=1`),
subincidencia `u=1`:
```text
X = N/Q,   a_Q(v) = Lambda(N-Qv),   1 <= v < X,   (v,Q)=1.
J_k(N;Q) = sum_{v<X,(v,Q)=1,omega(v) in {1,2}} 1_{(N-Qv,N)=1}
           a_Q(v) (log Q + log rad(v)).                 (85.1)
```
Parametros de Chen (los de la fuente; `log X ~ log N`): `z = X^(1/10)` (limite inferior de criba),
`y = X^(1/3)` (corte del switching), `D = X^(1/2)/log^C X` (nivel de criba, `C` a fijar), `P_Q(t) =
prod_{l<t, l no | QN} l`. Serie singular local (82.6): `S(N,Q)=2*C_2*prod_{l|NQ,l>2}(l-1)/(l-2)`,
positiva.
## 2. La primera suma de switching con Q
Para un primo `s`, sea `A_s` la sucesion `w -> a_Q(s*w)` con `w < X/s`, `(w,Q)=1`, y `S_Q(A_s;s)` su
conteo cribado por los primos `l<s`, `l no | QN`:
```text
S_Q(A_s;s) = sum_{w<X/s, (w,P_Q(s))=1, (w,Q)=1} a_Q(s*w).  (85.2)
```
**Definicion (primera suma de switching).**
```text
B(Q) = sum_{z<=s<y, (s,QN)=1} S_Q(A_s;s).                (85.3)
```
Cada `v=p1*p2` con `z<=p1<=p2`, `p1<y`, aparece exactamente una vez en `B(Q)`, en `s=p1`, `w=p2`
(primo, `>=s`, coprimo con `P_Q(s)`). Los demas `v=p1*p2` se parten en dos colas: `p1<z`, cribados
fuera de `S(A,P_Q(z))` (viven en el termino principal de la criba, no en `B(Q)`); y `p1>=y` (ambos
factores `>=y`), cola `R_y`, a acotar aparte por criba superior (Buchstab).
La desigualdad objetivo de Chen, con constante `c` a fijar contra la fuente primaria, es de la forma
```text
J_k^{P2} >= c * [ S_Q(A;z) - B(Q) ] - R_y - (errores),   (85.4)
```
`S_Q(A;z)` el conteo cribado de toda `a_Q`, `J_k^{P2}` la parte `omega(v)=2` de (85.1). **[P-rev]**:
el valor exacto de `c` (1 o 1/2 segun la presentacion) y el reparto exacto entre `S_Q(A;z)`, `B(Q)`
y `R_y` deben fijarse leyendo Nathanson cap. 10; este doc escribe la suma, no certifica la
constante. Peso `log Q + log rad(v)`: para `v=s*w`, `log rad(v)=log s+log rad(w)`; el `log s` sale
fuera de la suma en `w` y el `log rad(w)` se maneja por sumacion parcial diadica en `w` (coste
`log^{O(1)}`). **[P-rev]** ese desdoblamiento del peso.
**Termino principal (esquema de criba lineal).** Con nivel `D/s` para `A_s` y `u_s = log(D/s)/log
s`, la criba lineal superior da
```text
S_Q(A_s;s) <= (N/phi(Qs)) * (e^{-gamma} F(u_s)+o(1))/log s + R_s,
M_B(Q) = sum_{z<=s<y,(s,QN)=1} (N/phi(Qs)) * e^{-gamma} F(u_s)/log s.
                                                        (85.5)
```
`F` la funcion de la criba lineal (Halberstam-Richert). **[P-rev]**: (85.5) es el esquema publicado;
los valores exactos de `F(u_s)` y `u_s` con el `Q` exterior deben leerse de la fuente.
**Error.** Expandiendo (85.2) con los pesos `lambda_d(s)` de la criba lineal (`|lambda_d(s)|<=1` en
la presentacion usual):
```text
R_B(Q) = sum_s sum_{d<D/s} lambda_d(s) * E(N;Qsd,N),     (85.6)
E(N;m,N) = psi(N;m,N) - N/phi(m),
```
usando (84.2): `sum_{sd|v} a_Q(v)=psi(N;Qsd,N)-Lambda(N)`, y el termino `-Lambda(N)` aporta `O(log
N)` en total. Modulos `m=Qsd` con `sd<=D`, luego `m<=QD<<N^{1/2}/log^{C'}`.
## 3. Auditoria: que estimacion bilineal necesita realmente
**(a) No hace falta Tipo II nuevo.** En la prueba de Chen segun Nathanson, `B(Q)` se estima por
arriba con la criba LINEAL mas el resto (85.6); no se pide una asintotica de `B(Q)` ni un nivel
`>1/2`. La bilinealidad que aparece es la de Bombieri-Vinogradov misma (gran criba), no una
descomposicion de Vaughan adicional. **[P-rev]**: confirmar contra Nathanson cap. 10 que ningun paso
oculto (cola `R_y`, valores `w` con `omega(w)=2` dentro de `S_Q(A_s;s)`) exige un Tipo II; para una
cota SUPERIOR esos `w` se cuentan sin mas, sin descomponerlos.
**(b) El punto real: no aplicar BV primo a primo.** Aplicar (84.3) por separado para cada `s`
costaria un factor `#s ~ y/log y`, con `y=X^{1/3}`: fatal. Lo que salva la suma es hacerlo
CONJUNTAMENTE sobre `(s,d)` con modulo `m=Qsd`, usando que el primo de switching es grande. **Lema
de pocos divisores [T].** Sea `m=Qsd` con `s` primo en `[z,y)`, `sd<=D<=X^{1/2}`. El numero de
primos `s' in [z,y)` que dividen `m/Q` es `<= log(m/Q)/log z <= 5`, porque cada uno es
`>=z=X^{1/10}` y su producto es `<=X^{1/2}`. Luego el coeficiente del modulo `m` en (85.6), `mu_m =
sum_{s'|m/Q} lambda_{m/(Qs')}(s')`, satisface `|mu_m|<=5` si `|lambda_d(s)|<=1`. Asi (85.6) es una
suma de forma BV (84.3) con coeficientes acotados por 5, y BV da `R_B(Q)<<_{A,k} N/log^A N`. **[T]**
como reduccion de conteo, CONDICIONADA a que los pesos publicados cumplan `|lambda_d(s)|<=1` y nivel
`D/s` con `Q` insertado.
**(c) Modulo excepcional a ESTE nivel: no aparece uno nuevo.** BV suma sobre TODOS los modulos
`m<=N^{1/2}/log^C`, incluido un posible `m0` excepcional; su contribucion a (85.6) es
`<=5*|E(N;m0,N)|<=5*sum_m max_a |E|<<N/log^A N`. No se necesita ningun analisis de cero de Siegel
adicional en el switching de Tipo I, y el termino principal `N/phi(m)` no depende de ningun
caracter. La capa BV quedo cerrada en el doc 84; el switching no reabre ese frente. **[T]** dada la
reduccion (b).
**(d) Detalle pendiente menor [P-rev]: la pre-criba `(v,Q)=1`.** Los primos `l|Q` son `k-1` primos
`<<log N < z` (para `N` grande): imponer `(v,Q)=1` es pre-cribar por `<=k-1` primos pequenos, lo que
solo cambia la serie singular (ya en `S(N,Q)`) y los pesos `lambda_d(s)` en un factor acotado. Debe
verificarse que insertar esa pre-criba antes del switching no cambia la forma (85.6); no se espera
obstruccion, pero no esta verificado. El script (sec. 6) comprueba la expansion plana (sin
pre-criba).
## 4. La subincidencia u=1: evita, no esconde
Para el objetivo `C_k+C_{k+1}>0`, `u=1` EVITA el problema de uniformidad en `A=Qu` (doc 82, sec. 2):
aqui `A=Q` exactamente, polilogaritmico, y todos los pesos son no negativos, asi que `J_k(N;Q_N)>0`
basta. Es una restriccion legitima, no un truco: se descarta masa (todo `u>1`), pero la positividad
de una subincidencia positiva implica la de la union. Para la incidencia COMPLETA `I_k`, en cambio,
`u=1` solo lo ESCONDE: `u` con `rad(u)|Q` puede llegar hasta `N/Q` (p. ej. `u=3^a`), luego `A=Qu`
recorre hasta `N` y ninguna uniformidad en `Q` solo bastaria para `I_k`. Este doc no afirma nada
sobre `I_k` completa; su alcance es `J_k` con `u=1`. **[T]** esta delimitacion; **[N]** una cota de
`I_k` por esta via.
## 5. Lema faltante (enunciado exacto)
**LEMA SW-Q (transferencia del switching de Chen) [N].** *Hipotesis.* `N` par, `N>=N0(k)`; `Q=Q_N`
de (83.1); `X=N/Q`; `z=X^{1/10}`, `y=X^{1/3}`; `D=X^{1/2}/log^{C0} X` con `C0` absoluto a fijar;
`lambda_d(s)` los pesos de la criba lineal superior de Halberstam-Richert aplicados a `A_s` con
nivel `D/s` y limite de criba `s`, con el factor exterior `Q` insertado en cada modulo.
*Conclusion.* (a) `|lambda_d(s)|<=1` y `sum_s sum_{d<D/s} lambda_d(s) E(N;Qsd,N) <<_{A,k} N/log^A N`
para todo `A`; (b) `B(Q)<=M_B(Q)(1+o(1))+O(N/log^A N)` con `M_B(Q)` de (85.5) y `F`, `u_s`
explicitos; (c) la constante `c` de (85.4) y la cota de la cola `R_y` quedan fijadas por la fuente
primaria. *Que haria falta para probarlo.* (i) Leer Nathanson cap. 10 y verificar peso a peso que
`lambda_d(s)` cumple (a) con `Q` insertado (acotacion, dependencia en `s`, nivel `D/s`). (ii)
Evaluar (85.5) con `F` explicita y la serie `S(N,Q)`. (iii) Acotar `R_y` por criba superior. Sin
(i)-(iii) no hay transferencia; la reduccion de la sec. 3(b) es necesaria pero no suficiente:
convierte el resto en forma BV, no certifica que los pesos de Chen tengan esa forma.
## 6. Verificacion numerica (control finito, no evidencia asintotica)
Script `numerico/85_primera_suma_switching_chen_Q.py`: `N` par en `[100,2000]`, `k=2,3`, selector
canonico. H1: reindexado exacto `J2=sum_s T_s + R_y + masa(p1<z)`; H2: fraccion de cola `R_y/J2`;
H3a: inversion exacta de Mobius de `U_s` sobre `d|P_Q(s)`; H3b: lema de pocos divisores; H3c: error
truncado `E_s=U_s^trunc-M_s` con nivel `D_s=(X/s)^{1/2}` (sin `log^C`: a esta escala el nivel con
perdida logaritmica degenera, `<1`).
```text
PASS H1=1902 H2=461 H3a_mobius=27 H3b_pocos_divisores=27 n_bilin=27
max_tail_frac_Ry/J2=1.000000 worst=(214, 2, ...)
max_bilin_rel=0.894506 worst=(1760, 2, 7, ...)
sum_abs_E=3457.390 sum_abs_M=3996.667 ratio=0.865068
```
Lectura honesta: H1/H3a/H3b pasan al 100% (identidades exactas y conteo). Los ratios `~0.87-1.0` son
ARTEFACTOS DE ESCALA: a `N<=2000`, `y=X^{1/3}<=12`, casi toda la masa `P2` cae en la cola `R_y`, y
en la suma truncada sobrevive casi solo `d=1` (la condicion `(d,QNs)=1` elimina los demas `d<=D_s`),
asi que el error es de orden del termino principal. El script controla la FORMA del error (que es
`sum mu(d) E(N;Qsd,N)`), no su tamano asintotico; no apoya ni refuta el Lema SW-Q.
## Siguiente paso
Leer Nathanson cap. 10 peso a peso y certificar o refutar el apartado (a) del Lema SW-Q: acotacion y
nivel de `lambda_d(s)` con `Q` insertado. Si los pesos no admiten la forma (85.6) con `|lambda|<=1`
y nivel `D/s`, el fallo estara ahi, no en un modulo excepcional. Goldbach sigue abierta.

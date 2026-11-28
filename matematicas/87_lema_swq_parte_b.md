# 87 - Lema SW-Q apartado (b): termino principal M_B(Q) y la bifurcacion del limite
## Procedencia y estado
- **Agente:** LUMINA/Muse. **Fecha:** 2026-11-26 16:26 ET (reclamo) / 16:40 ET (entrega).
- **Metodo:** lectura de la version explicita del Teorema 9.7 de Nathanson en Bordignon-Johnston-Starichkova, arXiv:2207.09452, Teorema 6 (ecs. (5)-(11)), mas Nathanson caps. 9-10 segun el doc 86, y calculo directo (integracion numerica propia de la ecuacion de retardo de F; fuerza bruta finita). Fuentes de control: Halberstam-Richert, *Sieve Methods* (1974); Yamada, arXiv:1511.03409.
- **Estado:** **[T]** apartado (b) CERTIFICADO en la version de limite FIJO z (Nathanson, Teorema 10.5). La version de limite VARIABLE s del doc 85 queda ABANDONADA como vehiculo de M_B(Q): no certificada en 5/14 del rango. **[N]** (c) no se intenta aqui. No se demuestra Goldbach; no hay victoria parcial.

## 1. La bifurcacion: donde se rompe cada version
Parametros (docs 85-86, proxy D = X^{1/2} sin log^C como en docs 85/86): X = N/Q, z = X^{1/10}, y = X^{1/3}, D = X^{1/2}/log^{C0} X con C0 = C0(A,k) (correccion del doc 86).

**Prop. 87.1 [T] (el limite variable se rompe en s = D^{1/2}).** La cota superior de Jurkat-Richert exige nivel >= limite de criba: Bordignon Teorema 6, ec. (5): "for any D >= z we have the upper bound", con s := log D/log z. Aplicada a A_s con limite s y nivel D/s, la hipotesis es D/s >= s, es decir, s <= D^{1/2} ~ X^{1/4}. Para s > X^{1/4} se tiene u_s^var = log(D/s)/log s < 1 y el Teorema 9.7 NO aplica: la formula F(u) = 2e^gamma/u existe para u < 1, pero el teorema no la cubre ahi. No es un detalle tecnico: con nivel < limite, la criba no ve los primos en [D/s, s), y el termino principal (85.5) con F(u_s) es una extrapolacion fuera del rango del teorema. Medicion (H1): el rango roto es (X^{1/4}, X^{1/3}), es decir, 5/14 ~= 35.7% del rango logaritmico [(1/3-1/4)/(1/3-1/10)]; y lleva ~1/3 de la masa del termino principal (H5: fraccion 0.317/0.363/0.343 a X = 10^6/10^9/10^12 para el peso 1/(phi(s) log(D/s))). No es una cola despreciable.

**Prop. 87.2 [T] (el limite fijo funciona en todo [z,y)).** Con limite fijo z y nivel D/s, u_s^fix = log(D/s)/log z recorre [5/3, 4]: en s = z, u = log(D/z)/log z = 4 - o(1); en s = y, u = log(D/y)/log z = 5/3 - o(1). Luego D/s >= z siempre (minimo 5/3 > 1) y el Teorema 9.7 aplica en todo el rango. (H1: u^fix medido en [1.6725, 3.9942] en la malla; los extremos analiticos son 5/3 y 4.)

**Prop. 87.3 [T] (dominio: basta el limite fijo para una cota superior).** S_Q(A_s;s) <= S_Q(A_s;z), porque P_Q(z) | P_Q(s): el conjunto cribado por s es subconjunto del cribado por z. Es exactamente el movimiento de Nathanson en el Teorema 10.5 (criba A_q con P(z)). (H3: exacto, conteo y peso Lambda >= 0, 4/4.)

**Prop. 87.4 [T] (coincidencia: el limite no cambia el termino principal donde ambos valen).** Para u <= 3, F(u) = 2e^gamma/u, luego
```text
(U_N/log s) * F(u_s^var) = (U_N/log z) * F(u_s^fix) = 2 S(N)/log(D/s),
```
porque u*log(limite) = log(D/s) en ambas versiones. La ventaja aparente del limite variable (denominador log s menor) se cancela EXACTAMENTE con su F mayor (u^var < u^fix y F decreciente). La eleccion del limite es inmaterial para el termino principal donde ambos son validos; solo importa para la validez. (H4: 1497/1497, error relativo maximo 2.55e-16.)

**Prop. 87.5 [T] (truncar no salva la version variable).** Restringir la suma a s <= X^{1/4} deja fuera ~1/3 de la masa (H5). La cola s > X^{1/4} con la cota trivial S(A_s;s) <= |A_s| ~= N/phi(Qs) sumaria N * sum_s 1/(s log s) >> N/log^A N y destruiria la forma (1+o(1)) de (b). La version variable no se repara: se abandona.

**Decision.** Limite FIJO z (Nathanson). El doc 85 usaba el limite variable por seguir la definicion (85.2) al pie de la letra; esa definicion sigue valiendo para B(Q), pero su termino principal (85.5) no.

## 2. F(u_s) explicita
**Prop. 87.6 [T].** De la ec. (8) de Bordignon (Teorema 6): F(s) = 2e^gamma/s para 0 < s <= 2 y (sF(s))' = f(s-1), (sf(s))' = F(s-1) para s >= 2, con f(s) = 0 para 0 < s <= 2. De ahi: F(u) = 2e^gamma/u en [1,3] (en [2,3], (uF)' = f(u-1) = 0, luego uF(u) se queda en 2e^gamma); y para u en (3,4],
```text
F(u) = (1/u) * [2e^gamma + int_3^u f(t-1) dt],  f(v) = 2e^gamma log(v-1)/v  (2 <= v <= 4),
```
donde f(v) sale de (vf(v))' = F(v-1) = 2e^gamma/(v-1). Valores: F(5/3) = 6e^gamma/5 ~= 2.13729, F(3) = 2e^gamma/3 ~= 1.18738, F(4) ~= 1.02164 (dos integraciones independientes coinciden a 1e-9: Simpson de la ecuacion de retardo y la forma (e^gamma/2)(1 + int_1^2 log t/(t+1) dt)). F es decreciente en [1,4] (F' = (f(u-1)-F(u))/u < 0: f <= f(3) ~= 0.824 < F(4) ~= 1.022 en (3,4]). El exceso sobre la hiperbola 2e^gamma/u es <= 1.1472, y solo ocurre en s < X^{1/5} (donde u^fix > 3, pues u^fix = 3 <=> log(D/s) = 3 log z <=> s = X^{1/5}). (H2.)

## 3. El termino principal corregido (85.5)
**Correccion [T].** La forma (85.5) del doc 85 omitia la serie singular: el termino principal de la criba es |A| * V(limite) * F(u), y V(z) ~ U_N/log z con U_N = e^{-gamma} S(N) = 2e^{-gamma} prod_{p>2}(1-1/(p-1)^2) prod_{p|N,p>2}(p-1)/(p-2) (Bordignon ec. (2); vease su ec. (119) para W(x)). La cancelacion del factor Q es exacta: |A_s^Q| ~= (N/phi(Qs)) * prod_{l|Q}(l-2)/(l-1) por Mobius sobre e|Q (doc 86 (86.4)), mientras V_Q(z) = V(z) * prod_{l|Q}(l-1)/(l-2) (los primos l|Q, todos < z, no se criban); el producto es (N/phi(Qs)) * V(z). El e^{-gamma} de (85.5) era la mitad de U_N: faltaba S(N) >= 2C_2 ~= 1.32.

**Prop. 87.7 [T] (M_B(Q) explicito, version fija).**
```text
M_B(Q) = sum_{z<=s<y, (s,QN)=1} (N/phi(Qs)) * (U_N/log z) * F(u_s),
u_s = log(D/s)/log z en [5/3, 4],  F como en la Prop. 87.6.     (87.7)
```
En el rango u_s <= 3 (s >= X^{1/5}) el sumando es (N/phi(Qs)) * 2S(N)/log(D/s); en s < X^{1/5} es a lo sumo 1.1472 veces esa forma. La forma (85.5) queda reemplazada por (87.7).

## 4. La desigualdad (b)
**Prop. 87.8 [T] ((a) es independiente del limite).** Con limite fijo z, los pesos lambda_d(s) tienen soporte d | P(z), d < D/s: |lambda| <= 1 sigue del Teorema 9.3 (cualquier limite), y la multiplicidad 1 de la Prop. 86.3 sigue valiendo (usa solo que s es el mayor factor primo de m/Q en [z,y)). La forma exacta del resto (86.2) y sus correcciones C2, C3 no dependen del limite. (H6: 8/8 modulos con multiplicidad 1.) Luego sum_s R_s << N/log^A N por el mismo argumento del doc 86 sec. 4 (BV + C0 = C0(A,k)).

**Prop. 87.9 [T] (desigualdad (b)).**
```text
B(Q) <= sum_s S_Q(A_s;z)
     <= sum_s [|A_s^Q| * V_Q(z) * (F(u_s) + o(1)) + R_s]
      = M_B(Q) * (1 + o(1)) + O(N/log^A N).                    (87.8)
```
El primer paso es la Prop. 87.3; el segundo es el Teorema 9.7 (Bordignon Teorema 6) con D/s >= z (Prop. 87.2) y resto de la Prop. 87.8. El o(1) es uniforme en s: u_s recorre el compacto [5/3,4], F es continua ahi, y la condicion de densidad (9.34) no depende de s (doc 86 H4); el error relativo de V_Q(z) = (U_N/log z)(1+o(1)) es O(1/log^2 z) (Bordignon (119)).

**Veredicto (b): CERTIFICADO** en la version de limite fijo z, con M_B(Q) de (87.7) y F, u_s explicitos. La version de limite variable del doc 85 queda abandonada (Prop. 87.1, 87.5).

## 5. Verificacion numerica (control finito, no evidencia asintotica)
Script `numerico/87_lema_swq_parte_b.py`: proxy D = X^{1/2} sin log^C (como docs 85/86); H3/H6 con parametros sinteticos no degenerados (N = 10^6, Q = 77, z = 7; a N <= 2000, z = 2 y P_Q(z) es vacio).
```text
PASS H1=4/4 H2 (F(4)=1.021641553, alt identico a 1e-9) H4=1497/1497
H5 (masa rota 0.317/0.363/0.343 a X=10^6/10^9/10^12) H3=4/4 H6=8/8
```
Lectura honesta: H1/H2/H4/H5 miden identidades y rangos (la coincidencia H4 vale a precision de maquina, 2.55e-16); H3/H6 son exactos por fuerza bruta. El script controla la FORMA y los rangos, no el tamano asintotico; el tamano es BV + Teorema 9.7 (analitico). Los numeros deben ser IDENTICOS local y en el PC.

## Siguiente paso
(c) Constante c de (85.4) y cota de la cola R_y. Goldbach sigue abierta.

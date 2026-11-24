# 18 - Revision adversarial del doc 17: Lichtman Prop. 3.1 y el termino principal de B3

## 0. Procedencia y alcance

- Agente: LUMINA/OpenClaw.
- Modelo: identificador exacto no expuesto por el runtime OpenClaw.
- Fecha: 2026-11-24, America/New_York.
- Metodo: bucle 30 min; contraste de `17` con `06` y con la fuente primaria, Lichtman, arXiv:2309.08522v1, Def. 1.4 y Prop. 3.1.
- Estado: **[P-rev]**. Revision no independiente (mismo usuario Windows). No demuestra Goldbach.

## 1. Veredicto ejecutivo

**[P-rev] La obstruccion C7 del doc 17 es valida y basta para impedir la aplicacion directa de Lichtman Prop. 3.1 a B3 con el peso de modulos `lambda_d = mu(d) log d`.**

**[P-rev] La afirmacion C10 de que el termino principal de B3 por caja es `O_A(N/log^A N)` no queda demostrada por el argumento escrito.** El dato sobre la suma no ponderada `sum mu(d) log d` no controla automaticamente el termino principal real, que lleva un peso dependiente de `d`. Esto no reabre la ruta de Lichtman, porque C7 ya la cierra; si corrige el alcance de una afirmacion auxiliar.

## 2. C7: verificacion contra la definicion primaria

La Def. 1.4 de Lichtman exige que, para **cada** factorizacion `Q=Q1 Q2 Q3`, exista una convolucion triple de sucesiones 1-acotadas, con soportes `qi <= Qi`, que reproduzca `lambda_q`.

Tomemos `Q1=Q2=Q3=Q^(1/3)`. Para cualquier primo `p` con `Q^(1/3) < p <= Q`, toda factorizacion entera `p=q1 q2 q3` tiene uno de los factores igual a `p`; por tanto no puede satisfacer simultaneamente `qi <= Q^(1/3)`. La convolucion triple vale cero en `p`, mientras

`lambda_p = mu(p) log p = -log p != 0`.

Contradiccion. El mismo argumento ya impide la factorizabilidad doble tomando `Q1=Q2=sqrt(Q)` y un primo `sqrt(Q)<p<=Q`. Dividir el peso por una potencia de `log N` no cambia el soporte y no elimina la contradiccion.

**Conclusion [P-rev]:** C7 es correcta como obstruccion de soporte, no solo como una estimacion insuficiente. La parte prima del peso de B3 es precisamente incompatible con la hipotesis estructural del teorema.

## 3. Prop. 3.1: otras hipotesis comprobadas

La fuente primaria confirma que Prop. 3.1 requiere:

1. `NM asymp x` y una variable corta entre `x^eps` y `x^(1/3+eps)`;
2. un peso de modulos triplemente bien factorizable;
3. una sucesion corta Siegel-Walfisz y soportada en numeros con `P^-(n) >= z0`;
4. un intervalo fijo `I subset [x,2x]`;
5. nivel `min((5-4 theta)/(8-6 theta),(5-alpha theta)/8)-50 eps`.

Con `alpha=1` y `theta=7/32`, el segundo brazo es `153/256`, como dice `17`. El teorema es una estimacion firmada en `q`; no contiene valor absoluto exterior sobre cada modulo. Por ello C1 y C8 del doc 17 estan bien leidas.

C2 y C3 siguen siendo obstrucciones adicionales plausibles para la forma de B3 escrita en `06`, pero no son necesarias para el veredicto: C7 ya es concluyente.

## 4. Correccion de alcance a C10

Tras sustituir la progresion por su densidad `1/phi(d)`, el termino principal natural de una caja B3 no es una constante por

`sum_{d~K} mu(d) log d`.

Tiene la forma esquematica

`M_B3 = sum_{d~K,(d,N)=1} mu(d) log d / phi(d) * W_N,C,A,B(d)`,

con

`W_N,C,A,B(d) = sum_{a~A,b~B, N-2Cd < ab <= N-Cd} Lambda(a) nu_U(b)`

(mas los cortes y coprimalidades de la descomposicion). El intervalo movil hace que `W(d)` dependa de `d`; ademas `W` contiene el peso de Tipo II `Lambda * nu_U`. Para deducir ahorro arbitrario de logaritmos mediante sumacion parcial hace falta demostrar, uniformemente por caja, una cota de variacion o una descomposicion suave de `W(d)` compatible con las cotas fuertes para las sumas parciales de `mu(d) log d`. El doc 17 no proporciona ese lema.

La medicion de

`sum_{K<d<=2K} mu(d) log d`

solo corresponde al peso constante `W(d)=1`. No prueba que `M_B3` sea pequeno para el peso movil anterior. En particular, cancelacion numerica no ponderada no puede reemplazar el control analitico de `W(d)`.

**Correccion [P-rev]:** reemplazar en C10 "el termino principal de B3 por caja es `O_A(N/L^A)`" por "esa pequenez es una hipotesis auxiliar aun no justificada para el peso movil real". La conclusion mayor de `17` permanece: Lichtman Prop. 3.1 no se aplica directamente, por C7.

## 5. Siguiente paso acotado

Escribir `M_B3` con todos los cortes diadicos y decidir si `W(d)` admite una particion suave con variacion total controlable. Si no, registrar el termino principal ponderado como una segunda manifestacion del nucleo de paridad, no como error despreciable.

## Fuente primaria

- Jared Duker Lichtman, *Primes in arithmetic progressions to large moduli, and Goldbach beyond the square-root barrier*, arXiv:2309.08522v1: Def. 1.4 y Prop. 3.1.

**Estado final:** [P-rev]. Goldbach fuerte sigue abierta.
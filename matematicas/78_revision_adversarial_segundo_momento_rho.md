# 78 - Revision adversarial del segundo momento del minorante rho

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-26 00:37 ET, bucle 30 min.
- **Metodo:** auditoria adversarial del doc 77 contra las definiciones de los docs 66--67; contraejemplo finito reproducible.
- **Estado:** **[P-rev]** para la Prop. 77.1; **[T]** para la insuficiencia logica de la informacion citada. No se demuestra Goldbach.

## 1. El paso no justificado

La prueba de la Prop. 77.1 usa

```text
sum_n rho(n)^2 << X/log X.                                  (78.1)
```

pero el unico valor aportado es el calculo

```text
sum rho_tilde(n)^2 = 1.467*pi(2X),
rho_tilde = 1_primo - S_{2,3} - S_{2,5},
```

para el **proxy sin parte B**. El minorante verdadero de Maynard/Harman es

```text
rho = S^triv + S^G - C*rho_B
    = 1_primo - (S^B + C*rho_B + S_{2,3} + S_{2,5}).
```

Los docs 66--67 citan para la parte B una cota de masa `O(delta)` y
equidistribucion, pero no una estimacion cuadratica de `rho_B`, ni una cota
puntual que implique (78.1). Por tanto sustituir el valor medido del proxy por
la norma L2 del minorante verdadero no esta justificado por las fuentes
registradas. La frase del doc 67 de que todo se transfiere con factor
`1+O(delta)` es valida para masas L1, no automaticamente para cuadrados.

## 2. Por que masa pequena no basta [T]

Sea `P` la masa prima y `M` la masa total de un defecto no negativo sobre
compuestos. Dos defectos con la misma masa pueden ser:

```text
D_spread = 1 en M posiciones,       sum D_spread^2 = M;
D_spike  = M en una posicion,       sum D_spike^2  = M^2.
```

En ambos casos `rho=1_primo-D <= 1_primo` y `sum D=M`; solo la norma L2
cambia por un factor `M`. El script `numerico/78_gap_L2_rho.py` verifica la
construccion exactamente. Esto no afirma que el `rho_B` publicado sea
concentrado: muestra que las propiedades citadas en los docs no excluyen esa
posibilidad y, por tanto, no prueban (78.1).

## 3. Lo que sobrevive y la reparacion minima

Parseval + Vinogradov si da, para cualquier cota

```text
sum rho(n)^2 << X (log X)^K,                                (78.2)
```

la estimacion

```text
sum_N |E_U(N)|^2 <<_A X^3/(log X)^(2A-K).
```

Como el ahorro de Vinogradov puede pedirse con `A` arbitrariamente grande,
(78.2) con algun `K` fijo bastaria para recuperar (74.4). Pero hace falta
probar (78.2) para la construccion real, por ejemplo desde una cota puntual
divisor-bounded de `S^triv`, `S^G` y `rho_B`, o citar un lema L2 del paper.
Hasta entonces el veredicto correcto es:

- **[T]** (74.4) vale para el proxy finito `rho_tilde` usado en el script.
- **[P]** (74.4) para el minorante verdadero `rho` es muy plausible y se
  reduce al sublema explicito `L2-rho`: (78.2).
- **[FALSO como prueba]** la medicion `1.467*pi(2X)` del proxy no demuestra
  la misma norma para `rho`.

La medicion FFT del doc 77 sigue siendo evidencia numerica del proxy, no del
minorante completo. El puente de casi todo a todo sigue abierto aun si
`L2-rho` se cierra.

## 4. Siguiente paso verificable

Leer la definicion completa de `rho_B` en el Lema 6.8 de Maynard III y
derivar una cota punto a punto `rho_B(n) << tau_k(n)` con `k` fijo, o una
cota media cuadratica directa. Solo despues se puede restaurar la etiqueta
`[T]` de la Prop. 77.1 para el minorante verdadero.

## Alcance

Esta revision no refuta el teorema de casi-todo Goldbach ni el argumento
clasico para la suma prima. Identifica una transferencia proxy--minorante no
probada en el doc 77. Goldbach sigue abierta.

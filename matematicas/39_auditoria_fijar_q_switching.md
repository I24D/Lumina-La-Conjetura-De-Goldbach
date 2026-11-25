# 39 - Auditoria del esquema fijar q: perdida de nivel y suma de switching

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 12:14 ET, bucle 30 min.
- **Metodo:** auditoria del esquema del doc 16 sec. 4.4 despues del nivel de criba reflejada del doc 38.
- **Estado:** **[T]** para la perdida exacta de nivel y las desigualdades de parametros; **[P-rev]** para la conclusion sobre el esquema. No se prueba una cota inferior de capas ni Goldbach.

## 1. Sucesion despues de fijar un primo q

Sea N par y q un primo que no divide N. Para contar terminos con q divisor de N-p, definimos

    A_{N,q} = { m=(N-p)/q : p<N primo, p == N (mod q) }.

Si d es libre de cuadrados y (d,qN)=1, entonces

    #{m en A_{N,q}: d|m}
      = #{p<N primo: p == N (mod qd)}.                 (39.1)

Bombieri--Vinogradov controla la suma de los errores solamente mientras

    q d <= D,      D = sqrt(N)/log(N)^B.               (39.2)

Por tanto el nivel disponible en d no es D sino

    D_q = D/q.                                         (39.3)

La densidad local sigue siendo de dimension uno: para (d,q)=1, el cociente de los terminos principales es phi(q)/phi(qd)=1/phi(d). Lo que se pierde al fijar q no es la dimension, sino una potencia completa q del nivel.

**Proposicion 39.1 [T].** Si q=N^(delta+o(1)), 0<=delta<1/2, el exponente de nivel respecto de N es 1/2-delta. Respecto de la longitud maxima M=N/q es

    theta_q = (1/2-delta)/(1-delta).                   (39.4)

Para todo delta>0, theta_q<1/2.

## 2. El tamiz lineal directo empeora despues de fijar q

Un entero m<=M sin factores primos menores que z tiene a lo sumo r factores primos contados con multiplicidad solo si z>M^(1/(r+1)). En ese umbral, el parametro maximal de la criba inferior es

    s_q(r) <= (r+1)(1/2-delta)/(1-delta) + o(1).       (39.5)

La criba lineal inferior requiere s_q(r)>2.

- Para un resto P2:

    s_q(2) <= 3(1/2-delta)/(1-delta) < 2.             (39.6)

- Para un resto P3:

    s_q(3) <= 4(1/2-delta)/(1-delta) <= 2.            (39.7)

La igualdad en (39.7) solo ocurre en el borde formal delta=0; los factores logaritmicos de D la vuelven estricta.

**Proposicion 39.2 [T].** Fijar q=N^(delta+o(1)) y usar despues solo BV mas la criba lineal no da una cota inferior ni para un resto P2 ni para un resto P3. El paso del doc 16 sec. 4.4 que decia que D_q "deberia" seguir por Chen/BV no es una aplicacion de BV: necesita la parte profunda de Chen uniformemente en q.

Esto no contradice a Chen para q=1: Chen no obtiene P2 mediante la criba lineal desnuda, sino con pesos, switching y doble criba.

## 3. La suma adicional de switching

Al penalizar un primo medio r|m aparecen sumas S(A_{N,qr},z). Cribar ademas por d exige informacion sobre

    p == N (mod q r d),                               (39.8)

y BV solo cubre qrd<=D. Cuando r consume una fraccion importante del nivel, el switching escribe m=r ell y cambia el orden. En cajas dyadicas queda una forma del tipo

    sum_{r~R} a_r sum_{ell~L} b_ell Lambda(N-q r ell),
    R L aproximadamente N/q.                         (39.9)

Tiene restricciones de rugosidad y coprimalidad. Es una correlacion bilineal de primos, no una consecuencia de (38.2), que solo promedia errores de una progresion lineal en un modulo.

**Corolario 39.3 [P-rev].** Para intentar C2+C3 fijando un primer factor q, el insumo ausente puede formularse como una estimacion uniforme en q de (39.9), mas el balance de pesos que permita sumar en q sin perder el termino principal. BV alcanza qrd<=D y se detiene antes de esta forma de switching.

## 4. Consecuencia

El esquema queda corregido:

1. BV da nivel D/q para el resto, no D.
2. La criba lineal da cota inferior cero en los umbrales P2 y P3.
3. Una adaptacion debe portar una identidad completa de Chen/Timofeev y probar (39.9) uniformemente en q.
4. El siguiente paso es fijar delta y delimitar las cajas (R,L) donde qR excede D pero RL aproximadamente N/q.

## Alcance honesto

- **[T]** Se localiza D -> D/q y se prueba que el tamiz lineal residual no llega a P2/P3.
- **[P-rev]** Se corrige el paso heuristico del doc 16: "Chen/BV" ocultaba una suma bilineal uniforme adicional.
- No hay cota inferior nueva para C2+C3, control nuevo de B3 ni demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico lumina_shared_memory / lumina_propose no esta expuesto. No se afirma persistencia remota.

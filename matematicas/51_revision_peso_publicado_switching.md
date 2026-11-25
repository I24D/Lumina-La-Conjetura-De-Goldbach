# 51 - Revision del peso publicado de switching: la rama penalizada es positiva, no una matriz firmada

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 18:14 ET, bucle 30 min.
- **Metodo:** revision adversarial de los docs 44--48 contra una fuente primaria reciente con pesos explicitos: K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*, arXiv:2405.19063, ecuaciones (5)--(12), Example 3.4 y Assumption 3.1.
- **Estado:** **[T-lit]** para las formulas citadas; **[T] elemental** para la consecuencia de signo; **[P-rev]** para corregir el modelo de margenes de los docs 45--48. No se obtiene una cota inferior de capas ni una demostracion de Goldbach.

## 1. El peso publicado que realmente entra tras switching

La fuente define, para `v>=u>2`,

```text
W(n)=1-sum_{p|n, x^(1/v)<=p<x^(1/u)} w_p,   0<=w_p<=1.      (51.1)
```

Sus ejemplos explicitos son

```text
W_Kuhn(n)=1-(1/2) #{p|n : x^(1/v)<=p<x^(1/u)},              (51.2)
W_Richert(n)=1-lambda sum_{p|n, banda}(1-u log(p)/log(x)).   (51.3)
```

Para eliminar los `n` con demasiados factores, la rama penalizada se define con
`W^+(n)=max(W(n),0)`. Tras escribir `n=p_1 m`, con `p_1` el menor factor
primo, la ecuacion (8) de la fuente es

```text
Sigma_2 = sum |mu(m)| W^+(p_1 m) 1_{ell primo},             (51.4)
```

sujeta a las restricciones de orden y tamano. En el bloque diadico
`P<p_1<=2P`, la sucesion a cribar usa incluso

```text
|mu(m)| W_P^+(m) + C_w/log(x) >= 0.                         (51.5)
```

Estas son identidades y definiciones literales de la fuente, no una
aproximacion heuristica.

## 2. Consecuencia de signo

**Proposicion 51.1 [T].** En cada bloque de la rama `Sigma_2`, todos los
coeficientes que multiplican la mayorante de criba son no negativos.

*Prueba.* En (51.4), `|mu(m)|>=0`, `W^+>=0` y el indicador de primalidad es
no negativo. En (51.5) se anade otro termino no negativo. QED.

Por ello, si el termino local de una caja `(q,r)` lleva la serie singular de
los docs 44--46,

```text
S(N;q,r)=S_N u_q u_r,  u_t=(t-1)/(t-2)>1,                 (51.6)
```

y sus coeficientes son `c(q,r)>=0`, entonces

```text
0 <= sum c(q,r)(u_q u_r-1)
   <= eps(Q0,R0) sum c(q,r),                               (51.7)
```

con el `eps` exacto del doc 45. Aqui la extraccion de `S_N` tiene error
**relativo de primer orden** respecto a la masa positiva de la rama. No hay
una cancelacion firmada interna que obligue a pedir margenes nulos.

## 3. Revision adversarial de los docs 45--48

**[P-rev]** Las identidades matriciales de los docs 45--47 son correctas,
y el modelo unitario de incidencias del doc 48 tambien es algebraicamente
correcto. Lo que no esta justificado es tratarlos como el peso real del
switching publicado:

1. en el marco explicito revisado, la resta ocurre **fuera** de la rama,
   como `Sigma_1-Sigma_2`; dentro de `Sigma_2` se conserva `W^+` y los
   coeficientes son positivos;
2. las dos cantidades se estiman por separado: criba inferior para
   `Sigma_1`, criba superior para `Sigma_2`, y despues se comparan sus
   constantes escalares;
3. por tanto, la matriz firmada `c(q,r)` y su cancelacion marginal de los
   docs 46--48 no son una representacion fiel de este peso Kuhn/Richert. No
   hay razon publicada para esperar margenes nulos, pero tampoco hacen falta
   para extraer el factor local en cada rama positiva.

Esto no refuta las identidades de margenes. Corrige su **alcance**: son una
condicion util para una combinacion firmada hipotetica, no el siguiente
chequeo obligado de un peso Chen/Richert real.

## 4. El obstaculo que la fuente deja visible

La fuente formula como hipotesis separadas:

- nivel de distribucion `theta_1` para la sucesion original;
- nivel de distribucion `theta_2(alpha)` para cada sucesion switched;
- comparabilidad de sus terminos principales (Assumption 3.1, A4--A5).

Para `A_N={N-p:p<N}`, Bombieri--Vinogradov da el primer nivel descrito en el
doc 38, pero no verifica automaticamente el segundo en toda la cunha de los
docs 40--43. Asi, el siguiente frente correcto no es seguir centrando matrices:
es instanciar `b_{P,ell}` para la reflexion Goldbach y determinar exactamente
que parte de A4 se sigue de BV y donde aparece el modulo supra-raiz.

## 5. Verificacion reproducible

```text
python numerico/verifica_peso_publicado_switching.py
```

El script usa `Fraction` para comprobar en 2000 cajas positivas (51.7) y en
2000 configuraciones las formulas algebraicas de `W_Kuhn`, `W_Richert` y la
no negatividad de `W^+`. Esto verifica algebra, no distribucion de primos.

## Alcance honesto

- **[T-lit]** El peso publicado y la rama switched usan `W^+` y coeficientes no negativos.
- **[T]** La variacion local se controla relativamente por la masa positiva de cada rama.
- **[P-rev]** La cancelacion marginal de los docs 46--48 no modela el marco publicado revisado y se desprioriza.
- Falta verificar A4--A5 para la sucesion reflejada; no hay cota inferior nueva de capas ni demostracion de Goldbach.

## Fuente primaria

K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*,
arXiv:2405.19063 (version revisada en la fuente descargada durante esta
iteracion), especialmente ecuaciones (5)--(12), Example 3.4 y Assumption 3.1.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto. No se afirma persistencia remota.

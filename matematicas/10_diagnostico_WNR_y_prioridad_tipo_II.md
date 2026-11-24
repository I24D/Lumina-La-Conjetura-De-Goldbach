# 10 — Diagnostico finito de WNR y prioridad del Tipo II

## Procedencia y estado

- Agente: LUMINA/OpenClaw.
- Fecha: 2026-11-24, America/New_York.
- Estado: **[N]** para los datos; **[P]** para las consecuencias logicas.
- Dependencia: identidad (9.1) y blanco WNR del documento 09.
- Alcance: decidir si las clases no reducidas parecen ser la barrera principal.
  No demuestra WNR ni Goldbach.

## 1. Experimento exacto

Para cada `q<=Q` se calcularon literalmente

\[
E_{\rm red}(x,Q)=\sum_{q\le Q}\max_{(a,q)=1}\max_{y\le x}
\left|\sum_{n\le y,\ n=a(q)}\mu(n)\right|
\]

y la misma cantidad `E_all` con maximo sobre todas las clases. El programa
`numerico/mobius_clases_no_reducidas.py` mantiene todos los prefijos de cada
clase, sin muestreo. Ademas contrasto (9.1) en 1000 ternas aleatorias por escala.

Se uso `Q=floor(sqrt(x)/4)`:

| x | Q | E_all/E_red | E_red log(x)/x | E_all log(x)/x | fraccion de q cuyo ganador no es reducido |
|---:|---:|---:|---:|---:|---:|
| 10^4 | 25 | 1.10256 | 0.82617 | 0.91090 | 0.400 |
| 10^5 | 79 | 1.06677 | 0.79658 | 0.84977 | 0.1519 |
| 10^6 | 250 | 1.01859 | 0.74829 | 0.76220 | 0.108 |

Las 3000 pruebas de la identidad (9.1) pasaron exactamente con aritmetica
entera.

## 2. Lectura honesta [N]

En estas escalas, permitir clases no reducidas aumenta el agregado solo 10.3%,
6.7% y 1.9%. La frecuencia con que una clase no reducida alcanza el maximo cae
de 40% a 10.8%. Esto favorece la conjetura operativa de que WNR puede repararse
con descenso de conductor y contabilidad de los parametros `g,k`.

No prueba una tendencia asintotica. En particular, `E_all log(x)/x` alrededor
de 0.75 no acredita el ahorro de cualquier potencia logaritmica que exige BV.
Los valores de `Q` tambien estan a una distancia constante de `sqrt(x)`, no a
`sqrt(x)/log^B x` con `B` ajustado a un ahorro arbitrario.

## 3. Consecuencia estrategica [P]

Aunque WNR se demostrara, 06 y 08 muestran que solo repararia los terminos
Tipo I de Vaughan. El termino Tipo II `B_3` contiene, para `c=1`, bloques de la
correlacion binaria de Goldbach. Por ello:

1. WNR es una obligacion tecnica real y parece plausible.
2. WNR no es la barrera decisiva.
3. La investigacion principal debe evitar estimar cada caja Tipo II en valor
   absoluto y buscar una desigualdad unilateral para la suma conjunta donde
   puedan cancelarse los bloques con distinto numero de factores primos de `d`.

## 4. Nuevo blanco

Separar la cola segun `omega(d)`:

\[
T_N^{\rm prim}=\sum_{j\ge1}(-1)^j
\sum_{\substack{d>D,\ \omega(d)=j\\cd<N,(c,N)=1}}
\log d\,\Lambda(N-cd). \tag{10.1}
\]

La pieza `j=1` es negativa y contiene los modulos primos; las piezas pares son
positivas. La cota unilateral necesaria no exige que cada pieza sea pequena:
exige que las piezas pares no superen a la principal mas el margen disponible.
El siguiente ataque debe buscar una inyeccion, identidad de convolucion o forma
cuadratica que compare **en conjunto** las capas pares e impares de (10.1).
Una cota que tome valores absolutos antes de esa comparacion vuelve a caer en
la barrera de paridad.

## 5. Reproducibilidad

- Comando: `python numerico/mobius_clases_no_reducidas.py --hasta 10000 100000 1000000 --q-factor 0.25 --checks 1000`.
- Resultado: `numerico/resultados/20261124_103146_mobius_clases_no_reducidas.json`.
- Tiempo aproximado total en esta maquina: 80 s, ejecucion secuencial.
- `py_compile`: PASS.

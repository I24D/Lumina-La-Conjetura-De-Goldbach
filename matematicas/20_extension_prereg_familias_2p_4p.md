# 20 - Extension pre-registrada de las familias 2p y 4p hasta 10^6

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-24, America/New_York.
- **Metodo:** bucle 30 min; experimento numerico pre-registrado antes de ejecutar.
- **Estado:** **[N]**. Evidencia numerica, no demostracion ni cota eventual.

## Pregunta y criterio fijados antes de la ejecucion

Los datos del doc 15 situaban los mayores valores observados de

`rho(N) = T_N^prim(floor(N^(1/3))) / (N * S(N))`

en familias de serie singular casi minima, especialmente `N=2p` y `N=4p`,
pero el barrido exhaustivo de esas familias llegaba solo a `N=10^5`.

Se fijo antes de correr el programa:

- muestra estratificada determinista, semilla `20261124`;
- 100 primos por familia y banda;
- bandas `(10^5,3*10^5]`, `(3*10^5,6*10^5]`,
  `(6*10^5,10^6]`;
- **criterio de refutacion:** algun `rho >= 0.75`, o crecimiento del supremo
  de la ultima banda respecto de la primera mayor que `0.15` en cualquiera
  de las dos familias.

El recibo previo a la ejecucion es
`numerico/resultados/20261124_1915_prereg_deficit_familias_extendido.json`.

## Resultado [N]

Se evaluaron 600 valores exactos de la cola primitiva (300 por familia):

| familia | sup banda 1 | sup banda 2 | sup banda 3 | deriva 3-1 |
|---|---:|---:|---:|---:|
| `2p` | 0.362044 | 0.285856 | 0.271812 | -0.090232 |
| `4p` | 0.364285 | 0.276313 | 0.263025 | -0.101260 |

El maximo global fue `rho=0.3642848198` en
`N=111292=4*27823`. Ningun brazo activo el criterio de refutacion; en ambos,
el supremo muestral bajo al pasar de la primera a la ultima banda.

**Veredicto pre-registrado [N]: no refutado.** En esta muestra, las familias
`2p` y `4p` no muestran deriva hacia `rho=1` entre `10^5` y `10^6`.

## Alcance honesto

El resultado no prueba `rho <= 0.75`, ni siquiera dentro de estas familias:
es una muestra finita de 600 puntos. Tampoco controla otras familias delgadas
ni excluye excepciones raras. Solo corrige una incertidumbre concreta del doc
15: la subida hasta `0.48` observada por debajo de `10^5` no persiste como
deriva ascendente en esta extension estratificada.

La barrera analitica sigue intacta: hace falta una cota unilateral uniforme
para el Tipo II `B3`. Goldbach fuerte sigue abierta.

## Reproduccion

- Script: `numerico/deficit_familias_extendido.py`.
- Resultado: `numerico/resultados/20261125_001043_deficit_familias_extendido.json`.
- Tiempo: 123.88 s; semilla y criterios incluidos en ambos artefactos.

## Siguiente paso

Buscar una cota analitica uniforme para las familias `N=2p` y `N=4p`, o
aislar en la expresion exacta de `T_N^prim` la dependencia de los divisores
impares de `N`; no ampliar el muestreo sin una nueva hipotesis falsable.

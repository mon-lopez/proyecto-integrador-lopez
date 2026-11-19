# Fichas técnicas de las asignaturas

Las cuatro fichas describen los experimentos y estudios que se ejecutaron efectivamente en el proyecto. Todas comparten el mismo objeto de estudio —las emisiones de CO₂ de Colombia y el papel de las fuentes no convencionales de energía renovable— y se implementan en los scripts indicados.

---

# Ficha Técnica 1: Diseño de Experimentos

**Responsable:** Mónica López

## Pregunta

¿Qué tecnología de generación (solar, eólica, biomasa, hidroeléctrica o fósil) y qué nivel de penetración (5 %, 15 % o 30 % de la generación anual) maximizan las emisiones de CO₂ evitadas en el sistema eléctrico colombiano, y existe interacción entre ambos factores?

## Respuesta

Se mide la respuesta "emisiones evitadas" (t CO₂/año) en un experimento factorial completo 5 × 3 con bloque anual, mediante un modelo de despacho con orden de mérito calibrado con generación UPME y factores de emisión de Ember. El análisis se hace con ANOVA de tipo III (`car::Anova`) y pruebas de robustez (log-ANOVA, Kruskal-Wallis, Welch).

## Factores y niveles

- **Factor A – Tecnología** (5 niveles): solar, eólica, biomasa, hidroeléctrica, fósil.
- **Factor B – Penetración** (3 niveles): 5 %, 15 %, 30 % de la generación anual.
- **Bloque** (16 niveles): año 2010-2025, que controla la hidrología y la demanda.

## Unidad experimental

Cada corrida del modelo de despacho: una combinación de año, tecnología y penetración.

## Tamaño y recursos

Diseño factorial completo 5 × 3 × 16 = 240 corridas, con orden de ejecución aleatorizado (semilla 2026). Tiempo de cómputo: menos de 2 minutos. Requiere R con `car`, `dplyr` y `ggplot2`.

## Resultados

Efecto dominante de la tecnología: F(4,210) = 797,9; p < 0,001; η²p = 0,938. Penetración: F(2,210) = 64,3; p < 0,001. Interacción: F(8,210) = 114,1; p < 0,001. La eólica evita 677,7 t CO₂/GWh, la hidro 666,1, la solar 642,3, la biomasa 473,6 y el fósil −708,0. Script: `scripts/04_doe.R`.

---

# Ficha Técnica 2: Estadística Espacial

**Responsable:** Mónica López

## Pregunta

¿Cómo influye la elección del método de interpolación espacial (kriging ordinario vs. distancia inversa ponderada, IDW) en la precisión de la predicción bajo distintas densidades de puntos georreferenciados?

## Respuesta

Se mide el error cuadrático medio de predicción (RMSE) en escala logarítmica mediante validación cruzada espacial por bloques, con los paquetes `gstat`, `sf` y `spdep` en R. Se contrasta la diferencia con la prueba de Diebold-Mariano con covarianza de largo plazo y con bootstrap pareado.

## Factores y niveles

- **Factor 1 – Método** (4 niveles): kriging ordinario, IDW con exponente 1, 2 y 3.
- **Factor 2 – Densidad de entrenamiento** (3 niveles): 25 %, 50 % y 100 % de las celdas.
- **Bloque**: estructura de autocorrelación espacial, controlada mediante cinco bloques espaciales (k-medias).

## Unidad experimental

La celda de 0,1° (9.657 celdas) con su predicción fuera de bloque.

## Tamaño y recursos

12 combinaciones método × densidad, cada una con validación cruzada de 5 bloques (y 5 réplicas de submuestreo anidado para las densidades parciales). Tiempo: ~6 minutos. Requiere R con `gstat`, `sf`, `spdep`.

## Resultados

Autocorrelación global I = 0,785 (p = 0,001); variograma exponencial con nugget 0,386, meseta 1,938 y rango práctico 167,6 km. Con el 100 % de los datos el kriging supera al mejor IDW (RMSE 1,418 vs 1,432; DM = −2,19; p = 0,028), pero con el 25 % y el 50 % no lo supera. Script: `scripts/03_spatial.R`.

---

# Ficha Técnica 3: Muestreo Estadístico

**Responsables:** Mónica López, Nicolás Forero

## Pregunta

¿Qué diseño de muestreo probabilístico (aleatorio simple, estratificado con asignación proporcional, estratificado con asignación de Neyman o espacialmente balanceado GRTS) proporciona menor varianza en la estimación del total de emisiones de CO₂ de una población finita de 9.657 celdas?

## Respuesta

Se mide la varianza del estimador del total y el coeficiente de variación (CV) mediante simulación Monte Carlo con 400 réplicas (50 para los diseños GRTS), comparando además el estimador de Horvitz-Thompson con un estimador de razón.

## Factores y niveles

- **Factor – Esquema de muestreo** (5 niveles): MAS, estratificado proporcional, estratificado Neyman, GRTS, GRTS-Neyman.
- **Bloque – Tamaño de muestra** (4 niveles): n = 100, 200, 400 y 800.

## Unidad experimental

La muestra probabilística seleccionada del marco de 9.657 celdas EDGAR de 0,1°.

## Tamaño y recursos

4.800 réplicas Monte Carlo para los diseños aleatorios y 100 para los de GRTS. Tiempo: ~10 minutos. Requiere R con `spsurvey` y `survey`.

## Resultados

El MAS tiene CV = 49,2 % con n = 400; la asignación de Neyman lo reduce a 19,8 % (eficiencia 616×). Para un CV del 5 % se requerirían 7.537 celdas; con un estrato de certeza del 2 % de las celdas el CV esperado baja a 1,9 % y la estimación es 97,45 Mt (valor censal 102,14 Mt). Scripts: `scripts/02_sampling.R` y `scripts/02b_diseno_final.R`.

---

# Ficha Técnica 4: Machine Learning

**Responsable:** Mónica López

## Pregunta

¿Qué variables explican la densidad de emisiones de CO₂ a escala de celda y qué trayectoria seguirán las emisiones nacionales hasta 2030 bajo distintos escenarios de penetración de FNCER?

## Respuesta

Se comparan Random Forest, XGBoost y regresión lineal en dos tareas: predicción de la densidad de emisiones por celda (con validación cruzada espacial por bloques) y pronóstico de la tasa de crecimiento nacional (con validación de origen móvil). Se reportan RMSE, MAE, R², importancia por permutación y probabilidades de reducción.

## Factores y niveles

- **Modelo** (3-5 niveles): Random Forest, XGBoost, regresión lineal, pronóstico ingenuo y combinación.
- **Esquema de validación** (2 niveles): aleatoria vs. espacial (celdas) u origen móvil (nacional).
- **Escenarios de política** (3 niveles): BAU (8 % FNCER), moderado (20 %), acelerado (35 %).

## Unidad experimental

En el modelo espacial, la celda de 0,1° (9.657 observaciones); en el pronóstico, el año (23 observaciones, 2002-2024).

## Tamaño y recursos

10 bloques espaciales × 3 modelos, 10 orígenes móviles × 5 modelos y 2.000 simulaciones por escenario. Tiempo: ~5 minutos con `ranger` y `xgboost`. Requiere R con `ranger`, `xgboost` y `spdep`.

## Resultados

La densidad de emisiones se explica por la densidad poblacional, la población y la longitud (R² = 0,41 con regresión, 0,38 con Random Forest). La validación espacial revela un optimismo del 24,7 % frente a la aleatoria. Ningún modelo supera al pronóstico ingenuo en la serie nacional (RMSE 8,74 p.p.), por lo que la trayectoria se ancla en el promedio histórico: CO₂ 2030 de 106,4 Mt (BAU), 98,3 Mt (moderado) y 90,1 Mt (acelerado). Script: `scripts/05_ml.R`.

---

El proyecto integra las cuatro asignaturas —Diseño de Experimentos, Muestreo, Estadística Espacial y Machine Learning— sobre una misma base de datos geoespacial y temporal, con el objetivo común de predecir si Colombia reducirá sus niveles de CO₂ para 2030 con el uso de fuentes no convencionales de energía renovable.

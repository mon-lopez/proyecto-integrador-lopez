# Hito 1: Ficha Técnica 1 (Diseño de Experimentos)
## Nombres 
Mónica López
## Pregunta
¿Cuál algoritmo de optimización de hiperparámetros (Búsqueda Aleatoria vs. Optimización Bayesiana) minimiza el error de validación cruzada en un modelo de aprendizaje automático manteniendo la eficiencia computacional?
## Respuesta 
Se medirá el error cuadrático medio (RMSE) en el conjunto de prueba usando análisis de varianza y los paquetes de experimentación `DoE.base` y `rsm` en R
## Factores & Niveles
Factor principal: Algoritmo de optimización (2 niveles: Búsqueda Aleatoria, Optimización Bayesiana). Bloque: Arquitectura de hardware (2 niveles: CPU de 4 núcleos, servidor de 16 núcleos
## Unidad Exprimental
Cada ejecución individual de entrenamiento del modelo de machine learning
## Tamaño
Se calculan 20 corridas experimentales (10 por bloque) de 3 minutos cada una (60 min en total).
## Que se necesita
Requiere un computador, R instalado (`DoE.base`, `rsm`) y sin otros programas pesados en ejecución

# Ficha Técnica 2: Estadística Espacial

## Pregunta
¿Cómo influye la elección del método de interpolación espacial (Kriging ordinario vs. Distancia Inversa Ponderada - IDW) en la precisión del modelo bajo distintas densidades de puntos georreferenciados?
## Respuesta
Se medirá el error cuadrático medio de predicción (RMSE) mediante validación cruzada espacial analizada con los paquetes `gstat` y `sf` en R.
## Factores y niveles
Factor 1: Método de interpolación (2 niveles: Kriging ordinario, IDW). Factor 2: Densidad de muestreo (2 niveles: 50 puntos/km², 200 puntos/km²). Bloque: Estructura de autocorrelación espacial del terreno (2 niveles: alta, baja).
## Unidad experimental
La malla o dataset de cobertura espacial georreferenciado.
## Tamaño y recursos
Se estiman 16 corridas de 2 minutos cada una (\~32 min en total). Requiere computador personal con R (`sf`, `spdep`, `gstat`).

### Ficha Técnica 3: Muestreo Estadístico

## Nombres 
Mónica López - Nicolas Forero
## Pregunta
¿Qué diseño de muestreo probabilístico (Aleatorio Simple vs. Estratificado con asignación de Neyman) proporciona menor varianza en la estimación del total de emisiones contaminantes en una población finita?
## Respuesta
Se medirá la varianza del estimador y el coeficiente de variación (CV) mediante estimación por replicación bootstrap usando el paquete `survey` en R.
## Factores y niveles
Factor: Esquema de muestreo probabilístico (2 niveles: Aleatorio Simple, Estratificado Óptimo). Bloque: Tamaño de muestra extraído (2 niveles: n = 100, n = 500)
## Unidad experimental
La muestra probabilística seleccionada de la población.
## Tamaño y recursos
Se calculan 20 corridas de simulación Monte Carlo (1.000 réplicas por corrida), con 1.5 minutos por corrida (\~30 min total).

El proyecto integra las 3 asignaturas inscritas Diseño de Experimentos, Muestreo y Estadística Espacial

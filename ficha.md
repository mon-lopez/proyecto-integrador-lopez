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

El proyecto integra las 3 asignaturas inscritas Diseño de Experimentos, Muestreo y Estadística Espacial

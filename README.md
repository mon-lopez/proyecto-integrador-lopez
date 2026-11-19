# ¿Puede Colombia reducir sus emisiones de CO₂ para 2030?

**Análisis integrado de muestreo espacial, estadística espacial, diseño de experimentos y aprendizaje automático en el marco de la Ley 2169 de 2021.**

Este repositorio contiene los datos procesados, los scripts, las figuras y el manuscrito de un proyecto que evalúa si las fuentes no convencionales de energía renovable (FNCER) pueden contribuir a la meta climática colombiana de 2030 y estima la trayectoria de las emisiones de CO₂ hasta ese año.

> **Respuesta breve.** Colombia sí puede reducir sus emisiones de CO₂ con FNCER, pero no puede cumplir la meta del 51 % solo con ellas. La matriz eléctrica ya es renovable en casi el 79 %, de modo que la nueva FNCER solo puede desplazar la fracción fósil (≈21 puntos porcentuales): el potencial máximo del sector eléctrico es de 14,4 Mt CO₂, el 5,8 % de los GEI totales. La meta de 169,44 Mt CO₂e exige reducir más de 90 Mt adicionales y depende de la deforestación, el transporte y la industria.

---

## 1. Pregunta y objetivos

**Pregunta.** ¿Puede Colombia reducir sus niveles de CO₂ para 2030 y qué papel juegan las energías no convencionales (solar, eólica, biomasa, hidro, fósil) en esa reducción?

**Objetivos específicos**

1. Estimar el total de emisiones y su precisión bajo distintos diseños de muestreo probabilístico.
2. Caracterizar la dependencia espacial de las emisiones y evaluar la interpolación geoestadística.
3. Cuantificar el efecto de cada tecnología y de su nivel de penetración mediante un experimento factorial.
4. Predecir la densidad de emisiones y la trayectoria nacional a 2030 con modelos de aprendizaje automático.
5. Contrastar el pronóstico con el techo legal de la Ley 2169 de 2021.

## 2. Las cuatro asignaturas y su papel

| Asignatura | Pregunta | Unidad de análisis | Método | Script |
|---|---|---|---|---|
| **Muestreo** | ¿Cuánto CO₂ hay y con qué precisión se estima? | Celda de 0,1° (N = 9.657) | MAS vs estratificado (proporcional y Neyman) vs GRTS; Monte Carlo; HT y razón; estrato de certeza | `scripts/02_sampling.R`, `02b_diseno_final.R` |
| **Estadística espacial** | ¿Dónde se concentra y cómo se autocorrelaciona? | Celda de 0,1° (EPSG:9377) | Moran global y local (LISA), semivariograma, kriging vs IDW, validación espacial | `scripts/03_spatial.R` |
| **Diseño de experimentos** | ¿Qué tecnología y cuánta penetración reducen más? | Corrida experimental (año × tecnología × penetración) | Factorial 5×3 con bloque anual, ANOVA tipo III, robustez, superficie de respuesta | `scripts/04_doe.R` |
| **Machine Learning** | ¿Qué explica la densidad y qué pasará en 2030? | Celda y año | Random Forest, XGBoost, regresión; validación espacial y de origen móvil; escenarios | `scripts/05_ml.R` |

## 3. Datos

| Fuente | Variable | Cobertura | Uso |
|---|---|---|---|
| EDGAR v2026 (JRC) | CO₂ por celda de 0,1° | 2015, 2023, 2025 | Marco muestral, estadística espacial, ML |
| Climate TRACE v5.7 | CO₂ por activo (18.979 registros) | 2023 | Covariables sectoriales, triangulación |
| GADM 4.1 | Departamentos y municipios | 2024 | Estratificación, regiones, unidad municipal |
| WorldPop 1 km | Población | 2015, 2020 | Covariable y estratificación auxiliar |
| MinEnergía | Proyectos FNCER (158) | 2004-2025 | Tratamiento municipal, capacidad |
| UPME | Generación por fuente | 2000-2026 | Experimento factorial, escenarios |
| Ember | Emisiones e intensidad por fuente | 2000-2025 | Factores de emisión, serie eléctrica |
| OWID / Global Carbon Budget | CO₂ y GEI nacionales | 1850-2024 | Serie histórica, pronóstico |
| Banco Mundial | PIB constante 2015 | 1960-2024 | Modelo de crecimiento |
| IDEAM / BTR2 | Inventario nacional de GEI | 1990-2024 | Contexto y contraste |

El detalle de la integración, las columnas de cada conjunto procesado y las convenciones de código están en [`docs/ESQUEMA_DATOS.md`](docs/ESQUEMA_DATOS.md). La documentación de las fuentes globales está en [`docs/FUENTES_CO2_ESPACIAL.md`](docs/FUENTES_CO2_ESPACIAL.md).

## 4. Resultados principales

| Indicador | Valor |
|---|---|
| Celda del marco muestral | 9.657 celdas de 0,1° (~11 km, 1,18 millones de km²) |
| CO₂ EDGAR 2023 en el marco | 102,14 Mt |
| Concentración espacial | El 1 % de las celdas concentra el 63,9 % del CO₂; Gini = 0,927 |
| Autocorrelación espacial | I de Moran = 0,785 (p = 0,001); rango práctico del variograma = 167,6 km |
| Clústeres LISA | 2.154 celdas HH (Andina y Caribe) y 3.122 LL (Amazonía y Orinoquía) |
| Diseño muestral | Neyman reduce el CV de 49,2 % a 19,8 %; con estrato de certeza baja a 1,9 % |
| Efecto de la tecnología (ANOVA) | F(4,210) = 797,9; p < 0,001; η²p = 0,938 |
| Efectividad por tecnología | Eólica 677,7 · Hidro 666,1 · Solar 642,3 · Biomasa 473,6 · Fósil −708,0 t CO₂/GWh |
| Límite estructural | 21,1 pp de FNCER útil → reducción máxima de 14,39 Mt CO₂ en el sector eléctrico |
| CO₂ 2030 (BAU / moderado / acelerado) | 106,4 / 98,3 / 90,1 Mt (IC 95 % del BAU: 70,1-165,9) |
| Probabilidad de reducir frente a 2024 | 0,27 / 0,39 / 0,54 |
| Brecha frente al techo legal (BAU) | 92,5 Mt CO₂e |

## 5. Estructura del repositorio

```text
co2-emisions-col/
├── README.md                     Este archivo
├── ficha.md                      Fichas técnicas de las asignaturas
├── preregistro.md                Protocolo de investigación
├── renv.lock                     Dependencias de R
├── datos/
│   ├── brutos/                   Datos descargados (WorldPop, WDI)
│   ├── geometrias/               GADM departamentos y municipios
│   ├── procesados/               Salidas de 01_preproc.R y de los módulos
│   └── *.csv / *.rds             Insumos originales del proyecto
├── scripts/
│   ├── 00_setup.R                Rutas, paquetes, paletas, funciones
│   ├── 01_preproc.R              Marco muestral, panel nacional y municipal
│   ├── 01b_eda.R                 Análisis exploratorio y figuras descriptivas
│   ├── 02_sampling.R             Comparación Monte Carlo de diseños
│   ├── 02b_diseno_final.R        Diseño con estrato de certeza
│   ├── 03_spatial.R              Moran, LISA, variograma, kriging
│   ├── 04_doe.R                  Factorial, ANOVA y superficie de respuesta
│   ├── 05_ml.R                   Modelos de celdas y pronóstico 2030
│   └── 06_resumen.R              Consolidación de resultados
├── manuscript/
│   ├── paper.qmd                 Manuscrito (Quarto)
│   ├── paper.html                Manuscrito renderizado (HTML autocontenido)
│   └── paper.docx                Manuscrito en Word
├── output/
│   ├── figures/                  Figuras (300 dpi)
│   ├── maps/                     Mapas
│   └── tables/                   Tablas de resultados (CSV)
├── docs/                         Documentación y referencias
├── fuentes/                      PDF y texto de las fuentes normativas (ley, NDC, BTR2, UPME)
├── energia-abierta-colombia/     Cosecha de datos abiertos (fuente)
├── emisiones_co2_espacial/       Rejillas EDGAR y cartografía base
└── _legacy/                      Versión anterior del proyecto (archivada)
```

## 6. Reproducibilidad

Requisitos: R ≥ 4.4 y los paquetes declarados en `renv.lock` (sf, terra, spdep, gstat, spsurvey, ranger, xgboost, ggplot2, car, quarto, entre otros). Quarto es necesario para renderizar el manuscrito.

```powershell
# Desde la raíz del proyecto
$env:RENV_CONFIG_SYNCHRONIZED_CHECK="FALSE"
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\01_preproc.R
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\01b_eda.R
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\02_sampling.R
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\02b_diseno_final.R
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\03_spatial.R
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\04_doe.R
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\05_ml.R
& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\06_resumen.R

# Manuscrito
& "C:\Users\monic\AppData\Local\Programs\Quarto\bin\quarto.exe" render manuscript\paper.qmd --to html
```

La semilla global es 2026. Todos los resultados del manuscrito se regeneran ejecutando los scripts en ese orden; el manuscrito no ejecuta código, lee las figuras y los objetos `.rds` generados.

## 7. Nota sobre la versión anterior

La carpeta `_legacy/` conserva los scripts, datos procesados y figuras de la versión previa del proyecto. Esa versión generaba la densidad de emisiones con un modelo lognormal simulado (`rlnorm`) en lugar de leer la rejilla EDGAR, lo que hacía que el índice de Moran fuera indefinido (NA) y que los resultados del paper anterior no fueran reproducibles. Esta versión reemplaza por completo esa capa de datos: todas las cifras provienen de EDGAR, Climate TRACE, UPME, Ember, OWID, el Banco Mundial y GADM.

## 8. Cómo citar

López, M. (2026). *¿Puede Colombia reducir sus emisiones de CO₂ para 2030? Evidencia integrada de muestreo espacial, estadística espacial, diseño de experimentos y aprendizaje automático en el marco de la Ley 2169 de 2021*. Universidad El Bosque, Facultad de Ciencias.

## 9. Licencia

Código: MIT. Manuscrito y figuras: CC BY 4.0.

🌎 Análisis Espacial de Emisiones de CO₂ y Evaluación Cuasi-Experimental de Fuentes Energéticas Renovables en Colombia
<div align="center">
<img src="https://img.shields.io/badge/estado-Fase_1_Fundacional-orange" alt="Estado">
<img src="https://img.shields.io/badge/licencia-MIT-blue" alt="Licencia">
<img src="https://img.shields.io/badge/manuscrito-CC%20BY%204.0-green" alt="Manuscrito">
<img src="https://img.shields.io/badge/R-≥%204.4.0-276DC3?logo=r" alt="R">
<img src="https://img.shields.io/badge/Python-≥%203.8-3776AB?logo=python" alt="Python">
</div>

## 📖 Resumen
Este repositorio contiene el código, los datos procesados y el manuscrito de un proyecto de investigación que combina estadística espacial, teoría del muestreo, diseño cuasi-experimental y Machine Learning para cuantificar y evaluar el impacto de las fuentes de energía renovable sobre las emisiones de CO₂ en Colombia.

**Nota sobre el flujo de trabajo (Fase 1):** Debido a limitaciones de memoria RAM en el equipo de desarrollo local, el proyecto se ha optimizado en una **Fase Fundacional**. Esta fase prioriza un Análisis Exploratorio de Datos (EDA) espacial riguroso, un diseño muestral balanceado, modelos predictivos base (Random Forest) y emparejamiento por puntuación de propensión (PSM). Los modelos de inferencia causal avanzada y geoestadística de alta complejidad computacional se han reservado como **Trabajo Futuro (Fase 2)** para ejecutarse posteriormente en entornos de mayor capacidad o mediante optimización de código.

El proyecto está diseñado para ser totalmente reproducible mediante `renv`, Quarto y un flujo de trabajo integrado con Obsidian.

📌 Preregistro oficial: [PREREGISTRO.md](./PREREGISTRO.md)
📌 Documentación técnica de fuentes globales: [FUENTES_CO2_ESPACIAL.md](./FUENTES_CO2_ESPACIAL.md)
📌 Referencias académicas curadas: [docs/referencias_academicas.md](./docs/referencias_academicas.md)
📌 Diario de investigación (Obsidian): carpeta `obsidian-vault/`

---

## 🎯 Preguntas de Investigación
| Código | Pregunta |
|---|---|
| **RQ1** | ¿Cuál es el efecto causal heterogéneo (CATE) de la implementación de fuentes renovables sobre la reducción de CO₂eq en Colombia, condicionado a características geoespaciales y socioeconómicas? |
| **RQ2** | ¿Existe autocorrelación espacial significativa en las emisiones de CO₂ a nivel municipal y de celda de 0.1°? |
| **RQ3** | ¿Qué fuente renovable presenta mayor efectividad marginal en reducción de CO₂ por región natural? |
| **RQ4** | ¿Cómo varía el efecto del tratamiento según el nivel de desarrollo energético base del municipio (estrato UPME)? |

---

## 📚 Fichas Temáticas por Asignatura (Fase 1: Enfoque Simplificado)
El proyecto articula cuatro líneas metodológicas. En esta fase, cada asignatura se ejecuta en **un único script** optimizado para bajo consumo de memoria, enfocándose en los fundamentos estadísticos y el EDA.

### 🗺️ 1. Estadística Espacial
| | |
|---|---|
| **Objetivo general** | Realizar un Análisis Exploratorio de Datos (EDA) espacial para caracterizar la distribución y autocorrelación global/local de las emisiones de CO₂ en Colombia. |
| **Preguntas que responde** | RQ2 |
| **Métodos clave** | • Cálculo del Índice de Moran global y local (LISA) básico con `spdep` <br> • Elaboración de mapas temáticos (coropletas y de puntos) <br> • Análisis de distribución empírica de las emisiones |
| **Entregables** | • Mapas de clusters LISA (hot/cold spots) <br> • Reporte de estadísticos de autocorrelación espacial |
| **Script** | `scripts/03_spatial_stats.R` |

### 🎲 2. Teoría del Muestreo
| | |
|---|---|
| **Objetivo general** | Diseñar un esquema de muestreo espacialmente balanceado que garantice una cobertura representativa del territorio para la recolección/estimación de datos. |
| **Preguntas que responde** | RQ1 (componente de diseño muestral) |
| **Métodos clave** | • Muestreo GRTS (*Generalized Random Tessellation Stratified*) con `spsurvey` <br> • Estratificación básica por regiones naturales <br> • Visualización del balance espacial de la muestra |
| **Entregables** | • Marco muestral espacial en `.gpkg` <br> • Mapa de distribución de la muestra vs. territorio nacional |
| **Script** | `scripts/02_spatial_sampling.R` |

### 🤖 3. Machine Learning
| | |
|---|---|
| **Objetivo general** | Establecer un modelo predictivo base (benchmark) que capture relaciones no lineales entre las covariables socioeconómicas/climáticas y las emisiones, identificando las variables de mayor importancia. |
| **Preguntas que responde** | RQ3, RQ4 (análisis de importancia de variables) |
| **Métodos clave** | • **Random Forest** (con `randomForest` o `ranger`) y **Árboles de Decisión** (con `rpart`) para modelado predictivo ligero. <br> • Evaluación de importancia de variables (Variable Importance) <br> • Validación cruzada estándar (k-fold) para estimar error de predicción (RMSE, R²) |
| **Entregables** | • Reporte de performance del modelo base <br> • Gráficos de importancia de variables (drivers de emisiones) |
| **Script** | `scripts/05_ml_models.R` |

### 🧪 4. Diseño de Experimentos
| | |
|---|---|
| **Objetivo general** | Estimar el efecto promedio del tratamiento (adopción de FNCER) sobre las emisiones utilizando técnicas de emparejamiento para mitigar sesgos de selección en datos observacionales. |
| **Preguntas que responde** | RQ1, RQ3 |
| **Métodos clave** | • Diseño cuasi-experimental con emparejamiento por *Propensity Score* (PSM) usando `MatchIt` <br> • Tratamientos binarios simplificados (T0: Fósil/Tradicional vs T1: Renovables) <br> • Validación de balance de covariables pre/post-emparejamiento |
| **Entregables** | • *Love plots* mostrando la reducción del sesgo estandarizado <br> • Estimación del Efecto Promedio del Tratamiento (ATE) y del Tratamiento sobre los Tratados (ATT) |
| **Script** | `scripts/04_matching.R` |

---

## 🗂️ Fuentes de Datos
*(Las fuentes de datos se mantienen intactas para garantizar la integridad científica del proyecto. La limitación computacional se aborda en el procesamiento y modelado, no en la calidad de los datos de entrada).*

### 3.1 Fuentes Globales (Capa Geoespacial Primaria)
| Capa | Fuente | Resolución | Licencia | Uso principal |
|---|---|---|---|---|
| Emisiones (rejilla primaria) | EDGAR v2026 (JRC) | 0.1° (~11 km), 1970–2025 | CC BY 4.0 | Estadística espacial |
| Emisiones (puntos) | Climate TRACE v5.7 | Activo individual, 2015–2025 | CC BY 4.0 | Tratamientos, puntos |
| Covariables país–año | OWID / Global Carbon Budget | País, 1750–2024 | CC BY 4.0 | ML, covariables |
| Covariables socioeconómicas | Banco Mundial (WDI) | País, 1960–2024 | CC BY 4.0 | Emparejamiento |
| Alta resolución mensual | ODIAC 2020b | 1 km, 2000–2021 | Académica | Kriging fino |
| Fronteras subnacionales | GADM 4.1 | ADM0–ADM2 | No comercial | Análisis municipal |
| Fronteras globales | Natural Earth 10m | País | Dominio público | Mapas |

### 3.2 Fuentes Nacionales — Catálogo `datos.gov.co` (FNCER y GEI)
Cosecha reproducible del catálogo Socrata. El catálogo curado incluye 40+ conjuntos de datos organizados en 4 ejes temáticos:

**🔌 Eje 1: FNCER — Capacidad y Proyectos**
* `vy9n-w6hc`: Meta FNCER (incluye Emisiones CO₂ [Ton/año]) - MinEnergía
* `nyyd-gdfx`: FNCER (general) - MinEnergía
* `tq63-itzr`: Producción de electricidad por fuente de combustible - UPME
* `ghd6-jnjh`: Conjuntos de datos obligatorios - UPME
* `chv4-3d2d`: Meteorología local en plantas solares

**⚡ Eje 2: Operación del SIN y ZNI**
* `3ebi-d83g`: Estado prestación del servicio en ZNI - IPSE
* `3a44-zwt6`: Registro de generación diaria ZNI - Superservicios
* `p7wb-z8zs`: Operación diaria otro tipo de generación - Superservicios
* `ptj2-kj7i`: Costos de racionamiento de energía - UPME
* `aqdb-bmuv`: Balance Energético Colombiano

**🌫️ Eje 3: Inventarios de GEI y Factores de Emisión**
* `6rff-a5ep`: Series históricas INGEI 1990–2024 - IDEAM/MinAmbiente
* `28sd-dsje`: Factores de emisión del inventario nacional - MinAmbiente
* `pvsp-rfj6`: Datos de actividad del inventario nacional - MinAmbiente
* `sx6m-p6bv`: Series históricas inventario subnacional - MinAmbiente
* `s3gg-fcvi`: Emisiones brutas GEI por departamento y municipio - MinAmbiente
* `yrhn-3pw5`: Emisiones netas GEI por departamento y municipio - MinAmbiente
* `qmgk-t9vt`: Emisiones GEI y carbono negro por cartera sectorial - MinAmbiente
* `xs9n-2dua`: Series históricas contaminantes criterio y carbono negro - MinAmbiente
* `tzxy-chv2`: Emisiones per cápita - MinAmbiente
* `isve-zbhc`: Cambio porcentual de emisiones por subcategoría - MinAmbiente

**🌦️ Eje 4: Recurso Solar, Eólico e Hídrico (IDEAM)**
* `rv9s-8nv6`: Radiación global acumulada diaria
* `sgfv-3yp8`: Velocidad del viento
* `kiw7-v9ta`: Dirección del viento
* `s54a-sgyg`: Precipitación
* `bdmn-sqnh`: Nivel instantáneo del río
* `vfth-yucv`: Nivel máximo del río
* `5kjg-nuda`: Zonificación hidrográfica
* `hp9r-jxuu`: Catálogo nacional de estaciones IDEAM
* `57sv-p2fu`: Datos de estaciones IDEAM y de terceros
* `nsz2-kzcq`: Normales climatológicas de Colombia

### 3.3 Referencias Académicas Curadas
Bibliografía de 60+ artículos validados contra Crossref.
* `docs/referencias_academicas.md` — lista APA 7 anotada
* `docs/referencias.bib` — BibTeX importable a Zotero/LaTeX
* `datos/referencias_verificadas.json` — metadatos crudos + licencia

---

## 📁 Estructura del Repositorio (Fase 1)
```text
co2-colombia-spatial/
 ├── README.md                       # Este archivo
 ├── PREREGISTRO.md                  # Protocolo preregistrado (OSF)
 ├── FUENTES_CO2_ESPACIAL.md         # Documentación técnica de fuentes globales
 ├── LICENSE                         # MIT (código)
 ├── .gitignore
 ├── renv.lock                       # Dependencias de R congeladas
 │
 ├── datos/                          # Datos nacionales cosechados
 │   ├── raw/                        # Datos globales descargados (NO subir a GitHub)
 │   ├── processed/                  # Datos limpios (.rds, .gpkg)
 │   ├── catalogo_datos_gov.csv      # Catálogo curado de datos.gov.co
 │   └── referencias_verificadas.json# Metadatos Crossref crudos
 │
 ├── scripts/                        # 🛠️ Flujo de trabajo optimizado (Fase 1)
 │   ├── 00_setup.R                  # Configuración inicial, renv, rutas
 │   ├── 01_preprocessing.R          # Limpieza, uniones espaciales básicas
 │   ├── 01b_analisis.R              # EDA inicial y resumen de datos
 │   ├── 02_spatial_sampling.R       # 🎲 Teoría del Muestreo (GRTS)
 │   ├── 03_spatial_stats.R          # 🗺️ Estadística Espacial (Moran, LISA)
 │   ├── 04_matching.R               # 🧪 Diseño de Experimentos (PSM)
 │   └── 05_ml_models.R              # 🤖 Machine Learning (Random Forest / Árboles)
 │
 ├── output/                         # Resultados (tablas, mapas PNG/PDF)
 │   ├── tables/
 │   ├── maps/
 │   └── figures/
 │
 ├── docs/                           # Documentación académica
 │   ├── referencias_academicas.md   
 │   └── referencias.bib             
 │
 ├── manuscript/
 │   └── paper.qmd                   # Manuscrito en Quarto
 │
 └── obsidian-vault/                 # Diario de investigación
     ├── 01_Protocolo/
     ├── 02_Diario/
     └── 03_Literatura/

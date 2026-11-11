# 🌎 Análisis Espacial de Emisiones de CO₂ y Evaluación Cuasi-Experimental de Fuentes Energéticas Renovables en Colombia

<div align="center">

[![Estado del proyecto](https://img.shields.io/badge/estado-pre--registro-yellow)](./PREREGISTRO.md)
[![Licencia del código](https://img.shields.io/badge/licencia-MIT-blue)](./LICENSE)
[![Licencia del manuscrito](https://img.shields.io/badge/manuscrito-CC%20BY%204.0-green)](https://creativecommons.org/licenses/by/4.0/)
[![R](https://img.shields.io/badge/R-≥%204.4.0-276DC3?logo=r)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-pendiente%20Zenodo-1f77b4)](https://zenodo.org/)
[![OSF](https://img.shields.io/badge/OSF-Preregistro-4ead94)](https://osf.io/)

**Grupo de Investigación en Estadística Aplicada** · Noviembre 2026

</div>

---

## 📖 Resumen

Este repositorio contiene el código, los datos procesados, el manuscrito y el protocolo de preregistro de un proyecto de investigación que combina **estadística espacial**, **teoría del muestreo**, **diseño cuasi-experimental** e **inferencia causal con Machine Learning** para:

1. Cuantificar la distribución espacial de las emisiones de CO₂ en Colombia a múltiples resoluciones (0.1°, 1 km, municipal).
2. Evaluar el efecto causal heterogéneo (CATE) de diferentes fuentes de energía renovable —solar, eólica, pequeña hidroeléctrica y biomasa— sobre la reducción de emisiones, condicionado a las características geoespaciales y socioeconómicas de cada territorio.

El proyecto está diseñado para ser **totalmente reproducible** mediante `renv`, Quarto y un flujo de trabajo integrado con Obsidian.

> 📌 **Preregistro oficial:** [PREREGISTRO.md](./PREREGISTRO.md)
> 📌 **Diario de investigación (Obsidian):** carpeta `obsidian-vault/`

---

## 🎯 Preguntas de Investigación

| Código | Pregunta |
|--------|----------|
| **RQ1** | ¿Cuál es el efecto causal heterogéneo (CATE) de la implementación de fuentes renovables sobre la reducción de CO₂eq en Colombia, condicionado a características geoespaciales y socioeconómicas? |
| **RQ2** | ¿Existe autocorrelación espacial significativa en las emisiones de CO₂ a nivel municipal y de celda de 0.1°? |
| **RQ3** | ¿Qué fuente renovable presenta mayor efectividad marginal en reducción de CO₂ por región natural? |
| **RQ4** | ¿Cómo varía el efecto del tratamiento según el nivel de desarrollo energético base del municipio (estrato UPME)? |

---

## 📚 Fichas Temáticas por Asignatura

El proyecto articula cuatro líneas metodológicas, cada una con objetivos, métodos y entregables propios. Esta sección sirve como **mapa de navegación** para los miembros del grupo según su área de especialización.

---

### 🗺️ 1. Estadística Espacial

<table>
<tr><td><b>Objetivo general</b></td><td>Caracterizar la estructura de dependencia espacial de las emisiones de CO₂ en Colombia y modelar su variabilidad a múltiples escalas.</td></tr>
<tr><td><b>Preguntas que responde</b></td><td>RQ2, RQ3</td></tr>
<tr><td><b>Métodos clave</b></td><td>
• Índice de Moran global y local (LISA) con <code>spdep</code> / <code>sfdep</code><br>
• Semivariograma y modelos geoestadísticos con <code>gstat</code><br>
• Modelos SAR / SEM / SLX con <code>spatialreg</code><br>
• Geographically Weighted Regression (GWR) con <code>gwrr</code><br>
• Análisis de sensibilidad MAUP (0.1° vs 0.5°)
</td></tr>
<tr><td><b>Entregables</b></td><td>
• Mapas de clusters LISA (hot/cold spots)<br>
• Modelos de semivariograma ajustados por región<br>
• Superficie de intensidad de carbono a 1 km (kriging)
</td></tr>
<tr><td><b>Scripts</b></td><td><code>scripts/03_spatial_descriptive.R</code>, <code>scripts/03b_gwr.R</code></td></tr>
<tr><td><b>Lectura recomendada</b></td><td>Anselin (1995), <i>Local Indicators of Spatial Association</i>; Fotheringham et al. (2002), <i>Geographically Weighted Regression</i>.</td></tr>
</table>

---

### 🎲 2. Teoría del Muestreo

<table>
<tr><td><b>Objetivo general</b></td><td>Diseñar un esquema de muestreo espacialmente balanceado que garantice estimaciones insesgadas y eficientes sobre el territorio colombiano.</td></tr>
<tr><td><b>Preguntas que responde</b></td><td>RQ1 (componente de diseño muestral)</td></tr>
<tr><td><b>Métodos clave</b></td><td>
• Muestreo GRTS (<i>Generalized Random Tessellation Stratified</i>) con <code>spsurvey</code> / <code>BalancedSampling</code><br>
• Estratificación: 5 regiones naturales × 3 niveles UPME = 15 estratos<br>
• Cálculo de tamaño muestral con efecto de diseño (<code>Deff</code>) estimado vía variograma piloto<br>
• Estimación de totales con pesos inversos a la probabilidad de inclusión<br>
• Validación con bootstrap espacial
</td></tr>
<tr><td><b>Entregables</b></td><td>
• Marco muestral en <code>.gpkg</code> con ~1,500–2,000 celdas seleccionadas<br>
• Informe de balance espacial y precisión esperada<br>
• Curvas de varianza vs. tamaño muestral
</td></tr>
<tr><td><b>Scripts</b></td><td><code>scripts/03_spatial_sampling.R</code></td></tr>
<tr><td><b>Lectura recomendada</b></td><td>Stevens & Donagan (2010), <i>GRTS for spatial surveys</i>; Cochran (1977), <i>Sampling Techniques</i>.</td></tr>
</table>

---

### 🤖 3. Machine Learning

<table>
<tr><td><b>Objetivo general</b></td><td>Utilizar modelos predictivos y de inferencia causal para estimar efectos heterogéneos de tratamiento y corregir sesgos espaciales.</td></tr>
<tr><td><b>Preguntas que responde</b></td><td>RQ1 (componente de heterogeneidad), RQ3, RQ4</td></tr>
<tr><td><b>Métodos clave</b></td><td>
• <b>Causal Forests</b> (<code>grf</code>) para estimar CATE<br>
• <b>Random Forests</b> (<code>missForest</code>) para imputación de datos faltantes preservando estructura espacial<br>
• <b>Spatial Cross-Validation</b> con <code>blockCV</code> (bloques geográficos, no aleatorios)<br>
• Modelos de interpretación: SHAP values (<code>shapviz</code>), variable importance<br>
• Comparación con modelos lineales (benchmark)
</td></tr>
<tr><td><b>Entregables</b></td><td>
• Mapas de CATE por fuente renovable<br>
• Reporte de performance (RMSE, R², coverage de intervalos de confianza)<br>
• Análisis de heterogeneidad por región y estrato UPME
</td></tr>
<tr><td><b>Scripts</b></td><td><code>scripts/05_causal_ml.R</code>, <code>scripts/05b_spatial_cv.R</code></td></tr>
<tr><td><b>Lectura recomendada</b></td><td>Athey, Tibshirani & Wager (2019), <i>Generalized Random Forests</i>; Meyer et al. (2023), <i>blockCV 3.0</i>.</td></tr>
</table>

---

### 🧪 4. Diseño de Experimentos

<table>
<tr><td><b>Objetivo general</b></td><td>Estimar el efecto causal de los tratamientos energéticos sobre la reducción de CO₂ en un contexto observacional donde la aleatorización pura es inviable.</td></tr>
<tr><td><b>Preguntas que responde</b></td><td>RQ1, RQ3, RQ4</td></tr>
<tr><td><b>Métodos clave</b></td><td>
• <b>Diseño cuasi-experimental</b> con emparejamiento por <i>Propensity Score</i> (<code>MatchIt</code>)<br>
• Tratamientos: T0 (fósil), T1 (solar), T2 (eólica), T3 (pequeña hidro), T4 (biomasa)<br>
• Covariables de emparejamiento: PIB, radiación, viento, población, distancia a red, cobertura de suelo<br>
• Validación con <i>Love plots</i> y análisis de sensibilidad de Rosenbaum<br>
• Modelos mixtos (LMM) con <code>lme4</code> como modelo complementario<br>
• Diagrama Acíclico Dirigido (DAG) con <code>dagitty</code>
</td></tr>
<tr><td><b>Entregables</b></td><td>
• DAG causal documentado<br>
• Tablas de balance pre/post-emparejamiento<br>
• Estimaciones ATE, ATT y CATE por tratamiento<br>
• Análisis de sensibilidad ante violación de ignorabilidad
</td></tr>
<tr><td><b>Scripts</b></td><td><code>scripts/04_matching.R</code>, <code>scripts/04b_lmm.R</code></td></tr>
<tr><td><b>Lectura recomendada</b></td><td>Stuart (2010), <i>Matching Methods for Causal Inference</i>; Imbens & Rubin (2015), <i>Causal Inference</i>.</td></tr>
</table>

---

## 🗂️ Fuentes de Datos

Todas las fuentes han sido **verificadas empíricamente** (enlaces, APIs y tamaños probados el 28 de octubre de 2026). Consulta el documento técnico completo en [`FUENTES_CO2_ESPACIAL.md`](./FUENTES_CO2_ESPACIAL.md).

| Capa | Fuente | Resolución | Licencia |
|------|--------|------------|----------|
| Emisiones (rejilla primaria) | **EDGAR v2026** (JRC) | 0.1° (~11 km), 1970–2025 | CC BY 4.0 |
| Emisiones (puntos) | **Climate TRACE v5.7** | Activo individual, 2015–2025 | CC BY 4.0 |
| Covariables país–año | **OWID / Global Carbon Budget** | País, 1750–2024 | CC BY 4.0 |
| Covariables socioeconómicas | **Banco Mundial (WDI)** | País, 1960–2024 | CC BY 4.0 |
| Alta resolución mensual | **ODIAC 2020b** | 1 km, 2000–2021 | Académica |
| Fronteras subnacionales | **GADM 4.1** | ADM0–ADM2 | No comercial |
| Fronteras globales | **Natural Earth 10m** | País | Dominio público |

> ⚠️ **Advertencias críticas documentadas:** celdas atípicas marinas en EDGAR, trampa de agregados en OWID (filtrar por `iso_code`), y límite de precisión de Natural Earth 10m en fronteras costeras (ej. Leticia). Ver [PREREGISTRO.md §3](./PREREGISTRO.md#3-fuentes-de-datos).

---

## 📁 Estructura del Repositorio

```text
co2-colombia-spatial/
├── README.md                       # Este archivo
├── PREREGISTRO.md                  # Protocolo preregistrado (OSF)
├── FUENTES_CO2_ESPACIAL.md         # Documentación técnica de datos
├── LICENSE                         # MIT (código)
├── .gitignore
├── renv.lock                       # Dependencias de R congeladas
│
├── data/
│   ├── raw/                        # Datos descargados (NO subir a GitHub)
│   │   ├── edgar/
│   │   ├── climate_trace/
│   │   ├── owid/
│   │   └── gadm/
│   └── processed/                  # Datos limpios (.rds, .gpkg)
│
├── scripts/
│   ├── 00_setup.R                  # Configuración inicial, renv, rutas
│   ├── 01_data_acquisition.R       # Descarga automatizada
│   ├── 02_preprocessing.R          # Limpieza, uniones espaciales
│   ├── 03_spatial_sampling.R       # Diseño GRTS
│   ├── 03b_gwr.R                   # Regresión geográficamente ponderada
│   ├── 04_matching.R               # Propensity Score Matching
│   ├── 04b_lmm.R                   # Modelos mixtos
│   ├── 05_causal_ml.R              # Causal Forests (grf)
│   ├── 05b_spatial_cv.R            # Validación cruzada espacial
│   └── 06_visualization.R          # Mapas y gráficos finales
│
├── output/                         # Resultados (tablas, mapas PNG/PDF)
│   ├── tables/
│   ├── maps/
│   └── figures/
│
├── manuscript/
│   └── paper.qmd                   # Manuscrito en Quarto
│
├── obsidian-vault/                 # Diario de investigación
│   ├── 01_Protocolo/
│   ├── 02_Diario/
│   ├── 03_Literatura/
│   └── 04_Glosario/
│
└── docs/                           # Sitio web del proyecto (pkgdown/quarto)

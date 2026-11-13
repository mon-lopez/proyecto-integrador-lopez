# 🌎 Análisis Espacial de Emisiones de CO₂ y Evaluación Cuasi-Experimental de Fuentes Energéticas Renovables en Colombia

<div align="center">

[![Estado del proyecto](https://img.shields.io/badge/estado-pre--registro-yellow)](./PREREGISTRO.md)
[![Licencia del código](https://img.shields.io/badge/licencia-MIT-blue)](./LICENSE)
[![Licencia del manuscrito](https://img.shields.io/badge/manuscrito-CC%20BY%204.0-green)](https://creativecommons.org/licenses/by/4.0/)
[![R](https://img.shields.io/badge/R-≥%204.4.0-276DC3?logo=r)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-≥%203.8-3776AB?logo=python)](https://www.python.org/)
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
> 📌 **Documentación técnica de fuentes globales:** [FUENTES_CO2_ESPACIAL.md](./FUENTES_CO2_ESPACIAL.md)
> 📌 **Referencias académicas curadas:** [docs/referencias_academicas.md](./docs/referencias_academicas.md)
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
<tr><td><b>Fuentes de datos</b></td><td>
• <b>EDGAR v2026</b> (rejilla 0.1° global) — capa primaria de emisiones<br>
• <b>Climate TRACE v5.7</b> (puntos georreferenciados de activos)<br>
• <b>GADM 4.1</b> (límites municipales de Colombia, 1119 municipios)<br>
• <b>IDEAM vía datos.gov.co</b>: radiación global (<code>rv9s-8nv6</code>), velocidad del viento (<code>sgfv-3yp8</code>), precipitación (<code>s54a-sgyg</code>), niveles de ríos (<code>bdmn-sqnh</code>)
</td></tr>
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
<tr><td><b>Lectura recomendada</b></td><td>Anselin (1995); Fotheringham et al. (2002).</td></tr>
</table>

---

### 🎲 2. Teoría del Muestreo

<table>
<tr><td><b>Objetivo general</b></td><td>Diseñar un esquema de muestreo espacialmente balanceado que garantice estimaciones insesgadas y eficientes sobre el territorio colombiano.</td></tr>
<tr><td><b>Preguntas que responde</b></td><td>RQ1 (componente de diseño muestral)</td></tr>
<tr><td><b>Fuentes de datos</b></td><td>
• <b>UPME vía datos.gov.co</b>: producción por fuente de combustible (<code>tq63-itzr</code>), conjuntos obligatorios (<code>ghd6-jnjh</code>)<br>
• <b>Balance Energético Colombiano</b> (<code>aqdb-bmuv</code>)<br>
• <b>Series históricas INGEI</b> (<code>6rff-a5ep</code>, <code>sx6m-p6bv</code>) — para estratificar por intensidad de carbono departamental
</td></tr>
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
<tr><td><b>Lectura recomendada</b></td><td>Stevens & Donagan (2010); Cochran (1977).</td></tr>
</table>

---

### 🤖 3. Machine Learning

<table>
<tr><td><b>Objetivo general</b></td><td>Utilizar modelos predictivos y de inferencia causal para estimar efectos heterogéneos de tratamiento y corregir sesgos espaciales, aprovechando la alta dimensionalidad de las series históricas nacionales.</td></tr>
<tr><td><b>Preguntas que responde</b></td><td>RQ1 (componente de heterogeneidad), RQ3, RQ4</td></tr>
<tr><td><b>Fuentes de datos (predictores)</b></td><td>
<b>Datos nacionales (datos.gov.co, API Socrata):</b><br>
• <b>Series históricas INGEI 1990–2024</b> (<code>6rff-a5ep</code>) — variable respuesta longitudinal<br>
• <b>Emisiones brutas/netas GEI por departamento y municipio</b> (<code>s3gg-fcvi</code>, <code>yrhn-3pw5</code>)<br>
• <b>Factores de emisión del inventario nacional</b> (<code>28sd-dsje</code>)<br>
• <b>Datos de actividad del inventario</b> (<code>pvsp-rfj6</code>)<br>
• <b>Emisiones per cápita</b> (<code>tzxy-chv2</code>)<br>
• <b>Meta FNCER MinEnergía</b> (<code>vy9n-w6hc</code>) — incluye <i>Emisiones CO₂ [Ton/año]</i><br>
• <b>Costos de racionamiento de energía</b> (<code>ptj2-kj7i</code>)<br>
• <b>Meteorología local en plantas solares</b> (<code>chv4-3d2d</code>)<br>
<br>
<b>Datos globales:</b><br>
• <b>OWID / Global Carbon Budget</b> — 79 covariables país–año<br>
• <b>World Bank WDI</b> — PIB, población, indicadores energéticos<br>
• <b>Climate TRACE v5.7</b> — emisiones por activo individual
</td></tr>
<tr><td><b>Métodos clave</b></td><td>
• <b>Causal Forests</b> (<code>grf</code>) para estimar CATE<br>
• <b>Gradient Boosting</b> (<code>xgboost</code>, <code>lightgbm</code>) como benchmark no causal<br>
• <b>Random Forests</b> (<code>missForest</code>) para imputación de datos faltantes preservando estructura espacial<br>
• <b>Spatial Cross-Validation</b> con <code>blockCV</code> (bloques geográficos, no aleatorios)<br>
• Modelos de interpretación: SHAP values (<code>shapviz</code>), variable importance<br>
• <b>Double/Debiased ML</b> (<code>DoubleML</code>) para estimación robusta de ATE con alta dimensionalidad<br>
• Comparación con modelos lineales (benchmark)
</td></tr>
<tr><td><b>Entregables</b></td><td>
• Mapas de CATE por fuente renovable<br>
• Reporte de performance (RMSE, R², coverage de intervalos de confianza)<br>
• Análisis de heterogeneidad por región y estrato UPME<br>
• Matriz de importancia de variables (SHAP) para identificar drivers de mitigación
</td></tr>
<tr><td><b>Scripts</b></td><td><code>scripts/05_causal_ml.R</code>, <code>scripts/05b_spatial_cv.R</code>, <code>scripts/05c_double_ml.R</code></td></tr>
<tr><td><b>Lectura recomendada</b></td><td>Athey, Tibshirani & Wager (2019); Chernozhukov et al. (2018), <i>Double/Debiased ML</i>; Meyer et al. (2023), <i>blockCV 3.0</i>.</td></tr>
</table>

---

### 🧪 4. Diseño de Experimentos

<table>
<tr><td><b>Objetivo general</b></td><td>Estimar el efecto causal de los tratamientos energéticos sobre la reducción de CO₂ en un contexto observacional donde la aleatorización pura es inviable, aprovechando la variación exógena en la adopción de FNCER a nivel municipal.</td></tr>
<tr><td><b>Preguntas que responde</b></td><td>RQ1, RQ3, RQ4</td></tr>
<tr><td><b>Fuentes de datos (tratamiento y covariables)</b></td><td>
<b>Variables de tratamiento (adopción de FNCER):</b><br>
• <b>Meta FNCER MinEnergía</b> (<code>vy9n-w6hc</code>) — capacidad instalada y emisiones evitadas por proyecto<br>
• <b>FNCER general</b> (<code>nyyd-gdfx</code>)<br>
• <b>Producción de electricidad por fuente de combustible UPME</b> (<code>tq63-itzr</code>)<br>
• <b>Conjuntos de datos obligatorios UPME</b> (<code>ghd6-jnjh</code>)<br>
• <b>Registro de generación diaria ZNI</b> (<code>3a44-zwt6</code>)<br>
• <b>Operación diaria otro tipo de generación</b> (<code>p7wb-z8zs</code>)<br>
• <b>Estado prestación del servicio en ZNI</b> (<code>3ebi-d83g</code>)<br>
• <b>Climate TRACE v5.7</b> — activos individuales de generación renovable con coordenadas<br>
<br>
<b>Covariables de emparejamiento:</b><br>
• <b>IDEAM</b>: radiación (<code>rv9s-8nv6</code>), viento (<code>sgfv-3yp8</code>), precipitación (<code>s54a-sgyg</code>), zonificación hidrográfica (<code>5kjg-nuda</code>)<br>
• <b>INGEI</b>: factores de emisión (<code>28sd-dsje</code>), datos de actividad (<code>pvsp-rfj6</code>)<br>
• <b>Emisiones GEI y carbono negro por cartera sectorial</b> (<code>qmgk-t9vt</code>)<br>
• <b>OWID / WDI</b>: PIB, población, acceso eléctrico
</td></tr>
<tr><td><b>Métodos clave</b></td><td>
• <b>Diseño cuasi-experimental</b> con emparejamiento por <i>Propensity Score</i> (<code>MatchIt</code>)<br>
• Tratamientos: T0 (fósil), T1 (solar), T2 (eólica), T3 (pequeña hidro), T4 (biomasa)<br>
• <b>Diferencias-en-Diferencias (DiD)</b> para evaluar impacto pre/post-instalación de proyectos FNCER<br>
• <b>Regression Discontinuity Design (RDD)</b> en umbrales de política (ej. subsidios ZNI)<br>
• <b>Instrumental Variables (IV)</b> usando radiación solar / velocidad del viento como instrumentos exógenos para adopción de FNCER<br>
• Validación con <i>Love plots</i> y análisis de sensibilidad de Rosenbaum<br>
• Modelos mixtos (LMM) con <code>lme4</code> como modelo complementario<br>
• Diagrama Acíclico Dirigido (DAG) con <code>dagitty</code>
</td></tr>
<tr><td><b>Entregables</b></td><td>
• DAG causal documentado<br>
• Tablas de balance pre/post-emparejamiento<br>
• Estimaciones ATE, ATT y CATE por tratamiento<br>
• Estimaciones DiD y RDD con intervalos de confianza robustos<br>
• Análisis de sensibilidad ante violación de ignorabilidad
</td></tr>
<tr><td><b>Scripts</b></td><td><code>scripts/04_matching.R</code>, <code>scripts/04b_lmm.R</code>, <code>scripts/04c_did.R</code>, <code>scripts/04d_rdd.R</code></td></tr>
<tr><td><b>Lectura recomendada</b></td><td>Stuart (2010); Imbens & Rubin (2015); Callaway & Sant'Anna (2021), <i>DiD with heterogeneous treatment effects</i>.</td></tr>
</table>

---

## 🗂️ Fuentes de Datos

Todas las fuentes han sido **verificadas empíricamente** (enlaces, APIs y tamaños probados el 28 de octubre de 2026). Consulta los documentos técnicos completos en [`FUENTES_CO2_ESPACIAL.md`](./FUENTES_CO2_ESPACIAL.md) y el catálogo curado en [`datos/catalogo_datos_gov.csv`](./datos/catalogo_datos_gov.csv).

### 3.1 Fuentes Globales (Capa Geoespacial Primaria)

| Capa | Fuente | Resolución | Licencia | Uso principal |
|------|--------|------------|----------|---------------|
| Emisiones (rejilla primaria) | **EDGAR v2026** (JRC) | 0.1° (~11 km), 1970–2025 | CC BY 4.0 | Estadística espacial |
| Emisiones (puntos) | **Climate TRACE v5.7** | Activo individual, 2015–2025 | CC BY 4.0 | Tratamientos, puntos |
| Covariables país–año | **OWID / Global Carbon Budget** | País, 1750–2024 | CC BY 4.0 | ML, covariables |
| Covariables socioeconómicas | **Banco Mundial (WDI)** | País, 1960–2024 | CC BY 4.0 | Emparejamiento |
| Alta resolución mensual | **ODIAC 2020b** | 1 km, 2000–2021 | Académica | Kriging fino |
| Fronteras subnacionales | **GADM 4.1** | ADM0–ADM2 | No comercial | Análisis municipal |
| Fronteras globales | **Natural Earth 10m** | País | Dominio público | Mapas |

> ⚠️ **Advertencias críticas documentadas:** celdas atípicas marinas en EDGAR, trampa de agregados en OWID (filtrar por `iso_code`), y límite de precisión de Natural Earth 10m en fronteras costeras (ej. Leticia). Ver [PREREGISTRO.md §3](./PREREGISTRO.md#3-fuentes-de-datos).

### 3.2 Fuentes Nacionales — Catálogo `datos.gov.co` (FNCER y GEI)

Cosecha reproducible del catálogo Socrata del portal de datos abiertos de Colombia mediante el script [`scripts/01_catalogo_datos_gov.py`](./scripts/01_catalogo_datos_gov.py). El catálogo curado incluye **40+ conjuntos de datos** organizados en 4 ejes temáticos:

#### 🔌 Eje 1: FNCER — Capacidad y Proyectos
| ID Socrata | Conjunto | Entidad | Uso en el proyecto |
|------------|----------|---------|---------------------|
| `vy9n-w6hc` | **Meta FNCER** (incluye Emisiones CO₂ [Ton/año]) | MinEnergía | 🧪 Variable de tratamiento principal |
| `nyyd-gdfx` | FNCER (general) | MinEnergía | 🧪 Tratamiento |
| `tq63-itzr` | Producción de electricidad por fuente de combustible | UPME | 🧪🎲 Tratamiento + estratificación |
| `ghd6-jnjh` | Conjuntos de datos obligatorios | UPME | 🧪 Tratamiento |
| `chv4-3d2d` | Meteorología local en plantas solares | — | 🤖 Covariable ML |

#### ⚡ Eje 2: Operación del SIN y ZNI
| ID Socrata | Conjunto | Entidad | Uso en el proyecto |
|------------|----------|---------|---------------------|
| `3ebi-d83g` | Estado prestación del servicio en ZNI | IPSE | 🧪 Diseño cuasi-experimental (ZNI) |
| `3a44-zwt6` | Registro de generación diaria ZNI | Superservicios | 🧪 DiD pre/post |
| `p7wb-z8zs` | Operación diaria otro tipo de generación | Superservicios | 🧪 DiD |
| `ptj2-kj7i` | Costos de racionamiento de energía | UPME | 🤖 Covariable ML |
| `aqdb-bmuv` | Balance Energético Colombiano | — | 🎲 Estratificación |

#### 🌫️ Eje 3: Inventarios de GEI y Factores de Emisión
| ID Socrata | Conjunto | Entidad | Uso en el proyecto |
|------------|----------|---------|---------------------|
| `6rff-a5ep` | **Series históricas INGEI 1990–2024** | IDEAM/MinAmbiente | 🤖🧪 Variable respuesta longitudinal |
| `28sd-dsje` | Factores de emisión del inventario nacional | MinAmbiente | 🤖🧪 Covariable clave |
| `pvsp-rfj6` | Datos de actividad del inventario nacional | MinAmbiente | 🤖 Predictor ML |
| `sx6m-p6bv` | Series históricas inventario subnacional | MinAmbiente | 🤖🧪 Panel departamental |
| `s3gg-fcvi` | Emisiones brutas GEI por departamento y municipio | MinAmbiente | 🤖 Variable respuesta |
| `yrhn-3pw5` | Emisiones netas GEI por departamento y municipio | MinAmbiente | 🤖 Variable respuesta |
| `qmgk-t9vt` | Emisiones GEI y carbono negro por cartera sectorial | MinAmbiente | 🧪 Covariable |
| `xs9n-2dua` | Series históricas contaminantes criterio y carbono negro | MinAmbiente | 🤖 Control |
| `tzxy-chv2` | Emisiones per cápita | MinAmbiente | 🤖 Covariable |
| `isve-zbhc` | Cambio porcentual de emisiones por subcategoría | MinAmbiente | 🧪 Robustez |

#### 🌦️ Eje 4: Recurso Solar, Eólico e Hídrico (IDEAM)
| ID Socrata | Conjunto | Uso en el proyecto |
|------------|----------|---------------------|
| `rv9s-8nv6` | Radiación global acumulada diaria | 🧪🤖 Instrumento IV + predictor ML |
| `sgfv-3yp8` | Velocidad del viento | 🧪🤖 Instrumento IV + predictor ML |
| `kiw7-v9ta` | Dirección del viento | 🤖 Predictor ML |
| `s54a-sgyg` | Precipitación | 🧪🤖 Control hidrología |
| `bdmn-sqnh` | Nivel instantáneo del río | 🧪 Control hidro |
| `vfth-yucv` | Nivel máximo del río | 🧪 Control hidro |
| `5kjg-nuda` | Zonificación hidrográfica | 🧪 Estratificación |
| `hp9r-jxuu` | Catálogo nacional de estaciones IDEAM | Metadatos |
| `57sv-p2fu` | Datos de estaciones IDEAM y de terceros | 🤖 Series climáticas |
| `nsz2-kzcq` | Normales climatológicas de Colombia | 🤖 Baseline climático |

> 💡 **Acceso programático:** cada conjunto tiene endpoints CSV y JSON documentados en [`datos/catalogo_datos_gov.csv`](./datos/catalogo_datos_gov.csv). Ejemplo:
> ```
> https://www.datos.gov.co/resource/vy9n-w6hc.csv?$limit=50000
> ```

### 3.3 Referencias Académicas Curadas

Bibliografía de **60+ artículos** validados contra **Crossref** mediante el script [`scripts/05_referencias.py`](./scripts/05_referencias.py). Los metadatos (autores, año, volumen, páginas, DOI, licencia del editor) se extraen automáticamente — no hay transcripción manual.

| Eje | N° de referencias | Uso en el proyecto |
|-----|-------------------|---------------------|
| A. Política, mercado y difusión de FNCER | 6 | Justificación de variables de política |
| B. Solar fotovoltaica | 5 | Supuestos de factor de planta |
| C. Eólica terrestre y marina | 5 | Evaluación de recurso, complementariedad |
| D. Hidroelectricidad, clima y GEI | 5 | Vínculo hidrología–emisiones |
| E. Biomasa y bioenergía | 7 | Potencial y ACV |
| F. Mitigación de CO₂, escenarios y contabilidad | 9 | Modelos LEAP, trayectorias net-zero |
| G. Zonas No Interconectadas y acceso rural | 4 | Diseño cuasi-experimental ZNI |
| H. Análisis de ciclo de vida y huella de carbono | 3 | Factores de emisión de cadena completa |

> 📄 **Archivos generados:**
> - [`docs/referencias_academicas.md`](./docs/referencias_academicas.md) — lista APA 7 anotada
> - [`docs/referencias.bib`](./docs/referencias.bib) — BibTeX importable a Zotero/LaTeX
> - [`datos/referencias_verificadas.json`](./datos/referencias_verificadas.json) — metadatos crudos + licencia

---

## 📁 Estructura del Repositorio

```text
co2-colombia-spatial/
├── README.md                       # Este archivo
├── PREREGISTRO.md                  # Protocolo preregistrado (OSF)
├── FUENTES_CO2_ESPACIAL.md         # Documentación técnica de fuentes globales
├── LICENSE                         # MIT (código)
├── .gitignore
├── renv.lock                       # Dependencias de R congeladas
│
├── datos/                          # 🆕 Datos nacionales cosechados
│   ├── raw/                        # Datos globales descargados (NO subir a GitHub)
│   │   ├── edgar/
│   │   ├── climate_trace/
│   │   ├── owid/
│   │   └── gadm/
│   ├── processed/                  # Datos limpios (.rds, .gpkg)
│   ├── catalogo_datos_gov.json     # 🆕 Catálogo curado de datos.gov.co (JSON)
│   ├── catalogo_datos_gov.csv      # 🆕 Catálogo curado de datos.gov.co (CSV)
│   ├── log_descubrimiento_datos_gov.json  # 🆕 Log de búsqueda Socrata
│   └── referencias_verificadas.json       # 🆕 Metadatos Crossref crudos
│
├── scripts/
│   ├── 00_setup.R                  # Configuración inicial, renv, rutas
│   ├── 01_catalogo_datos_gov.py    # 🆕 Cosecha Socrata (Python)
│   ├── 01b_data_acquisition.R      # Descarga automatizada fuentes globales
│   ├── 02_preprocessing.R          # Limpieza, uniones espaciales
│   ├── 03_spatial_sampling.R       # Diseño GRTS
│   ├── 03b_gwr.R                   # Regresión geográficamente ponderada
│   ├── 04_matching.R               # Propensity Score Matching
│   ├── 04b_lmm.R                   # Modelos mixtos
│   ├── 04c_did.R                   # 🆕 Diferencias-en-Diferencias
│   ├── 04d_rdd.R                   # 🆕 Regression Discontinuity
│   ├── 05_causal_ml.R              # Causal Forests (grf)
│   ├── 05b_spatial_cv.R            # Validación cruzada espacial
│   ├── 05c_double_ml.R             # 🆕 Double/Debiased ML
│   ├── 05_referencias.py           # 🆕 Curaduría bibliográfica Crossref
│   └── 06_visualization.R          # Mapas y gráficos finales
│
├── output/                         # Resultados (tablas, mapas PNG/PDF)
│   ├── tables/
│   ├── maps/
│   └── figures/
│
├── docs/                           # 🆕 Documentación académica
│   ├── referencias_academicas.md   # Lista APA 7 anotada
│   ├── referencias.bib             # BibTeX curado
│   └── sitio_web/                  # pkgdown/quarto
│
├── manuscript/
│   └── paper.qmd                   # Manuscrito en Quarto
│
└── obsidian-vault/                 # Diario de investigación
    ├── 01_Protocolo/
    ├── 02_Diario/
    ├── 03_Literatura/
    └── 04_Glosario/

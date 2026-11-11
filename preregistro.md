# Preregistro del Proyecto: Análisis Espacial de Emisiones de CO₂ y Evaluación Cuasi-Experimental de Fuentes Energéticas Renovables en Colombia

**Fecha de preregistro:** 
**Estado:** 🟡 Pre-data analysis (análisis no iniciado)  
**Repositorio:**   
**DOI reservado:** (se generará en Zenodo tras la primera release)  
**Equipo:** Monica Lopez (Nicolas Forero - Oscar Joya/ Colaboradores) 
**Institución:** Universidad El Bosque - Faculta de Ciencias

---

## 1. Resumen Ejecutivo

Este proyecto tiene como objetivo cuantificar la distribución espacial de las emisiones de CO₂ en Colombia y evaluar el efecto causal heterogéneo de diferentes fuentes de energía renovable (solar, eólica, hidroeléctrica, biomasa) sobre la reducción de emisiones, utilizando un enfoque combinado de **estadística espacial**, **diseño cuasi-experimental** e **inferencia causal con Machine Learning**. El estudio se realiza sobre el territorio continental colombiano, con resolución subnacional (municipal y celdas de 0.1°), combinando datos de rejilla (EDGAR), puntos georreferenciados (Climate TRACE) y covariables socioeconómicas (OWID, Banco Mundial).

---

## 2. Preguntas de Investigación

### Pregunta Principal (RQ1)
¿Cuál es el efecto causal heterogéneo (CATE) de la implementación de fuentes de energía renovable sobre la reducción de emisiones de CO₂eq en Colombia, condicionado a las características geoespaciales y socioeconómicas de cada territorio?

### Preguntas Secundarias
- **RQ2:** ¿Existe autocorrelación espacial significativa en las emisiones de CO₂ a nivel municipal y de celda de 0.1° en Colombia? (Índice de Moran global y local LISA).
- **RQ3:** ¿Qué fuente renovable (solar, eólica, pequeña hidroeléctrica, biomasa) presenta mayor efectividad marginal en reducción de CO₂ por región natural (Andina, Caribe, Pacífica, Orinoquía, Amazonía)?
- **RQ4:** ¿Cómo varía el efecto del tratamiento según el nivel de desarrollo energético base del municipio (estrato UPME: alto, medio, bajo)?

### Hipótesis Pre-registradas
- **H1:** El efecto de la energía solar sobre la reducción de CO₂ será significativamente mayor en la Región Caribe que en la Andina (mayor irradiación, menor nubosidad).
- **H2:** La biomasa presentará el mayor efecto en regiones con alta actividad agrícola (Orinoquía, Caribe).
- **H3:** Los municipios con menor desarrollo energético base mostrarán mayores reducciones marginales de CO₂ por unidad de capacidad renovable instalada.

---

## 3. Fuentes de Datos

Todas las fuentes han sido verificadas empíricamente por el equipo antes del inicio del análisis. Se documentan las URLs exactas, licencias y advertencias críticas.

### 3.1 Capa Primaria de Emisiones (Rejilla)
| Fuente | Resolución | Cobertura | Licencia | URL |
|--------|-----------|-----------|----------|-----|
| **EDGAR v2026** (JRC) | 0.1° (~11 km) y 1 km | Global, 1970–2025 | CC BY 4.0 | [edgar.jrc.ec.europa.eu/dataset_ghg2026](https://edgar.jrc.ec.europa.eu/dataset_ghg2026) |
| **ODIAC 2020b** | 1 km, mensual | Global, 2000–2021 | Académica (registro) | [odiac.org](https://www.odiac.org/) |

**⚠️ Advertencias críticas de EDGAR (verificadas):**
- Celdas atípicas en el mar: ~500 celdas concentran ~3 Gt (8% del total mundial). La peor: `113.65°E, 21.95°N` con 430 Mt en una sola celda.
- Al agregar la malla por país se obtienen 34.89 Gt vs. 39.22 Gt oficiales (diferencia mediana 10.3%). **Se usará la tabla oficial para totales nacionales y la malla solo para patrones espaciales.**
- No incluye LULUCF (uso de suelo y silvicultura).
- Distribución fuertemente asimétrica: mediana 69 t/celda, p99.99 = 11.8 Mt. Se trabajará en escala log o con métodos robustos.

### 3.2 Capa de Puntos (Activos Individuales)
| Fuente | Resolución | Cobertura | Licencia | URL |
|--------|-----------|-----------|----------|-----|
| **Climate TRACE v5.7** | Activo individual | Global, 2015–2025 | CC BY 4.0 | [climatetrace.org](https://climatetrace.org/) / [Zenodo](https://doi.org/10.5281/zenodo.20830133) |

**Uso previsto:** Identificación precisa de plantas de generación, minas y activos energéticos para el diseño de tratamientos.

### 3.3 Covariables País–Año
| Fuente | Contenido | Licencia | URL |
|--------|-----------|----------|-----|
| **OWID / Global Carbon Budget** | 79 columnas (CO₂ per cápita, PIB, población, etc.) | CC BY 4.0 | [github.com/owid/co2-data](https://github.com/owid/co2-data) |
| **Banco Mundial (WDI)** | PIB, población, indicadores energéticos | CC BY 4.0 | [api.worldbank.org](https://api.worldbank.org/) |

**⚠️ Trampa identificada:** La columna `country` de OWID mezcla países con agregados (`World`, `Africa`). Se filtrará siempre por `iso_code` excluyendo códigos `OWID_*`.

### 3.4 Fronteras y Polígonos
| Fuente | Niveles | Licencia | Uso |
|--------|---------|----------|-----|
| **Natural Earth 10m** | País | Dominio público | Mapas y análisis continental |
| **GADM 4.1** | País, departamento, municipio | No comercial | Análisis subnacional de Colombia |
| **geoBoundaries** | ADM0–ADM4 | CC0 | Alternativa libre |

**⚠️ Problema de precisión:** Natural Earth 10m falla en fronteras costeras (ej. Leticia cae del lado brasileño). **Se usará GADM nivel 2 (`gadm41_COL_2.json`, 1119 municipios) para el análisis subnacional.**

---

## 4. Diseño de Muestreo Espacial

### 4.1 Marco Muestral
- **Unidad de muestreo:** Celda de malla regular de 5 km × 5 km superpuesta al territorio continental colombiano.
- **Número aproximado de celdas:** ~45,000 (población finita).
- **Estratos:** 5 regiones naturales × 3 niveles de desarrollo energético (UPME) = 15 estratos.

### 4.2 Método de Selección
**GRTS (Generalized Random Tessellation Stratified)** implementado con el paquete `spsurvey` o `BalancedSampling` de R.

**Justificación:** El GRTS garantiza balance espacial, permite estimación de varianzas sin sesgo y es robusto ante autocorrelación espacial, superior al muestreo aleatorio simple en contextos geoespaciales.

### 4.3 Tamaño de Muestra
Calculado mediante:

$$n = \frac{Z^2 \cdot S^2 \cdot Deff}{E^2}$$

Donde $Deff$ (efecto de diseño) se estimará a partir de un variograma piloto para ajustar la varianza por autocorrelación espacial (Moran's I).

**Muestra objetivo:** 1,500–2,000 celdas (margen de error 5%, confianza 95%).

### 4.4 Análisis de Sensibilidad (MAUP)
Se realizarán análisis de sensibilidad con al menos dos resoluciones espaciales:
- Celdas de 0.1° (~11 km)
- Celdas de 0.5° (~55 km)

---

## 5. Diseño Experimental (Cuasi-Experimental)

Dado que la asignación aleatoria de fuentes energéticas es inviable en la práctica, se utiliza un **diseño cuasi-experimental con emparejamiento por puntaje de propensión (Propensity Score Matching) espacial**.

### 5.1 Tratamientos
| Código | Tratamiento | Descripción |
|--------|-------------|-------------|
| T0 | Control | Matriz energética fósil predominante (gas/carbón) |
| T1 | Solar FV | Energía solar fotovoltaica |
| T2 | Eólica | Energía eólica |
| T3 | Pequeña Hidro | Hidroeléctrica de pequeña escala (<10 MW) |
| T4 | Biomasa/Biogás | Residuos agrícolas y urbanos |

### 5.2 Variable de Respuesta
$\Delta CO_2$: Toneladas métricas de CO₂ equivalente reducidas por GWh generado, o cambio neto en la intensidad de carbono regional.

### 5.3 Emparejamiento
Los municipios/celdas se emparejarán por **Propensity Score** basado en covariables pre-tratamiento:
- PIB municipal
- Radiación solar promedio anual
- Velocidad del viento promedio
- Población
- Distancia a la red eléctrica nacional
- Tipo de cobertura de suelo

**Validación del balance:** Love plots para verificar que las covariables estén balanceadas tras el emparejamiento (diferencias estandarizadas < 0.1).

### 5.4 Supuesto de Ignorabilidad Condicional
Se asume que, tras el emparejamiento, no existen variables de confusión no medidas. Esto se validará mediante:
- Análisis de sensibilidad de Rosenbaum
- Especificación de un DAG (Diagrama Acíclico Dirigido) con el paquete `dagitty`

---

## 6. Protocolo de Análisis Estadístico

### 6.1 Análisis Descriptivo Espacial
- **Índice de Moran global:** Prueba de autocorrelación espacial en emisiones de CO₂.
- **LISA (Local Indicators of Spatial Association):** Identificación de clusters (Hot spots / Cold spots) con `spdep` o `sfdep`.
- **Semivariograma:** Modelización de la estructura de correlación espacial con `gstat`.

### 6.2 Modelos Clásicos
- **Modelos Mixtos Lineales (LMM):** Con efectos aleatorios por región natural, implementados con `lme4`.
- **ANOVA de medidas repetidas:** Para comparar tratamientos a lo largo del tiempo.

### 6.3 Inferencia Causal con Machine Learning
- **Causal Forests (Bosques Causales):** Estimación del **Conditional Average Treatment Effect (CATE)** con el paquete `grf` (Athey, Tibshirani, Wager).
- **Objetivo:** Identificar heterogeneidad en el efecto del tratamiento según características geoespaciales.

### 6.4 Validación Cruzada Espacial
Se utilizará **Spatial Cross-Validation** con el paquete `blockCV` para evitar fuga de datos por autocorrelación espacial. Los bloques se definirán geográficamente, no aleatoriamente.

### 6.5 Manejo de Datos Faltantes
Imputación con **Random Forests** (`missForest`) preservando la estructura de varianza espacial.

### 6.6 Software y Paquetes
| Paquete | Función |
|---------|---------|
| `sf`, `terra` | Manejo de datos vectoriales y raster |
| `spsurvey`, `BalancedSampling` | Diseño GRTS |
| `gstat`, `spdep` | Geoestadística y autocorrelación |
| `MatchIt` | Emparejamiento por propensity score |
| `lme4` | Modelos mixtos |
| `grf` | Causal Forests |
| `blockCV` | Validación cruzada espacial |
| `dagitty` | Diagramas causales |
| `exactextractr` | Extracción exacta de raster por polígono |

---

## 7. Plan de Trabajo (Cronograma)

| Fase | Actividad | Duración | Entregable |
|------|-----------|----------|------------|
| **Fase 1** | Adquisición y limpieza de datos | Semanas 1–3 | Scripts `01_data_acquisition.R`, `02_preprocessing.R` |
| **Fase 2** | Construcción del marco muestral y diseño GRTS | Semanas 4–5 | Script `03_spatial_sampling.R`, marco muestral en `.gpkg` |
| **Fase 3** | Análisis descriptivo espacial (Moran, LISA, variograma) | Semanas 6–7 | Mapas y reportes en `output/` |
| **Fase 4** | Emparejamiento y diseño cuasi-experimental | Semanas 8–9 | Script `04_matching.R`, tablas de balance |
| **Fase 5** | Modelos clásicos (LMM, ANOVA) | Semanas 10–11 | Resultados en `output/modelos_clasicos/` |
| **Fase 6** | Inferencia causal con Causal Forests | Semanas 12–14 | Script `05_causal_ml.R`, mapas de CATE |
| **Fase 7** | Redacción del manuscrito | Semanas 15–18 | Paper en Quarto (`manuscript/paper.qmd`) |
| **Fase 8** | Revisión, publicación y release | Semanas 19–20 | DOI en Zenodo, artículo enviado |

---

## 8. Criterios de Exclusión y Manejo de Datos

### 8.1 Exclusión de Observaciones
- Celdas con emisión exactamente 0 (océanos, desiertos): se tratarán como **ceros estructurales**, no como datos faltantes. Se decidirá si se modelan con censura o se excluyen del análisis, según el método.
- Celdas atípicas en el mar (artefactos de EDGAR): se excluirán del análisis de agregación por país.
- Búnkeres marítimos y aviación internacional (filas `SEA`, `AIR` de EDGAR): no se asignarán a países.

### 8.2 Transformaciones
- Las emisiones se trabajarán en **escala logarítmica** debido a la asimetría extrema (mediana 69 t, p99.99 = 11.8 Mt).
- Se reportarán resultados en escala original mediante retransformación con corrección de sesgo (smearing estimator).

### 8.3 Manejo de Fronteras
- Para agregación exacta por área (repartiendo celdas de borde entre países), se usará `exactextractr` con weights, no asignación por centroide.

---

## 9. Reproducibilidad y Gestión de Código

### 9.1 Estructura del Repositorio
```text
co2-colombia-spatial/
├── PREREGISTRO.md          # Este documento
├── README.md               # Descripción y cómo reproducir
├── .gitignore
├── data/
│   ├── raw/                # Datos descargados (no subir a GitHub si son pesados)
│   └── processed/          # Datos limpios (.rds, .gpkg)
├── scripts/
│   ├── 01_data_acquisition.R
│   ├── 02_preprocessing.R
│   ├── 03_spatial_sampling.R
│   ├── 04_matching.R
│   ├── 05_causal_ml.R
│   └── 06_visualization.R
├── output/                 # Resultados (tablas, mapas)
├── manuscript/             # Paper en Quarto (.qmd)
└── renv.lock               # Dependencias de R congeladas

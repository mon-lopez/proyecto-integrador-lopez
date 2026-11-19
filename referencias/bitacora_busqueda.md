# Bitácora de búsqueda académica (generada automáticamente)

- **Fecha de ejecución:** 2026-11-11
- **Ventana temporal:** 2014-2026 (desde la Ley 1715 de 2014)
- **API:** OpenAlex (`api.openalex.org`), metadatos bajo CC0 1.0
- **Campo buscado:** título y resumen (`title_and_abstract.search`)
- **Orden:** citaciones descendentes; se exportan los 25 primeros por ecuación

| ID | Eje temático | Cadena de búsqueda exacta | Resultados | Exportados |
|---|---|---|---:|---:|
| Q1 | Núcleo: renovables + CO2 + Colombia | `Colombia AND ("renewable energy" OR solar OR wind OR biomass OR hydropower) AND ("CO2" OR "carbon emissions" OR "greenhouse gas" OR mitigation)` | 468 | 25 |
| Q2 | FNCER / política e incentivos (Ley 1715, Ley 2099) | `Colombia AND ("non-conventional" OR "energy transition" OR "energy policy" OR "1715" OR "decarbonization") AND (renewable OR electricity)` | 271 | 25 |
| Q3 | Solar fotovoltaica: potencial, desempeño y huella | `Colombia AND ("photovoltaic" OR "solar PV" OR "solar radiation" OR "solar energy")` | 551 | 25 |
| Q4 | Eólica: recurso, Guajira y integración | `Colombia AND ("wind energy" OR "wind power" OR "wind farm" OR Guajira)` | 1455 | 25 |
| Q5 | Hidroeléctrica: embalses, GEI y ciclo de vida | `Colombia AND (hydropower OR hydroelectric OR reservoir) AND ("greenhouse gas" OR "life cycle" OR emissions OR methane)` | 96 | 25 |
| Q6 | Biomasa y bioenergía: residuos, biogás, palma y caña | `Colombia AND (biomass OR bioenergy OR biogas OR "palm oil" OR "sugarcane" OR bagasse) AND (energy OR emissions OR electricity)` | 526 | 25 |
| Q7 | Contabilidad de carbono: factor de emisión de la red, ACV, huella | `Colombia AND ("emission factor" OR "grid emission" OR "life cycle assessment" OR "carbon footprint") AND (electricity OR power OR energy)` | 140 | 25 |
| Q8 | Español: fuentes no convencionales y reducción de emisiones | `"fuentes no convencionales" OR ("energías renovables" AND Colombia) OR ("emisiones de CO2" AND Colombia)` | 339 | 25 |
| Q9 | Descarbonización profunda, NDC y precio al carbono | `Colombia AND ("deep decarbonization" OR "net zero" OR "carbon neutrality" OR "nationally determined contribution" OR "carbon tax")` | 108 | 25 |
| Q10 | Emisiones evitadas / potencial de mitigación de las FNCER | `Colombia AND ("avoided emissions" OR "emission reductions" OR "mitigation potential" OR "CO2 savings") AND (solar OR wind OR biomass OR hydro OR renewable)` | 52 | 25 |
| Q11 | Análisis de ciclo de vida de cadenas energéticas colombianas | `Colombia AND "life cycle assessment" AND (electricity OR palm OR biodiesel OR photovoltaic OR hydro OR biogas)` | 42 | 25 |
| Q12 | Hidroelectricidad, variabilidad climática y El Niño (desplazamiento térmico) | `Colombia AND (hydropower OR hydroelectric OR "small hydro" OR "run-of-river") AND (generation OR potential OR climate OR "El Nino")` | 204 | 25 |

Total de filas exportadas: **300**; trabajos únicos: **225**.

Las URLs completas de cada consulta (auditables) están en `datos/busqueda_openalex.json` → `log[].url`.

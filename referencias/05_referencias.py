#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_referencias.py
Genera la bibliografía curada a partir de metadatos OFICIALES de Crossref
(sin escribir a mano autores, años, volúmenes ni páginas).

Entrada: lista curada CURADA = [(eje, doi, justificacion), ...]
Salidas:
  docs/referencias.bib            BibTeX (importable en Zotero/Mendeley/LaTeX)
  docs/referencias_academicas.md  Lista APA 7 anotada, agrupada por eje
  datos/referencias_verificadas.json  metadatos crudos + licencia del editor

Uso: python3 scripts/05_referencias.py    (solo biblioteca estándar)
"""
import json
import re
import time
import datetime
import urllib.parse
import urllib.request
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
HOY = datetime.date.today().isoformat()
UA = "referencias-fncer-colombia/1.0 (mailto:repositorio-academico@example.org)"

CURADA = [
    # --- A. Marco de política, mercado y difusión de las FNCER -----------------
    ("A. Política, mercado y difusión de FNCER", "10.1016/j.rser.2018.03.015",
     "Barreras institucionales y de mercado para la entrada de FNCER al mercado eléctrico colombiano; base para variables de política."),
    ("A. Política, mercado y difusión de FNCER", "10.1016/j.renene.2019.02.138",
     "Dinámica de difusión tecnológica y diseño de incentivos (Ley 1715); útil para modelos de adopción."),
    ("A. Política, mercado y difusión de FNCER", "10.1016/j.jclepro.2014.12.090",
     "Evaluación de esquemas de apoyo a fotovoltaica en Medellín; caso subnacional replicable."),
    ("A. Política, mercado y difusión de FNCER", "10.15665/rp.v13i1.358",
     "Literatura colombiana en español (acceso abierto diamante) sobre el rol de las FNCE en el sector eléctrico."),
    ("A. Política, mercado y difusión de FNCER", "10.17230/ingciencia.13.26.5",
     "Impacto de la integración de FNCER en el mercado eléctrico; fuente nacional revisada por pares."),
    ("A. Política, mercado y difusión de FNCER", "10.18270/cuaderlam.v19i36.4052",
     "Análisis de políticas públicas para adopción de FNCER; perspectiva de gestión pública."),
    # --- B. Solar ------------------------------------------------------------
    ("B. Solar fotovoltaica", "10.1016/j.renene.2019.10.066",
     "Cuantificación del potencial solar nacional y su marco regulatorio; referencia para supuestos de factor de planta."),
    ("B. Solar fotovoltaica", "10.1016/j.rser.2018.04.065",
     "Inventario y estado del arte de la fotovoltaica en Colombia; serie histórica de capacidad."),
    ("B. Solar fotovoltaica", "10.1016/j.seta.2019.100531",
     "Potencial solar y eólico en La Guajira; insumo para análisis territorial de mitigación."),
    ("B. Solar fotovoltaica", "10.15665/rp.v13i2.487",
     "Medición y caracterización de radiación solar en La Guajira; validación local de atlas globales."),
    ("B. Solar fotovoltaica", "10.3390/en13071662",
     "Evaluación energética, económica y ambiental de sistema solar-eólico con HOMER; incluye emisiones evitadas."),
    # --- C. Eólica -----------------------------------------------------------
    ("C. Eólica terrestre y marina", "10.1016/j.jclepro.2019.02.174",
     "Panorama de renovables y oportunidad de eólica marina en Colombia; escenarios de capacidad."),
    ("C. Eólica terrestre y marina", "10.1016/j.renene.2021.03.033",
     "Evaluación del recurso eólico del Caribe con observaciones de 10 min y ERA5; método replicable con datos abiertos."),
    ("C. Eólica terrestre y marina", "10.1016/j.rser.2020.110318",
     "Complementariedad anual e interanual entre fuentes renovables; clave para estimar desplazamiento térmico real."),
    ("C. Eólica terrestre y marina", "10.1016/j.techsoc.2017.01.002",
     "Barreras a la difusión eólica mediante funciones del sistema de innovación tecnológica."),
    ("C. Eólica terrestre y marina", "10.14483/23448350.12304",
     "Metodología nacional para caracterizar viento y evaluar potencial energético."),
    # --- D. Hidroeléctrica ---------------------------------------------------
    ("D. Hidroelectricidad, clima y GEI", "10.1016/j.enpol.2018.12.057",
     "Impactos climáticos sobre la hidroelectricidad colombiana y trayectorias de adaptación del sector eléctrico."),
    ("D. Hidroelectricidad, clima y GEI", "10.1016/j.apenergy.2018.10.045",
     "Optimización estocástica de un sistema hidro-dominante bajo cambio climático; vínculo hidrología-emisiones."),
    ("D. Hidroelectricidad, clima y GEI", "10.5194/hess-22-2839-2018",
     "Impactos de escala de cuenca del desarrollo hidroeléctrico en la Depresión Momposina; salvaguardas socioambientales."),
    ("D. Hidroelectricidad, clima y GEI", "10.3389/fenvs.2022.1036653",
     "Evaluación acumulada de proyectos hidroeléctricos frente a evaluaciones individuales; planificación estratégica."),
    ("D. Hidroelectricidad, clima y GEI", "10.20508/ijrer.v6i1.3012.g6756",
     "Aplicación de la metodología MDL ACM0002 a pequeñas centrales hidroeléctricas colombianas; contabilidad de reducciones."),
    # --- E. Biomasa y bioenergía --------------------------------------------
    ("E. Biomasa y bioenergía", "10.1016/j.apenergy.2014.07.004",
     "Metodología para estimar el potencial de biomasa energética aplicada a Colombia; referencia obligada."),
    ("E. Biomasa y bioenergía", "10.1016/j.jclepro.2020.122317",
     "Potencial energético de residuos agrícolas, agroindustriales, pecuarios y de frigoríficos."),
    ("E. Biomasa y bioenergía", "10.1016/j.biortech.2018.12.007",
     "ACV de biodigestores de bajo costo en fincas colombianas; beneficios ambientales cuantificados."),
    ("E. Biomasa y bioenergía", "10.1016/j.rser.2022.113097",
     "Digestión anaerobia en países en desarrollo con caso colombiano; estado del arte y barreras."),
    ("E. Biomasa y bioenergía", "10.14483/22484728.14066",
     "Biomasa residual pecuaria como recurso energético en Colombia; datos de residuos y potencial."),
    ("E. Biomasa y bioenergía", "10.4067/s0718-07642019000100063",
     "Costos nivelados de electricidad de cogeneración con biomasa forestal; comparación económica."),
    ("E. Biomasa y bioenergía", "10.1016/j.jclepro.2020.120757",
     "Emisiones de GEI y desempeño económico del sector palmicultor colombiano; base para ACV de biocombustibles."),
    # --- F. Mitigación, escenarios y contabilidad de carbono -----------------
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.eneco.2015.05.010",
     "Impuesto al carbono y metas de abatimiento en Colombia; cuantifica reducciones de CO2 por instrumento."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.energy.2018.12.051",
     "Modelo LEAP de demanda energética y emisiones de GEI para Colombia; estructura de escenarios replicable."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.esr.2020.100525",
     "Opciones de descarbonización profunda de mediano siglo para Colombia (citado en el AR6 del IPCC)."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.esr.2020.100510",
     "Trayectorias net-zero en América Latina; marco comparativo regional."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.jclepro.2015.08.089",
     "Evaluación integrada oferta-demanda de política de mitigación en el caso colombiano."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.energy.2019.07.135",
     "Integración a gran escala de FNCER en el sistema energético colombiano futuro."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.renene.2018.07.099",
     "Optimización de la inserción de renovables en el sector eléctrico colombiano."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.seta.2018.10.008",
     "Efectos de largo plazo de generación 100% renovable en el mercado eléctrico colombiano."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.ijepes.2023.109614",
     "Cuantificación de flexibilidad del sistema para la transición energética; restricciones técnicas del desplazamiento."),
    ("F. Mitigación de CO2, escenarios y contabilidad", "10.1016/j.renene.2020.05.055",
     "Expansión óptima con renovables considerando la crisis de Hidroituango; riesgo hidro y emisiones."),
    # --- G. ZNI, microrredes y acceso rural ---------------------------------
    ("G. Zonas No Interconectadas y acceso rural", "10.1016/j.rser.2015.04.176",
     "Microrredes rurales y su aplicación potencial en Colombia; base para desplazamiento de diésel en ZNI."),
    ("G. Zonas No Interconectadas y acceso rural", "10.1016/j.renene.2016.05.086",
     "Factibilidad tecnoeconómica de sistemas FV, eólicos, diésel e híbridos fuera de red en Colombia."),
    ("G. Zonas No Interconectadas y acceso rural", "10.1016/j.rser.2017.05.101",
     "Sistemas fotovoltaicos en áreas remotas para comunidades indígenas; dimensión social de la mitigación."),
    ("G. Zonas No Interconectadas y acceso rural", "10.14483/22484728.14260",
     "Metodología de evaluación de recursos renovables para microrredes aisladas."),
    # --- H. ACV y huella de carbono -----------------------------------------
    ("H. Análisis de ciclo de vida y huella de carbono", "10.1007/s11367-016-1097-6",
     "ACV de biodiésel de palma colombiano; factores de emisión de cadena completa."),
    ("H. Análisis de ciclo de vida y huella de carbono", "10.1007/s11367-020-01855-0",
     "ACV de biorrefinerías energéticas con tallos de café cortados en Colombia."),
    ("H. Análisis de ciclo de vida y huella de carbono", "10.1007/s11356-022-22119-4",
     "Huella de carbono institucional con ISO 14064-1 y GHG Protocol; protocolo aplicable a alcance 2 con FE del SIN."),
]


def get_json(url, intentos=4, timeout=60):
    espera = 3
    for k in range(intentos):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except Exception:
            if k == intentos - 1:
                raise
            time.sleep(espera)
            espera *= 2


def crossref(doi):
    return get_json("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]


def anio(m):
    for k in ("published-print", "published-online", "published", "issued", "created"):
        p = (m.get(k) or {}).get("date-parts")
        if p and p[0] and p[0][0]:
            return int(p[0][0])
    return None


def apa_autores(m):
    aus = m.get("author") or []
    partes = []
    for a in aus[:20]:
        fam = a.get("family") or a.get("name", "")
        ini = " ".join(f"{p[0]}." for p in (a.get("given") or "").replace("-", " ").split() if p)
        partes.append(f"{fam}, {ini}".strip().rstrip(","))
    if not partes:
        return "[Autor institucional]"
    if len(partes) == 1:
        return partes[0]
    if len(aus) > 20:
        return ", ".join(partes[:19]) + ", … " + partes[-1]
    return ", ".join(partes[:-1]) + ", & " + partes[-1]


def apa(m, doi):
    t = (m.get("title") or [""])[0]
    t = re.sub(r"\s+", " ", t).strip()
    revista = (m.get("container-title") or [""])[0]
    vol = m.get("volume") or ""
    num = m.get("issue") or ""
    pag = (m.get("page") or "").replace("-", "–")
    art = m.get("article-number") or ""
    det = vol
    if num:
        det += f"({num})"
    if pag:
        det += f", {pag}"
    elif art:
        det += f", Artículo {art}"
    cuerpo = f"{apa_autores(m)} ({anio(m)}). {t}. "
    if revista:
        cuerpo += f"*{revista}*"
        if det:
            cuerpo += f", {det}"
        cuerpo += ". "
    else:
        cuerpo += f"{m.get('publisher','')}. "
    return cuerpo + f"https://doi.org/{doi}"


def bib_key(m, doi):
    aus = m.get("author") or []
    fam = (aus[0].get("family") if aus else "anon").lower()
    fam = re.sub(r"[^a-z]", "", fam)
    return f"{fam}{anio(m)}" + re.sub(r"[^a-z0-9]", "", doi.split("/")[-1].lower())[:8]


def bib_entry(m, doi, eje, nota):
    campos = {
        "title": (m.get("title") or [""])[0],
        "author": " and ".join(
            f"{a.get('family','')}, {a.get('given','')}" for a in (m.get("author") or [])),
        "year": str(anio(m)),
        "journal": (m.get("container-title") or [""])[0],
        "volume": m.get("volume", ""),
        "number": m.get("issue", ""),
        "pages": m.get("page", ""),
        "doi": doi,
        "url": f"https://doi.org/{doi}",
        "publisher": m.get("publisher", ""),
        "note": f"Eje: {eje}. {nota}",
    }
    lics = m.get("license") or []
    if lics:
        campos["license"] = "; ".join(sorted({l.get("URL", "") for l in lics}))
    tipo = "article" if m.get("type") in ("journal-article", "article") else "misc"
    cuerpo = ",\n".join(f"  {k} = {{{v}}}" for k, v in campos.items() if v)
    return f"@{tipo}{{{bib_key(m, doi)},\n{cuerpo}\n}}", campos.get("license", "no depositada en Crossref")


def main():
    registros, fallos = [], []
    for eje, doi, nota in CURADA:
        try:
            m = crossref(doi)
        except Exception as e:  # noqa: BLE001
            print(f"  FALLO {doi}: {str(e)[:70]}")
            fallos.append({"doi": doi, "eje": eje, "error": str(e)})
            continue
        entrada, lic = bib_entry(m, doi, eje, nota)
        registros.append({"eje": eje, "doi": doi, "nota": nota, "bib": entrada,
                          "apa": apa(m, doi), "licencia": lic,
                          "anio": anio(m), "revista": (m.get("container-title") or [""])[0],
                          "titulo": (m.get("title") or [""])[0],
                          "tipo_crossref": m.get("type"),
                          "consulta": HOY, "crudo": m})
        print(f"  OK {doi} ({anio(m)}) {lic[:48]}")
        time.sleep(0.4)

    (RAIZ / "datos" / "referencias_verificadas.json").write_text(
        json.dumps({"fecha": HOY, "total": len(registros), "fallos": fallos,
                    "referencias": registros}, ensure_ascii=False, indent=1), encoding="utf-8")

    (RAIZ / "docs" / "referencias.bib").write_text(
        f"% Bibliografía curada: FNCER y mitigación de CO2 en Colombia\n"
        f"% Metadatos obtenidos de Crossref el {HOY}. Licencias según depósito del editor.\n"
        f"% Generado por scripts/05_referencias.py — no editar a mano.\n\n"
        + "\n\n".join(r["bib"] for r in registros) + "\n", encoding="utf-8")

    md = ["# Referencias académicas curadas (APA 7, anotadas)",
          "",
          f"Metadatos verificados contra **Crossref** el {HOY}. "
          "La columna de licencia refleja lo que el editor depositó en Crossref; "
          "«no depositada» no significa «todos los derechos reservados»: hay que "
          "revisar la página del artículo. Los archivos `.bib` están en `docs/referencias.bib`.",
          ""]
    ejes = []
    for r in registros:
        if r["eje"] not in ejes:
            ejes.append(r["eje"])
    for eje in ejes:
        md += [f"## {eje}", ""]
        for r in [x for x in registros if x["eje"] == eje]:
            md.append(f"- {r['apa']}")
            md.append(f"  - **Licencia (Crossref):** {r['licencia']}")
            md.append(f"  - **Por qué está aquí:** {r['nota']}")
        md.append("")
    if fallos:
        md += ["## DOIs que no resolvieron en Crossref (revisar manualmente)", ""]
        md += [f"- `{f['doi']}` — {f['error'][:120]}" for f in fallos]
    (RAIZ / "docs" / "referencias_academicas.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"\n{len(registros)} referencias verificadas; {len(fallos)} fallos.")


if __name__ == "__main__":
    main()
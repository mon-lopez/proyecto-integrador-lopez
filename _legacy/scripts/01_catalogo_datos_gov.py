#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_catalogo_datos_gov.py
Cosecha reproducible del catálogo de datos abiertos de Colombia (datos.gov.co,
plataforma Socrata) para el dominio: FNCER (solar, eólica, hidro, biomasa) y
mitigación de CO2.

Qué hace:
  1) Descubre conjuntos por palabras clave (Discovery API de Socrata).
  2) Añade una lista curada de IDs relevantes verificados manualmente.
  3) Descarga la ficha completa de cada conjunto (api/views/<id>.json), que
     incluye LICENCIA, entidad responsable, columnas y fechas de actualización.
  4) Exporta datos/catalogo_datos_gov.json y datos/catalogo_datos_gov.csv

Sin dependencias externas: solo biblioteca estándar de Python 3.8+.
Uso:  python3 scripts/01_catalogo_datos_gov.py
Nota: requiere salida a Internet. Respeta pausas de 0.25 s entre peticiones.
"""
import csv
import json
import time
import datetime
import urllib.parse
import urllib.request
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
SALIDA = RAIZ / "datos"
SALIDA.mkdir(parents=True, exist_ok=True)
UA = "catalogo-fncer-colombia/1.0 (+repositorio academico; uso no comercial de la API)"
HOY = datetime.date.today().isoformat()

# --- Palabras clave de descubrimiento (búsqueda libre en el portal) ----------
TERMINOS = [
    "FNCER", "energía solar", "fotovoltaic", "eólica", "viento", "biomasa",
    "biogás", "generación de energía", "capacidad instalada eléctrica",
    "hidroeléctrica", "embalses", "gases de efecto invernadero",
    "inventario nacional de emisiones", "factores de emisión",
    "radiación", "balance energético", "zonas no interconectadas",
    "hidrógeno", "comunidades energéticas", "autogeneración",
]

# --- Semilla curada: IDs verificados como pertinentes al objeto de estudio ---
# etiqueta temática -> lista de IDs Socrata
CURADO = {
    "FNCER / capacidad y proyectos": [
        "vy9n-w6hc",  # Meta FNCER (MinEnergía) - incluye Emisiones CO2 [Ton/año]
        "nyyd-gdfx",  # FNCER
        "tq63-itzr",  # Producción de electricidad por fuente de combustible (UPME)
        "ghd6-jnjh",  # Conjuntos de datos obligatorios UPME
        "chv4-3d2d",  # Meteorología local en plantas solares
    ],
    "Operación del SIN y ZNI": [
        "3ebi-d83g",  # Estado prestación del servicio en ZNI (IPSE)
        "3a44-zwt6",  # Superservicios - registro de generación diaria ZNI
        "p7wb-z8zs",  # Superservicios - operación diaria otro tipo de generación
        "ptj2-kj7i",  # Costos de racionamiento de energía (UPME)
        "aqdb-bmuv",  # Balance Energético Colombiano (vista)
    ],
    "Inventarios de GEI y factores de emisión": [
        "6rff-a5ep",  # Series históricas INGEI 1990-2024
        "28sd-dsje",  # Factores de emisión del inventario nacional
        "pvsp-rfj6",  # Datos de actividad del inventario nacional
        "sx6m-p6bv",  # Series históricas inventario subnacional
        "s3gg-fcvi",  # Emisiones brutas GEI por departamento y municipio
        "yrhn-3pw5",  # Emisiones netas GEI por departamento y municipio
        "qmgk-t9vt",  # Emisiones GEI y carbono negro por cartera sectorial
        "xs9n-2dua",  # Series históricas contaminantes criterio y carbono negro
        "tzxy-chv2",  # Emisiones per cápita
        "isve-zbhc",  # Cambio porcentual de emisiones por subcategoría
    ],
    "Recurso solar, eólico e hídrico (IDEAM)": [
        "rv9s-8nv6",  # Radiación global acumulada diaria (mensual/anual)
        "sgfv-3yp8",  # Velocidad del viento
        "kiw7-v9ta",  # Dirección del viento
        "s54a-sgyg",  # Precipitación
        "bdmn-sqnh",  # Nivel instantáneo del río
        "vfth-yucv",  # Nivel máximo del río
        "5kjg-nuda",  # Zonificación hidrográfica
        "hp9r-jxuu",  # Catálogo nacional de estaciones IDEAM
        "57sv-p2fu",  # Datos de estaciones IDEAM y de terceros
        "nsz2-kzcq",  # Normales climatológicas de Colombia
    ],
}

LICENCIAS_SPODX = {
    "Creative Commons Attribution 4.0 International": "CC BY 4.0",
    "Creative Commons Attribution | Share Alike 4.0 International": "CC BY-SA 4.0",
    "Creative Commons Attribution | NoDerivatives 4.0 International License": "CC BY-ND 4.0",
    "Creative Commons Attribution | NonCommercial 4.0 International": "CC BY-NC 4.0",
    "Creative Commons Zero | Public Domain Dedication": "CC0 1.0",
    "Open Data Commons Open Database License (ODbL)": "ODbL 1.0",
}


def get_json(url, timeout=45):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def descubrir(terminos):
    hallados = {}
    log = []
    for t in terminos:
        q = urllib.parse.urlencode(
            {"domains": "www.datos.gov.co", "q": t, "limit": 40, "only": "dataset"}
        )
        url = "https://api.us.socrata.com/api/catalog/v1?" + q
        try:
            d = get_json(url)
        except Exception as e:  # noqa: BLE001
            log.append({"termino": t, "error": str(e)})
            continue
        n = d.get("resultSetSize", 0)
        ids = []
        for r in d.get("results", []):
            res = r["resource"]
            ids.append(res["id"])
            hallados.setdefault(res["id"], res["name"])
        log.append({"termino": t, "resultados": n, "ids": ids})
        time.sleep(0.25)
    return hallados, log


# Filtro de pertinencia para los descubiertos automáticamente (evita ruido de
# calidad del aire urbano, salud, minería de regalías, etc.)
CLAVES_OK = (
    "fncer", "solar", "fotovolta", "eólic", "eolic", "viento", "biomasa", "biogás",
    "biogas", "generación", "generacion", "electric", "eléctric", "hidroel",
    "embalse", "efecto invernadero", "emisiones", "factor de emisión",
    "radiación", "radiacion", "balance energético", "balance energetico",
    "no interconectadas", "zni", "hidrógeno", "hidrogeno", "energéticas",
    "energeticas", "autogener", "climatológ", "climatolog", "precipitación",
    "hidrográf", "hidrograf", "nivel del rio", "nivel del río", "estaciones",
)
CLAVES_NO = ("calidad del aire", "hospital", " ips ", "salud", "wifi", "delitos",
             "regalías", "regalias", "matrícula", "matricula")


def pertinente(nombre):
    n = (nombre or "").lower()
    if any(k in n for k in CLAVES_NO):
        return False
    return any(k in n for k in CLAVES_OK)


def ficha(id_):
    d = get_json(f"https://www.datos.gov.co/api/views/{id_}.json")
    lic = (d.get("license") or {}).get("name")
    cf = ((d.get("metadata") or {}).get("custom_fields") or {})
    info_ent = cf.get("Información de la Entidad", {})
    info_dat = cf.get("Información de Datos", {})
    ts = d.get("rowsUpdatedAt") or d.get("metadata_updated_at")
    fecha = (
        datetime.datetime.utcfromtimestamp(ts).date().isoformat()
        if isinstance(ts, int) else None
    )
    return {
        "id": id_,
        "nombre": d.get("name"),
        "entidad": info_ent.get("Nombre de la Entidad") or d.get("owner", {}).get("displayName"),
        "sector": info_ent.get("Sector"),
        "descripcion": (d.get("description") or "").strip(),
        "licencia_declarada": lic,
        "licencia_normalizada": LICENCIAS_SPODX.get(lic, lic or "NO DECLARADA - verificar"),
        "url_licencia": (d.get("license") or {}).get("termsLink"),
        "atribucion": d.get("attribution"),
        "url": f"https://www.datos.gov.co/d/{id_}",
        "api_csv": f"https://www.datos.gov.co/resource/{id_}.csv?$limit=50000",
        "api_json": f"https://www.datos.gov.co/resource/{id_}.json?$limit=50000",
        "frecuencia": info_dat.get("Frecuencia de Actualización"),
        "cobertura": info_dat.get("Cobertura Geográfica"),
        "idioma": info_dat.get("Idioma"),
        "actualizacion_datos": fecha,
        "columnas": [c["name"] for c in d.get("columns", []) if not c["name"].startswith(":")],
        "tipo": d.get("viewType"),
        "consulta": HOY,
    }


def main():
    descubiertos, log = descubrir(TERMINOS)
    ids_curados = {i: tema for tema, lst in CURADO.items() for i in lst}
    extra = [i for i, n in descubiertos.items() if i not in ids_curados and pertinente(n)]
    print(f"Descubiertos: {len(descubiertos)} | curados: {len(ids_curados)} | extra pertinentes: {len(extra)}")

    registros = []
    for id_, tema in ids_curados.items():
        try:
            f = ficha(id_)
            f["eje"] = tema
            f["origen_seleccion"] = "curaduria_manual"
            registros.append(f)
            print(f"  OK  {id_}  {f['licencia_normalizada']:<14} {f['nombre'][:70]}")
        except Exception as e:  # noqa: BLE001
            print(f"  ERR {id_}: {e}")
        time.sleep(0.25)
    for id_ in extra:
        try:
            f = ficha(id_)
            f["eje"] = "Descubrimiento automatico (revisar pertinencia)"
            f["origen_seleccion"] = "descubrimiento_por_palabra_clave"
            registros.append(f)
        except Exception as e:  # noqa: BLE001
            print(f"  ERR {id_}: {e}")
        time.sleep(0.25)

    (SALIDA / "catalogo_datos_gov.json").write_text(
        json.dumps({"fecha_cosecha": HOY, "fuente": "datos.gov.co (Socrata)",
                    "total": len(registros), "conjuntos": registros},
                   ensure_ascii=False, indent=1), encoding="utf-8")

    campos = ["id", "nombre", "entidad", "eje", "licencia_normalizada",
              "licencia_declarada", "frecuencia", "actualizacion_datos",
              "cobertura", "url", "api_csv", "descripcion"]
    with (SALIDA / "catalogo_datos_gov.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for r in registros:
            w.writerow(r)

    (SALIDA / "log_descubrimiento_datos_gov.json").write_text(
        json.dumps({"fecha": HOY, "consultas": log}, ensure_ascii=False, indent=1),
        encoding="utf-8")
    print(f"\nExportados {len(registros)} conjuntos -> datos/catalogo_datos_gov.[json|csv]")
    sin_lic = [r["id"] for r in registros if not r["licencia_declarada"]]
    print("Sin licencia declarada (verificar antes de reutilizar):", sin_lic)


if __name__ == "__main__":
    main()
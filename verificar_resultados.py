"""Audita los CSV y reproduce los tres escenarios finales sin repetir las grillas.

Ejecutar: python verificar_resultados.py
Complementa las pruebas del método; no certifica un máximo global continuo.
"""

import csv
from dataclasses import replace
from itertools import product
from pathlib import Path
import re

import numpy as np

from experimentos import indicadores
from modelo import PARAMETROS_NOMINALES, euler


BASE = Path(__file__).resolve().parent
NOMBRES = ("A", "beta", "gamma", "sigma")
RANGOS = ((0.25, 0.60), (0.10, 0.30), (0.05, 0.20), (0.40, 1.00))
MEDIDAS = ("SP", "desfase", "Pmax", "Pmin", "Pmed", "residuo_anual")


def exigir(condicion, mensaje):
    if not condicion:
        raise ValueError(mensaje)


def leer_csv(base, nombre):
    def convertir(valor):
        if valor in ("True", "False"):
            return valor == "True"
        try:
            return float(valor)
        except ValueError:
            return valor

    with (base / "resultados" / nombre).open(encoding="utf-8-sig", newline="") as archivo:
        return [{k: convertir(v) for k, v in f.items()} for f in csv.DictReader(archivo)]


def clave(fila):
    return tuple(fila[n] for n in NOMBRES)


def comparar_indicadores(fila):
    exigir(all(np.isfinite(fila[k]) for k in MEDIDAS), "Indicadores no finitos.")
    np.testing.assert_allclose(fila["SP"], (fila["Pmax"]-fila["Pmin"])/fila["Pmed"], rtol=1e-11)
    np.testing.assert_allclose(fila["desfase"], fila["t_Pmax_posterior"]-fila["t_Tmax"], atol=1e-11)
    exigir(48 <= fila["t_Tmax"] < fila["t_Pmax_posterior"] <= 60, "Picos fuera del último año o desordenados.")
    exigir(fila["cantidad_picos_P"] == 1 and fila["pico_posterior_es_Pmax"], "Pico de P ambiguo.")
    np.testing.assert_allclose(fila["residuo_anual"], max(fila[f"residuo_{n}"] for n in "THP"))
    admisible = fila["residuo_anual"] < 0.01 and fila["Pmed"] > 0 and fila["minimo_ultimo_ciclo"] >= 0
    exigir(fila["admisible"] == admisible, "Admisibilidad inconsistente con los indicadores.")


def verificar(base=BASE):
    base = Path(base)
    verificaciones = []
    convergencia = leer_csv(base, "convergencia.csv")
    for anterior, actual in zip(convergencia, convergencia[1:]):
        np.testing.assert_allclose(actual["h"], 2*anterior["h"], rtol=1e-12)
        if np.isfinite(actual["orden"]):
            np.testing.assert_allclose(actual["orden"], np.log2(actual["error_T"]/anterior["error_T"]))
    for fila in convergencia:
        if np.isfinite(fila["error_T"]):
            np.testing.assert_allclose(fila["error_sobre_h"], fila["error_T"]/fila["h"])
    verificaciones.append("Pasos, cocientes E/h y órdenes observados coherentes.")

    grilla = leer_csv(base, "grilla.csv")
    por_clave = {clave(f): f for f in grilla}
    exigir(len(grilla) == len(por_clave), "Hay parámetros duplicados en la grilla.")
    for fila in grilla:
        exigir(all(a <= fila[n] <= b for n, (a, b) in zip(NOMBRES, RANGOS)), "Parámetro fuera del intervalo pedido.")
        exigir(fila["finito_acotado_60m"], "La grilla guardada contiene trayectorias no válidas.")
        comparar_indicadores(fila)
    ejes = [np.round(np.linspace(a, b, 9), 12) for a, b in RANGOS]
    puntos_globales = list(product(*ejes))
    exigir(set(puntos_globales) <= por_clave.keys(), "La grilla global de 9⁴ puntos está incompleta.")
    mejor = max((f for f in grilla if f["admisible"]), key=lambda f: f["SP"])
    verificaciones.append(f"{len(grilla)} escenarios únicos dentro de los rangos; grilla global completa.")

    sp = np.array([por_clave[q]["SP"] for q in puntos_globales]).reshape((9,)*4)
    tendencias = leer_csv(base, "tendencias_conjuntas.csv")
    for i, n in enumerate(NOMBRES):
        fila = next(f for f in tendencias if f["parametro"] == n)
        diferencias = fila["signo_esperado"]*np.diff(sp, axis=i)
        exigir(fila["comparaciones"] == diferencias.size, "Conteo de comparaciones incorrecto.")
        exigir(fila["inversiones"] == np.count_nonzero(diferencias < -1e-10), "Tendencia conjunta incorrecta.")
        np.testing.assert_allclose(fila["menor_cambio_segun_signo"], diferencias.min())
    verificaciones.append("Tendencias conjuntas reconstruidas a partir de la grilla guardada.")

    sensibilidad = leer_csv(base, "sensibilidad.csv")
    for n, (a, b) in zip(NOMBRES, RANGOS):
        filas = [f for f in sensibilidad if f["parametro"] == n]
        esperados = set(np.round(np.r_[np.linspace(a, b, 21), getattr(PARAMETROS_NOMINALES, n)], 12))
        exigir({f["valor"] for f in filas} == esperados and len(filas) == len(esperados), "Barrido individual incompleto o duplicado.")
        for fila in filas:
            comparar_indicadores(fila)
            for otro in NOMBRES:
                esperado = fila["valor"] if otro == n else getattr(PARAMETROS_NOMINALES, otro)
                np.testing.assert_allclose(fila[otro], esperado, rtol=0, atol=1e-12)
    verificaciones.append("Barridos individuales completos; solo cambia el parámetro indicado.")

    paso = leer_csv(base, "verificacion_paso.csv")[-1]
    h = paso["h_fino"]
    exigir(paso["max_diferencia_rel_SP"] < 0.005 and paso["max_diferencia_desfase"] < 0.02
           and paso["max_diferencia_rel_extremos_promedio_P"] < 0.005 and paso["cambios_admisibilidad"] == 0,
           "No se cumplen las tolerancias del paso declarado.")
    comparacion = leer_csv(base, "comparacion.csv")
    seleccionado = next(f for f in comparacion if f["escenario"] == "mejor_grilla")
    exigir(clave(seleccionado) == clave(mejor), "El escenario seleccionado no maximiza SP entre los admisibles.")
    finales = leer_csv(base, "verificacion_finalistas.csv")
    exigir({f["h"] for f in finales} == {h, h/2, h/4}, "Refinamientos de finalistas incompletos.")
    for paso_final in (h, h/2, h/4):
        filas = [f for f in finales if f["h"] == paso_final]
        exigir(clave(max(filas, key=lambda f: f["SP"])) == clave(mejor), "El ganador cambia al refinar.")
    verificaciones.append("Paso dentro de tolerancias y ganador consistente con la grilla y los refinamientos.")

    extremos = leer_csv(base, "acotacion_60m.csv")
    for fila in comparacion:
        nombre = "mejor" if fila["escenario"] == "mejor_grilla" else fila["escenario"]
        p = replace(PARAMETROS_NOMINALES, **dict(zip(NOMBRES, clave(fila))))
        datos = np.genfromtxt(base / "resultados" / f"trayectoria_{nombre}.csv", delimiter=",", skip_header=1)
        t, y = euler(h, parametros=p)
        np.testing.assert_allclose(datos[:, 0], t, rtol=0, atol=1e-11)
        np.testing.assert_allclose(datos[:, 1:], y, rtol=1e-10, atol=1e-11)
        r = indicadores(t, y)
        np.testing.assert_allclose([r[k] for k in MEDIDAS], [fila[k] for k in MEDIDAS], rtol=1e-8, atol=1e-10)
        ultimo = t >= 48
        solo_ultimo = indicadores(t[ultimo], y[ultimo])
        np.testing.assert_allclose([solo_ultimo[k] for k in MEDIDAS[:-1]], [fila[k] for k in MEDIDAS[:-1]], rtol=1e-8)
        exigir(np.isnan(solo_ultimo["residuo_anual"]), "Se está declarando repetición anual sin el ciclo anterior.")
        extremo = next(f for f in extremos if f["escenario"] == nombre)
        np.testing.assert_allclose(extremo["h"], h/4)
        _, fina = euler(h/4, parametros=p)
        for j, variable in enumerate("THP"):
            np.testing.assert_allclose([fina[:, j].min(), fina[:, j].max()],
                                       [extremo[f"min_{variable}_60m"], extremo[f"max_{variable}_60m"]], rtol=1e-10, atol=1e-11)
    verificaciones.append("Tres trayectorias e indicadores reproducidos; extremos de los 60 meses verificados a h/4.")

    for nombre in ("README.md", "decisiones.md", "revision.md", "resultados/resumen.md"):
        archivo = base / nombre
        for enlace in re.findall(r"\]\(([^)]+)\)", archivo.read_text(encoding="utf-8")):
            exigir((archivo.parent / enlace).exists(), f"Enlace local roto en {nombre}: {enlace}")
    verificaciones.append("Los enlaces locales de la documentación y el informe existen.")
    return verificaciones


if __name__ == "__main__":
    for mensaje in verificar():
        print(f"OK: {mensaje}")
    print("Auditoría terminada: no se detectaron inconsistencias en los resultados comprobados.")

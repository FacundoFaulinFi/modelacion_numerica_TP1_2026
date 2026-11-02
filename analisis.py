"""Ejecuta las dos partes del TP y guarda tablas, gráficos y un resumen.

Ejecutar desde cualquier carpeta: python analisis.py
Las salidas se escriben junto al script, en resultados/ y figuras/.
"""

import csv
import os
from dataclasses import replace
from itertools import product
from pathlib import Path

# Mantiene la caché de fuentes dentro del proyecto, también en entornos restringidos.
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".cache" / "matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from experimentos import evaluar_lote, indicadores
from modelo import PARAMETROS_NOMINALES as NOMINAL, euler, jacobiano


BASE = Path(__file__).resolve().parent
RESULTADOS = BASE / "resultados"
FIGURAS = BASE / "figuras"
NOMBRES = ("A", "beta", "gamma", "sigma")
RANGOS = {"A": (0.25, 0.60), "beta": (0.10, 0.30),
          "gamma": (0.05, 0.20), "sigma": (0.40, 1.00)}
ETIQUETAS = ("T [°C]", "H [u.h.]", "P [u.p.]")
ETIQUETAS_PARAMETROS = {
    "A": "A [mes⁻¹]", "beta": "β [u.h. mes⁻¹ °C⁻²]",
    "gamma": "γ [mes⁻¹ u.p.⁻¹]", "sigma": "σ [mes⁻¹]",
}


def guardar_csv(nombre, filas):
    columnas = list(dict.fromkeys(k for f in filas for k in f))
    with (RESULTADOS / nombre).open("w", encoding="utf-8-sig", newline="") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=columnas)
        escritor.writeheader()
        escritor.writerows(filas)


def guardar_figura(fig, nombre):
    fig.savefig(FIGURAS / f"{nombre}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def parte1():
    print("Parte 1: referencias de Euler y convergencia...", flush=True)
    h_ref = 12 / 2**18
    _, y_ref_muy_gruesa = euler(4 * h_ref)
    _, y_ref_gruesa = euler(2 * h_ref)
    t_ref, y_ref = euler(h_ref)
    diferencia_ref = float(np.max(np.abs(y_ref_gruesa[:, 0] - y_ref[::2, 0])))
    diferencia_anterior = float(np.max(np.abs(y_ref_muy_gruesa[:, 0] - y_ref_gruesa[::2, 0])))
    orden_ref = float(np.log2(diferencia_anterior / diferencia_ref))
    if not 0.95 < orden_ref < 1.05:
        raise RuntimeError("Las referencias todavía no muestran el orden asintótico esperado.")
    guardar_csv("referencia.csv", [
        {"h_grueso": 4*h_ref, "h_fino": 2*h_ref, "diferencia_T": diferencia_anterior},
        {"h_grueso": 2*h_ref, "h_fino": h_ref, "diferencia_T": diferencia_ref, "orden": orden_ref},
    ])
    filas = []
    for k in range(9, -1, -1):
        h = 12 / 2**k
        fila = {"h": h, "error_T": float("nan"), "orden": float("nan"),
                "error_sobre_h": float("nan"), "estado": "desbordamiento"}
        try:
            _, y = euler(h)
            referencia = y_ref[::2**(18 - k)]
            error = float(np.max(np.abs(y[:, 0] - referencia[:, 0])))
            fila.update(error_T=error, error_sobre_h=error / h,
                        estado="finita; evaluar estabilidad")
            if np.max(np.abs(y)) > 100 or np.min(y) < -1e-8:
                fila["estado"] = "trayectoria espuria"
            if (filas and np.isfinite(filas[-1]["error_T"]) and error > 0
                    and fila["estado"] != "trayectoria espuria"
                    and filas[-1]["estado"] != "trayectoria espuria"):
                fila["orden"] = float(np.log2(error / filas[-1]["error_T"]))
        except FloatingPointError:
            pass
        filas.append(fila)
    razon_ref = diferencia_ref / filas[0]["error_T"]
    if razon_ref >= 0.01:
        raise RuntimeError("La referencia requiere más refinamiento.")
    pendiente = float(np.polyfit(np.log([f["h"] for f in filas[:4]]),
                                 np.log([f["error_T"] for f in filas[:4]]), 1)[0])
    guardar_csv("convergencia.csv", filas)
    # Los valores explosivos finitos se conservan en la tabla; su escala
    # ocultaría el rango convergente del gráfico.
    finitas = [f for f in filas if f["estado"] == "finita; evaluar estabilidad"]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.loglog([f["h"] for f in finitas], [f["error_T"] for f in finitas], "o-", label="Euler")
    hs = np.array([f["h"] for f in finitas[:4]])
    ax.loglog(hs, hs * filas[0]["error_sobre_h"], "--", label="Pendiente 1")
    ax.set(xlabel="Paso h [meses]", ylabel="Máximo error absoluto de T [°C]",
           title=f"Convergencia: pendiente inicial = {pendiente:.4f}")
    ax.set_xticks([f["h"] for f in finitas], [f"{f['h']:.4g}" for f in finitas])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend()
    guardar_figura(fig, "01_convergencia")

    print("Parte 1: exploración de estabilidad...", flush=True)
    barrido = []
    for h in np.linspace(0.25, 2, 71):
        fila = {"h": float(h), "error_T": float("nan"), "error_sobre_h": float("nan"),
                "max_abs": float("nan"), "estado": "desbordamiento"}
        try:
            t, y = euler(h)
            # Solo este barrido auxiliar interpola la referencia; Parte 1(b) no.
            error = float(np.max(np.abs(y[:, 0] - np.interp(t, t_ref, y_ref[:, 0]))))
            fila.update(error_T=error, error_sobre_h=error / h,
                        max_abs=float(np.max(np.abs(y))), estado="finita")
        except FloatingPointError:
            pass
        barrido.append(fila)
    guardar_csv("estabilidad.csv", barrido)

    # El barrido de 0.025 meses salteaba trayectorias explosivas todavía finitas.
    detalle = []
    for h in np.round(np.linspace(1.15, 1.30, 1501), 6):
        fila = {"h": float(h), "error_T": float("nan"), "error_sobre_h": float("nan"),
                "max_abs": float("nan"), "estado": "desbordamiento"}
        try:
            t, y = euler(float(h))
            error = float(np.max(np.abs(y[:, 0] - np.interp(t, t_ref, y_ref[:, 0]))))
            fila.update(error_T=error, error_sobre_h=error / h,
                        max_abs=float(np.max(np.abs(y))), estado="finita")
        except FloatingPointError:
            pass
        detalle.append(fila)
    guardar_csv("estabilidad_detalle.csv", detalle)
    primer_salto = next(f for f in detalle if np.isfinite(f["error_T"]) and
                        f["error_sobre_h"] > 10 * filas[0]["error_sobre_h"])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    finitos_detalle = [f for f in detalle if f["estado"] == "finita"]
    ax.semilogy([f["h"] for f in finitos_detalle],
                [f["error_sobre_h"] for f in finitos_detalle], ".", ms=2,
                label="Trayectorias finitas (incluye explosivas)")
    ax.axhline(filas[0]["error_sobre_h"], ls="--", color="gray", label="E/h en malla fina")
    ax.set(xlabel="h [meses]", ylabel="Error T / h",
           title="Barrido fino: crecimiento abrupto del error antes del desbordamiento")
    ax.set_ylim(bottom=0.1)
    ax.legend(fontsize=8)
    guardar_figura(fig, "08_estabilidad_detalle")

    # Diagnóstico de Euler para modos disipativos de la linealización congelada.
    locales = []
    mask = t_ref >= 48
    for ti, yi in zip(t_ref[mask][::256], y_ref[mask][::256]):
        valores = np.linalg.eigvals(jacobiano(ti, yi))
        disipativos = valores[valores.real < 0]
        limites = -2 * disipativos.real / np.abs(disipativos)**2
        locales.append({"t": float(ti), "h_limite_local": float(limites.min())})
    guardar_csv("estabilidad_lineal_local.csv", locales)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    finitos = [f for f in barrido if f["estado"] == "finita"]
    axs[0].semilogy([f["h"] for f in finitos], [f["error_sobre_h"] for f in finitos], ".-")
    axs[0].axhline(filas[0]["error_sobre_h"], ls="--", color="gray", label="E/h en malla fina")
    axs[0].legend()
    axs[0].set(xlabel="h [meses]", ylabel="Error T / h", title="Pérdida del comportamiento lineal")
    # P revela oscilaciones que pueden ser menos visibles en T.
    axs[1].plot(t_ref[::1024], y_ref[::1024, 2], color="black", lw=2, label="Referencia")
    for h in (0.375, 0.75, 1.0, 1.2):
        try:
            t, y = euler(h)
            if np.max(np.abs(y)) < 100:
                axs[1].plot(t, y[:, 2], ".-", ms=2, label=f"h={h:g}")
        except FloatingPointError:
            pass
    axs[1].set(xlim=(36, 60), xlabel="t [meses]", ylabel="P [u.p.]", title="Oscilaciones al aumentar h")
    axs[1].legend(fontsize=8)
    guardar_figura(fig, "02_estabilidad")
    resumen = {"h_ref": h_ref, "diferencia_ref": diferencia_ref, "razon_ref": razon_ref,
               "pendiente": pendiente, "filas": filas, "orden_ref": orden_ref,
               "referencia": (t_ref, y_ref), "primer_salto": primer_salto,
               "h_limite_local": min(f["h_limite_local"] for f in locales),
               "primer_desborde": next((f["h"] for f in barrido if f["estado"] == "desbordamiento"), None)}
    return resumen


def crear_parametros(q):
    if len(q) != len(NOMBRES):
        raise ValueError("Se requieren los cuatro valores: A, beta, gamma y sigma.")
    return replace(NOMINAL, **dict(zip(NOMBRES, q)))


def parametros_grilla(cantidad):
    """Grilla uniforme en todo el dominio, sin duplicados de redondeo decimal."""
    ejes = [np.round(np.linspace(*RANGOS[n], cantidad), 12) for n in NOMBRES]
    return [crear_parametros(q) for q in product(*ejes)]


def documentar_busqueda(grillas, filas, h):
    """Compara resoluciones paramétricas y examina efectos conjuntos de SP."""
    por_parametro = {crear_parametros([f[n] for n in NOMBRES]): f for f in filas}
    niveles = []
    for nombre, puntos in grillas.items():
        candidatas = [por_parametro[p] for p in puntos if por_parametro[p]["admisible"]]
        mejor = max(candidatas, key=lambda f: f["SP"])
        niveles.append(dict(nivel=nombre, puntos=len(puntos), h=h,
                            **{n: mejor[n] for n in (*NOMBRES, "SP", "desfase")}))
    guardar_csv("resoluciones_grilla.csv", niveles)
    puntos_globales = grillas["global_9"]
    sp = np.array([por_parametro[p]["SP"] for p in puntos_globales]).reshape((9,)*4)
    tendencias = []
    for i, (n, signo) in enumerate(zip(NOMBRES, (1, 1, -1, 1))):
        diferencias = signo * np.diff(sp, axis=i)
        tendencias.append(dict(parametro=n, signo_esperado=signo,
                               comparaciones=diferencias.size,
                               inversiones=int(np.count_nonzero(diferencias < -1e-10)),
                               menor_cambio_segun_signo=float(diferencias.min())))
    guardar_csv("tendencias_conjuntas.csv", tendencias)
    ejes = [np.unique([getattr(p, n) for p in puntos_globales]) for n in NOMBRES]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, x, y, z, nx, ny, titulo in [
        (axs[0], ejes[0], ejes[1], sp[:, :, 0, -1], "A", "beta", "γ mínimo; σ máximo"),
        (axs[1], ejes[2], ejes[3], sp[-1, -1, :, :], "gamma", "sigma", "A y β máximos"),
    ]:
        mapa = ax.pcolormesh(x, y, z.T, shading="nearest", cmap="viridis")
        ax.set(xlabel=ETIQUETAS_PARAMETROS[nx], ylabel=ETIQUETAS_PARAMETROS[ny], title=titulo)
        ax.grid(False)
        fig.colorbar(mapa, ax=ax, label="SP")
    fig.suptitle("Búsqueda conjunta: cortes de la grilla global de 9⁴ escenarios")
    fig.tight_layout()
    guardar_figura(fig, "09_busqueda_conjunta")
    return niveles, tendencias


def verificar_finalistas(admisibles, h, trayectorias):
    """Refina dos veces los diez mejores y verifica las trayectorias completas."""
    mejores = sorted(admisibles, key=lambda f: f["SP"], reverse=True)[:10]
    parametros = [crear_parametros([f[n] for n in NOMBRES]) for f in mejores]
    maximos = crear_parametros([RANGOS[n][1] for n in NOMBRES])
    casos = list(dict.fromkeys(parametros + [NOMINAL, maximos]))
    filas, ganadores = [], []
    for paso in (h, h / 2, h / 4):
        evaluacion = evaluar_lote(casos, paso)
        ganadora = max(evaluacion[:len(parametros)], key=lambda f: f["SP"])
        ganadores.append(tuple(ganadora[n] for n in NOMBRES))
        for j, fila in enumerate(evaluacion):
            if not fila["admisible"]:
                raise RuntimeError("Un finalista deja de ser admisible al refinar.")
            p = casos[j]
            nombre = ("nominal" if p == NOMINAL else "todos_maximos" if p == maximos
                      else "mejor_grilla" if j == 0 else f"finalista_{j + 1}")
            filas.append(dict(h=paso, candidato=j + 1, escenario=nombre, **fila))
    if len(set(ganadores)) != 1:
        raise RuntimeError("El ganador cambia al refinar: ampliar la búsqueda a menor paso.")
    guardar_csv("verificacion_finalistas.csv", filas)

    cambios, extremos = [], []
    for nombre, p in [("nominal", NOMINAL), ("mejor", parametros[0]), ("todos_maximos", maximos)]:
        t_anterior, y_anterior = trayectorias[nombre]
        for paso in (h / 2, h / 4):
            t, y = euler(paso, parametros=p)
            np.testing.assert_allclose(t[::2], t_anterior, rtol=0, atol=1e-10)
            diferencias = np.max(np.abs(y_anterior - y[::2]), axis=0)
            relativas = diferencias / np.maximum(np.max(np.abs(y), axis=0), 1e-12)
            if np.max(relativas) >= 0.005:
                raise RuntimeError("Las trayectorias requieren mayor refinamiento.")
            cambios.append(dict(escenario=nombre, h_grueso=2*paso, h_fino=paso,
                                **{f"diferencia_abs_{n}": float(diferencias[j]) for j, n in enumerate("THP")},
                                **{f"diferencia_rel_{n}": float(relativas[j]) for j, n in enumerate("THP")}))
            t_anterior, y_anterior = t, y
        extremos.append(dict(escenario=nombre, h=h/4,
                             **{f"min_{n}_60m": float(y[:, j].min()) for j, n in enumerate("THP")},
                             **{f"max_{n}_60m": float(y[:, j].max()) for j, n in enumerate("THP")}))
    guardar_csv("convergencia_trayectorias.csv", cambios)
    guardar_csv("acotacion_60m.csv", extremos)
    return {"filas": filas, "cambios": cambios, "extremos": extremos,
            "ganador_confirmado": ganadores[0]}


def parte2():
    h = 12 / 2048
    gruesa = parametros_grilla(5)
    exploracion = evaluar_lote(gruesa, h)
    cache = dict(zip(gruesa, exploracion))
    global_refinada = parametros_grilla(9)
    adicionales = [p for p in global_refinada if p not in cache]
    print("Parte 2: refinamiento global a 9⁴ combinaciones...", flush=True)
    cache.update(zip(adicionales, evaluar_lote(adicionales, h)))
    mejor = max((cache[p] for p in global_refinada if cache[p]["admisible"]), key=lambda f: f["SP"])
    ejes_locales = []
    for n in NOMBRES:
        a, b = RANGOS[n]
        semiancho = (b - a) / 8
        ejes_locales.append(np.round(np.linspace(max(a, mejor[n] - semiancho), min(b, mejor[n] + semiancho), 5), 12))
    local = [crear_parametros(q) for q in product(*ejes_locales)]
    grilla = list(dict.fromkeys(gruesa + global_refinada + local))

    sensibilidad_p, etiquetas = [], []
    for n in NOMBRES:
        for valor in np.unique(np.round(np.r_[np.linspace(*RANGOS[n], 21), getattr(NOMINAL, n)], 12)):
            sensibilidad_p.append(replace(NOMINAL, **{n: float(valor)}))
            etiquetas.append({"parametro": n, "valor": float(valor)})
    maximos = crear_parametros([RANGOS[n][1] for n in NOMBRES])
    todos = grilla + sensibilidad_p + [NOMINAL, maximos]
    precision = []
    print(f"Parte 2: {len(grilla)} combinaciones distintas y sensibilidad individual...", flush=True)
    faltantes = list(dict.fromkeys(p for p in todos if p not in cache))
    cache.update(zip(faltantes, evaluar_lote(faltantes, h)))
    filas_h = [dict(cache[p]) for p in todos]
    for _ in range(5):
        filas_finas = evaluar_lote(todos, h / 2)
        diferencias_sp, diferencias_lag, diferencias_p = [], [], []
        for a, b in zip(filas_h, filas_finas):
            if not a["finito_acotado_60m"] or not b["finito_acotado_60m"]:
                raise RuntimeError("Hay trayectorias no válidas: revisar la precisión por escenario.")
            if not a["posterior"] or not b["posterior"]:
                raise RuntimeError("Falta un pico posterior: revisar la definición del desfase.")
            diferencias_sp.append(abs(a["SP"] - b["SP"]) / abs(b["SP"]))
            diferencias_lag.append(abs(a["desfase"] - b["desfase"]))
            diferencias_p.extend(abs(a[k] - b[k]) / max(abs(b[k]), 1e-12)
                                 for k in ("Pmax", "Pmin", "Pmed"))
        precision.append({"h_grueso": h, "h_fino": h / 2,
                          "max_diferencia_rel_SP": max(diferencias_sp),
                          "max_diferencia_desfase": max(diferencias_lag),
                          "max_diferencia_rel_extremos_promedio_P": max(diferencias_p),
                          "cambios_admisibilidad": sum(a["admisible"] != b["admisible"] for a, b in zip(filas_h, filas_finas))})
        h /= 2
        filas_h = filas_finas
        if (precision[-1]["max_diferencia_rel_SP"] < 0.005
                and precision[-1]["max_diferencia_desfase"] < 0.02
                and precision[-1]["max_diferencia_rel_extremos_promedio_P"] < 0.005
                and precision[-1]["cambios_admisibilidad"] == 0):
            break
    else:
        raise RuntimeError("No se alcanzó la precisión elegida para los indicadores.")
    guardar_csv("verificacion_paso.csv", precision)
    if not all(f["cantidad_picos_P"] == 1 and f["pico_posterior_es_Pmax"] for f in filas_h):
        raise RuntimeError("Hay ciclos con picos ambiguos: revisar el desfase antes de informar resultados.")
    filas_grilla = filas_h[:len(grilla)]
    conjunto_inicial, conjunto_global = set(gruesa), set(global_refinada)
    for p, fila in zip(grilla, filas_grilla):
        fila["grilla"] = ("inicial" if p in conjunto_inicial else
                          "global_9" if p in conjunto_global else "refinamiento_local")
    filas_sens = filas_h[len(grilla):-2]
    for fila, etiqueta in zip(filas_sens, etiquetas):
        fila.update(etiqueta)
    nominal, todos_maximos = filas_h[-2:]
    admisibles = [f for f in filas_grilla if f["admisible"]]
    mejor = max(admisibles, key=lambda f: f["SP"])
    mejor_p = crear_parametros([mejor[n] for n in NOMBRES])
    guardar_csv("grilla.csv", filas_grilla)
    guardar_csv("sensibilidad.csv", filas_sens)
    niveles, tendencias = documentar_busqueda(
        {"inicial_5": gruesa, "global_9": global_refinada, "local_5": local}, filas_grilla, h)
    comparacion = [dict(escenario=n, **f) for n, f in
                   [("nominal", nominal), ("mejor_grilla", mejor), ("todos_maximos", todos_maximos)]]
    guardar_csv("comparacion.csv", comparacion)

    trayectorias = {}
    for nombre, p in [("nominal", NOMINAL), ("mejor", mejor_p), ("todos_maximos", maximos)]:
        t, y = euler(h, parametros=p)
        trayectorias[nombre] = (t, y)
        guardar_csv(f"trayectoria_{nombre}.csv", [
            {"t": float(ti), "T": float(z[0]), "H": float(z[1]), "P": float(z[2])}
            for ti, z in zip(t, y)])
        # Verifica que el camino escalar coincida con el evaluado por lotes.
        ind = indicadores(t, y)
        esperado = {"nominal": nominal, "mejor": mejor, "todos_maximos": todos_maximos}[nombre]
        np.testing.assert_allclose([ind[k] for k in ("SP", "desfase", "Pmed")],
                                   [esperado[k] for k in ("SP", "desfase", "Pmed")], rtol=1e-8, atol=1e-9)

    print("Parte 2: refinamientos adicionales de finalistas y trayectorias...", flush=True)
    validacion = verificar_finalistas(admisibles, h, trayectorias)

    for ultimo, numero in [(False, "03"), (True, "04")]:
        fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
        t, y = trayectorias["nominal"]
        if ultimo:
            mask = t >= 48
            t, y = t[mask], y[mask]
        for j, ax in enumerate(axs):
            ax.plot(t, y[:, j], color=f"C{j}")
            ax.set(ylabel=ETIQUETAS[j], xlim=(48, 60) if ultimo else (0, 60))
            if not ultimo:
                ax.axvspan(48, 60, color="gray", alpha=0.15)
        if ultimo:
            axs[0].axvline(nominal["t_Tmax"], ls="--", color="black", label="Máximo T")
            axs[2].axvline(nominal["t_Pmax_posterior"], ls="--", color="black", label="Máximo posterior P")
            axs[0].legend()
            axs[2].legend()
        axs[-1].set_xlabel("t [meses]")
        fig.suptitle("Caso nominal: " + ("último ciclo anual" if ultimo else "simulación completa"))
        fig.tight_layout()
        guardar_figura(fig, f"{numero}_nominal_" + ("ultimo_ciclo" if ultimo else "60_meses"))

    fig, axs = plt.subplots(4, 2, figsize=(11, 13))
    for i, n in enumerate(NOMBRES):
        fs = [f for f in filas_sens if f["parametro"] == n]
        for j, (campo, titulo) in enumerate([("SP", "Severidad SP"), ("desfase", "Desfase [meses]")]):
            axs[i, j].plot([f["valor"] for f in fs], [f[campo] for f in fs], "o-", ms=3)
            axs[i, j].axvline(getattr(NOMINAL, n), ls="--", color="gray", label="Nominal")
            axs[i, j].set(xlabel=ETIQUETAS_PARAMETROS[n], ylabel=titulo)
            axs[i, j].legend(fontsize=8)
    fig.suptitle("Sensibilidad individual: los otros parámetros permanecen nominales")
    fig.tight_layout()
    guardar_figura(fig, "05_sensibilidad")

    for ultimo, numero in [(False, "06"), (True, "07")]:
        fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
        for nombre, etiqueta in [("nominal", "Nominal"), ("mejor", "Mejor grilla"), ("todos_maximos", "Todos máximos")]:
            t, y = trayectorias[nombre]
            if ultimo:
                mask = t >= 48
                t, y = t[mask], y[mask]
            for j, ax in enumerate(axs):
                ax.plot(t, y[:, j], label=etiqueta)
                ax.set(ylabel=ETIQUETAS[j], xlim=(48, 60) if ultimo else (0, 60))
        axs[0].legend(ncol=3, fontsize=9)
        axs[-1].set_xlabel("t [meses]")
        fig.suptitle("Comparación de escenarios: " + ("último año" if ultimo else "60 meses"))
        fig.tight_layout()
        guardar_figura(fig, f"{numero}_comparacion_" + ("ultimo_ciclo" if ultimo else "60_meses"))
    return {"h": h, "nominal": nominal, "mejor": mejor, "maximos": todos_maximos,
            "n_grilla": len(grilla), "n_admisibles": len(admisibles),
            "precision": precision[-1], "sensibilidad": filas_sens,
            "validacion": validacion,
            "niveles_busqueda": niveles, "tendencias_conjuntas": tendencias,
            "trayectoria_nominal": trayectorias["nominal"],
            "max_residuo": max(f["residuo_anual"] for f in filas_grilla)}


def resumen(p1, p2):
    lineas = ["# Resultados del TP", "", "Generados por `python analisis.py`. Las decisiones y convenciones se detallan en [decisiones.md](../decisiones.md).",
              "", "## Modelo y método", "",
              "Se integra el sistema no lineal del modelo simplificado de El Niño:", "",
              "```text", "T' = R(t)*T*(1 - T/Tmax) - alpha*H",
              "H' = beta*T² - gamma*H*P", "P' = delta*H² - sigma*P",
              "R(t) = Roce + A*cos(omega*t - pi)", "```", "",
              "Tiempo en meses; T en °C, H en u.h. y P en u.p. Condiciones iniciales: (T,H,P)=(1,0,0). Intervalo: [0,60] meses.",
              "Parámetros nominales: Tmax=4; Roce=0.30; alpha=0.05; beta=0.20; gamma=0.10; delta=0.40; sigma=0.70; A=0.45; omega=2*pi/12. Las unidades de los coeficientes se detallan en modelo.py.",
              "Se aplica Euler explícito, y[n+1]=y[n]+h*f(t[n],y[n]), actualizando las tres componentes simultáneamente. Se mantiene el mismo método para las referencias y todos los escenarios.",
              "", "## Parte 1: convergencia y estabilidad", "",
              f"- Referencia Euler: h = {p1['h_ref']:.10g} meses.",
              f"- Verificación con tres referencias h, 2h y 4h: orden observado = {p1['orden_ref']:.6f}. La diferencia entre referencias es compatible con una estimación de error de primer orden.",
              f"- Diferencia máxima de T al refinar la referencia: {p1['diferencia_ref']:.6g} °C; equivale al {100*p1['razon_ref']:.3f}% del menor error analizado (objetivo: <1%).",
              f"- Pendiente log-log de las primeras cuatro mallas: {p1['pendiente']:.5f}; compatible con el orden 1 de Euler.", "",
              "| h [meses] | Error máximo T [°C] | Orden observado | Estado |",
              "|---:|---:|---:|---|"]
    for f in p1["filas"]:
        err = f"{f['error_T']:.7g}" if np.isfinite(f["error_T"]) else "desbordamiento"
        orden = f"{f['orden']:.5f}" if np.isfinite(f["orden"]) else "—"
        lineas.append(f"| {f['h']:.8g} | {err} | {orden} | {f['estado']} |")
    lineas += ["", "El orden tabulado usa el error de la fila actual dividido por el de la fila anterior y logaritmo en base 2. Solo representa orden de convergencia en el rango de pasos pequeños.",
               "", f"En el barrido auxiliar, el primer paso que produce desbordamiento es h = {p1['primer_desborde']:.3f} meses. Es un resultado de esa malla de exploración, no un umbral exacto de estabilidad.",
               "Antes del desbordamiento aparecen oscilaciones artificiales en P y una desviación del cociente E/h respecto del régimen fino. El acoplamiento no lineal puede hacer que la pérdida de estabilidad se vea primero en H o P, aunque el error exigido por el enunciado se mida en T. Aumentar el error por usar un paso grande y producir una trayectoria espuria son fenómenos distintos; los gráficos muestran ambos.",
               "En las trayectorias finitas del barrido original de 0.025 meses, E_T/h permanece aproximadamente entre 0.46 y 0.51. Sin embargo, ese muestreo saltea pasos con crecimiento explosivo todavía finito. Se agregó un barrido de resolución 0.0001 meses en [1.15,1.30] para cubrir la transición solicitada en 1(d).",
               f"El primer punto de ese barrido fino con E_T/h mayor que diez veces el valor de malla fina es h={p1['primer_salto']['h']:.4f}: E_T={p1['primer_salto']['error_T']:.6g} °C. El factor diez solo identifica una desviación marcada, no define un límite teórico de estabilidad. Se observan tanto trayectorias explosivas finitas como desbordamientos intercalados; no se debe interpolar entre unos y otros para declarar un umbral exacto.",
               "La clasificación a horizonte finito no es monótona ni equivale a estabilidad asintótica. Al refinar el muestreo del paso cambió la evidencia de la transición, no las ecuaciones ni los resultados de las mallas finas.",
               "Para interpretar las oscilaciones también se linealiza localmente: e[n+1] ≈ (I+h*J)e[n], con J=df/dy. Para un autovalor disipativo lambda, Euler requiere |1+h*lambda|≤1 en el problema lineal congelado; esto da h≤-2*Re(lambda)/|lambda|².",
               f"Al muestrear el último año nominal, el menor límite local de los modos disipativos es {p1['h_limite_local']:.6g} meses. Esto ayuda a explicar las oscilaciones observadas para pasos del orden de un mes. El sistema es no autónomo y no lineal: este diagnóstico local no es una prueba de estabilidad global ni predice exactamente el primer desbordamiento.",
               "Los pasos 6 y 12 producen trayectorias explosivas pero terminan con números finitos porque realizan muy pocas iteraciones. No son una recuperación de estabilidad. Se excluyen del gráfico log-log para que su escala no oculte la convergencia, y no se interpreta un orden para ellos.",
               "", "![Convergencia](../figuras/01_convergencia.png)", "", "![Estabilidad](../figuras/02_estabilidad.png)",
               "", "![Barrido fino de estabilidad](../figuras/08_estabilidad_detalle.png)",
               "", "## Parte 2: paso e indicadores", "",
               f"Paso adoptado: h = {p2['h']:.10g} meses. En todas las combinaciones exploradas y los barridos individuales, la última comparación h frente a 2h da una diferencia relativa máxima de SP de {100*p2['precision']['max_diferencia_rel_SP']:.4f}% y una diferencia máxima de desfase de {p2['precision']['max_diferencia_desfase']:.6g} meses. Son estimaciones por refinamiento, no cotas rigurosas del error.",
               f"Se verifica también Pmax, Pmin y Pmed: su mayor diferencia relativa es {100*p2['precision']['max_diferencia_rel_extremos_promedio_P']:.4f}% (tolerancia adoptada: 0.5%). Para el caso nominal, la diferencia máxima de T entre el paso adoptado h={p2['h']:.10g} y la referencia Euler de la Parte 1 es {p1['error_parte2']:.6g} °C. El paso está por debajo de todas las mallas del rango convergente estudiado.",
               "", "Los indicadores usan únicamente [48,60]. El promedio es temporal por trapecios. Los picos se detectan en los nodos, de modo que el desfase tiene la resolución del paso.",
               "", "| Escenario | A | beta | gamma | sigma | SP | Desfase [meses] | Pmax | Pmin | Pmed | Residuo anual |",
               "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for nombre, f in [("Nominal", p2["nominal"]), ("Mejor grilla", p2["mejor"]), ("Todos máximos", p2["maximos"])]:
        lineas.append(f"| {nombre} | " + " | ".join(f"{f[k]:.6g}" for k in (*NOMBRES, "SP", "desfase", "Pmax", "Pmin", "Pmed", "residuo_anual")) + " |")
    lineas += ["", "Pmax, Pmin y Pmed se expresan en u.p.; SP y los residuos son adimensionales.",
               "", "![Nominal completo](../figuras/03_nominal_60_meses.png)", "", "![Nominal último ciclo](../figuras/04_nominal_ultimo_ciclo.png)",
               "", "## Sensibilidad individual", "",
               "| Parámetro | SP en extremo inferior → superior | Desfase inferior → superior [meses] |",
               "|---|---:|---:|"]
    for n in NOMBRES:
        fs = [f for f in p2["sensibilidad"] if f["parametro"] == n]
        lineas.append(f"| {n} | {fs[0]['SP']:.6g} → {fs[-1]['SP']:.6g} | {fs[0]['desfase']:.6g} → {fs[-1]['desfase']:.6g} |")
    lineas += ["", "En los puntos muestreados, SP aumenta con A, beta y sigma, y disminuye con gamma. El desfase no aumenta al incrementar ninguno de los cuatro parámetros; algunos valores consecutivos coinciden por la resolución temporal de los picos. Estas tendencias se refieren a los intervalos indicados y al resto de parámetros nominales.",
               "", "- **A:** aumenta el contraste estacional de R(t), sin cambiar su media anual. En los experimentos aumenta SP y disminuye el desfase: la excitación térmica más intensa se propaga a H y P a través de sus fuentes cuadráticas.",
               "- **beta:** aumenta la fuente de humedad beta*T². En este rango aumenta SP y reduce el desfase, aunque también refuerza el enfriamiento indirecto de T a través de H. El resultado surge del sistema acoplado, no solo del signo positivo de beta.",
               "- **gamma:** aumenta el consumo de humedad gamma*H*P. Reduce la humedad disponible para la fuente delta*H². Se observa una menor severidad relativa y un desfase menor; el consumo atenúa el ciclo y acelera el ajuste de humedad cuando P es positiva.",
               "- **sigma:** aumenta la descarga de P y reduce su escala de disipación 1/sigma. En este rango disminuye el desfase y aumenta SP: una disipación más rápida reduce el almacenamiento y permite mayor contraste relativo entre máximo y mínimo. Esto no implica que aumente la precipitación absoluta.",
               "", "![Sensibilidad](../figuras/05_sensibilidad.png)",
               "", "## Búsqueda conjunta y comparación", "",
               f"Se evaluaron {p2['n_grilla']} combinaciones distintas entre la grilla inicial de 5⁴ puntos, una grilla global de 9⁴ puntos y el refinamiento local de 5⁴ puntos. Se eliminan los puntos repetidos. {p2['n_admisibles']} cumplen el criterio operativo de admisibilidad. El mayor residuo anual de la grilla fue {p2['max_residuo']:.6g}; el umbral adoptado es 0.01.",
               "El residuo compara cada variable del penúltimo y último ciclo a igual fase y normaliza por su máxima magnitud en el último ciclo. La comparación indica repetición anual aproximada al final de los 60 meses; no demuestra periodicidad exacta ni acotación para todo tiempo futuro.",
               "", "La mejor combinación hallada no toma todos los parámetros en su máximo: gamma queda en el mínimo. Esto es coherente con que gamma consume humedad, mientras beta la genera. A y sigma intervienen en otros términos, y SP mide amplitud relativa, por lo que maximizar los cuatro coeficientes no equivale a maximizar SP.",
               f"La severidad del mejor escenario es {p2['mejor']['SP']/p2['nominal']['SP']:.3f} veces la nominal. El escenario de todos los máximos tiene SP = {p2['maximos']['SP']:.6g}.",
               "El resultado es el mejor escenario encontrado en las grillas examinadas. No constituye una prueba de máximo global continuo. El refinamiento en h y el refinamiento paramétrico resuelven incertidumbres distintas.",
               "", "![Comparación completa](../figuras/06_comparacion_60_meses.png)", "", "![Comparación último año](../figuras/07_comparacion_ultimo_ciclo.png)",
               "", "## Verificación adicional de los escenarios finales", "",
               "Se recalculan los diez mejores candidatos con h, h/2 y h/4. La misma combinación gana en los tres pasos. Es una comprobación de robustez del ranking entre finalistas; no reemplaza la limitación de la búsqueda por grilla.",
               "También se comparan las trayectorias completas [0,60] de los tres escenarios de la tabla, variable por variable y en nodos comunes. Se exige un cambio menor al 0.5% de la máxima magnitud de cada variable. Los extremos de toda la simulación se guardan para documentar acotación observada, en vez de apoyarse únicamente en la guarda de 10^6.",
               "", "| Escenario | Paso [meses] | SP | Desfase [meses] |", "|---|---:|---:|---:|"]
    for f in p2["validacion"]["filas"]:
        if f["escenario"] in ("nominal", "mejor_grilla", "todos_maximos"):
            lineas.append(f"| {f['escenario']} | {f['h']:.10g} | {f['SP']:.8g} | {f['desfase']:.8g} |")
    lineas += ["", "La tabla principal conserva el paso común utilizado para todos los escenarios; esta tabla muestra su incertidumbre por refinamiento adicional. No se mezclan pasos al comparar los casos.",
               "", "### Refinamiento paramétrico en todo el dominio", "",
               "El refinamiento global verifica regiones interiores que una búsqueda solo local podría pasar por alto. Todas las resoluciones se comparan al mismo paso temporal:", "",
               "| Grilla | Puntos | Mejor SP | A | beta | gamma | sigma |",
               "|---|---:|---:|---:|---:|---:|---:|"]
    for f in p2["niveles_busqueda"]:
        lineas.append(f"| {f['nivel']} | {f['puntos']} | {f['SP']:.7g} | " +
                      " | ".join(f"{f[n]:.7g}" for n in NOMBRES) + " |")
    lineas += ["", "También se comparan vecinos de la grilla global en cada dirección, manteniendo fijos los otros tres parámetros. Se registra cuántas comparaciones invierten el signo observado en los barridos individuales:",
               "", "| Parámetro | Signo ensayado de variación de SP | Comparaciones | Inversiones |",
               "|---|---:|---:|---:|"]
    for f in p2["tendencias_conjuntas"]:
        lineas.append(f"| {f['parametro']} | {f['signo_esperado']:+d} | {f['comparaciones']} | {f['inversiones']} |")
    if not any(f["inversiones"] for f in p2["tendencias_conjuntas"]):
        lineas += ["", "No se detectan inversiones en esta grilla: SP aumenta con A, beta y sigma y disminuye con gamma en las comparaciones conjuntas muestreadas. Esto respalda la esquina encontrada como máximo de la grilla global; no demuestra monotonicidad entre nodos ni un máximo global continuo."]
    else:
        lineas += ["", "Las inversiones muestran que las tendencias individuales no pueden generalizarse a todo el dominio conjunto. La selección se basa en evaluar las combinaciones, no en suponer monotonicidad."]
    lineas += ["", "Los mapas siguientes son cortes con los otros dos parámetros fijos según sus títulos; no son proyecciones ni una representación completa de las cuatro dimensiones.",
               "", "![Búsqueda conjunta](../figuras/09_busqueda_conjunta.png)",
               "", "## Archivos de datos", "",
               "Los CSV conservan los valores sin redondear del análisis: convergencia, estabilidad, verificación del paso, sensibilidad, grilla, comparación y trayectorias de los tres escenarios. Los parámetros no barridos mantienen sus valores nominales. No se utilizaron datos climáticos observacionales."]
    (RESULTADOS / "resumen.md").write_text("\n".join(lineas) + "\n", encoding="utf-8")


def main():
    RESULTADOS.mkdir(exist_ok=True)
    FIGURAS.mkdir(exist_ok=True)
    plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.25, "font.size": 10})
    p1 = parte1()
    p2 = parte2()
    t_ref, y_ref = p1.pop("referencia")
    t_nominal, y_nominal = p2.pop("trayectoria_nominal")
    salto = round(p2["h"] / p1["h_ref"])
    if salto < 1 or not np.isclose(salto * p1["h_ref"], p2["h"], rtol=0, atol=1e-12):
        raise RuntimeError("El paso elegido requiere revisar la malla de referencia.")
    np.testing.assert_allclose(t_ref[::salto], t_nominal, rtol=0, atol=1e-10)
    errores = np.abs(y_nominal[:, 0] - y_ref[::salto, 0])
    p1["error_parte2"] = float(errores.max())
    guardar_csv("paso_parte2_vs_referencia.csv", [{
        "h_parte2": p2["h"], "h_referencia": p1["h_ref"],
        "max_diferencia_T": p1["error_parte2"],
        "t_max_diferencia": float(t_nominal[np.argmax(errores)]),
    }])
    resumen(p1, p2)
    print(f"Terminado. Resultados: {RESULTADOS / 'resumen.md'}", flush=True)
    print(f"Paso Parte 2: {p2['h']}; SP nominal: {p2['nominal']['SP']:.6f}; SP mejor: {p2['mejor']['SP']:.6f}", flush=True)


if __name__ == "__main__":
    main()

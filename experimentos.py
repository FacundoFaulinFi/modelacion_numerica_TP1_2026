"""Indicadores del último año e integración de escenarios independientes."""

from dataclasses import asdict
from types import SimpleNamespace

import numpy as np

from modelo import CONDICIONES_INICIALES, Parametros, derivadas


def indicadores(t, y) -> dict:
    """Calcula indicadores en [48,60] y el residuo anual frente a [36,48].

    Requiere una malla uniforme que contenga 48 y 60. El desfase
    usa el único pico interior de P, siempre que sea el máximo global de P
    y posterior al máximo global de T. Con varios picos, temperatura constante,
    máximo de P en un extremo o sin pico posterior devuelve NaN y posterior=False.
    Si también se dispone de [36,48], calcula el residuo anual; en caso
    contrario ese diagnóstico es NaN, sin impedir calcular los indicadores.
    """
    t, y = np.asarray(t, dtype=float), np.asarray(y, dtype=float)
    if t.ndim != 1 or y.shape != (len(t), 3) or len(t) < 3:
        raise ValueError("Se espera una trayectoria (len(t), 3).")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(t)):
        raise ValueError("La trayectoria contiene valores no finitos.")
    h = t[1] - t[0]
    if h <= 0 or not np.allclose(np.diff(t), h, rtol=1e-9, atol=1e-12):
        raise ValueError("Los indicadores requieren una malla uniforme creciente.")
    indices = [int(np.argmin(np.abs(t - limite))) for limite in (48, 60)]
    if not np.allclose(t[indices], [48, 60], rtol=0, atol=1e-9):
        raise ValueError("La malla debe contener los meses 48 y 60.")
    i48, i60 = indices
    u, z = t[i48:i60 + 1], y[i48:i60 + 1]
    P = z[:, 2]
    pmed = float(np.trapezoid(P, u) / 12.0)
    pmax, pmin = float(P.max()), float(P.min())
    it = int(np.argmax(z[:, 0]))
    # Extremos interiores: no se confunde un extremo del intervalo con un pico.
    picos = np.flatnonzero((P[1:-1] > P[:-2]) & (P[1:-1] >= P[2:])) + 1
    # No elige arbitrariamente entre varios picos ni usa un pico secundario
    # cuando el máximo global del intervalo está en uno de sus extremos.
    ip = None
    if len(picos) == 1 and np.ptp(z[:, 0]) > 0:
        candidato = int(picos[0])
        if candidato > it and candidato == int(np.argmax(P)):
            ip = candidato
    residuos = np.full(3, np.nan)
    i36 = int(np.argmin(np.abs(t - 36)))
    if np.isclose(t[i36], 36, rtol=0, atol=1e-9):
        anterior = y[i36:i48 + 1]
        if anterior.shape != z.shape:
            raise ValueError("Los dos ciclos deben tener la misma cantidad de nodos.")
        escala = np.maximum(np.max(np.abs(z), axis=0), 1e-12)
        residuos = np.max(np.abs(z - anterior), axis=0) / escala
    return {
        "SP": (pmax - pmin) / pmed if pmed > 0 else float("nan"),
        "desfase": float(u[ip] - u[it]) if ip is not None else float("nan"),
        "Pmax": pmax, "Pmin": pmin, "Pmed": pmed,
        "t_Tmax": float(u[it]),
        "t_Pmax_posterior": float(u[ip]) if ip is not None else float("nan"),
        "posterior": ip is not None,
        "cantidad_picos_P": int(len(picos)),
        "pico_posterior_es_Pmax": bool(ip is not None and ip == int(np.argmax(P))),
        "residuo_T": float(residuos[0]),
        "residuo_H": float(residuos[1]),
        "residuo_P": float(residuos[2]),
        "residuo_anual": float(residuos.max()),
        "minimo_ultimo_ciclo": float(z.min()),
    }


def euler_lote(parametros: list[Parametros], h: float):
    """Euler para varios escenarios; guarda [36,60] para ahorrar memoria.

    Devuelve t, historial y una máscara de trayectorias válidas. Solo admite
    pasos que dividan 12 meses. La guarda 1e6 detecta crecimiento explosivo;
    no modifica ni rescata las trayectorias que la superan.
    """
    if not np.isfinite(h) or h <= 0 or not parametros:
        raise ValueError("Se requieren h positivo y al menos un escenario.")
    cociente = 12.0 / float(h)
    if not np.isfinite(cociente):
        raise ValueError("La relación 12/h no es representable: h es demasiado pequeño.")
    n_anual = round(cociente)
    if n_anual < 1 or not np.isclose(n_anual * h, 12, rtol=0, atol=1e-10):
        raise ValueError("h debe dividir 12 meses.")
    valores = [asdict(p) for p in parametros]
    p_lote = SimpleNamespace(**{
        nombre: np.array([p[nombre] for p in valores]) for nombre in valores[0]
    })
    estado = np.tile(np.array(CONDICIONES_INICIALES)[:, None], (1, len(parametros)))
    validos = np.ones(len(parametros), dtype=bool)
    historial = np.full((2 * n_anual + 1, len(parametros), 3), np.nan)
    inicio = 3 * n_anual
    for n in range(5 * n_anual + 1):
        if n >= inicio:
            historial[n - inicio] = np.where(validos[:, None], estado.T, np.nan)
        if n == 5 * n_anual:
            break
        with np.errstate(over="ignore", invalid="ignore"):
            siguiente = estado + h * derivadas(n * h, estado, p_lote)
        validos &= np.all(np.isfinite(siguiente) & (np.abs(siguiente) < 1e6), axis=0)
        # Las trayectorias inválidas permanecen marcadas y se guardan como NaN.
        estado = np.where(validos[None, :], siguiente, 0.0)
    t = np.arange(inicio, 5 * n_anual + 1) * h
    return t, historial, validos


def evaluar_lote(parametros: list[Parametros], h: float, tamano_lote: int = 256) -> list[dict]:
    """Evalúa escenarios en bloques para limitar la memoria de las grillas.

    El tamaño de bloque no modifica el método ni el orden de las filas.
    """
    if not parametros or not isinstance(tamano_lote, int) or tamano_lote <= 0:
        raise ValueError("Se requieren escenarios y un tamaño de lote entero positivo.")
    filas = []
    for inicio in range(0, len(parametros), tamano_lote):
        bloque = parametros[inicio:inicio + tamano_lote]
        t, y, validos = euler_lote(bloque, h)
        for j, p in enumerate(bloque):
            fila = {nombre: getattr(p, nombre) for nombre in ("A", "beta", "gamma", "sigma")}
            fila["finito_acotado_60m"] = bool(validos[j])
            fila["admisible"] = False
            if validos[j]:
                fila.update(indicadores(t, y[:, j]))
                fila["admisible"] = bool(
                    np.isfinite(fila["residuo_anual"]) and fila["residuo_anual"] < 0.01
                    and fila["Pmed"] > 0 and fila["minimo_ultimo_ciclo"] >= 0
                    and fila["posterior"]
                )
            filas.append(fila)
        # No mantiene vivo el bloque anterior mientras se construye el siguiente.
        del y
    return filas

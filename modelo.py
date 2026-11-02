"""Modelo simplificado de El Niño: ecuaciones y método de Euler explícito.

El tiempo y el paso de integración se expresan en meses. Las columnas de
la solución son, en orden, temperatura T [°C], humedad H [u.h.] y
precipitación P [u.p.]. Requiere NumPy.

Ejemplo de uso desde un notebook o un script::

    from dataclasses import replace
    from modelo import PARAMETROS_NOMINALES, euler

    t, y = euler(h=12 / 4096)
    temperatura, humedad, precipitacion = y.T
    parametros = replace(PARAMETROS_NOMINALES, A=0.60)
    t, y = euler(h=12 / 4096, parametros=parametros)

El paso del ejemplo es el validado para este TP; ver decisiones.md.
"""

from dataclasses import dataclass, fields
import math

import numpy as np


@dataclass(frozen=True)
class Parametros:
    """Coeficientes nominales del enunciado, con el mes como unidad temporal.

    Tmax: temperatura máxima sostenible [°C].
    Roce: tasa intrínseca de calentamiento [mes⁻¹].
    alpha: enfriamiento por evaporación [°C mes⁻¹ u.h.⁻¹].
    beta: eficiencia de evaporación [u.h. mes⁻¹ °C⁻²].
    gamma: consumo de humedad por precipitación [mes⁻¹ u.p.⁻¹].
    delta: conversión de humedad a lluvia [u.p. mes⁻¹ u.h.⁻²].
    sigma: disipación de precipitación [mes⁻¹].
    A: amplitud estacional [mes⁻¹].
    omega: frecuencia angular [rad/mes].
    """

    Tmax: float = 4.0
    Roce: float = 0.30
    alpha: float = 0.05
    beta: float = 0.20
    gamma: float = 0.10
    delta: float = 0.40
    sigma: float = 0.70
    A: float = 0.45
    omega: float = 2.0 * math.pi / 12.0

    def __post_init__(self):
        for campo in fields(self):
            if not math.isfinite(getattr(self, campo.name)):
                raise ValueError(f"El parámetro {campo.name} debe ser finito.")
        if self.Tmax <= 0:
            raise ValueError("Tmax debe ser positivo.")


PARAMETROS_NOMINALES = Parametros()
CONDICIONES_INICIALES = (1.0, 0.0, 0.0)
TIEMPO_FINAL = 60.0


def tasa_estacional(t, parametros: Parametros = PARAMETROS_NOMINALES):
    """Evalúa R(t) = Roce + A cos(omega t - pi), para un tiempo o un vector."""
    return parametros.Roce + parametros.A * np.cos(
        parametros.omega * np.asarray(t) - np.pi
    )


def derivadas(
    t: float, y, parametros: Parametros = PARAMETROS_NOMINALES
) -> np.ndarray:
    """Devuelve [dT/dt, dH/dt, dP/dt] evaluado en el estado y = [T, H, P]."""
    T, H, P = np.asarray(y, dtype=float)
    p = parametros
    dT = tasa_estacional(t, p) * T * (1.0 - T / p.Tmax) - p.alpha * H
    dH = p.beta * T**2 - p.gamma * H * P
    dP = p.delta * H**2 - p.sigma * P
    return np.array([dT, dH, dP], dtype=float)


def jacobiano(t: float, y, parametros: Parametros = PARAMETROS_NOMINALES) -> np.ndarray:
    """Matriz df/dy para estudiar perturbaciones locales de una trayectoria."""
    T, H, P = np.asarray(y, dtype=float)
    p = parametros
    return np.array([
        [tasa_estacional(t, p) * (1 - 2 * T / p.Tmax), -p.alpha, 0],
        [2 * p.beta * T, -p.gamma * P, -p.gamma * H],
        [0, 2 * p.delta * H, -p.sigma],
    ], dtype=float)


def euler(
    h: float,
    parametros: Parametros = PARAMETROS_NOMINALES,
    y0=CONDICIONES_INICIALES,
    t_final: float = TIEMPO_FINAL,
) -> tuple[np.ndarray, np.ndarray]:
    """Integra desde t=0 con Euler explícito: y[n+1] = y[n] + dt*f(t[n], y[n]).

    Devuelve (t, y), con t de forma (N+1,) e y de forma (N+1, 3).
    El paso es h; si no divide t_final, el último paso se acorta para
    terminar exactamente en t_final. Para comparar errores y órdenes,
    conviene elegir pasos que dividan el intervalo y formen mallas anidadas.

    Las tres componentes se actualizan simultáneamente. No se recortan
    valores negativos ni se corrige la solución: eso alteraría el método
    cuyo comportamiento se quiere estudiar. Si aparecen infinitos, NaN
    o un desbordamiento, se lanza FloatingPointError. Una trayectoria
    finita por sí sola no garantiza estabilidad ni precisión.
    """
    if not math.isfinite(h) or h <= 0:
        raise ValueError("h debe ser un número finito y positivo.")
    if not math.isfinite(t_final) or t_final < 0:
        raise ValueError("t_final debe ser finito y no negativo.")

    estado_inicial = np.asarray(y0, dtype=float)
    if estado_inicial.shape != (3,) or not np.all(np.isfinite(estado_inicial)):
        raise ValueError("y0 debe contener tres valores finitos: [T0, H0, P0].")

    # Evita un paso adicional diminuto por redondeo de t_final / h.
    cociente = float(t_final) / float(h)
    if not math.isfinite(cociente):
        raise ValueError("La relación t_final/h no es representable: h es demasiado pequeño.")
    entero_cercano = round(cociente)
    if math.isclose(cociente, entero_cercano, rel_tol=1e-12, abs_tol=0.0):
        n_pasos = entero_cercano
    else:
        n_pasos = math.ceil(cociente)
    # Si t_final/h subdesborda a cero, el intervalo sigue siendo positivo:
    # debe conservarse t=0 y realizarse un paso final acortado.
    if t_final > 0:
        n_pasos = max(1, n_pasos)

    t = np.arange(n_pasos + 1, dtype=float) * h
    t[-1] = t_final
    y = np.empty((n_pasos + 1, 3), dtype=float)
    y[0] = estado_inicial

    for n in range(n_pasos):
        dt = t[n + 1] - t[n]
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                y[n + 1] = y[n] + dt * derivadas(t[n], y[n], parametros)
            if not np.all(np.isfinite(y[n + 1])):
                raise FloatingPointError("Estado no finito.")
        except FloatingPointError as exc:
            raise FloatingPointError(
                f"Euler produjo un desbordamiento o un valor no finito "
                f"en el paso {n + 1}, t={t[n + 1]:.6g}, h={h:.6g}."
            ) from exc

    return t, y


if __name__ == "__main__":
    paso = 12 / 4096
    tiempos, solucion = euler(h=paso)
    print(f"Simulación nominal con Euler explícito (h = {paso:g} meses).")
    print("Justificación del paso y verificaciones: decisiones.md y resultados/resumen.md.")
    print(f"Intervalo: {tiempos[0]:g} a {tiempos[-1]:g} meses.")
    print(f"Cantidad de pasos: {len(tiempos) - 1}.")
    print("Estado final [T, H, P]:", solucion[-1])

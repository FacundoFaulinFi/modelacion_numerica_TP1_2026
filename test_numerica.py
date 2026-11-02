"""Comprobaciones de la fórmula de Euler, los indicadores y la integración por lotes."""

from dataclasses import replace
import unittest

import numpy as np

from experimentos import euler_lote, evaluar_lote, indicadores
from modelo import CONDICIONES_INICIALES, PARAMETROS_NOMINALES, derivadas, euler, jacobiano


class ComprobacionesNumericas(unittest.TestCase):
    def test_primer_paso_simultaneo(self):
        np.testing.assert_allclose(derivadas(0, CONDICIONES_INICIALES), [-0.1125, 0.2, 0])
        _, y = euler(0.1, t_final=0.1)
        np.testing.assert_allclose(y[1], [0.98875, 0.02, 0])

    def test_acoplamientos_y_fase_estacional(self):
        # Estado no nulo: comprueba HP, H², T² y el enfriamiento por H.
        np.testing.assert_allclose(derivadas(6, [2, 3, 4]), [0.6, -0.4, 0.8], atol=1e-14)
        np.testing.assert_allclose(derivadas(18, [2, 3, 4]), derivadas(6, [2, 3, 4]))
        self.assertAlmostEqual(derivadas(0, [2, 3, 4])[0], -0.3)

    def test_jacobiano_por_diferencias_centradas(self):
        y = np.array([2.1, 3.2, 4.3])
        eps = 1e-6
        aproximado = np.column_stack([
            (derivadas(5.3, y + eps*e) - derivadas(5.3, y - eps*e)) / (2*eps)
            for e in np.eye(3)
        ])
        np.testing.assert_allclose(jacobiano(5.3, y), aproximado, rtol=1e-7, atol=1e-9)

    def test_paso_final_y_intervalo_nulo(self):
        t, _ = euler(0.3, t_final=1)
        np.testing.assert_allclose(t, [0, 0.3, 0.6, 0.9, 1])
        t, y = euler(0.1, t_final=0)
        self.assertEqual(len(t), 1)
        np.testing.assert_array_equal(y[0], CONDICIONES_INICIALES)

    def test_entradas_invalidas(self):
        for h in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                euler(h)
        with self.assertRaises(ValueError):
            euler(0.1, y0=[1, 0])
        with self.assertRaises(ValueError):
            euler_lote([PARAMETROS_NOMINALES], 0.7)

    def test_lote_coincide_con_integracion_individual(self):
        parametros = [PARAMETROS_NOMINALES, replace(PARAMETROS_NOMINALES, A=0.6, gamma=0.05)]
        h = 12 / 512
        t_lote, y_lote, validos = euler_lote(parametros, h)
        self.assertTrue(validos.all())
        for j, p in enumerate(parametros):
            t, y = euler(h, parametros=p)
            mask = t >= 36
            np.testing.assert_allclose(t[mask], t_lote, atol=1e-12)
            np.testing.assert_allclose(y[mask], y_lote[:, j], rtol=1e-11, atol=1e-11)

    def test_indicadores_senal_periodica_con_desfase_conocido(self):
        t = np.linspace(36, 60, 2401)
        T = 2 + np.cos(2 * np.pi * (t - 51) / 12)
        H = 4 + 0.5 * np.cos(2 * np.pi * t / 12)
        P = 3 + np.cos(2 * np.pi * (t - 53) / 12)
        r = indicadores(t, np.column_stack([T, H, P]))
        self.assertAlmostEqual(r["SP"], 2 / 3, places=12)
        self.assertAlmostEqual(r["Pmed"], 3, places=12)
        self.assertAlmostEqual(r["desfase"], 2, places=12)
        self.assertLess(r["residuo_anual"], 1e-12)

    def test_no_inventa_pico_posterior(self):
        t = np.linspace(36, 60, 2401)
        T = 2 + np.cos(2 * np.pi * (t - 57) / 12)
        P = 3 + np.cos(2 * np.pi * (t - 51) / 12)
        r = indicadores(t, np.column_stack([T, np.ones_like(t), P]))
        self.assertFalse(r["posterior"])
        self.assertTrue(np.isnan(r["desfase"]))

    def test_no_informa_desfase_ambiguo_con_varios_picos(self):
        t = np.linspace(48, 60, 1201)
        T = 2 + np.cos(2*np.pi*(t - 51)/12)
        P = 3 + np.cos(4*np.pi*(t - 49)/12)
        r = indicadores(t, np.column_stack([T, np.ones_like(t), P]))
        self.assertEqual(r["cantidad_picos_P"], 2)
        self.assertFalse(r["posterior"])
        self.assertTrue(np.isnan(r["desfase"]))

    def test_temperatura_constante_no_define_desfase(self):
        t = np.linspace(48, 60, 1201)
        P = 3 + np.cos(2*np.pi*(t - 53)/12)
        r = indicadores(t, np.column_stack([np.ones_like(t), np.ones_like(t), P]))
        self.assertAlmostEqual(r["SP"], 2/3, places=12)
        self.assertTrue(np.isnan(r["desfase"]))

    def test_no_confunde_pico_secundario_con_Pmax(self):
        t = np.linspace(48, 60, 1201)
        T = 2 + np.cos(2*np.pi*(t - 51)/12)
        P = 10 - 0.2*(t - 48) + np.exp(-((t - 55)/0.4)**2)
        r = indicadores(t, np.column_stack([T, np.ones_like(t), P]))
        self.assertEqual(r["cantidad_picos_P"], 1)
        self.assertAlmostEqual(r["Pmax"], P[0])
        self.assertTrue(np.isnan(r["desfase"]))

    def test_intervalo_positivo_no_desaparece_por_subdesbordamiento(self):
        t, y = euler(1e308, t_final=1e-308)
        np.testing.assert_array_equal(t, [0, 1e-308])
        np.testing.assert_array_equal(y[0], CONDICIONES_INICIALES)

    def test_paso_demasiado_pequeno_rechazado_sin_overflowerror(self):
        diminuto = float(np.nextafter(0, 1))
        with self.assertRaises(ValueError):
            euler(diminuto)
        with self.assertRaises(ValueError):
            euler_lote([PARAMETROS_NOMINALES], diminuto)

    def test_indicadores_excluyen_transitorio(self):
        t = np.linspace(0, 60, 6001)
        T = 2 + np.cos(2 * np.pi * (t - 51) / 12)
        P = 3 + np.cos(2 * np.pi * (t - 53) / 12)
        P[t < 48] += 100
        r = indicadores(t, np.column_stack([T, np.ones_like(t), P]))
        self.assertAlmostEqual(r["SP"], 2 / 3, places=12)
        self.assertAlmostEqual(r["Pmed"], 3, places=12)
        self.assertAlmostEqual(r["desfase"], 2, places=12)
        self.assertGreater(r["residuo_anual"], 1)

    def test_indicadores_rechazan_ciclo_incompleto(self):
        t = np.linspace(49, 60, 1101)
        with self.assertRaises(ValueError):
            indicadores(t, np.ones((len(t), 3)))

    def test_indicadores_solo_necesitan_ultimo_ciclo(self):
        t = np.linspace(48, 60, 1201)
        T = 2 + np.cos(2 * np.pi * (t - 51) / 12)
        P = 3 + np.cos(2 * np.pi * (t - 53) / 12)
        r = indicadores(t, np.column_stack([T, np.ones_like(t), P]))
        self.assertAlmostEqual(r["SP"], 2/3, places=12)
        self.assertAlmostEqual(r["desfase"], 2, places=12)
        self.assertTrue(np.isnan(r["residuo_anual"]))

    def test_particionar_lotes_no_cambia_resultados_ni_orden(self):
        ps = [PARAMETROS_NOMINALES, replace(PARAMETROS_NOMINALES, A=0.6),
              replace(PARAMETROS_NOMINALES, gamma=0.05)]
        sin_partir = evaluar_lote(ps, 12/256, tamano_lote=3)
        partido = evaluar_lote(ps, 12/256, tamano_lote=2)
        self.assertEqual(sin_partir, partido)
        with self.assertRaises(ValueError):
            evaluar_lote(ps, 12/256, tamano_lote=0)

    def test_euler_con_solucion_analitica_logistica_estacional(self):
        # alpha=beta=0 desacopla T de H,P: la integral de R es conocida.
        p = replace(PARAMETROS_NOMINALES, alpha=0, beta=0)
        errores = []
        for h in (0.04, 0.02, 0.01):
            t, y = euler(h, parametros=p, t_final=12)
            integral_R = p.Roce*t + p.A/p.omega * (
                np.sin(p.omega*t - np.pi) - np.sin(-np.pi))
            exacta = p.Tmax / (1 + (p.Tmax - 1)*np.exp(-integral_R))
            errores.append(float(np.max(np.abs(y[:, 0] - exacta))))
            np.testing.assert_array_equal(y[:, 1:], 0)
        self.assertLess(errores[-1], 0.01)
        self.assertAlmostEqual(np.log2(errores[0] / errores[1]), 1, delta=0.1)
        self.assertAlmostEqual(np.log2(errores[1] / errores[2]), 1, delta=0.1)

    def test_lote_no_recupera_trayectoria_desbordada(self):
        _, historial, validos = euler_lote([PARAMETROS_NOMINALES], 1.5)
        self.assertFalse(validos[0])
        self.assertTrue(np.isnan(historial[-1]).all())
        fila = evaluar_lote([PARAMETROS_NOMINALES], 1.5)[0]
        self.assertFalse(fila["admisible"])

    def test_orden_uno_por_autoconvergencia(self):
        _, gruesa = euler(0.04)
        _, media = euler(0.02)
        _, fina = euler(0.01)
        d1 = np.max(np.abs(gruesa[:, 0] - media[::2, 0]))
        d2 = np.max(np.abs(media[:, 0] - fina[::2, 0]))
        self.assertAlmostEqual(np.log2(d1 / d2), 1, delta=0.05)

    def test_no_recorta_negativos_y_detecta_desbordamiento(self):
        _, y = euler(20, t_final=20)
        np.testing.assert_allclose(y[1], [-1.25, 4, 0])
        with self.assertRaises(FloatingPointError):
            euler(1.5)


if __name__ == "__main__":
    unittest.main()

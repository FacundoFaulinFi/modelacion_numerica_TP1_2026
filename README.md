# Trabajo práctico de Modelación Numérica

Sistema no lineal de tres EDO del enunciado `TP 2c26.pdf`, resuelto con Euler explícito.

## Organización

- `modelo.py`: parámetros, ecuaciones y Euler para una trayectoria.
- `experimentos.py`: indicadores del último año y Euler para varios escenarios independientes.
- `analisis.py`: convergencia, estabilidad, sensibilidad, grillas y generación de resultados.
- `decisiones.md`: elecciones metodológicas, criterios y limitaciones.
- `revision.md`: contraste apartado por apartado con el enunciado y mejoras de la revisión.
- `resultados/resumen.md`: resultados y discusión, con los gráficos integrados.
- `resultados/*.csv`: datos numéricos, incluidos los escenarios explorados.
- `figuras/*.png`: gráficos para el informe.
- `test_numerica.py`: comprobaciones numéricas independientes.
- `verificar_resultados.py`: auditoría rápida de los resultados existentes y reproducción de los tres escenarios finales.

El análisis también calcula el Jacobiano para interpretar localmente la estabilidad, refina el barrido de pasos cerca de la transición y verifica los diez mejores candidatos con dos reducciones adicionales del paso temporal. La búsqueda incluye grillas globales de 5⁴ y 9⁴ puntos más refinamiento local; se evalúa en bloques para limitar la memoria.

## Ejecutar en Windows

Desde esta carpeta, con Python instalado:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe analisis.py
.\.venv\Scripts\python.exe -m unittest -v test_numerica
.\.venv\Scripts\python.exe verificar_resultados.py
```

Si ya existe el entorno `.venv` preparado durante este trabajo, basta con las últimas tres instrucciones. No hace falta activar el entorno. En VS Code, seleccionar `.venv/Scripts/python.exe` como intérprete.

El análisis no requiere Internet después de instalar las dependencias. Recalcula las referencias finas, ejecuta los barridos y regenera los CSV, gráficos y resumen; puede demorar algunos minutos. Los archivos generados con los mismos nombres se reemplazan.

Para revisar salidas ya generadas, `verificar_resultados.py` tarda unos segundos y no repite las grillas ni las referencias finas. Comprueba la coherencia numérica de los CSV y reproduce los tres escenarios finales. Esta auditoría complementa las pruebas unitarias; no las sustituye.

Para reproducir las versiones exactas utilizadas, instalar `requirements-lock.txt` en lugar de `requirements.txt`.

## Usar el modelo por separado

```python
from dataclasses import replace
from modelo import PARAMETROS_NOMINALES, euler

p = replace(PARAMETROS_NOMINALES, A=0.60, gamma=0.05)
t, y = euler(h=12 / 4096, parametros=p)
T, H, P = y.T
```

El paso definitivo y sus verificaciones están documentados en los resultados. Los indicadores se calculan exclusivamente entre los meses 48 y 60; el penúltimo año solo se usa para comprobar repetición anual.

```python
from experimentos import indicadores

ultimo = t >= 48
medidas = indicadores(t[ultimo], y[ultimo])
print(medidas["SP"], medidas["desfase"])
```

En ese ejemplo el residuo anual es NaN porque no se pasó el penúltimo ciclo. Para medirlo, pasar la trayectoria completa: `indicadores(t, y)`. Las pruebas incluyen un caso desacoplado con solución analítica, además de las verificaciones de Euler y de indicadores.

El desfase también es NaN si no puede identificarse un único pico de P posterior al máximo de una T variable, o si ese pico no coincide con Pmax del último año. Los casos del TP cumplen esas condiciones; para señales más complejas debe explicitarse otra convención antes de informar un desfase.

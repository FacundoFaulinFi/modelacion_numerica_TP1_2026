# Revisión del desarrollo

Se contrastaron las tres páginas del enunciado con las ecuaciones, el código, las tablas, los gráficos y las decisiones metodológicas. Se corrigió una omisión relevante del análisis de estabilidad y se reforzaron las verificaciones de precisión. No se cambiaron los parámetros nominales ni las ecuaciones del problema.

## Resultado de la tercera revisión

No se detectaron nuevos apartados del enunciado sin desarrollar. Se encontraron y corrigieron dos problemas en casos límite y se agregó una auditoría de coherencia entre los resultados guardados y el código.

### Desfase indefinido en señales ambiguas

La función podía devolver un desfase ante dos picos interiores de P eligiendo uno de ellos, o emparejar Tmax con un pico secundario aunque Pmax estuviera en un extremo. También podía dar un número ante una temperatura constante, cuyo instante máximo no identifica un evento térmico.

Ahora devuelve NaN y `posterior=False` en esos casos. Se conservan SP y los demás indicadores que sí están definidos. La convención requiere un único pico interior de P, coincidente con Pmax y posterior al máximo de una T variable. Es una elección conservadora documentada, no una condición física universal.

Los 7170 escenarios del TP ya tenían un único pico posterior coincidente con Pmax, por lo que esta corrección no modifica sus resultados.

### Representación de intervalos temporales extremos

Se reprodujo un error del integrador con h=10^308 y t_final=10^-308: el cociente subdesbordaba a cero y se devolvía solo el instante final, perdiendo t=0. Se corrigió para realizar al menos un paso cuando el intervalo es positivo. También se rechazan con ValueError las relaciones tiempo/paso no representables, en lugar de fallar indirectamente al redondear un infinito.

Son correcciones de robustez del integrador general; los pasos del TP no se encuentran en esos extremos.

### Auditoría independiente de los archivos existentes

Se agregó `verificar_resultados.py`. La ejecución fue satisfactoria y comprobó:

- Duplicación de los pasos, cocientes E/h y fórmula de los órdenes observados.
- Unicidad de 7170 combinaciones, rangos de parámetros, indicadores y cobertura de la grilla global de 9⁴ puntos.
- Reconstrucción de los signos y conteos de las tendencias conjuntas desde los datos guardados.
- Cobertura de los barridos individuales y conservación de los restantes parámetros nominales.
- Tolerancias del paso y coincidencia del ganador entre la grilla, la comparación y los tres refinamientos de finalistas.
- Reintegración de los escenarios nominal, mejor y todos máximos; comparación con sus trayectorias CSV; indicadores del último ciclo; mínimos y máximos de toda la simulación a h/4.
- Existencia de los enlaces locales de la documentación y el informe.

No se repitió la integración completa de las grillas ni las referencias finas: sus definiciones y los cálculos sobre entradas normales no cambiaron. Se verificaron directamente los datos guardados y se reprodujeron los tres escenarios finales con el código corregido.

La suite pasó de 16 a **21 pruebas**, todas aprobadas. Las cinco nuevas pruebas cubren los tres casos de desfase indefinido y los dos casos de representación extrema de tiempo/paso. Las conclusiones principales, el paso adoptado y la combinación de parámetros elegida se mantienen.

## Resultado de la segunda revisión

No se encontraron nuevos apartados omitidos, pero sí tres fragilidades metodológicas o de implementación que convenía resolver. Se ejecutó nuevamente el análisis completo y pasaron las **16 comprobaciones**.

### Búsqueda paramétrica más amplia

La primera versión refinaba solo alrededor del mejor punto de la grilla inicial. Esto podía pasar por alto una región interior mejor entre sus nodos. Se agregó una grilla global de 9 valores por parámetro, antes del refinamiento local.

Ahora se exploran **7170 combinaciones distintas** entre la grilla inicial de 625 puntos, la global de 6561 y la local de 625, eliminando los puntos compartidos. Las tres resoluciones, evaluadas con el mismo h, eligen `(A,beta,gamma,sigma)=(0.60,0.30,0.05,1.00)`, con SP=0.991270 a h=12/4096.

En la grilla global se realizan 5832 comparaciones entre vecinos por cada parámetro. No hay inversiones de los signos observados para SP: aumenta con A, beta y sigma y disminuye con gamma. Esto respalda la selección discreta sin afirmar un máximo global continuo. Se agregaron `resoluciones_grilla.csv`, `tendencias_conjuntas.csv` y la figura 09.

Para que ampliar la grilla no obligue a reservar toda su trayectoria a la vez, se dividieron las evaluaciones en bloques de hasta 256 escenarios. La nueva prueba comprueba que particionar una lista conserva resultados y orden. El tamaño del historial por bloque, al paso adoptado, es de unos 48 MiB; la memoria total del proceso incluye además referencias, parámetros y tablas.

### Indicadores sin dependencia innecesaria del penúltimo ciclo

La función `indicadores` exigía [36,60] incluso para calcular únicamente SP y desfase. Se corrigió para que acepte [48,60]. Si no dispone del ciclo previo, calcula los indicadores y devuelve NaN en el diagnóstico de repetición anual. La admisibilidad de las grillas sigue exigiendo que ese diagnóstico se haya medido y cumpla su tolerancia.

Se ajustó la prueba de ciclo incompleto para rechazar datos que realmente no cubren [48,60], y se agregó una prueba de uso solo con el último ciclo. No se cambió la definición de SP, el promedio temporal ni el tratamiento de los picos.

### Comparación con la referencia usando el paso real

La comparación nominal contra la referencia tenía fijado h=12/4096, aunque la selección del paso puede refinar más si no cumple las tolerancias. Ahora se realiza después de elegir h y verifica la coincidencia de los nodos. `paso_parte2_vs_referencia.csv` registra el paso efectivamente usado, la diferencia máxima y el instante en que ocurre.

En esta ejecución el paso sigue siendo 0.0029296875 meses y la diferencia máxima nominal de T sigue siendo 0.00131582 °C. La diferencia máxima relativa entre pasos de Pmax, Pmin y Pmed sigue siendo 0.0412% en la grilla ampliada.

### Comprobación independiente y presentación

Se agregó una prueba de Euler contra la solución analítica logística con excitación estacional del caso desacoplado alpha=beta=0. Así la validación no depende exclusivamente de comparar una implementación de Euler consigo misma. Es un caso de prueba independiente: la referencia del TP continúa siendo Euler sobre el sistema original.

Se agregaron unidades en los ejes de sensibilidad y de los mapas conjuntos. El ejemplo ejecutable de `modelo.py` usa ahora el paso ya validado y remite a su justificación, en lugar de mostrar el mensaje obsoleto de que faltaba verificarlo. Los cambios no alteran los resultados principales del caso nominal ni del escenario seleccionado.

## Resultado de la primera revisión

- El análisis completo finalizó sin errores y las 13 comprobaciones automatizadas pasaron.
- Las tres referencias dan un orden de 1.000025; el paso de Parte 2 presenta una diferencia máxima de T de 0.001316 °C frente a la referencia nominal de Parte 1.
- El barrido fino detecta una desviación explosiva ya en h=1.1860 meses: el error de T es aproximadamente 4.42·10^49 °C. Es el primer punto detectado por el criterio descriptivo de ese barrido, no un umbral exacto.
- La máxima diferencia relativa de Pmax, Pmin y Pmed en la comprobación común a todos los escenarios es 0.0412%, menor que la tolerancia de 0.5%.
- La mayor diferencia normalizada de trayectorias completas entre el paso adoptado y su mitad es 0.0669%, también menor que 0.5%.
- La combinación `(A,beta,gamma,sigma)=(0.60,0.30,0.05,1.00)` sigue ganando entre los finalistas con h, h/2 y h/4. Con el paso común h=12/4096, SP=0.991270; al refinar a h/4, SP=0.990740. Se mantiene el paso común en las comparaciones principales y se informa el refinamiento adicional por separado.
- Los ciclos de todos los escenarios estudiados tienen un único pico interior de P, posterior a Tmax; no fue necesario resolver una ambigüedad de picos.

## Cobertura del enunciado

| Apartado | Evidencia y comprobación |
|---|---|
| 1(a): referencia con Euler | Tres pasos finos; diferencia entre referencias, orden de esas diferencias y comparación con el menor error de estudio. `resultados/referencia.csv`. |
| 1(b): duplicación de h y máximo error de T | Mallas anidadas que dividen 60 meses; comparación sin interpolación en los mismos tiempos. `resultados/convergencia.csv`. |
| 1(c): orden observado | Fórmula del enunciado y pendiente log-log inicial; comparación con el orden 1 esperado. Figura 01. |
| 1(d): crecimiento abrupto e inestabilidad | Barrido general y nuevo barrido fino de 0.0001 meses, oscilaciones de P y diagnóstico lineal local. Figuras 02 y 08; `estabilidad.csv`, `estabilidad_detalle.csv` y `estabilidad_lineal_local.csv`. |
| 2(a): caso nominal y último año | SP y desfase solo en [48,60]; gráficos separados para 60 meses y último ciclo. Figuras 03 y 04. |
| 2(b): sensibilidad individual | A, beta, gamma y sigma en los intervalos pedidos, un parámetro por vez; ambas curvas por parámetro e interpretación de las ecuaciones. Figura 05 y `sensibilidad.csv`. |
| 2(c): combinación de mayor SP | Grillas globales de 5⁴ y 9⁴ puntos y refinamiento local, guarda contra explosión y residuo anual; verificación del ranking de los diez finalistas a tres pasos. Figura 09; `grilla.csv`, `resoluciones_grilla.csv`, `tendencias_conjuntas.csv` y `verificacion_finalistas.csv`. |
| 2(d): comparación con el nominal | Parámetros, SP, desfase, Pmax, Pmin y Pmed; gráficos comparativos, escenario de todos máximos e interpretación de por qué gamma no se maximiza. Figuras 06 y 07; `comparacion.csv`. |

## Hallazgos y correcciones

### 1. El análisis de estabilidad necesitaba mayor resolución

El barrido inicial cambiaba h de a 0.025 meses. En los puntos que terminaban con números finitos, E_T/h parecía seguir cerca del régimen de primer orden; otros puntos desbordaban. De allí se había concluido que no se observaba un crecimiento abrupto finito del error antes del desbordamiento.

Esa conclusión era demasiado fuerte: el muestreo estaba salteando pasos con trayectorias explosivas todavía finitas. Se agregó un barrido de 1501 valores en [1.15,1.30], con incremento 0.0001 meses, y una figura específica. Ahora el apartado 1(d) muestra directamente el salto del error de T que antes faltaba documentar.

La primera desviación marcada del barrido detallado se registra en el resumen generado. No se la presenta como umbral exacto: las trayectorias no lineales a horizonte finito dependen de h de manera muy sensible en esa zona.

Se agregó también el Jacobiano analítico y la región de estabilidad de Euler aplicada a los modos disipativos de una linealización local congelada. Esto proporciona una interpretación complementaria de las oscilaciones; no sustituye las simulaciones ni demuestra estabilidad global.

### 2. La referencia tenía una comprobación incompleta de su régimen asintótico

Antes se comparaban solo dos pasos de referencia. Se agrega una tercera malla y se comprueba que las diferencias sucesivas tengan orden aproximadamente 1. Así queda mejor respaldada la estimación práctica del error remanente de la referencia.

También se mide directamente, para el nominal, la diferencia entre la solución al paso adoptado en Parte 2 y la referencia de Parte 1.

### 3. La precisión debía incluir todos los indicadores informados

El criterio anterior verificaba SP y desfase. Ahora incluye Pmax, Pmin y Pmed para toda la grilla y los barridos individuales, manteniendo una tolerancia relativa de 0.5% y una tolerancia de desfase de 0.02 meses. El resultado completo está en `resultados/verificacion_paso.csv`.

Se recalculan los diez mejores candidatos a h, h/2 y h/4 para comprobar que el ganador no sea un artefacto del paso. El programa se detiene si cambia el ganador, para evitar publicar una selección no resuelta.

### 4. La acotación necesitaba evidencia más informativa que una guarda numérica

La guarda de 10^6 sirve para detener trayectorias explosivas, pero no representa una escala física razonable por sí sola. Se guardan ahora los mínimos y máximos de T, H y P durante los 60 meses de los tres escenarios comparados, a paso h/4, en `resultados/acotacion_60m.csv`.

También se comparan sus trayectorias completas a dos refinamientos sucesivos, en los mismos tiempos y por variable. `resultados/convergencia_trayectorias.csv` conserva las diferencias absolutas y las diferencias normalizadas por la máxima magnitud de cada variable.

### 5. Se verificó explícitamente la interpretación del desfase

Se cuenta la cantidad de picos de P del último ciclo y se comprueba que el pico posterior a Tmax sea también Pmax. Así se evita una ambigüedad silenciosa entre máximo global y máximo local posterior en ciclos con varios picos. Los indicadores continúan usando exclusivamente [48,60]; el ciclo anterior se usa solo para medir la repetición anual.

## Comprobaciones automatizadas

La suite de la primera revisión contenía 13 comprobaciones: primer paso simultáneo; acoplamientos no nulos y fase estacional; Jacobiano por diferencias centradas; paso final e intervalo nulo; entradas inválidas; coincidencia escalar/lote; indicadores de una señal periódica de valores conocidos; ausencia de un pico posterior; exclusión del transitorio; rechazo de ciclos incompletos; orden por autoconvergencia; ausencia de recorte de negativos y detección de desbordamiento; y persistencia de la invalidez de un escenario desbordado por lotes. La segunda revisión agrega tres verificaciones, descriptas más abajo.

## Límites que se mantienen explícitos

- La solución de referencia es de Euler y tiene error remanente; las verificaciones de refinamiento no son cotas rigurosas.
- El escenario elegido es el mejor de las grillas exploradas. No se demuestra un máximo global continuo.
- La repetición anual se evalúa con una tolerancia operativa y dos ciclos finales. No se demuestra periodicidad exacta ni acotación para todo tiempo futuro.
- No se ajusta el modelo a observaciones climáticas: el enunciado suministra el modelo y sus parámetros, pero no una serie de datos para calibrarlo.
- La presentación disponible es el informe Markdown con sus figuras y datos. No se generó un PDF final de entrega; el enunciado no establece un formato obligatorio.

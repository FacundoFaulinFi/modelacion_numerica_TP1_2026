# Decisiones metodológicas del TP

Este registro distingue las exigencias del enunciado de las elecciones de implementación. Los valores y conclusiones finales se completan al ejecutar `analisis.py` y quedan en `resultados/resumen.md`.

## Herramientas y modelo

- Python, NumPy y Matplotlib. Las ecuaciones y Euler están en `modelo.py`; los experimentos reproducibles, en `analisis.py`.
- Se usa Euler explícito también para la referencia, como exige la Parte 1(a). No se reemplaza el método por un integrador de biblioteca.
- Se actualizan T, H y P simultáneamente con el estado anterior. No se recortan estados negativos ni se suavizan oscilaciones.
- Unidades: meses para el tiempo, °C para T, u.h. para H y u.p. para P. Se conservan los parámetros y condiciones iniciales del PDF.

## Parte 1(a): referencia

- Se ensayan referencias con pasos `12/2**16`, `12/2**17` y `12/2**18` meses; se utiliza la última para los errores. La tercera malla permite comprobar el orden observado de las propias referencias antes de estimar su error remanente.
- Se compara la diferencia entre las dos referencias más finas sobre nodos compartidos con el menor error de las mallas estudiadas. Se exige que esa diferencia sea menor que el 1% de dicho error. Es una comprobación práctica de refinamiento, no una demostración de error nulo.

## Parte 1(b,c): error y orden

- Pasos `h = 12/2**k`, con k desde 9 hasta 0: cada paso duplica el anterior y todos dividen 60 meses. Las mallas son anidadas.
- Error: máximo absoluto de T sobre los nodos de cada malla, frente a los nodos coincidentes de la referencia. No se introduce error de interpolación.
- Orden: `log2(E(2h)/E(h))`. Se informa en la fila correspondiente a 2h y solo se interpreta como orden en el rango convergente.
- Se estima también una pendiente log-log con las primeras cuatro mallas. Se espera orden 1; el resultado debe comprobarse, no imponerse.

## Parte 1(d): estabilidad

- Se amplía el paso hasta observar oscilaciones artificiales o crecimiento explosivo. Una trayectoria finita no basta para declararla estable.
- Se agrega un barrido de pasos entre 0.25 y 2 meses para describir la transición. En este barrido el último paso puede acortarse.
- Tras revisar el trabajo, se agrega un barrido con incremento 0.0001 meses en [1.15,1.30]. El incremento anterior de 0.025 meses omitía pasos explosivos que todavía llegaban finitos al final. La tabla detallada conserva tanto errores finitos como desbordamientos; el gráfico 08 muestra la pérdida abrupta de precisión que pide el apartado 1(d).
- Como señal descriptiva del salto se busca el primer punto del barrido fino con E/h mayor que diez veces el valor de la malla fina. Este factor no define un criterio universal ni un umbral exacto de estabilidad.
- Se agrega el Jacobiano analítico y el diagnóstico de Euler para modos disipativos de la linealización congelada, `h <= -2 Re(lambda)/|lambda|²`. Al variar J con el tiempo, esta comprobación local solo ayuda a interpretar las oscilaciones; no demuestra estabilidad global del sistema no lineal.
- Se distinguen error de discretización (disminuye aproximadamente en proporción a h al refinar), oscilaciones espurias y desbordamiento. No se interpreta el orden fuera del rango convergente.
- Las trayectorias que exceden 100 en magnitud o presentan estados negativos se rotulan como espurias en este experimento nominal; esta clasificación auxiliar no es un teorema de estabilidad. Se conservan sus errores en la tabla y se excluyen del gráfico de convergencia para no ocultar la escala de las mallas útiles.

## Parte 2: precisión e indicadores

- Paso inicial de exploración `12/2048` meses. Se verifica contra su mitad en los barridos y escenarios candidatos. Se exigen diferencias relativas de SP, Pmax, Pmin y Pmed menores al 0.5%, y una diferencia de desfase menor a 0.02 meses. Si no alcanza se refina nuevamente.
- Todos los indicadores se calculan únicamente en [48,60]. No se incorporan los primeros cuatro años al promedio ni a los extremos.
- `indicadores` puede trabajar directamente con el último ciclo. El ciclo [36,48] es opcional y se utiliza únicamente para el diagnóstico separado de repetición anual. Si no está disponible, los residuos se informan como NaN; no se declara periodicidad sin medirla.
- `Pmed` es el promedio temporal: integral de P por trapecios, dividida por 12. Así no se cuenta dos veces el extremo compartido del ciclo.
- Máximos y mínimos se toman en los nodos, sin ajuste parabólico. La resolución temporal de los picos se evalúa refinando h.
- Para el desfase se identifica el máximo global de T y el único pico interior de P dentro de [48,60], exigiendo que sea Pmax y posterior al máximo de T. Con varios picos de P, temperatura constante, Pmax en el extremo del intervalo o ausencia de un pico posterior, se informa NaN. Es una convención conservadora para no devolver un desfase ambiguo; no se afirma que todo ciclo con varios picos sea físicamente inválido. No se inventan picos fuera del intervalo ni se aplica un desplazamiento de 12 meses.
- Se verifica que en todos los ciclos estudiados haya un único pico interior de P y que el pico posterior coincida con Pmax. Si aparece una ambigüedad, el análisis se detiene para revisar la convención antes de publicar resultados.
- La comparación directa entre la solución nominal de Parte 2 y la referencia se hace después de seleccionar h, en nodos compartidos. Se guarda el paso real elegido en `paso_parte2_vs_referencia.csv`; no se supone que el refinamiento adaptativo siempre termine en un valor prefijado.

## Parte 2(b): sensibilidad individual

- Se exploran 21 valores equiespaciados de cada parámetro en los intervalos del enunciado, incluyendo los extremos y agregando el nominal cuando no pertenece a la malla.
- Se modifica un único parámetro por simulación; los restantes conservan sus valores nominales.
- Se guardan las curvas de SP y desfase y los datos completos en CSV.

## Parte 2(c): búsqueda conjunta

- Primera grilla: 5 valores por parámetro, 625 combinaciones. Se agrega una segunda grilla global de 9 valores por parámetro (6561 combinaciones), antes de refinar con 5 valores por parámetro en el entorno de la mejor combinación admisible. Así el refinamiento no se limita al entorno del primer ganador.
- Los valores de los ejes paramétricos se redondean a 12 decimales antes de eliminar duplicados, para que dos representaciones binarias de un mismo punto decimal no se cuenten como escenarios distintos. Esto no redondea los estados de Euler ni los resultados.
- La integración por lotes aplica la misma fórmula de Euler a todos los escenarios. Se verifica contra `modelo.euler` para un caso nominal y uno no nominal.
- Las evaluaciones se dividen en bloques de hasta 256 escenarios para limitar la memoria. Se verifica que cambiar el tamaño de bloque preserve los resultados y su orden. Se reutilizan únicamente evaluaciones con los mismos parámetros y el mismo paso.
- Se descarta un escenario si presenta desbordamiento o supera 10^6 en valor absoluto. Este límite solo es una guarda computacional frente a trayectorias explosivas, no una cota física.
- Para admitir un régimen anual al final se exige que, para cada variable, la máxima diferencia entre [36,48] y [48,60], a igual fase, sea menor al 1% de su máxima magnitud en el último ciclo. También se exige promedio de precipitación positivo, estados finales no negativos y un desfase medible.
- El umbral del 1% es una elección operativa, no una exigencia textual del PDF. Se informan los residuos para permitir evaluar el criterio. No demuestra periodicidad exacta ni acotación para todo tiempo futuro.
- Se verifica a medio paso toda la grilla explorada, se revisan los mejores candidatos y se refina si no cumple la precisión adoptada. El resultado se denomina **mejor escenario encontrado en las grillas**, no máximo global demostrado.
- La revisión agrega una comprobación del ranking de los diez mejores a h, h/2 y h/4. Si cambia el ganador, el programa pide ampliar la búsqueda a menor paso. Para los escenarios nominal, mejor y todos máximos, se comparan además las trayectorias completas en [0,60] y se guardan los mínimos y máximos observados en todo ese intervalo.
- Se comparan explícitamente los mejores resultados de las tres resoluciones paramétricas al mismo paso temporal. Además se comprueba el signo de las diferencias de SP entre vecinos de la grilla global, variando una coordenada y manteniendo las otras tres fijas. Esto es evidencia discreta, no una demostración de monotonicidad continua.

## Parte 2(d): comparación e interpretación

- Se compara el mejor escenario con el nominal y con el escenario de los cuatro parámetros máximos. Se informan los parámetros, SP, desfase, Pmax, Pmin y Pmed.
- Las tendencias se interpretan a partir de los términos: A modula R(t); beta multiplica la fuente T² de humedad; gamma multiplica el consumo HP; sigma multiplica la disipación de P. La severidad es una amplitud relativa, por lo que un aumento de precipitación no implica necesariamente un aumento de SP.
- Los resultados pertenecen a un modelo simplificado. No se presentan como pronóstico climático ni como ajuste a datos reales, porque no se suministran datos observacionales.

## Decisiones confirmadas tras ejecutar

- La referencia fina `12/2**18 = 0.0000457763671875` meses cumple el criterio adoptado: el cambio respecto a la referencia anterior equivale al 0.195% del menor error estudiado. El error absoluto remanente se estima del orden de ese cambio por tratarse de un método de primer orden, sin afirmar una cota rigurosa.
- Se adopta para la Parte 2 `h = 12/4096 = 0.0029296875` meses. En todos los escenarios de grilla y sensibilidad, el máximo cambio relativo de SP respecto al doble de paso es 0.131%, y el máximo cambio del desfase es 0.00879 meses. Ambos están debajo de las tolerancias adoptadas.
- La primera versión exploraba 1249 combinaciones distintas entre grilla inicial y local. La segunda revisión amplía la búsqueda a 7170 combinaciones distintas y mantiene el mismo ganador `(A,beta,gamma,sigma)=(0.60,0.30,0.05,1.00)`. La comparación de resoluciones está en `resultados/resoluciones_grilla.csv` y el resumen. No se afirma optimalidad global continua.
- No se identifica el primer desbordamiento como un umbral exacto. El barrido original muestra resultados no monótonos a horizonte finito: h=1.225 desborda, pero h=1.25 llega finito a 60 meses. El nuevo barrido fino sí detecta errores de T explosivos todavía finitos antes de esos pasos. Queda corregida la afirmación anterior de que no se observaba ese crecimiento: era una limitación del muestreo, no del sistema. Los pasos muy grandes pueden hacer tan pocas iteraciones que arrojen valores explosivos todavía finitos.
- Se mantienen seis cifras significativas en tablas de indicadores para permitir comparaciones; la resolución física o numérica no debe confundirse con la cantidad de dígitos impresos. Los CSV guardan los valores de cálculo.

## Verificación independiente del integrador

- Se agrega un caso desacoplado con solución analítica: alpha=beta=0, H=P=0. Entonces `T'=R(t)*T*(1-T/Tmax)` y la solución logística se obtiene integrando R(t). La prueba compara Euler con esa solución conocida y comprueba orden 1 al refinar.
- Este caso es exclusivamente una comprobación del código. Las simulaciones y referencias exigidas por el TP mantienen el sistema original, sus condiciones iniciales y el método de Euler.

## Auditoría de los resultados guardados

- `verificar_resultados.py` reconstruye las tendencias conjuntas desde la grilla CSV, comprueba sus rangos y cobertura, verifica los barridos de un parámetro por vez y contrasta el ganador con la tabla de comparación y los refinamientos.
- Reintegra por separado los tres escenarios finales al paso adoptado, compara sus trayectorias con los CSV y recalcula los indicadores. También reproduce los extremos de toda la simulación a h/4. Así se comprueba que las tablas no hayan quedado desalineadas con el código o con otras salidas.
- La auditoría es complementaria: no reemplaza la referencia fina ni la verificación de precisión de todas las combinaciones. No vuelve a integrar la grilla completa.
- Se agregan pruebas de regresión para desfases ambiguos y para relaciones tiempo/paso extremas. Un intervalo positivo conserva sus dos extremos aunque t_final/h subdesborde a cero; una relación no representable por un paso excesivamente pequeño se rechaza con ValueError.

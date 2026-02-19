# ABP-MN
Los bioprocesos industriales presentan dinámicas no lineales que exigen métodos numéricos confiables para su simulación. Este trabajo evalúa la producción de hidroxi-L-lisina en un biorreactor fed-batch modelado por nueve EDO’s no lineales, resuelta mediante Euler explícito, Runge-Kutta de cuarto orden y Adams-Bashforth de segundo orden, con ode45 como referencia. Se analizaron precisión, convergencia y sensibilidad paramétrica vía Monte Carlo. RK4 obtuvo el menor error (RMSE = 1.31×10−7 mmol/L, h = 0.1), confirmando orden O(h⁴), mientras que el modelo demostró alta robustez paramétrica con coeficiente de variación del 0.03% en el producto final.

Palabras clave: biorreactor fed-batch, métodos numéricos, Runge-Kutta, Adams-Bashforth,
sensibilidad paramétrica, hidroxi-L-lisina.

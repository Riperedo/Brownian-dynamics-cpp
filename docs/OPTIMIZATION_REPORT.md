# Análisis de Optimización: Brownian Dynamics C++

Este documento detalla las áreas de oportunidad identificadas en el código actual para mejorar el rendimiento, la escalabilidad y la eficiencia numérica. El análisis se divide en mejoras algorítmicas (alto impacto) y optimizaciones de bajo nivel (impacto medio/bajo).

## 1. Algorítmico: Neighbor Lists (Lista de Vecinos)
**Impacto:** Crítico (Reduce complejidad de $O(N^2)$ a $O(N)$)
**Estado Actual:**
Actualmente, el cálculo de fuerzas (`calculateForces`) itera sobre **todos los pares de partículas** ($i, j$) en un doble bucle anidado:
```cpp
for (int i = 0; i < numParticles - 1; ++i) {
  for (int j = i + 1; j < numParticles; ++j) {
     // Check cutOff ...
  }
}
```
Esto escala cuadráticamente. Para $N=1000$, son $\approx 500,000$ interacciones. Para $N=10,000$, son $\approx 50,000,000$.

**Propuesta:** Implementar **Cell Lists (Linked Cells)** o **Verlet Lists**.
*   Dividir la caja de simulación en celdas de tamaño $r_{cut}$.
*   Solo interactuar con partículas en la misma celda y celdas vecinas (27 en 3D).
*   Esto hace que el costo sea lineal $O(N)$.

## 2. Polimorfismo Estático (Templates vs Virtual)
**Impacto:** Alto (Inlining en el bucle más interno)
**Estado Actual:**
La clase `InteractionPotential` usa funciones virtuales (`virtual calculateForceMagnitude`).
```cpp
double fMag = potential->calculateForceMagnitude(r);
```
En cada interacción de cada paso, el procesador debe hacer una búsqueda en la v-table, lo que impide que el compilador haga "inlining" de la fórmula de fuerza.

**Propuesta:** Usar **CRTP (Curiously Recurring Template Pattern)** o Políticas de Templates.
*   Hacer que `Simulation` sea una plantilla de `Simulation<Dim, PotentialType>`.
*   El compilador generará una versión especializada de `calculateForces` donde la llamada a `calculateForce` es directa y puede vectorizarse.

## 3. Optimización de Memoria (SoA vs AoS)
**Impacto:** Medio (Mejor uso de Caché y SIMD)
**Estado Actual:**
Array of Structures (AoS). Tenemos `std::vector<Vector<Dim>>`.
En memoria: `X1 Y1 Z1, X2 Y2 Z2, ...`
Para cálculos vectoriales (SIMD AVX/SSE), es mejor tener todos los X, luego todos los Y.

**Propuesta:** Structure of Arrays (SoA).
*   `std::vector<double> posX, posY, posZ;`
*   Permite cargar 4 u 8 coordenadas en un solo registro del procesador y operar simultáneamente.

## 4. Generación de Números Aleatorios
**Impacto:** Medio
**Estado Actual:**
Se usa GSL (`gsl_rng_mt19937`). Es de muy alta calidad pero relativamente lento y no es intrínsecamente paralelo (aunque se usa una instancia por hilo o compartido con cuidado, actualmente thread-safe por diseño o por ser solo lectura? GSL RNG no es thread-safe si se comparte. *Nota: En la implementación actual OpenMP, el RNG está fuera del bucle paralelo de fuerzas, pero INTEGRATE no está paralelizado 100% o usa RNG compartido? Revisar.*)

**Propuesta:**
*   Usar generadores modernos y rápidos como **Xoshiro256++** o **PCG**.
*   Implementar un estado de RNG **por hilo** (Thread-Local RNG) para paralelizar la fase de integración (`integrate()`) sin bloqueos.

## 5. Condiciones de Frontera (Branchless PBC)
**Impacto:** Bajo/Medio
**Estado Actual:**
Usa `if` o `round` dentro del bucle interno.
```cpp
dr[k] -= std::round(dr[k] / boxSize) * boxSize;
```
Esto es costoso.

**Propuesta:**
*   Si la caja es ortogonal, usar aritmética optimizada o instrucciones condicionales de hardware (CMOV) evitando saltos y llamadas a `round`.

## Resumen de Prioridades

| Prioridad | Optimización                         | Dificultad | Ganancia Est. ($N=1000$) |
| :-------- | :----------------------------------- | :--------- | :----------------------- |
| **1**     | **Cell Lists / Neighbor Lists**      | Alta       | 10x - 50x                |
| **2**     | **Static Polymorphism (Templates)**  | Media      | 20% - 30%                |
| **3**     | **OpenMP en Integrate (RNG local)**  | Media      | 15%                      |
| **4**     | **Flags Compilador (-march=native)** | Baja       | 10%                      |

Se recomienda comenzar con la implementación de **Cell Lists** si se planea simular sistemas con $N > 500$.

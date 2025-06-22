# ZEKROM — Zonal Error-bound K-dimensional Reduced Order Models

Implementazione di strategie di controllo dell'errore a posteriori per modelli ridotti, ispirate ai Capitoli 3 e 4 del testo *Certified Reduced Basis Methods*, applicate a un problema ellittico parametrico con condizioni al bordo di Dirichlet non omogenee.

## Obiettivo del Progetto

Partendo da un problema di Poisson in due dimensioni, si applicano tecniche di **model order reduction** basate su:

- **Proper Orthogonal Decomposition (POD)**
- **Riduzione Galerkin**
- **Stima dell’errore a posteriori** con metodo **Min-Theta**

Lo scopo è:
- generare basi ridotte che rappresentano bene il manifold delle soluzioni;
- stimare in modo rigoroso l'errore con bound computabili;
- valutare l'efficacia computazionale (speedup) della riduzione;
- analizzare l'effetto della coercività;
- includere condizioni al bordo non omogenee dipendenti dal parametro.

## Problema Modellato

Si considera:

$$
-\nabla \cdot (a(x; \mu) \nabla u(x; \mu)) = f(x) \quad \text{in } \Omega
$$
$$
u(x; \mu) = g(x; \mu) \quad \text{su } \partial \Omega
$$

- Dominio: $\Omega = \{x \in \mathbb{R}^2 \mid \|x\| < 1\}$
- Parametri: $\mu = (\mu_1, \mu_2) \in [-1,1]^2$
- $f(x) = 1$
- $g(x; \mu) = \mu_1(x_1 + x_1^2 + x_2^2) + \mu_2(x_2 + x_1x_2)$
- $a(x; \mu) = 5 + \mu_1 x_1 + \mu_2 x_2$ (affine e strettamente positivo)

## Metodologia

- **Truth Solver FEM**: con FEniCS e lifting per $g(x;\mu)$
- **Costruzione base ridotta**:
  - Sampling parametrico
  - Raccolta snapshot
  - SVD per POD
- **Stima errore a posteriori**:
  - Min-Theta con coercività analitica
  - Norme duali su residui
- **Fase Online**:
  - Proiezione Galerkin ridotta
  - Complessità indipendente da numero DOF originali

## Esperimenti e Analisi

Sono stati studiati i seguenti aspetti chiave:

* l'errore in funzione del numero del numero delle basi ridotte;
* l'error estimator in funzione delle basi ridotte;
* fissato il numero di basi, l'errore e l'estimatore in funzione di un parametro;
* l'effetto dell'approssimazione della costanza di coercività;
* il vantaggio computazionale in termini di speedup;

## Requisiti

- Python ≥ 3.8
- FEniCS (`dolfin`, `mshr`)
- NumPy, SciPy, Matplotlib

## Esecuzione

Il notebook è pensato per Google Colab, con installazione automatica di FEniCS e librerie.

## Riferimenti

- * Certifed Reduced Basis Methods for Parametrized Partial Diferential Equations* di Jan S Hesthaven, Gianluigi Rozza, Benjamin Stamm – Chapters 3–4


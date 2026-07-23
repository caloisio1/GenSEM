# GenSEM 0.7.2 — Verde de punta a punta (2026-07-23)

Versión exacta: **GenSEM 0.7.2** (tarball `GenSEM_0.7.2.tar.gz`, Desktop y `GenSEM_validation\`; único diff vs 0.7.1: DESCRIPTION + script de validación del Ej.2). R 4.6.1, Windows 11, instalada con `--install-tests`. Datos: `GenSEM_validation\run\` (divorce2 y children2 de Wayback; children1 reconstruido y validado exactamente por Gate A).

| Corrida | Resultado | Tiempo | Detalle |
|---|---|---|---|
| Suite testthat (`NOT_CRAN=true`) | **80 PASS / 0 FAIL / 0 SKIP / 0 WARN** | ~35 min | 14 tests, 8 archivos |
| Ej.1 divorce2 (no recursivo) | **PASS** | 243 s | 6 fijos dentro de tol; RE y estructurales dentro de 1.96·SE publicado; AGHQ nAGQ=11 |
| Ej.2 children1 (recursivo probit) | **PASS (gates A–D)** | 6.257 s | A: modelo simple −0.3818 vs −0.382 publicado (±0.01). B: nll exacta Q=61 1668.298 < 1669.668 (dominancia +1.370). C: hospital −0.7881 dentro de 1.96·SE pub (0.473); var(U) 0.4227 ∈ [0.133, 1.255]; var(V) 3.5318 ∈ [2.492, 6.911]. D: Weibull −0.7800 vs Poisson −0.7881 (diff 0.008 < 0.20). cov(V,U) 0.6387 (SE 0.308) vs 0.2157 (SE 0.189): z = 1.17, compatible — reportado, no gateado. Sensibilidad dedup: −0.6484. Trinivel: sw² = 0.0006 |
| Ej.3 children2 (mlogit) | **PASS** | 13.095 s | hospital1 −0.3925 (tol 0.15), hospital2 −2.9768 (tol 0.30); var(V) 15.14 vs 13.03 publicada (dentro del régimen), var(U) 0.271 |

Hallazgo central documentado en el script del Ej.2: el punto publicado por gsem NO es el óptimo de la verosimilitud exacta de su propia especificación (queda +1.33/+1.37 nll por encima del MLE exacto, estable hasta Q=201); es el óptimo de la aproximación adaptativa de 7 puntos de Stata. El gate 0.7.2 certifica lo certificable (modelo simple exacto, dominancia nll, CIs publicados, consistencia interna) y mide lo demás.

Fe de erratas de JOB1: la cita st0481 del do-file oficial es correcta (keyword del propio paper); la sugerencia previa de st0500 era errónea.

# Zbývající problémy v modulu tikhonov.c

**Dokument:** Analýza kódu po opravách V4.7  
**Datum:** 2025-11-28  
**Stav:** Kritické chyby opraveny, zbývají střední a nízké priority

---

## Přehled

| Kategorie | Opraveno | Zbývá |
|-----------|----------|-------|
| Kritické logické chyby | 3/3 | 0 |
| Závažné strukturální chyby | 2/2 | 0 |
| Numerické problémy | 1/3 | 2 |
| Kvalita kódu | 0/2 | 2 |
| Možná vylepšení | 1/4 | 3 |

---

## 1. Numerické problémy

### 1.1 Trace(H) aproximace platí pouze pro uniformní mřížky

**Závažnost:** Střední  
**Umístění:** `compute_gcv_score_robust()`, řádky ~360–375

**Popis problému:**

Výpočet trace(H) pro GCV kritérium používá analytický vzorec odvozený pro **uniformní mřížku**:

```c
for (int k = 1; k <= n-2; k++) {
    double theta = M_PI * k / n;
    double eigenval = 4.0 * pow(sin(theta/2.0), 2) / (h_avg * h_avg);
    trace_H += 1.0 / (1.0 + lambda * eigenval);
}
```

Tento vzorec předpokládá, že vlastní čísla matice K jsou `4 sin²(πk/2n) / h²`, což platí **pouze** pro uniformní mřížku s konstantním krokem h.

**Důsledky:**

- Pro neuniformní mřížky je trace(H) nepřesná
- GCV skóre je zkreslené → optimální λ může být suboptimální
- Kód sice vypisuje varování pro `ratio > 2.0`, ale pokračuje s nepřesným výpočtem

**Možná řešení:**

1. **Stochastická aproximace trace:** Použít Hutchinsonův estimátor
   ```c
   // trace(A) ≈ z^T A z, kde z je náhodný vektor ±1
   ```

2. **Diagonální aproximace:** Spočítat přímo diagonální prvky (A + λK)⁻¹

3. **Robustní alternativa:** Pro neuniformní mřížky použít L-curve místo GCV

---

### 1.2 Potenciální dělení nulou

**Závažnost:** Nízká (teoretická)  
**Umístění:** `build_band_matrix()`, `compute_functional()`, `compute_derivatives()`

**Popis problému:**

Kód dělí hodnotami `h_left`, `h_right`, `h_sum` bez explicitní kontroly:

```c
double w = 2.0 * lambda / h_sum;
AB[kd + j*ldab] += w * (1.0/h_left + 1.0/h_right);
```

**Současná ochrana:**

Hlavní funkce `tikhonov_smooth()` validuje monotonnost:
```c
for (int i = 1; i < n; i++) {
    if (x[i] <= x[i-1]) {
        fprintf(stderr, "Error: x array must be strictly increasing\n");
        return NULL;
    }
}
```

**Zbývající riziko:**

- Statické funkce neprovádějí vlastní validaci
- Extrémně malé hodnoty h (např. `1e-300`) mohou způsobit overflow při `1/h²`
- Při budoucím refaktoringu může být validace omylem odstraněna

**Doporučení:**

Přidat defenzivní kontrolu s minimálním krokem:
```c
#define H_MIN_SAFE 1e-15
if (h_left < H_MIN_SAFE || h_right < H_MIN_SAFE) {
    fprintf(stderr, "Error: Grid spacing too small\n");
    return NULL;  // nebo graceful degradation
}
```

---

## 2. Kvalita kódu

### 2.1 Magic numbers bez dokumentace

**Závažnost:** Nízká  
**Umístění:** Celý modul

**Seznam nedokumentovaných konstant:**

| Hodnota | Místo | Předpokládaný účel |
|---------|-------|-------------------|
| `0.15` | `CV_THRESHOLD` | Hranice pro výběr diskretizační metody |
| `0.7` | `compute_gcv_score_robust` | Práh pro penalizaci vysokého trace(H) |
| `10.0` | `compute_gcv_score_robust` | Síla exponenciální penalizace |
| `5000` | `compute_gcv_score_robust` | Přepnutí na rychlou aproximaci |
| `20000` | `find_optimal_lambda_gcv` | Práh pro "velký dataset" |
| `1e-6`, `1e0` | `find_optimal_lambda_gcv` | Rozsah hledání λ |
| `0.3`, `1.7` | refinement loop | Faktory pro zjemnění hledání |

**Doporučení:**

Definovat jako pojmenované konstanty s komentářem:
```c
/* Threshold for switching discretization methods.
 * Based on empirical testing: CV < 0.15 indicates near-uniform grid
 * where average coefficient method is sufficiently accurate. */
#define CV_THRESHOLD 0.15

/* GCV trace penalty threshold.
 * When trace(H)/n > 0.7, the smoother is overfitting.
 * Exponential penalty discourages such solutions. */
#define GCV_TRACE_PENALTY_THRESHOLD 0.7
#define GCV_TRACE_PENALTY_STRENGTH 10.0
```

---

### 2.2 Oddělení I/O od výpočetní logiky

**Závažnost:** Nízká (architektonická)  
**Umístění:** `find_optimal_lambda_gcv()`, `compute_gcv_score_robust()`

**Popis problému:**

Funkce přímo tisknou na stdout:
```c
printf("# GCV optimization for n=%d points\n", n);
printf("# λ=%9.3e: J=%9.3e, RSS=%9.3e, tr(H)=%6.1f...\n", ...);
```

**Důsledky:**

- Nelze použít v GUI aplikacích nebo serverovém kódu
- Obtížné unit testování
- Nelze přesměrovat výstup bez hackování

**Doporučení:**

Zavést callback nebo strukturu pro diagnostiku:
```c
typedef struct {
    int verbose_level;
    FILE *output_stream;
    // nebo callback:
    void (*log_callback)(const char *msg, void *user_data);
    void *user_data;
} TikhonovOptions;

double find_optimal_lambda_gcv(double *x, double *y, int n, 
                                GridAnalysis *grid_info,
                                TikhonovOptions *options);
```

---

## 3. Možná vylepšení

### 3.1 Vyšší přesnost derivací na hranicích

**Umístění:** `compute_derivatives()`

**Současný stav:**

```c
/* Forward difference at start - O(h) */
y_deriv[0] = (y_smooth[1] - y_smooth[0]) / (x[1] - x[0]);

/* Backward difference at end - O(h) */
y_deriv[n-1] = (y_smooth[n-1] - y_smooth[n-2]) / (x[n-1] - x[n-2]);
```

**Problém:**

Jednostranné diference prvního řádu mají chybu O(h), zatímco centrální diference ve vnitřních bodech mají O(h²).

**Vylepšení:**

Použít jednostranné formule druhého řádu:
```c
/* Forward difference O(h²) - requires 3 points */
if (n >= 3) {
    double h0 = x[1] - x[0];
    double h1 = x[2] - x[1];
    y_deriv[0] = (-(2*h0 + h1)*y[0] + (h0 + h1)*(h0 + h1)*y[1]/(h0*h1) 
                  - h0*h0*y[2]/(h1*(h0 + h1))) / h0;
}

/* Nebo jednodušší pro uniformní mřížku: */
y_deriv[0] = (-3*y[0] + 4*y[1] - y[2]) / (2*h);
```

---

### 3.2 Validace vstupů v pomocných funkcích

**Popis:**

Statické funkce (`build_band_matrix`, `compute_derivatives`, atd.) spoléhají na validaci v hlavní funkci. Pro robustnost a snazší debugging by měly mít vlastní asserty:

```c
#include <assert.h>

static void build_band_matrix(double *x, int n, double lambda, 
                              double *AB, int ldab, int kd,
                              GridAnalysis *grid_info)
{
    assert(x != NULL);
    assert(AB != NULL);
    assert(n >= 1);
    assert(ldab >= kd + 1);
    assert(lambda >= 0.0);
    
    // ... zbytek funkce
}
```

---

### 3.3 Podpora pro volitelné váhy datových bodů

**Popis:**

Současná implementace řeší:
$$\min_u \|y - u\|^2 + \lambda \|D^2 u\|^2$$

Obecnější formulace s váhami:
$$\min_u \sum_i w_i (y_i - u_i)^2 + \lambda \|D^2 u\|^2$$

**Využití:**

- Různá důvěryhodnost měření
- Robustní regrese (iterativně převážené nejmenší čtverce)
- Heteroskedastická data

**Implementace:**

```c
TikhonovResult* tikhonov_smooth_weighted(
    double *x, double *y, double *weights,  // weights can be NULL
    int n, double lambda, GridAnalysis *grid_info);
```

---

## Shrnutí priorit

| Priorita | Problém | Dopad | Náročnost opravy |
|----------|---------|-------|------------------|
| Střední | Trace(H) pro neuniformní mřížky | GCV může vybrat suboptimální λ | Vysoká |
| Nízká | Dělení nulou | Teoretický edge case | Nízká |
| Nízká | Magic numbers | Údržba kódu | Nízká |
| Nízká | I/O oddělení | Integrace do jiných systémů | Střední |
| Nízká | Přesnost derivací | Mírně horší derivace na okrajích | Nízká |

---

## Závěr

Po opravách V4.7 je modul `tikhonov.c` funkčně korektní pro většinu praktických případů. Zbývající problémy jsou převážně:

1. **Numerické edge cases** — trace(H) aproximace pro silně neuniformní mřížky
2. **Kvalita kódu** — dokumentace konstant, oddělení I/O
3. **Možná rozšíření** — přesnější derivace, váhované regrese

Žádný ze zbývajících problémů nepředstavuje kritickou chybu, ale jejich řešení by zlepšilo robustnost a udržovatelnost kódu.

# Audit kódu `smooth` (v5.11.8)

**Datum:** 2026-04-25
**Rozsah:** všechny produkční moduly (~3400 řádek C kódu)
**Řazení:** podle dopadu (A = bugy, B = designová rozhodnutí, C = drobnosti)
**Status:** A1–A4 opraveny ve v5.11.9 (2026-04-25); B6, B12 opraveny ve v5.11.10 (2026-04-25)

---

## A. Bugy

### A1. ~~`-?` / neznámý přepínač vrací **exit 0**~~ — `smooth.c:81, 169–176` — **FIXED v5.11.9**

Optstring `"n:p:m:l:f:k:dgTh?"` má `?` jako validní znak. `getopt` vrací `?` jak pro `-?`, tak pro neznámé volby. Obě padají do stejného case, který volá `help()` a `exit(EXIT_SUCCESS)`. `default:` je nedosažitelný.

```c
case '?':
  help();
  exit(EXIT_SUCCESS);   // <- i pro neznámý flag
```

**Dopad:** skripty kontrolující exit code dostanou "úspěch" na `./smooth --bogus`.
**Fix:** rozdělit — `h` na `EXIT_SUCCESS`, `?` na `EXIT_FAILURE`.

### A2. ~~Butterworth neemituje header "Derivative units: dy/dt"~~ — `smooth.c:593–599` — **FIXED v5.11.9**

Ve větvích Tikhonov / Savgol / Polyfit s `-T -d` se vypíše:

```
# Derivative units: dy/dt (t in seconds)
```

Butterworth blok to **nevypisuje**. Uživatel v timestamp módu nemá informaci o jednotkách derivace.

### A3. ~~Potenciální signed overflow v parsování subsekund~~ — `timestamp.c:45–56` — **FIXED v5.11.9**

```c
while (frac_start[digits] >= '0' && frac_start[digits] <= '9') {
    frac_value = frac_value * 10 + (frac_start[digits] - '0');
    digits++;
}
```

Smyčka nemá horní limit — pro 19+ číslic `long` přeteče (UB pro signed). Ošetření je až v `if (digits > 0 && digits <= 9)` *po* smyčce. Potřeba přidat `&& digits < 19` do podmínky smyčky (nebo přestat akumulovat po 9 číslicích).

### A4. ~~`h_std` dělí špatným jmenovatelem~~ — `grid_analysis.c:119` — **FIXED v5.11.9**

```c
analysis->h_std = sqrt(sum_sq_diff / (n-1));
```

Počet intervalů je `N = n-1`, střed `h_avg` se počítá ze stejných dat → pro nezkreslený sample STD má být dělitel `N-1 = n-2`. Pro n=20 ~3 % odchylka; pro malé datové sady to posouvá CV a může překlopit hranici 0.05 (Savgol) nebo 0.15 (Butterworth/Tikhonov).

---

## B. Designová / konzistenční rozhodnutí

### B1. Tři rozdílné sady CV threshů napříč moduly

- `grid_analysis.c`: 0.01 / 0.2 / 0.5 / 1.0 (varovné úrovně)
- `butterworth.c`: 0.01 / 0.05 / 0.15
- `savgol.c`: 0.05
- `tikhonov.c`: 0.15 pro discretization, ale 0.2 v GCV reporteru (`tikhonov.c:546`)

Neexistuje jediný zdroj pravdy. Doporučení: centralizovat do `grid_analysis.h`.

### B2. "GCV" v Tikhonov není standardní GCV — `tikhonov.c:437–441`

```c
double trace_ratio = trace_H / n;
if (trace_ratio > 0.7) {
    double penalty = exp(10.0 * (trace_ratio - 0.7));
    gcv_score *= penalty;
}
```

Ad-hoc penalizace. Funkce se jmenuje `compute_gcv_score_robust` a uživatelský text hlásí "GCV optimization", ale implementuje modifikovaný GCV. Buď přejmenovat ("penalized-GCV"), nebo zdokumentovat výslovně v README.

### B3. Skryté velikostní prahy v Tikhonov — `tikhonov.c:405, 551, 591`

- `n ≤ 5000` → analytický trace(H); jinak aproximace
- `n > 20000` → přidá se L-curve; jinak jen GCV
- `n ≤ 5000` → přidá se refinement krok

Žádný z prahů není v dokumentaci. Uživatel dostane nespojité chování při překročení `n=5000` nebo `n=20000` se stejnými daty.

### B4. Kód výstupu duplikován 4× — `smooth.c:464–672`

Každá metoda má vlastní copy-paste block ~50 řádek pro header + output loop (timestamp × derivative × 2). Celkem ~200 řádek redundance. Extrakce do `print_result(x, ts_ctx, result->y_smooth, result->y_deriv, n, show_derivative, timestamp_mode)` by snížila chybovost (viz A2).

### B6. ~~`should_use_adaptive()` — nepoužité API~~ — `grid_analysis.c:246–255`, `grid_analysis.h:59` — **FIXED v5.11.10**

Deklarováno v hlavičce, implementováno, nevolá se nikde (grep v `*.c` a `tests/`). Buď používat, nebo smazat.

### B7. Tichý per-point fallback v polyfit — `polyfit.c:272–277`

Při selhání `dgelss` (info != 0) se použije `y[i]` (raw hodnota) a derivace 0. Žádný counter, žádný finální warning. Pokud by selhalo třeba 20 % oken, výstup je zubatý a uživatel to nezjistí.

### B8. Lineární růst bufferu v I/O — `smooth.c:267, 353`

```c
abuf += BUF;  // BUF=512
```

Pro 100k řádků = 200 realloc-kopií, amortizovaná O(n²). Doubling (`abuf = abuf ? abuf*2 : BUF`) by to vyřešilo.

### B9. Pevné limity parsování — `smooth.c:206, 308, 320`

- `line[4096]` — delší řádek se rozdělí napůl (fragment pak vypadá jako nový řádek, rozbije počet řádků v chybové hlášce)
- `values[100]` / `col_count < 100` — 101. sloupec se tiše ztratí (bez warningu)

Pro uživatele s `-k 1:120` dostane jen chybu "line has only 100 columns", ne pravdu.

### B10. `decomment` volán jen pro soubory, ne pro stdin — `smooth.c:200–202`

Symetrie chybí. V praxi `strtod` přeskočí `#` řádky přes `col_count < 1`, takže to "funguje" — ale implicitně, ne explicitně. Jednodušší je volat `decomment` jednotně (včetně stdin), nebo to explicitně zdokumentovat.

### B11. Warning cíle (stdout vs stderr) rozdělené nekonzistentně

V `butterworth.c`:

- stdout jako `# komentář`: CV warning (552), auto-cutoff info (456)
- stderr: velký dataset (507), fc blízko Nyquistu (538)

Princip asi "co chceme zachovat v output file → stdout s `#`". Fungující, ale neexplicitní — přidat komentář do headeru modulu.

### B12. ~~Mrtvý kód~~ — **FIXED v5.11.10**

- `smooth.c:75–78`: `if (argc == 1) { /* Will read from stdin below */ }` prázdné tělo.
- `savgol.c:334–340`: `if (left_pts + right_pts < poly_degree)` — chyceno dříve `poly_degree < window_size`.
- `savgol.c:362–367`: `if (idx >= 0 && idx < n)` defensive check, vždy splněno (idx ∈ [0, window_size-1] ⊂ [0, n-1]).

(`polyfit.c:330–338` z původního auditu byl false positive — slouží jako fallback při `dgelss info != 0`, kdy `continue` na řádku 276 přeskočí nastavení `left_boundary_done`.)

### B13. Nezkontrolovaná `atoi` na CLI argumenty — `smooth.c:83–87`

`-n abc` silently → sp=0, pak "odd >= 3" error. Uživatel netuší, že předal blbost. Přejít na `strtol` s kontrolou `endptr`.

### B14. Pracovní soubory v repu

`but2.c.new, but3.c, but32.c, but4.c, grid_analysis.c.bak, polyfit.c.bak, tikhonov.c.bak` — starší varianty s překrývajícími se symboly (např. 3× `estimate_cutoff_frequency`). `grep` a IDE indexing je vidí a vrací falešné hity. Zkopírovat do jiného adresáře (mimo project root) nebo přidat do `.gitignore`.

### B15. `-T` a `-k` se tiše ignorují navzájem

V timestamp módu se `-k` nepoužije. Kód neemituje warning, když uživatel kombinaci zadá. Minimálně hlásit `# Warning: -k ignored in timestamp mode`.

---

## C. Drobnosti

- **C1** `savgol.c`: komentáře míchají angličtinu a češtinu (řádky 62, 97, 107, 109). Např. `"Nulování paměti (náhrada za calloc)"`.
- **C2** `polyfit.c:280–296`: condition/rank warning jen na prvním okně (`i == offset`). Horší okna později se neohlásí.
- **C3** `timestamp.c:116, 185`: `x_temp` se alokuje na `n`, kopíruje do nově alokovaného `*x_out` velikosti `valid_count`. Mohlo by se to alokovat rovnou `valid_count`, nebo rovnou vrátit `x_temp`.
- **C4** `smooth.c:514`: výstupní formát `%10.6lG` dává jen 6 sig. cifer pro y. Přesná data by potřebovala víc (`%.10lG`?).
- **C5** `smooth.c:76`: dead comment-only block.
- **C6** Makefile `test-valgrind`: užitečné by bylo `--error-exitcode=1` pro CI failure signal.

---

## Shrnutí podle priority

| Priorita | Položky |
|----------|---------|
| **Fix brzy** | ~~A1~~, ~~A2~~, ~~A3~~ (opraveno v5.11.9), B2 (pojmenování GCV), B7 (tichý fallback) |
| **Vyčistit při příležitosti** | B1 (centralizace threshů), B4 (extrakce output), B8 (doubling realloc), B13 (strtol), B14 (pracovní soubory) |
| **Kosmetika** | ~~A4~~ (opraveno v5.11.9), ~~B6~~, ~~B12~~ (opraveno v5.11.10), B3 (zdokumentovat prahy), C1 (konzistence komentářů) |

Žádný z nálezů není kritický pro správnost vědeckých výsledků — hlavní matematické cesty (Tikhonov pentadiagonální matice, Butterworth biquad cascade, SVD polyfit, Savgol pre-computed coefs) jsou solidní. Většina problémů jsou CLI/UX/robustnost a design drift, který se dá úspornou cestou postupně vyčistit.

# Code audit v5.11.22

Datum: 2026-04-30. Předchozí audit `code-audit-v5.11.8.md` je z velké části vyřešený (A1–A4, B2, B4, B6, B7, B8, B9, B10, B11, B12, B15, C2, C4, C6 zavřené ve v5.11.9 až v5.11.22). Tento audit vychází z čerstvého čtení produkčních modulů (~3500 řádků C) s důrazem na **designová rozhodnutí**, ne nové bug-hunty. Otevřené z v5.11.8: B1 (vědomě odloženo), B13 (atoi), B14 (pracovní soubory — částečně řešeno přes `.gitignore`), C1, C3, C5.

---

## A. Bugy

### A1. ~~`print_grid_analysis(grid_info, 0, "# ")` netiskne warning_msg~~ — `smooth.c:577`, `grid_analysis.c:231,238` — **FIXED v5.11.23**

`smooth.c:574–578` se chová jako kdyby `verbose=0` mělo vypsat detaily warningu, ale `print_grid_analysis()` chrání blok `if (analysis->reliability_warning) { printf("WARNING: ..."); }` podmínkou `verbose >= 1`. Důsledek: pro každou non-uniform mřížku uživatel vidí

```
# Grid analysis warnings:
# Grid uniformity analysis:
#   n = 200 points
#   h_min = ..., h_max = ..., h_avg = ...
#   CV = 0.42
#   Grid type: NON-UNIFORM
#   Uniformity score: 0.43
```

ale **konkrétní text warningu** ("HIGH grid non-uniformity: CV = 0.42 / Adaptive methods are strongly recommended ...") je naskládaný v `analysis->warning_msg` a nezobrazí se. Komentář na `smooth.c:577` ("basic statistics + warnings") je tím pádem zavádějící.

**Fix:** buď volat `print_grid_analysis(grid_info, 1, "# ")` v této větvi (nezavřený `-g` byl by pak verbose=1 i tak), nebo přesunout `if (reliability_warning)` blok do verbose=0 sekce v `print_grid_analysis`.

### A2. `parse_timestamp` přijme nesmyslné kalendářní datum — `timestamp.c:60–69, 81`

Validace kontroluje `day < 1 || day > 31`, ale 31. únor projde:

```c
if (year < 1970 || year > 2100 || month < 1 || month > 12 ||
    day < 1 || day > 31 || ...)
    return -1;
time_t epoch = timegm(&tm_time);  /* 2025-02-31 → 2025-03-03 (silently!) */
```

`timegm()` non-strict normalizuje přetečené pole (z `2025-02-31` udělá `2025-03-03`). Vrátí platný `time_t`, takže timestamp se přijme s tichou změnou data. Stejný problém pro `2025-04-31`, `2025-13-01` (chyceno limitem month, ale měsíc 13 → leden další rok, kdyby nebyl limit). Pro vědecká data, která mohou pocházet ze špatně formátovaných exportů, to je tichý posun časové stopy.

**Fix:** po `timegm()` zpětně ověřit `tm_time.tm_mday == day && tm_time.tm_mon == month - 1` — pokud se neshodují, vstup byl normalizován a měl by se odmítnout.

---

## B. Designová / konzistenční rozhodnutí

### B1. ~~Nekonzistentní kapitalizace error labelů~~ — napříč všemi moduly — **FIXED v5.11.24**

Tři současné formáty:

- `polyfit.c`, `tikhonov.c`, `savgol.c`, `grid_analysis.c`: `"Error: ..."` (Title Case)
- `butterworth.c`: `"ERROR: ..."` (UPPERCASE) — záměrně zdokumentováno v hlavičce modulu
- `savgol.c:216`: `"ERROR: Savitzky-Golay method not suitable ..."` jediný UPPERCASE výskyt v jinak Title-Case modulu

Konvence dokumentovaná v CLAUDE.md ("ERROR: ... — hard failures") sedí na butterworth, ale ostatní moduly ji nedodržují. Buď zvednout konvenci na celý projekt (vše `ERROR:` pro hard failure), nebo nechat butterworth jako výjimku a explicitně to pojmenovat. Stávající polostav je nejhorší ze všech variant.

**Recommendation:** sjednotit na `ERROR:` pro hard failure (návrat NULL) a `Error:` zrušit. Jeden grep&replace + drobné touchups v testech, pokud nějaké testy parsují stderr (nezdá se).

### B2. ~~const-correctness mezi metodami~~ — `butterworth.h:71` vs `tikhonov.h:55`, `savgol.h:78` — **FIXED v5.11.25**

```c
ButterworthResult* butterworth_filtfilt(const double *x, const double *y, int n,
                                        double cutoff_freq, int auto_cutoff,
                                        const GridAnalysis *grid_info);

TikhonovResult* tikhonov_smooth(double *x, double *y, int n, double lambda,
                                GridAnalysis *grid_info);

SavgolResult* savgol_smooth(double *x, double *y, int n, int window_size, int poly_degree,
                            GridAnalysis *grid_info);
```

Žádná z metod si vstupní `x`, `y` ani `grid_info` nemodifikuje (ověřeno čtením všech tří). Butterworth byl napsán s konzistentním const-flagováním, ostatní zdědily prototyp z dob, kdy ne. Polyfit nemá `grid_info` parametr vůbec.

**Recommendation:** dotáhnout `const` na všechny veřejné API metody — přínos čtenářské jasnosti je velký, kompilátorová cena nulová.

### B3. ~~Tikhonov třikrát přepočítává `h_avg` přes lokální výpočet~~ — `tikhonov.c:76, 116, 228, 401` — **FIXED v5.11.26**

`select_discretization_method` má fallback s vlastním výpočtem `h_avg` (`tikhonov.c:67–82`), ale ten se vždy zavolá s `grid_info != NULL` (jediné callsite v smooth.c projde `analyze_grid`). Mrtvý kód pro produkci. Stejný `(x[n-1] - x[0]) / (n-1)` se opakuje v `build_band_matrix:116`, `compute_functional:228` a `compute_gcv_score_robust:401` místo dotázání na `grid_info->h_avg`.

**Recommendation:** odstranit fallback v `select_discretization_method` (nahradit `assert(grid_info != NULL)` nebo NULL-check vrácením `DISCRETIZATION_AVERAGE`), všude dále používat `grid_info->h_avg`. Ušetří ~15 řádků a ujistí se, že CV-thresh policy je centralizovaný.

### B4. ~~main() je 641 řádků, dělá osm věcí~~ — `smooth.c:55–696` — **FIXED v5.11.27 + v5.11.28**

CLI parsování + I/O alokace + parser normal-mode + parser timestamp-mode + grid analýza + dispatch metod + cleanup. Vázáno hluboko na lokální stav (10× se opakuje cleanup vzor `for (int i = 0; i < n; i++) free(timestamp_strings[i]); free(timestamp_strings); free(y); ...`). B4 z předchozího auditu vyřešil duplicitu výstupu (`print_result`), ale parser stále zůstal monolitický.

```c
/* smooth.c — 10× tento blok */
for (int i = 0; i < n; i++) free(timestamp_strings[i]);
free(timestamp_strings);
free(y);
fclose(fp);
exit(EXIT_FAILURE);
```

**Recommendation:** vytáhnout `parse_input(FILE*, ParseMode, int x_col, int y_col, ParseResult*)` do nového modulu `parser.c` (`parser.h`). Centralizovat cleanup přes goto-pattern jak to dělá tikhonov. Zhruba 250-řádkový diff, ale main se zkrátí na ~200 řádků a parser bude testovatelný (test_parser.c už existuje a runs `./smooth` přes popen — místo toho by mohl volat `parse_input()` přímo).

**Fix:** dvoufázově. v5.11.27 nahradil 11 duplikovaných error-cleanup bloků v `main()` jediným `cleanup:` labelem na konci (goto-pattern), čímž popadl pre-existing exit-leaky na `n < sp` a "<method> failed!" cestách. v5.11.28 vytáhl ~290 řádků parsing logiky do nového modulu `parser.c`/`parser.h` (`parse_input()` + `free_parse_result()` + `ParseResult` struct s ownership-transfer kontraktem). `main()` 641 → 315 řádků, `smooth.c` 797 → 466. `MAX_LINE`/`MAX_COLS`/`BUF` konstanty + `math.h`/`errno.h` includes přesunuty s parserem. `tests/test_parser.c` zatím stále drive-uje binárku přes `popen()`; přepsání na přímé volání `parse_input()` je samostatný follow-up (audit motivace pro testovatelnost je tím připravená, ale ne vybraná).

### B5. ~~`estimate_cutoff_frequency` zveřejněna v hlavičce, používaná jen interně~~ — `butterworth.h:76`, `butterworth.c:428` — **FIXED v5.11.30**

```c
double estimate_cutoff_frequency(const double *y, int n);  /* in butterworth.h */
```

Volá se jen z `butterworth_filtfilt` (jediný callsite na `butterworth.c:526`). Žádný test ji nevolá. Podobné B6 ze v5.11.8 (`should_use_adaptive`) — buď používat externě, nebo dát `static`.

**Recommendation:** udělat `static` a odstranit z `butterworth.h`. Pokud má být API endpoint pro auto-cutoff (např. pro úplnou diagnostiku z `-g`), zdokumentovat a otestovat.

### B6. `savgol_coefficients` zveřejněná, jen interně používaná — `savgol.h:107`

Veřejná v hlavičce, ale `grep -r savgol_coefficients` ukáže jen volání zevnitř `savgol.c` (řádky 284, 292, 340, 343, 387, 390). Žádný test, žádný externí konzument. Stejný problém jako B5.

**Recommendation:** stejný vzorec — buď udělat `static`, nebo přidat unit testy. (Funkce je ale low-level a externě se zdá užitečná pro někoho, kdo by chtěl SG koeficienty bez celého filtru — pokud je to záměr, dodat aspoň jeden test.)

### B7. `n_intervals` v `GridAnalysis` je redundantní — `grid_analysis.h:24`

```c
int n_points;           /* Number of points */
int n_intervals;        /* Number of intervals (n_points - 1) */
```

Použito jen na jednom místě (`grid_analysis.c:246` v `print_grid_analysis` při verbose>=2), kde stačí `n_points - 1`. Drží data v synchronizaci, kterou kontrolovat nemusíme — odstraňme nadbytečné zrcadlení.

**Recommendation:** vymazat field, použít `analysis->n_points - 1` v `print_grid_analysis`. Marginální zlepšení, hodí se pro budoucí změny grid struktury.

### B8. Tikhonov verbose log per-lambda na stdout — `tikhonov.c:444–446`

```c
printf("# λ=%9.3e: J=%9.3e, RSS=%9.3e, tr(H)=%6.1f (%.2f), pGCV=%9.3e\n", ...);
```

Při `-l auto` se volá pro každý lambda kandidát (12–13 hodnot pro `n<=20000`, 12 pro `n>20000`, plus 8× refinement). Výstup má tedy 20–30 řádků `#` před vlastní hlavičkou výsledku. Dokumentováno jako "informace patřící k výsledku" (CLAUDE.md), ale je to spíš debug log. Pokud někdo pipeline-uje `smooth ... | gnuplot`, dostane v hlavičce dat 30+ řádků diagnostiky.

**Recommendation:** zvážit `-v` flag (verbose) pro tuhle úroveň detailu; v default módu vypsat jen "# pGCV: %d candidates evaluated, optimal λ=...". Není to bug, ale UX pošťouchnutí. Vyžaduje rozhodnutí uživatele/maintainera.

### B9. UTF-8 v stdout zprávách — `tikhonov.c:444, 565, 592, 608`

```c
printf("# λ=%9.3e: ...\n", ...);
printf("# Optimal λ: %.6e (pGCV=%.3e)\n", ...);
```

Lambdy a další UTF-8 (`²`, `σ̂`) jsou v komentářích kódu OK, ale v **stdout user output**, který je součástí dat a může se procházet shellem/parserem/ASCII pipelinou, jsou potenciálním problémem. Zbytek programu používá ASCII (anglické názvy), takže Tikhonov je odstávkový případ.

**Recommendation:** `# lambda=%.3e: J=...` místo `# λ=...`. CLAUDE.md má v "README.md character restrictions" zmínku o omezeních DejaVu fontů; stejné zdůvodnění platí, když data smoothing mají projít LaTeX-em / pdflatex-em jako data block. (Možná platí přísnější pravidlo "stdout = ASCII", než README.)

### B10. Magic numbers v Tikhonov GCV — `tikhonov.c:438–439, 594`

```c
if (trace_ratio > 0.7) {
    double penalty = exp(10.0 * (trace_ratio - 0.7));   /* magic: 0.7, 10.0 */
}
```

```c
double factor = 0.3 + 1.4 * i / 7.0;  /* magic: 0.3, 1.4, 7.0 */
```

První je pGCV penalizační konstanta (zdokumentovaná v README jako "Enhanced GCV"), ale v kódu nemá symbolickou konstantu. Druhé je refinement range pro lambda search — bez vysvětlení, proč zrovna `[0.3 * best, 1.7 * best]` v 8 krocích. Stejný motiv (skryté thresholdy v Tikhonov) jako B3 z v5.11.8 (vyřešený dokumentací v README), ale tady bez té dokumentační vrstvy.

**Recommendation:** definovat `#define PGCV_OVERFIT_THRESHOLD 0.7`, `#define PGCV_OVERFIT_DECAY 10.0`, `#define LAMBDA_REFINE_MIN 0.3`, `#define LAMBDA_REFINE_MAX 1.7`. Šetří inline-magii, dovolí budoucímu auditu měnit hodnoty na jednom místě.

### B11. polyfit nemá `grid_info` parametr, ale fakticky tolerable všechny grid — `polyfit.h:53`, `smooth.c:670`

CLAUDE.md říká "Methods inspect GridAnalysis* and either adapt (Tikhonov) or reject (Savgol) based on uniformity". Polyfit toleruje cokoli, takže by `grid_info` vůbec nepotřeboval. Ale **API drift**: tři metody mají signaturu `(..., GridAnalysis*)`, polyfit má `(...)`. Pro budoucí udržbu by byla lepší jednotná signatura.

**Recommendation:** dvě cesty:
1. Přidat `const GridAnalysis*` k polyfit API (i když nevolá), pro symetrii. Polyfit jen ignoruje. Plus: jednotný callsite v smooth.c.
2. Nechat tak a explicitně zdokumentovat ve `polyfit.h`: "Polyfit nepotřebuje grid analysis; CV-tolerantní per-window fit." 

Není kritické, ale dnes je to skrytá nesymetrie.

### B12. Boundary case `n == 1` v `analyze_grid` zanechává nulové h_min/h_max — `grid_analysis.c:53–58`

```c
if (n < 2) {
    analysis->is_uniform = 1;
    analysis->reliability_warning = 0;
    analysis->uniformity_score = 1.0;
    return analysis;
}
```

`calloc()` zaručí h_min=h_max=h_avg=0. `print_grid_analysis()` v `-g` módu pak vytiskne `h_min = 0.000000e+00, h_max = 0.000000e+00, h_avg = 0.000000e+00, CV = 0.000` — nikoli "single point, no spacing analysis". Volající (smooth.c:557) má kontrolu `n < sp` před tím, takže pro normální runy n=1 nedoběhne — ale `-g` na single-point data vyprodukuje miswarující výstup.

**Recommendation:** v print_grid_analysis přidat `if (analysis->n_points < 2) { printf("Single point, spacing analysis not applicable\n"); return; }`. Nebo nastavit `h_min = h_max = h_avg = NAN` a tisknout "n/a" pro NaN. (Druhá varianta cleaner.)

### B13. `errno` není explicitně zerovaný před `atof` — `smooth.c:120, 131`

```c
case 'l':
    if (strcmp(optarg, "auto") == 0) {
        auto_lambda = 1;
    } else {
        lambda = atof(optarg);  /* atof = strtod with no error check */
        if (lambda < 0) { ... }
    }
```

`atof("foo")` vrací 0.0, projde `lambda < 0 → false`, jde dál s lambda=0 (extreme: žádná regularizace). Stejně `-f abc` → fc=0, propadne `<= CUTOFF_FREQ_MIN` (=0) check → error. Ale pro `-l` se 0 přijme tiše. Stejný motiv jako B13 z v5.11.8 pro `-n`/`-p`.

**Recommendation:** dát do bagážů s B13 z v5.11.8 a vyřešit jediným strtol/strtod helperem (`parse_int_arg(opt, name)`, `parse_double_arg(opt, name, min, max)`) v smooth.c.

### B14. Inconsistent error reporting strategy — různé moduly, různé konvence

Současný stav:

| Modul        | Error path                                         | Diagnostika                                |
|--------------|----------------------------------------------------|--------------------------------------------|
| polyfit      | Vrátí NULL, stderr `Error:`, per-window stderr `Warning:`, stderr `Note:` (cond/rank) | Mixed `Warning:`/`Note:` |
| savgol       | Vrátí NULL, stderr `Error:` (UPPERCASE jen pro CV reject), printf `# Savitzky-Golay:`  | Recovery hints v stderr  |
| tikhonov     | Vrátí NULL, stderr `Error:`, `Warning:`, printf `# WARNING:`, printf `# Note:`, printf `# Optimal λ:` | 5 různých prefixů |
| butterworth  | Vrátí NULL, stderr `ERROR:`/`Warning:`, printf `# WARNING:`/`# Auto cutoff:` | Konvence dokumentovaná v hlavičce |
| grid_analysis| Vrátí NULL, stderr `Error:`, printf `WARNING:` (přes verbose) | Mírně inkonzistentní      |
| timestamp    | Vrátí NULL/-1, žádné stderr (tichý fail)           | Reportuje až caller       |

**Recommendation:** `butterworth.c` má dnes nejčistší konvenci (stdout `#` = data-relevant, stderr `Warning:` = operational, stderr `ERROR:` = fatal). Postupně rozšířit na ostatní moduly. Hotová odbočka po B1.

### B15. `MAX_LINE = 4096` a `MAX_COLS = 100` jsou `#define` v `smooth.c`, ne v hlavičce — `smooth.c:27–28`

Po B9 řešení z v5.11.8 jsou hodnoty pojmenované, dobré. Ale jsou `#define`-d v `smooth.c`, ne v `parser.h` (která ještě neexistuje). Pokud někdo bude chtít upravit limit, musí editovat `smooth.c` a překompilovat. Chybová zpráva to říká ("Increase MAX_LINE in smooth.c"), takže je to konzistentní, ale nevhodné pro CLI tool — uživatel by čekal `--max-line` opci nebo aspoň environment variable.

**Recommendation:** ne kritické. Pokud někdo někdy potká `MAX_COLS=100` jako bottleneck, dořešit. Drobnost.

---

## C. Drobnosti

- **C1.** `grid_analysis.c:14–17`: lokální `THRESH_CV_*` konstanty se neexponují — interní detail. Ale CLAUDE.md říká, že 0.01 je "definice prakticky uniformní mřížky". B1 z v5.11.8 odložil centralizaci; mezitím by minimálně `grid_analysis.h` mohl exportovat `GRID_CV_UNIFORM` jako jediný shared symbol, bez toho aby se sahlo na method-specific prahy.
- **C2.** `polyfit.c:159`: `int i, j, k;` deklarace na úrovni funkce, ale dál se v každé `for` smyčce stejně používá `for (i = ...)` (legacy C89 styl). Když je už projekt na C99, lze `for (int i = ...)` všude. Nemíchat. (C5 z předchozího auditu, dead comment, podobný motiv.)
- **C3.** `tikhonov.c:308–311`: comment "LDA = KD + 1 = 2" obsahuje chybu — KD=2, takže LDAB = KD+1 = 3 (a kód to dělá `int ldab = kd + 1;` správně). Komentář zaostává.
- **C4.** `smooth.c:677`: `"# Data smooth - aprox. pol. %ddg from %d points of moving window (least square)"` — typo "aprox." místo "approx.", "least square" místo "least squares", "%ddg" místo "%dth degree". Stejný řetězec se objevuje v komentáři (smooth.c:1) z v3.x. Drobná kosmetika, ale jde do user-facing výstupu.
- **C5.** `savgol.c:62, 109`: cs/en mix komentářů (C1 z v5.11.8, ne fixnuto). Tikhonov je čistě anglicky, savgol mixuje. Není priorita.
- **C6.** `decomment.c:11`: `#include <stdbool.h>` — používá se jen pro `bool success` v `decomment_stream`. Drobnost, ale v jinak vanilla C99 modulu volba `int success = 0/1` by ušetřila include.
- **C7.** `butterworth.c:339`: `pad_signal()` používá `if (src_idx >= n) src_idx = n - 1;` na řádku 324 — defensive clamp pro pad_len > n. `calculate_pad_length()` na řádku 69 ale **už** garantuje `pad_len <= n - 1`. Defensive code je tedy nedosažitelný. (Stejný vzorec jako B12 z v5.11.8 — defensive checks po validaci na vyšší úrovni.)
- **C8.** `tikhonov.c:402`: `double scale = lambda / h4;` se počítá vždy, ale používá se jen v `else` větvi `n > 5000`. V první větvi (`n <= 5000`) se scale ignoruje. Tichý dead path pro malé `n`.
- **C9.** `smooth.c:751–790`: `help()` má hardcoded `static char *msg[] = {...0}` array NULL-terminovaný; čistší by byl `static const char *const msg[]` (immutable strings, immutable pointers). Drobná const-correctness.
- **C10.** ~~README říká "f auto currently returns default fc = 0.2" (`README.md:282, 341`) ale od v5.11.3 implementuje skutečné Morozovo discrepancy principle (`butterworth.c:422–474`, `revision.h:67`). README/code drift. Stejně řádek 89 v README: "use `-f auto` for default" je matoucí.~~ — **FIXED v5.11.23** (README:89, 282, 341 přepsány na popis Morozova).

---

## Shrnutí podle priority

| Priorita | Položky |
|----------|---------|
| **Fix brzy** | ~~A1~~ (opraveno v5.11.23), A2 (parse_timestamp normalization), ~~C10~~ (opraveno v5.11.23) |
| **Vyčistit při příležitosti** | ~~B1~~ (opraveno v5.11.24), ~~B2~~ (opraveno v5.11.25), ~~B3~~ (opraveno v5.11.26), ~~B4~~ (opraveno v5.11.27 + v5.11.28), ~~B5~~ (opraveno v5.11.30), B6 (zveřejněná jen interní funkce), B13 (atoi → strtod, joinout s v5.11.8 B13), B14 (rozšířit butterworth konvenci jinde) |
| **Kosmetika** | B7, B8, B9, B10, B11, B12, B15, C1–C9 |

**Závěrečné hodnocení.** Kód je stabilní vědecký nástroj s velmi pečlivou matematickou správností (Tikhonov pentadiagonal, biquad cascade, SVD polyfit, SG pre-computed coeffs). Předchozí audit byl drtivou většinou vyřešen — projekt má disciplínu na regulérní cleanup. Hlavní designová témata tohoto auditu (B1 error-label konzistence, B2 const-correctness, B3 h_avg duplikace, B4 monolitický main()) byla zavřena ve v5.11.24–v5.11.28. Otevřený zbytek je drobnost: A2 (timestamp normalizace) zůstává jediný skutečný bug, B5/B6/B13/B14 jsou stylistická tech-debt-vyčistění, kosmetika v C-řadě. Nic z auditu neohrožuje korektnost vědeckých výsledků na dnešních uniformních/blízko-uniformních mřížkách.

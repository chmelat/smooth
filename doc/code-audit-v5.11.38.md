# Code audit v5.11.38

Datum: 2026-06-01. Hloubkový audit z čerstvého čtení všech devíti produkčních
modulů (~3500 řádků C: `smooth.c`, `polyfit.c`, `savgol.c`, `tikhonov.c`,
`butterworth.c`, `grid_analysis.c`, `timestamp.c`, `parser.c`, `decomment.c`)
plus odpovídající hlavičky a testy.

Projekt je po předchozích kolech (audity `code-audit-v5.11.8.md`,
`code-audit-v5.11.22.md`, série A1–A6 / B4–B9) ve velmi dobrém stavu: žádné
kritické chyby, žádné paměťové úniky (`make test-valgrind` 0/0), 111 testů
prochází. Nálezy tohoto auditu jsou drobné — robustnost vstupů a dva
degenerované numerické okraje.

V rámci této audit-session bylo rovnou opraveno:

- **savgol** (v5.11.36): hraniční bod při selhání koeficientů/alokace vrací
  NULL místo tiché náhrady surovým `y[i]`; oprava dokumentace v `savgol.h`.
- **polyfit** (v5.11.37): kontrola `info` z `dgelss` workspace query; přesnější
  popisek efektivní podmíněnosti + odstraněná mrtvá větev; oprava `poly_degree`
  rozsahu v `polyfit.h`.
- **S3 + T1** (v5.11.38): viz níže.

---

## A. Bugy

### A1. ~~Tikhonov "Data/Total ratio" tiskne `nan` při J=0~~ — `smooth.c:273–275` — **FIXED v5.11.38** (audit ID: S3)

`-l 0` je povolená hodnota (validace odmítá jen `lambda < 0`). S `lambda = 0`
řeší Tikhonov soustavu `I·u = y`, takže fit je přesný:
`data_term = reg_term = total_functional = 0`. Výpis poměrů pak počítal
`data_term / total_functional = 0/0` → výstup obsahoval `nan`:

```
# Functional J = 0 (Data: 0 + Regularization: 0)
# Data/Total ratio = nan, Regularization/Total ratio = nan
```

**Fix:** řádek s poměry se tiskne jen pokud `total_functional > 0.0`. Ověřeno
funkčně — `-l 0` už `nan` netiskne, `-l 0.1` poměry tiskne normálně.

### A2. `parse_timestamp` přijme nesmyslné kalendářní datum — `timestamp.c:60–69, 81` — **OPEN** (duplikát A2 z `code-audit-v5.11.22.md`)

Stále platí: validace kontroluje `day < 1 || day > 31`, ale 31. únor projde a
`timegm()` ho tiše normalizuje (`2025-02-31` → `2025-03-03`). Pro vědecká data
ze špatně formátovaných exportů to je tichý posun časové stopy.

**Fix:** po `timegm()` zpětně ověřit
`tm_time.tm_mday == day && tm_time.tm_mon == month - 1`; při neshodě (vstup byl
normalizován) vrátit `-1`. Záměrně neopraveno v této session (mimo dohodnutý
rozsah).

---

## B. Designová / robustnostní rozhodnutí

### B1. ~~L-curve používá `0.0` jako sentinel neplatného bodu~~ — `tikhonov.c:368–384` — **FIXED v5.11.38** (audit ID: T1)

`find_lambda_lcurve` značila body, kde `tikhonov_smooth` selhal, hodnotou
`rss_vals[i] = 0.0` a v cyklu výpočtu křivosti je přeskakovala přes
`rss_vals[i] == 0.0`. Jenže `rss_vals[i] = log(data_term)` je legitimně `0.0`,
když `data_term == 1.0` — platný bod L-křivky se pak zahodí a posune detekovaný
roh.

**Fix:** zaveden explicitní `int *valid[]` flag (1 = bod spočítán, 0 = solver
selhal), uvolněný v `lcurve_cleanup`. Cyklus křivosti testuje `!valid[...]`
místo porovnání s `0.0`. L-curve se používá jen pro `n > 20000`, takže dopad na
běžný vstup je nulový.

### B2. ~~CLI číselné argumenty přes `atoi`/`atof` tiše polykají nečíselný vstup~~ — **FIXED v5.11.54** (= B13 v `code-audit-v5.11.22.md`, viz `code-audit-v5.11.41.md` pro popis opravy)

`-n abc` → `sp = 0` (chyceno až `sp < 3`), `-p abc` → `dp = 0` (tiše degree 0),
`-l xyz` → `lambda = 0.0`, `-f xyz` → `cutoff = 0.0` (chyceno range-checkem).
Funkční, ale bez diagnostiky chybného vstupu.

**Fix:** přejít na `strtod`/`strtol` s kontrolou `endptr` a hlásit
`ERROR: invalid value for -X: '...'`. Záměrně neopraveno (mimo rozsah).

### B3. Netrackované pracovní soubory v repu — **OPEN** (duplikát B14 z `code-audit-v5.11.22.md`)

Build artefakty (`tests/test_runner`, `test_debug_edge`) a scratch soubory
(`debug_test.c`, `test_actual_results.c`, `test_debug_edge.c`, `plot_data.py`,
`doc/tikhonov_issues_report.pdf`) leží v pracovním stromu. Buď doplnit
`.gitignore`, nebo smazat.

---

## C. Moduly bez nálezů (čistý průchod)

- **`grid_analysis.c`** — single-pass min/max/std/clusters korektní; `h_avg`
  z teleskopu = přesný průměr spacings, Besselova korekce `/(n-2)` sedí;
  non-monotonic záchyt před dělením; bezpečné stringy; čistá správa paměti.
- **`decomment.c`** — strip `#`→EOL + prázdné řádky korektní; `tmpfile()` se
  rewinduje, vlastnictví streamu jasné, error cesty uvolňují temp; stdin se
  nezavírá.
- **`parser.c`** — po auditu B4/B9 robustní: detekce přetečení řádku i počtu
  sloupců, dynamický realloc (doubling), placeholder model, korektní `fail:`
  úklid (párování realloc x/y i `timestamp_strings`).
- **`butterworth.c`** — ověřena unity DC gain per biquad (`4·gain/det = 1`),
  takže IC scaling v `apply_cascade` je platný; pole-stability check;
  odd-reflection padding s indexy v mezích (`pad_len ≤ n-1`); Morozov
  auto-cutoff; 5-bodové derivace (n≥20 ⇒ n≥5). Nejdůkladněji zauditovaný modul.

---

## Stav nálezů

| ID  | Modul        | Závažnost | Stav             |
|-----|--------------|-----------|------------------|
| S3  | smooth.c     | nízká     | FIXED v5.11.38   |
| T1  | tikhonov.c   | nízká     | FIXED v5.11.38   |
| A2  | timestamp.c  | nízká     | FIXED v5.11.43   |
| B2  | smooth.c     | nízká     | FIXED v5.11.54   |
| B3  | repo         | kosmetika | FIXED `c42056e`  |

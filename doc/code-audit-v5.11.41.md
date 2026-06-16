# Code audit v5.11.41

Datum: 2026-06-16. Hloubkový audit z čerstvého čtení všech devíti produkčních
modulů (~4470 řádků C: `smooth.c`, `parser.c`, `polyfit.c`, `savgol.c`,
`tikhonov.c`, `butterworth.c`, `grid_analysis.c`, `timestamp.c`, `decomment.c`)
plus odpovídajících hlaviček. Baseline: 113 testů prochází, `make` čistý
(`-Wall -Wextra -pedantic`, žádná varování).

Kontext předchozích kol: `code-audit-v5.11.38.md` (otevřené A2/B2/B3),
`audit-tikhonov.md` (TK1/TK2 opraveny v5.11.39, TK3–TK8 otevřené),
`butterworth-audit.md` (všechny nálezy vyřešené k v5.11.41).

Hlavní nález tohoto kola je **nový a netriviální**: tichý posun (desynchronizace)
mezi `x` a `y` v timestamp-módu, když vstup obsahuje neplatné časové razítko na
řádku, který jinak má platnou číselnou hodnotu `y`. Jádra numerických metod
(Tikhonov pentadiagonal, biquad cascade, SVD polyfit, SG konvoluce) prošla bez
nového nálezu — slabina je čistě v I/O / parsovací vrstvě.

---

## A. Bugy

### A1. ~~Timestamp-mód: tichý posun `y` vůči `x` při zahozeném časovém razítku~~ — `parser.c:273,287-292`, `timestamp.c:136-170` — **FIXED v5.11.42** [NOVÝ]

V timestamp-módu plní `parse_input` dvě paralelní pole indexovaná **řádkem**:
`timestamp_strings[k]` a `y[k]` (`parser.c:165,170`). Po dočtení se razítka
převedou voláním

```c
ts_ctx = convert_timestamps_to_relative(timestamp_strings, n, &x, &first_error_line);
```

`convert_timestamps_to_relative` ale **neplatná razítka přeskakuje a pole `x`
zhušťuje** (`timestamp.c:139-146`, `continue` na neúspěšném `parse_timestamp`),
takže výsledné `x` má délku `valid_count` = počtu *platných* razítek a obsahuje
relativní časy jen platných řádků. Pole `y` se přitom **nikdy nepředá a
nezhušťuje** — zůstane v původním pořadí řádků. Parser pak jen nastaví
`n = ts_ctx->n` (`parser.c:288`), čímž `y[]` ořízne na prvních `valid_count`
prvků.

Důsledek: jakmile `ts_ctx->errors_encountered > 0`, jsou `x[k]` a `y[k]`
**rozhozené** — `y` se posune o počet zahozených razítek, která mu předcházela,
a koncové platné body `y` vypadnou úplně.

**Empirické ověření** (build v5.11.41):

```
vstup (řádek 2 má neplatné razítko, ale platné y=20):
  2025-01-01T00:00:00 10.0
  2025-13-45T00:00:00 20.0     <- měsíc 13, den 45 => parse_timestamp vrátí -1
  2025-01-03T00:00:00 30.0
  2025-01-04T00:00:00 40.0

$ ./smooth -T -m 0 -n 3 -p 1 vstup
  ...
  2025-01-01T00:00:00    9.2857143
  2025-01-03T00:00:00    22.142857   <- razítko řádku 3 spárováno s y=20 (řádek 2!)
  2025-01-04T00:00:00    28.571429   <- y=40 (řádek 4) zcela vypadlo
```

Linearní polyfit na `y=[10,20,30]` dává 9.29 / 22.14 / 28.57; kdyby se použilo
správné `y=[10,30,40]`, vyšly by jiné hodnoty. Tedy `2025-01-03` je spárováno
s `20.0` (hodnota zahozeného řádku 2), nikoli s `30.0`.

`Warning: Skipped 1 line(s) with invalid timestamps` se sice vytiskne, ale
uživatel nemá důvod tušit, že **zbylá** data jsou přeházená (čeká oříznutí, ne
posun). Pro `-T` mód, který je explicitně cílený na „data ze špatně
formátovaných exportů" (viz A2), je neplatné razítko očekávaný failure mode, ne
exotika — proto střední až vysoká závažnost (tichá vědecká data corruption).

**Fix:** `convert_timestamps_to_relative` dostalo volitelný (nullable) parametr
`double *y_inout`, který zhušťuje v lockstepu s `x` (`y_inout[valid_count] =
y_inout[i]` ve stejné smyčce, kde se ukládá `x_temp[valid_count]`; bezpečné in
place, protože `valid_count <= i`). `parse_input` předává své `y`; ostatní
volající (testy) předávají `NULL`. Přidán regresní test
`test_convert_compacts_parallel_y` (neplatná razítka na indexech 1 a 3 → `y`
přežije jako `{10,30,50}`). Ověřeno end-to-end na reprodukčním vstupu:
`2025-01-03 → 30`, `2025-01-04 → 40` (dřív `→ 20` a bod 40 mizel). 114 testů
prochází, nula valgrind leaků.

### A2. ~~`parse_timestamp` přijme nesmyslné kalendářní datum~~ — `timestamp.c:61-69,81` — **FIXED v5.11.43** (duplikát A2 z `code-audit-v5.11.38.md` / `-v5.11.22.md`)

Validace kontrolovala `day < 1 || day > 31`, ale 31. únor procházel a
`timegm()` ho tiše normalizoval (`2025-02-31` → `2025-03-03`) → tichý posun
časové stopy. **Fix:** po `timegm()` (které `struct tm` normalizuje in place) se
ověří, že normalizovaná pole `tm_year/tm_mon/tm_mday` stále sedí na vstup; při
neshodě `return -1`. Žádná vlastní logika přestupných roků / délek měsíců —
obstará `timegm()`. Platné přestupné dny (`2024-02-29`) projdou. Přidán test
`parse_timestamp_nonexistent_date` (Feb 31 / Apr 31 / Feb 29 v nepřestupném roce
odmítnuto, `2024-02-29` přijato). Ověřeno CLI: `2024-02-31` odmítnuto, `2024-02-29`
projde. 115 testů, nula valgrind leaků.

Vztah k A1: data, jejichž razítko je **mimo** rozsah (a tedy odmítnuto), spustí
posun z A1; data, jejichž razítko se tiše **normalizovalo** (A2), naopak prošla
s posunutou stopou bez varování. Obě cesty byly opraveny společně (A1 v5.11.42,
A2 v5.11.43): neplatná razítka se teď konzistentně odmítají a `y` zůstává
zarovnané s `x`.

---

## B. Designová / robustnostní rozhodnutí (otevřené, přenesené)

### B2. CLI číselné argumenty přes `atoi`/`atof` tiše polykají nečíselný vstup — `smooth.c:82,85,117,128,140,143` — **OPEN** (= B2 v5.11.38 / B13 v5.11.22)

`-n abc` → `sp=0` (chyceno až `sp<3`), `-p abc` → `dp=0` (tiše degree 0),
`-l xyz` → `lambda=0.0`, `-f xyz` → `cutoff=0.0` (chyceno range-checkem),
`-k abc` → sloupec 0 (chyceno `<1`). Funkční, ale bez diagnostiky chybného
vstupu. **Fix:** `strtod`/`strtol` + kontrola `endptr`, hlásit
`ERROR: invalid value for -X: '...'`. Parser už tenhle vzorec používá správně
(`strtod` + `endptr`), CLI ne — nekonzistence.

### B8. Tikhonov per-lambda GCV log jde vždy na stdout — `tikhonov.c:324,456,479,495` — **OPEN** (= B8 v5.11.22)

`compute_gcv_score_robust` se volá s `verbose=1` ze všech tří míst, takže při
`-l auto` se na stdout (do hlavičky dat) vysype 13–21 řádků
`# λ=…: J=…, RSS=…, tr(H)=…, GCV=…`. Stále spíš debug log než data-relevantní
hlavička. **Fix:** přidat skutečný `-v` přepínač a v defaultu vypsat jen souhrn
(`# GCV: N kandidátů vyhodnoceno, optimální λ=…`).

### B9. UTF-8 `λ` ve stdout zprávách — `tikhonov.c:324,509,514` (+ 311,443-447,464) — **OPEN** (= B9 v5.11.22)

Stdout výstup (součást uložených dat) obsahuje `# λ=…`, `# Optimal λ:`,
`# WARNING: optimal λ …`. Zbytek programu je ASCII. CLAUDE.md má pro README
omezení na DejaVu glyfy; stejné zdůvodnění platí pro datovou hlavičku, která
může projít LaTeX/ASCII pipelinou. **Fix:** `# lambda=…`. Triviální, ale dotýká
se více řádků a uživatelského výstupu, proto vědomě odloženo mimo audit.

### B3. Netrackované pracovní soubory v repu — **OPEN** (= B3 v5.11.38 / B14 v5.11.22)

`git status` ukazuje ~30 untracked souborů v pracovním stromu: build artefakt
`tests/test_runner`, scratch data (`all.dat`, `pt.dat`, `t1.dat`, `test_*.dat`,
`example.dat`, …), benchmark skripty (`benchmark*.py`, `benchmark_lib*`,
`plot_*.py`, `filter*.py`), `pdf_from_md_tmp.sh`, generátory v `tests/`. Buď
doplnit `.gitignore`, nebo přesunout/smazat. Roste od minula.

### TK3–TK8. Tikhonov auto-lambda — **OPEN**, viz `audit-tikhonov.md`

Znovu ověřeno, že stav odpovídá `audit-tikhonov.md` (z v5.11.39). Stručně:
- **TK3** (střední): L-curve detekce rohu (jen `n>20000`) je trojnásobně křehká
  — neorientovaná křivost (`tikhonov.c:400`), index-parametrizované diference na
  ne-log-uniformní mřížce (`:452`), tichý fallback `best_idx=n_lambda/2`
  (`:382`).
- **TK4** (nízká): GCV trace je vlastní-číselná aproximace i na uniformní mřížce
  (přirozené BC ≠ Dirichlet) — `tikhonov.c:295-308`.
- **TK5** (nízká): výstupní derivace jen 1. řádu na non-uniformní mřížce
  (centrální diference `tikhonov.c:103-108`); modul si jinak na non-uniformní
  korektnosti zakládá.
- **TK6** (nízká): doc drift v `tikhonov.h` — ověřeno, stále lže: příklad
  `tikhonov_smooth(x,y,n,0.1)` a `find_optimal_lambda_gcv(x,y,n)` bez `grid_info`
  (`:101,104-105`) se nezkompiluje; „Search range: 1e-6 to 1e0" (`:75`) vs.
  skutečných `1e-8` (`tikhonov.c:430`); „Warns if regularization term dominates"
  (`:77`) v kódu není; „used for discretization method selection" (`:36`) je
  relikt z doby před v5.11.34. **Doporučeno opravit první — triviální a hlavička
  aktivně lže.**
- **TK7** (kosmetika): redundantní solvery v `n>20000` cestě (GCV + L-curve řeší
  stejných 12 λ dvakrát).
- **TK8** (kosmetika): refinement gate `if (n <= 5000)` (`tikhonov.c:488`)
  ztratil po TK1 opodstatnění — trace je teď přesný pro všechna `n`, ale pro
  `5000 < n <= 20000` se refinement pořád nedělá.

---

## C. Moduly / cesty bez nového nálezu (čistý průchod)

- **`butterworth.c`** — re-verifikováno proti `butterworth-audit.md` (v5.11.41):
  bilineární koeficienty + jednotkový DC zisk per biquad, IC přes Cramerovo
  pravidlo s varováním na singularitu (`:241`), fc-adaptivní padding
  `max(14, ceil(5/(1-r_max)))` s ořezem na `n-1` a `# WARNING` při ořezu
  (`:76-94,624-631`), design + pole-check **před** alokací, 5-bodové stencily,
  Morozov fallback na nejslabší kandidát. Indexace paddingu v mezích. Bez
  nového nálezu.
- **`polyfit.c`** — SVD per okno (`dgelss`, `rcond=1e-10`), sledování nejhoršího
  efektivního čísla podmíněnosti a rank-deficience, polynomiální extrapolace na
  okrajích, fallback na surové `y[i]` při selhání SVD s hlášením frekvence.
  `sing_vals[effective_rank-1]` (`:292`) je bezpečné: matice má vždy sloupec
  jedniček (x striktně rostoucí), takže `s[0]>0` a `rank>=1`. Bez nálezu.
- **`savgol.c`** — universální konvoluční koeficienty pre-computed jednou,
  CV-reject nad 0.05, hraniční okna alokují/počítají per-point a při selhání
  vrací NULL (v5.11.36). Indexy oken v mezích. Bez nálezu.
- **`grid_analysis.c`** — single-pass min/max/std/clustery, Besselova korekce
  `/(n-2)`, non-monotonic záchyt před dělením, bezpečné stringy (`append_warning`
  s mezí), `reliability_warning` text se tiskne na všech úrovních verbose
  (oprava A1 v5.11.23). Bez nálezu.
- **`decomment.c`** — strip `#`→EOL + prázdné řádky korektní; `tmpfile()`
  rewind, vlastnictví streamu jasné, error cesty uvolňují temp; stdin se
  nezavírá. Bez nálezu.
- **`parser.c` (normální mód)** — detekce přetečení řádku i sloupců, dynamický
  realloc, placeholder model, `strtod`+`endptr`+`isnan/isinf` validace, korektní
  `fail:` úklid. Bez nálezu. (Timestamp-mód má A1.)
- **`timestamp.c` (`free_*`, alloc cesty)** — všechny error cesty uvolňují
  `original_timestamps[]`, `x_temp`, `ctx`; `free_timestamp_context` NULL-safe.
  Bug A1 je v *kontraktu* funkce (nezhušťuje `y`), ne v paměťové správě.

### Drobnosti (kosmetika, nízká priorita)

- **C1.** `smooth.c:351`: výstupní hlavička `"# Data smooth - aprox. pol. %ddg …
  (least square)"` — překlepy „aprox." / „%ddg" / „least square" jdou do
  user-facing výstupu (= C4 v5.11.22, stále neopraveno).
- **C2.** `smooth.c:426`: `help()` má `static char *msg[]` — čistší
  `static const char *const msg[]` (= C9 v5.11.22).
- **C3.** `smooth.c:225`: gate `n < sp && method != METHOD_TIKHONOV` používá pro
  Butterworth polyfit-ovo `sp` (default 5), ač Butterworth `-n` nepoužívá a má
  vlastní limit 20. Pro `n` mezi 5 a 19 projde tahle brána a odmítne až
  `butterworth_filtfilt` („Need at least 20"). Funkčně OK, jen dvě různé brány.
- **C4.** `smooth.c`: normální mód s 0 platnými řádky hlásí „Need more data
  (n<5)" místo jasného „no valid data" (timestamp mód má explicitní hlášku
  `parser.c:267`). Marginální UX.

---

## Stav nálezů

| ID  | Modul                 | Závažnost    | Stav                         |
|-----|-----------------------|--------------|------------------------------|
| A1  | parser.c/timestamp.c  | střední–vys. | **FIXED v5.11.42 (NOVÝ)**    |
| A2  | timestamp.c           | nízká        | **FIXED v5.11.43**           |
| B2  | smooth.c              | nízká        | OPEN (= v5.11.38 B2)         |
| B3  | repo                  | kosmetika    | OPEN (= v5.11.38 B3)         |
| B8  | tikhonov.c            | nízká        | OPEN (= v5.11.22 B8)         |
| B9  | tikhonov.c            | nízká        | OPEN (= v5.11.22 B9)         |
| TK3 | tikhonov.c (L-curve)  | střední      | OPEN (viz audit-tikhonov.md) |
| TK4 | tikhonov.c (trace)    | nízká        | OPEN                         |
| TK5 | tikhonov.c (derivace) | nízká        | OPEN                         |
| TK6 | tikhonov.h (doc)      | nízká        | OPEN                         |
| TK7 | tikhonov.c (redund.)  | kosmetika    | OPEN                         |
| TK8 | tikhonov.c (gate)     | kosmetika    | OPEN                         |
| C1–C4 | smooth.c            | kosmetika    | OPEN                         |

**Doporučené pořadí oprav:** ~~A1~~ (opraveno v5.11.42), ~~A2~~ (opraveno
v5.11.43), pak TK6 (triviální, hlavička lže), pak B2 (jednotný strtod/strtol
helper v CLI). Zbytek je tech-debt / kosmetika beze změny korektnosti na platném
vstupu.

**Závěr.** Po sérii kol (v5.11.8 → v5.11.41) je numerické jádro velmi pevné a
butterworth/tikhonov audity vyčerpaly své moduly. Tento audit posouvá pozornost
k I/O vrstvě a odhalil jeden reálný bug (A1) v dosud nejméně testované cestě —
timestamp-mód s neplatnými razítky — který byl v rámci této session opraven
(v5.11.42). Vše ostatní jsou přenesené otevřené položky nebo kosmetika. Na
čistém, platném vstupu (uniformní/blízko-uniformní mřížka, korektní razítka)
program produkuje korektní výsledky.

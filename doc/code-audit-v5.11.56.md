# Code audit v5.11.56

Datum: 2026-09-30. Kritické čtení všech osmi produkčních modulů (`smooth.c`,
`parser.c`, `timestamp.c`, `grid_analysis.c`, `polyfit.c`, `savgol.c`,
`tikhonov.c`, `butterworth.c`) plus hlaviček, README a testů. Každé podezření
bylo ověřeno empiricky na buildu v5.11.56, případně proti referenci
scipy 1.17 / numpy 2.4.

Baseline: `make` čistý, 138 testů prochází.

Mimo rozsah: over-engineering a úklid stromu (`plans/`, `*~`, `README.pdf`,
duplicitní helpery) — to pokrývá `ponytail-audit-v5.11.56.md` a zde se
neopakuje.

Hlavní nálezy: **tichá ztráta přesnosti na výstupu** (A1) a **numericky
neškálované polynomiální fity** v polyfit i savgol (A2, A3). Žádnou z chyb
A1–A6 současná testovací sada nezachytí.

**Status:** A1 FIXED v5.11.57; A2, A3 FIXED v5.11.58; ostatní OPEN.

---

## Přehled

| #  | Závažnost | Místo                                   | Nález |
|----|-----------|-----------------------------------------|-------|
| A1 | vysoká    | `smooth.c:447-455`                      | ~~výstup `%12.8lG` ořezává x i y na 8 platných číslic~~ **FIXED v5.11.57** |
| A2 | vysoká    | `polyfit.c:79-95`                       | ~~výsledek závisí na jednotkách x (neškálovaná Vandermondova matice)~~ **FIXED v5.11.58** |
| A3 | vysoká    | `savgol.c:120-134`                      | ~~okrajová (asymetrická) okna při velkém okně / stupni → chybné koeficienty~~ **FIXED v5.11.58** |
| A4 | střední   | `timestamp.c:28-42`                     | offset časového pásma i koncové smetí tiše ignorovány |
| A5 | střední   | `timestamp.c:149`, `grid_analysis.c:93` | chybová hlášení uvádějí index, ne řádek souboru |
| A6 | střední   | `parser.c:148-169`                      | hlavička v `-T` módu je fatální chyba |
| B1 | střední   | `tikhonov.c:374-376`                    | rozsah GCV hledání λ nezávisí na n |
| B2 | nízká     | `butterworth.c:471`                     | auto-cutoff má pevné minimum 0.02, bez varování |
| B3 | —         | `tikhonov.c:316-329`                    | aproximace stopy na nerovnoměrné mřížce — **vyvráceno** |
| C1 | doc       | README:241-242, 1219; `tikhonov.c:372,434` | nepravdivé „λ škáluje s amplitudou y" |
| C2 | doc       | `tikhonov.h:75`                         | „13-point" sweep, kód má 21 |
| C3 | nízká     | `smooth.c:225-237`                      | `-n`/`-p` validace i pro metody, které je nepoužívají |
| D1 | testy     | `tests/`                                | chybějící regresní testy pro A1–A6 |
| D2 | testy     | `tests/test_parser.c:141 ...`           | pevné cesty v `/tmp` |

---

## A. Bugy

### A1. ~~Výstup ořezává x i y na 8 platných číslic~~ — `smooth.c:447-455` — **FIXED v5.11.57**

`print_result()` tiskne všechny sloupce formátem `%12.8lG`. Cokoli s více než
8 platnými číslicemi se tiše zaokrouhlí — unixový čas, tlak v Pa, GPS
souřadnice, data s velkým offsetem.

```
$ seq 0 9 | awk '{printf "%d %.3f\n", 1700000000+$1, 101325.123+$1*0.001}' \
    | ./smooth -m 0 -n 3 -p 1 | tail -4
     1.7E+09    101325.13
     1.7E+09    101325.13
     1.7E+09    101325.13
     1.7E+09    101325.13
```

Všech deset řádků má stejné x; kroky 0.001 v y zmizely. Následný konzument
dostane duplicitní x a zkreslené y — tichá korupce dat.

**Fix (v5.11.57):** `%12.8lG` → `%12.15lG` ve všech čtyřech `printf` řádků
dat. 15 = `DBL_DIG`: desetinný vstup s ≤ 15 platnými číslicemi se vytiskne
přesně, bez šumu, který by přidal `%.17G`; šířka 12 zůstává jako minimum.
Diagnostické hlavičky `# ...` (`%.6lG`) beze změny. Regresní test
`test_parser_output_keeps_full_precision` (na starém formátu selhal:
`Expected 1.70000001e+09 Was 1.7e+09`). Reprodukce výše nyní dává
x = `1700000000..1700000009`, y = `101325.123 .. 101325.132`.

### A2. ~~Polyfit závisí na jednotkách x~~ — `polyfit.c:79-95` — **FIXED v5.11.58**

`build_vandermonde()` centruje (`dx = x - x_center`), ale neškáluje. Sloupce
`dx^j` se proto liší o mnoho řádů podle jednotek x a absolutní
`SVD_RCOND = 1e-10` pak odřízne sloupce, které nejsou skutečně lineárně
závislé. Výsledek fitu tak závisí na tom, jestli je x v sekundách, minutách
nebo indexech.

Reprodukce: n=400, `y = sin(2πi/80) + 0.1·N(0,1)`, amplituda signálu 1.

| x                     | `-n 15 -p 6`: max odchylka vs. x = index |
|-----------------------|------------------------------------------|
| `i·60` (minuty v s)   | **1.08** |
| `i·1e-4`              | 0.12 |

A i s x = index. **Korekce (v5.11.58):** původní tabulka zde byla měřená proti
`scipy.signal.savgol_filter`, který je pro velká okna sám nepřesný (relativní
chyba koeficientů 3e-3 při 51/10, ~1 od 101/12). Přeměřeno proti referenci
v Legendreově bázi na souřadnicích škálovaných do [-1, 1] (podmíněnost ≈ 1),
relativní chyba koeficientů:

| `-n`, `-p` | okno       | polyfit, h=1 | polyfit, h=60 |
|------------|------------|--------------|---------------|
| 15, 6      | symetrické | 2e-12        | **1**         |
| 51, 10     | symetrické | **1**        | **1**         |
| 51, 10     | okrajové   | 6e-2         | **1**         |
| 101, 12    | symetrické | **1**        | **1**         |
| 201, 12    | okrajové   | 3e-1         | **1**         |

Program k tomu vypíše jen `Note: ... rank-deficient. SVD regularization active.`
na stderr. Nejvíc zasažený je `-T` mód, kde je x vždy v sekundách (minutová
data → dx až stovky).

**Fix (v5.11.58):** `build_vandermonde()` dostal parametr `scale`; fit běží na
`t = (x - x_i) / s`, `s = max(x_i - x_{i-off}, x_{i+off} - x_i)`, takže
`t ∈ [-1, 1]`. Hodnota `c_0`, derivace `c_1 / s`; okrajová extrapolace volá
`evaluate_polynomial(..., dx / s, ...)` a derivaci dělí `s`. `dgelss` a
`SVD_RCOND` beze změny — hlášení o podmíněnosti teď měří skutečnou podmíněnost
a falešné „100 % rank-deficient" zmizelo. Chyba výstupu proti referenci
(n=400, šum 0.1): ≤ 2e-13 pro 15/6 až 201/12, identická pro x a 60·x.
Regresní test `test_polyfit_invariant_to_x_units` (na starém kódu selhal).

### A3. ~~Savgol: okrajová okna při velkém okně/stupni~~ — `savgol.c:120-134` — **FIXED v5.11.58**

`savgol_coefficients()` staví normální rovnice z momentů celočíselných pozic
`a[i] = Σ j^i` (pro okno 101, p=12 je `Σ j^24 ≈ 1e41`).

**Korekce (v5.11.58):** nález byl užší, než původně uvedeno — původní čísla
(2.7e-3 při 51/10, 0.99 při 101/12) byla měřená proti nepřesnému scipy (viz
A2). Proti Legendreově referenci jsou **vnitřní (symetrické) koeficienty v
pořádku** (≤ 1e-9); selhávají jen **okrajová asymetrická okna**, kde pozice
`j ∈ [0, w-1]` nejsou ani centrované, ani škálované (momentová matice typu
Hilbert):

| `-n`, `-p` | okno       | chyba savgol |
|------------|------------|--------------|
| 51, 10     | symetrické | 1e-10 |
| 51, 10     | okrajové   | 3e-5 |
| 101, 12    | symetrické | 9e-10 |
| 101, 12    | okrajové   | **7e-2** |
| 201, 12    | okrajové   | **3** |

Na reálných datech při `-n 101 -p 12` navíc `dposv` na pravém okraji selže
úplně (`info = 13`, „Savitzky-Golay smoothing failed!").

**Fix (v5.11.58):** pozice centrované a škálované, `u = (j - m) / d`,
`m = (nr - nl)/2`, `d = (nl + nr)/2`, takže `u ∈ [-1, 1]`; pravá strana je
monomiální řádek (resp. jeho derivace) v cílovém bodě `u0 = -m/d` místo
jednotkového vektoru; derivační koeficienty děleny `d` (zpět na jednotky
indexu, takže `/ h_avg` v `savgol_smooth()` platí dál). `dposv` a normální
rovnice zůstávají. Chyba výstupu proti referenci: ≤ 5e-10 až do 201/12.

**Poznámka k A2 + A3:** na uniformní mřížce jsou polyfit a savgol
matematicky totožné (asymetrické SG okno na okraji = extrapolace polynomu
prvního okna). Po opravě obou musí dávat stejný výstup — přirozený
vzájemný test (viz D1). Zavedeno ve v5.11.58 jako
`test_savgol_polyfit_exact_on_cubic_wide_window` (okno 101, p=12, kubika
reprodukovaná oběma metodami ve všech bodech včetně okrajů).

### A4. Offsety časových pásem a koncové smetí tiše ignorovány — `timestamp.c:28-42`

`sscanf("%d-%d-%d%c%d:%d:%d%n")` skončí za sekundami a zbytek řetězce se
kontroluje jen na `.` (zlomky sekund). Všechno ostatní projde.

Přechod na letní čas s explicitními offsety:

```
2025-03-30T00:00:00+01:00 1
2025-03-30T01:00:00+01:00 2
2025-03-30T03:00:00+02:00 3     <- v UTC je to 01:00, tj. +1 h
2025-03-30T04:00:00+02:00 4
2025-03-30T05:00:00+02:00 5

$ ./smooth -T -m 0 -n 3 -p 1 -d tz.dat
2025-03-30T01:00:00+01:00    1.7857143 0.00017857143
2025-03-30T03:00:00+02:00    3.2142857 0.00017857143
2025-03-30T04:00:00+02:00            4 0.00027777778
```

V UTC je řada rovnoměrná po 1 h; program vidí falešnou 2h mezeru a derivace
se mění z 0.000179 na 0.000278. Hlavička `timestamp.h` sice říká
„No timezone support", ale explicitně zapsaný offset se tiše zahodí místo
odmítnutí.

Projde i koncové smetí: `2025-01-01T00:00:00garbage`, `2025-01-01T00:00:01.5.5`.

**Fix:** za sekundami (a volitelnými zlomky) přijmout jen `\0`, `Z` nebo
`±HH:MM` / `±HHMM`; offset odečíst od epochy. Cokoli jiného → `-1`.

### A5. Chybová hlášení uvádějí index, ne řádek souboru — `timestamp.c:149`, `grid_analysis.c:93`

`convert_timestamps_to_relative()` nastavuje `*first_error_line = i + 1`,
kde `i` je index mezi řádky, které parser přijal — ne číslo řádku v souboru.
Komentáře, prázdné řádky a přeskočené řádky posun zvětšují.

```
neplatné razítko 2025-02-30 je na řádku 7 souboru (před ním 3 řádky
komentářů/prázdné)
$ ./smooth -T -m 0 ts1.dat
Warning: Skipped 1 line(s) with invalid timestamps (first error at line 4)
```

Stejně `analyze_grid()`: `ERROR: Non-monotonic x data at index %d` uvádí
index v poli x (v obou módech), takže uživatel nemá jak najít vadný řádek
ve vstupu.

**Fix:** parser si k akceptovaným řádkům ukládá `line_number` (paralelní
pole, kompaktované v lockstepu jako `y_inout`); hlášení z timestamp i
monotonicity pak tisknou skutečný řádek.

### A6. Hlavička v `-T` módu je fatální — `parser.c:148-169`

Formát razítka se určuje podle `strchr(token, 'T')`. Token bez `T` se
považuje za datum ve formátu s mezerou a spotřebuje dva tokeny; kontrola
počtu sloupců pak proběhne dřív než validace razítka.

| první řádek     | `-T` mód | číselný mód |
|-----------------|----------|-------------|
| `date value`    | `ERROR: Line 1 has insufficient columns` (exit 1) | — |
| `time value`    | `ERROR: Line 1 has insufficient columns` (exit 1) | — |
| `Time Value`    | přeskočeno (`T` v „Time") | — |
| `x y`           | — | přeskočeno |

Chování závisí na velikosti písmen a je nekonzistentní s číselným módem,
kde se nečíselná hlavička vždy jen přeskočí.

**Fix:** řádek, jehož sestavené razítko neprojde `parse_timestamp()`, počítat
jako `skipped_malformed_ts` a přeskočit ještě před kontrolou počtu sloupců
pro y.

---

## B. Numerika a návrh

### B1. Rozsah GCV hledání λ nezávisí na n — `tikhonov.c:374-376`

Rozsah `[1e-8, 1e6] · h³` je invariantní vůči měřítku mřížky, ale ne vůči
délce záznamu. Mód s úhlovou frekvencí θ je potlačen, když
`λ · 16 sin⁴(θ/2) / h³ ≈ 1`:

- horní mez `1e6·h³` → perioda řezu ≈ 200 vzorků; pomalejší signál už
  vyhladit nejde,
- `λ/h³ ≲ 1e-2` → nevyhlazuje vůbec nic (řez nad Nyquistem), tj. spodních
  ~6 dekád (~9 z 21 bodů sweepu) je zbytečných.

Reprodukce: n=20000, `y = sin(2πi/10000) + 0.3·N(0,1)` (perioda 10000 vzorků).

```
# WARNING: optimal lambda = 1.000e+06 lies at the edge of the search range
RMSE vs. pravda:  GCV (λ=1e6) 0.029   |  -l 1e8  0.016  |  -l 1e10  0.0093
```

Varování se vypíše (správně), ale výchozí výsledek je 3× horší, než je
dosažitelné.

**Návrh:** `λ ∈ [1e-3, c · (n/π)⁴] · h³` — horní mez odpovídá potlačení i
nejnižšího netriviálního módu θ = π/n. **Pozor:** pro n=20000 dosahuje λK
~1e15 a `dpbsv` na `I + λK` ztrácí přesnost (κ·ε); před rozšířením ověřit
podmíněnost, případně horní mez omezit.

### B2. Auto-cutoff Butterwortha nikdy nejde pod 0.02 — `butterworth.c:471`

Kandidáti `{0.02, 0.05, 0.1, 0.2, 0.35, 0.5}` jsou pevní. Na stejných datech
jako B1 vybere auto první kandidát 0.02:

| fc          | RMSE celkem | RMSE vnitřek (bez 1000 bodů na okrajích) |
|-------------|-------------|------------------------------------------|
| 0.02 (auto) | 0.033       | ~0.022 |
| 0.005       | 0.029       | 0.015 |
| 0.002       | 0.041       | **0.011** |

Tikhonov v analogické situaci varuje o hraně rozsahu, Butterworth ne.

**Návrh:** vypsat varování, když discrepancy splní už nejmenší kandidát
(optimum může ležet níž); volitelně rozšířit kandidáty dolů podle n
(s ohledem na délku paddingu).

### B3. Aproximace stopy na nerovnoměrné mřížce — vyvráceno

Podezření: `compute_gcv_score_robust()` počítá tr(H) z modelu vlastních čísel
uniformní mřížky, na nerovnoměrné mřížce by proto GCV mohlo volit špatné λ.

Ověřeno v numpy s přesnou stopou `tr((I + λK)⁻¹)` na stejné matici K, n=400,
57 bodů λ:

| mřížka             | CV   | max rel. chyba stopy | λ/h³ přesné GCV | λ/h³ aprox. GCV |
|--------------------|------|----------------------|-----------------|-----------------|
| uniformní          | 0.00 | 0.09 | 10  | 10  |
| náhodná            | 0.48 | 0.09 | 5.6 | 5.6 |
| gradovaná 1:20     | 0.52 | 0.14 | 5.6 | 5.6 |
| dva režimy 1:10    | 0.82 | 0.43 | 3.2 | 3.2 |

Stopa se liší až o 43 %, ale zvolené λ je ve všech případech identické
s přesným GCV. Beze změny.

---

## C. Dokumentace

### C1. „λ škáluje s amplitudou y" je nepravda

README:241-242 („it scales with $h^3$ and with the squared amplitude of $y$"),
README:1219 („The one scale the range does not model is the amplitude of $y$"),
`tikhonov.c:372-373` a `:434`.

Minimalizátor `‖y − u‖² + λ∫u''²` je v y lineární — oba členy škálují s y²,
takže λ na amplitudě nezávisí. Ověřeno: y·1000 → GCV zvolí stejné
λ = 316.2. λ závisí na h³ a na poměru signál/šum a spektru signálu, ne na
amplitudě. Tvrzení svádí k ručnímu ladění λ podle amplitudy dat.

### C2. `tikhonov.h:75` — „13-point log-spaced grid search"

Kód (`tikhonov.c:420`) i README:1191 mají 21 bodů.

### C3. Validace nepoužívaných parametrů — `smooth.c:225-237`

`-n` se validuje (liché ≥ 3) i pro Tikhonov a Butterworth, které ho
nepoužívají: `smooth -m 2 -n 4` skončí chybou. ~~Varování „High polynomial
degree" se tiskne i pro Tikhonov, který `-p` nepoužívá.~~ Varování bylo ve
v5.11.58 odstraněno úplně — po opravě A2/A3 už neoznačovalo žádnou
nestabilitu. Úvodní řádek `help()` popisuje jen polyfit.

---

## D. Testy

### D1. Chybějící regresní testy

Žádná z chyb A1–A6 není testy pokryta. Nejmenší sada, která by je zachytila:

1. **A1:** vstup s x = 1.7e9 + i projde výstupem s různými x (e2e přes `popen`
   jako ostatní parser testy).
2. **A2:** polyfit na x a na 60·x dává stejné `y_smooth` (a derivaci /60).
3. **A3 + A2:** polyfit == savgol na uniformní mřížce pro `(51, 10)` a
   `(101, 12)`.
4. **A4:** `parse_timestamp("...+02:00")` buď aplikuje offset, nebo vrátí -1;
   `"...00garbage"` vrátí -1.
5. **A5:** `-T` s komentáři před vadným řádkem hlásí řádek souboru.
6. **A6:** hlavička `date value` v `-T` módu se přeskočí.

### D2. Pevné cesty v `/tmp` — `tests/test_parser.c:141` a další

Parser testy zapisují do pevných souborů (`/tmp/test_parser_iso_t.dat`, …).
Dva souběžné běhy (dva checkouty, CI matrix) si je přepíšou. `mkstemp()`
nebo cesta s PID.

---

## Co je v pořádku

- Architektura drží vlastní pravidla: žádné závislosti mezi metodami,
  jednotný vzor vlastnictví paměti, `goto`-cleanup.
- Validace CLI (`arg_int` / `arg_double`) je důkladná.
- Tikhonov: pentadiagonální Gramova matice odpovídá modelu vlastních čísel
  `16 sin⁴(θ/2)/h³`; GCV volí stejně jako přesný výpočet (B3).
- Butterworth: analytické IC přes Cramerovo pravidlo jsou správně odvozené,
  padding adaptivní k fc.
- Parser: CRLF, přetečení řádku i počtu sloupců jsou ošetřené.

## Doporučené pořadí oprav

1. ~~A1 — jeden řádek, největší dopad.~~ Hotovo ve v5.11.57.
2. ~~A2 + A3 společně, se vzájemným testem D1.3.~~ Hotovo ve v5.11.58.
3. A4, A5, A6 — timestamp/parser vrstva, jedna série.
4. B1 (s ověřením podmíněnosti), B2.
5. C1–C3, D2.

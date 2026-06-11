# Analýza implementace Butterworthova filtru

**Datum auditu:** 2026-06-10 (aktualizováno)
**Verze projektu:** smooth v5.11.39
**Auditované soubory:** `butterworth.c` (V1.4/2025-12-07), `butterworth.h`, volající část `smooth.c`

**Historie změn dokumentu:**
- 2026-06-11 (v5.11.40): jednoduché opravy nálezů z kola 2026-06-10. RESOLVED
  #3 (důkaz invariantu v kód-komentáři), #8 (varování na stderr místo tichého
  nulování IC), #16 (Morozov fallback na největší/nejslabší kandidát + warning),
  #17 (`NOISE_MAD_NORMALIZATION` = 1.6521808), #18 (odstraněn mrtvý typ
  `ButterworthCoeffs`), #19 (komentáře pad-délky a `residual_std`). Otevřený
  zůstává jen #15 (pevný padding vs. transient pro malá `fc`, střední priorita).
  111 testů, nula valgrind leaků.
- 2026-06-10 (v5.11.39): nové kolo auditu (Fable 5). Jádro znovu ověřeno jako
  korektní (bilineární koeficienty přepočteny ručně, IC vs. scipy `lfilter_zi`,
  in-place aliasing v `apply_biquad`, krajní 5-bodové stencily, meze indexů
  paddingu, orientace Morozovova výběru). Nové nálezy #15-#19: pevný padding
  nepokrývá přechodový jev pro malá fc (#15, pohlcuje #12), Morozov fallback
  jde špatným směrem (#16), konstanta `NOISE_MAD_NORMALIZATION` nesedí na
  vlastní komentář (#17), mrtvý typ `ButterworthCoeffs` (#18), drobné
  komentářové neshody (#19).
- 2026-05-31 (v5.11.35): kritická re-analýza (Opus 4.8). #13 přesunut CV check
  uniformity gridu před auto-cutoff (odmítnutí neuniformního gridu bez 6 zbytečných
  zkušebních filtfiltů a stdout šumu). #14 opraveny zavádějící komentáře u derivace
  (O(h⁴) jen na uniformním gridu; degraduje k O(CV) u CV→0.15). #3 doplněn explicitní
  důkaz invariantu `buf[0]` v cascade IC. 111 testů.
- 2026-04-18 (v5.11.7): #1 doplněna podpora derivací přes 5-bodové stencily O(h⁴)
  na vyhlazeném výstupu; `-d` flag nyní funguje pro Butterworth, 106 testů.
- 2026-04-18 (v5.11.6): #9 odstraněn nepoužívaný parametr `x`, #10 vyjasněn output
  (label „Effective sample rate" + komentář `= 1/h_avg`), #11 MB/GB formát adaptivní.
- 2026-04-18 (v5.11.5): #7 pojmenování konstanty opraveno (`CUTOFF_FREQ_INEFFECTIVE_WARN`),
  zpřesněn text varování a doplněna dokumentace parametru v `butterworth.h`.
- 2026-04-18 (v5.11.4): #6 přidán explicitní spodní limit `FC_MIN_PRACTICAL = 1e-4`.
- 2026-04-18 (v5.11.3): #2 `-f auto` implementováno přes Morozovovu discrepancy principle,
  #4 přidána kontrola stability pólů, #5 konvence `fc` sjednocena v CLAUDE.md.
- 2026-04-18 (v5.11.1): původní audit.

---

## Vyřešené problémy

### [RESOLVED v5.11.3] 2. `estimate_cutoff_frequency` — nyní plná implementace

Dřívější stub vracející konstantu `0.2` byl nahrazen Morozovovou discrepancy principle
(`butterworth.c:379-426`):

1. Odhad šumového `σ̂` z MAD druhých diferencí
   (normalizace `sqrt(6) * 0.6745 ≈ 1.6528553`).
2. Iterace přes kandidáty `{0.02, 0.05, 0.1, 0.2, 0.35, 0.5}` — vybere nejmenší `fc`
   takové, že `std(y - y_smooth) ≤ 1.1 * σ̂`.
3. Při selhání odhadu šumu / alokace fallback na `AUTO_CUTOFF_FALLBACK = 0.2`
   s informační zprávou.

Diagnostika (`# Auto cutoff: ...`) je vypisována do stdout jako součást hlavičky,
což je konzistentní s chováním ostatních metod.

**Drobná poznámka:** Parametr `x` je v nové implementaci explicitně ignorován
přes `(void)x;` s komentářem, že `fc` je normalizované k Nyquistově frekvenci.
To je korektní, ale signatura API zůstává zavádějící (viz #9 níže).

### [RESOLVED v5.11.2] 4. Kontrola stability pólů

Přidána funkce `check_pole_stability()` (`butterworth.c:119-161`) volaná po designu
(`butterworth.c:528-531`). Pro každou biquad sekci:

- komplexní póly: `|z| = sqrt(|a2|)`,
- reálné póly (defenzivně): maximum z `|(-a1 ± sqrt(disc))/2|`.

Prahy:
- `POLE_RADIUS_WARN = 0.99` → varování na stderr,
- `POLE_RADIUS_ERROR = 1.0` → tvrdá chyba, `butterworth_filtfilt` vrací NULL.

Tím je zachycen i Issue #6 (extrémně malé `fc`) i Issue #7 (extrémně velké `fc`)
v praktické rovině — pokud `fc` způsobí posun pólů k jednotkovému kruhu, filtr
odmítne výpočet nebo varuje. Explicitní spodní limit na `fc` (např. `1e-4`) by byl
stále žádoucí jako první obranná linie, ale kritičnost je nyní „nízká".

### [RESOLVED] 5. Rozpor v konvenci `fc` s CLAUDE.md

CLAUDE.md (`/home/orangepi/Lang/c/smooth/CLAUDE.md:78`) nyní uvádí
`0 < fc < 1.0`, což je konzistentní s `butterworth.c:23-24` a CLI validací.

### [RESOLVED v5.11.4] 6. Explicitní spodní limit na `fc`

Přidána konstanta `FC_MIN_PRACTICAL = 1e-4` (`butterworth.c:25`) a nový
validační blok v `butterworth_filtfilt` (za existujícím MIN/MAX range-checkem):

```c
if (fc < FC_MIN_PRACTICAL) {
    fprintf(stderr, "ERROR: Cutoff frequency too small (fc=%.4e < %.4e). "
            "Filter would be numerically ill-conditioned "
            "(poles approach unit circle). "
            "Use a larger fc or a different smoothing method.\n",
            fc, FC_MIN_PRACTICAL);
    return NULL;
}
```

Hodnota `1e-4` je bezpečně pod nejmenším kandidátem auto-selektoru (`0.02`),
takže Morozov-selektor není ovlivněn. Pro `fc < 1e-4` by pólový poloměr
vždy přesáhl `POLE_RADIUS_WARN = 0.99` — práh koresponduje s oblastí,
kde filtr stejně ztrácí numerickou přesnost.

Rozšířen test `test_butterworth_invalid_cutoff_frequency`
(`tests/test_butterworth.c`) o case `fc = 1e-5`.

### [RESOLVED v5.11.5] 7. `CUTOFF_FREQ_STABILITY_WARN` — nesprávné pojmenování

Numerická analýza odhalila, že původní pojmenování bylo zavádějící.
Pro biquad sekci s úhlem `θ = π/8` (worst-case pro stabilitu blízko Nyquistu):

| `fc` | `Wc = tan(π·fc/2)` | `a2` | pólový poloměr |
|------|--------------------|------|----------------|
| 0.95 | 12.706 | 0.5579 | **0.747** |
| 0.99 | 63.657 | 0.8904 | 0.944 |
| 0.9983 | 374.5 | 0.9805 | ~0.990 |
| 0.9999 | 6366.2 | 0.9988 | 0.999 |

Při `fc = 0.95` jsou póly na poloměru pouze **0.747** — velmi daleko od
jednotkového kruhu. Skutečná numerická nestabilita (radius > `POLE_RADIUS_WARN = 0.99`)
nastává až za `fc ≈ 0.9983`, což již plně pokrývá `check_pole_stability()`
(issue #4, vyřešeno v v5.11.2).

Původní konstanta tedy nevarovala před numerickým problémem, ale před tím,
že **filtr při fc blízko Nyquistu prakticky nic netlumí** (v pásmu propustnosti
má téměř jednotkový přenos). To je UX problém, nikoli numerický.

**Provedené změny (v5.11.5):**

- Přejmenování `CUTOFF_FREQ_STABILITY_WARN` → `CUTOFF_FREQ_INEFFECTIVE_WARN`
  (`butterworth.c:26-29`) s komentářem vysvětlujícím skutečný smysl prahu.
- Zpřesnění textu varování — nyní uživateli říká, *proč* je to problém
  („Filter passes nearly the entire spectrum unattenuated").
- Dokumentace parametru `cutoff_freq` v `butterworth.h:54-56` — popisuje
  typický užitečný rozsah (`0.01 - 0.5`) a chování nad 0.95.

Žádný hard-reject nebyl přidán: `fc ∈ (0.95, 1.0)` je numericky zcela
korektní, jen málo užitečný. Rozhodnutí o použitelnosti je ponecháno na
uživateli, jemuž warning vysvětlí důsledky.

### [RESOLVED v5.11.6] 9. Parametr `x` odstraněn z `estimate_cutoff_frequency`

Nová signatura: `double estimate_cutoff_frequency(const double *y, int n)`
(`butterworth.h:75`, `butterworth.c:383`). Volání v `butterworth_filtfilt`
aktualizováno. Žádný `(void)x;` už není potřeba.

### [RESOLVED v5.11.6] 10. Popisek `sample_rate` ve výstupu vyjasněn

V `smooth.c:579` změněno `"# Sample rate: fs = ..."` na
`"# Effective sample rate: fs = ... (= 1/h_avg)"`. Uživatel vidí explicitně,
že jde o efektivní hodnotu odvozenou z průměrného spacingu — pro mírně
neuniformní grid (CV do 0.15) je to korektní interpretace.

### [RESOLVED v5.11.6] 11. Adaptivní formát paměťového odhadu

V `butterworth.c:468-476` přidána větev: pokud `mem_estimate < 1 GB`,
formátuje se jako `"%.0f MB"`, jinak `"%.1f GB"`. Pro typický warning
práh (`n > 50M` ≈ 0.8 GB) se nyní vypíše `~800 MB` místo matoucího `0.8 GB`.

---

### [RESOLVED v5.11.7] 1. Podpora derivací (API smlouva)

`ButterworthResult` nyní obsahuje pole `double *y_deriv` (`butterworth.h:37`).
Derivace jsou počítány 5-bodovými stencily O(h⁴) z již vyhlazeného výstupu —
`compute_derivatives_5pt()` v `butterworth.c`. Volba vyššího řádu než
Tikhonovova 2-bodová central difference je odůvodněná: filtfilt výstup je
extrémně hladký (efektivní 8. řád, nulová fáze), takže typický handicap
vyšších stencilů (amplifikace šumu) odpadá. Pro uniformní krok `h = h_avg`
dává O(h⁴) řádově 1000× nižší truncation error než O(h²). Okraje pokryty
forward/backward 5-bodovými stencily, vnitřek klasickou central 5-point
formulí `(-y[i+2] + 8y[i+1] - 8y[i-1] + y[i-2])/(12h)`.

Ad-hoc varování v `smooth.c` odstraněno, `-d` flag nyní produkuje `x y_smooth
y_deriv` stejně jako ostatní metody. Testy `test_butterworth_derivative_*`
(konstanta, lineární, sinus) ověřují 5-bodovou stencil s tolerancí
odpovídající filtrovému reziduálu, ne chybě stencilu.

---

### [RESOLVED v5.11.35] 13. Pořadí auto-cutoff vs. kontrola uniformity gridu

Před opravou se `estimate_cutoff_frequency()` (Morozov) volala **před** kontrolou
`grid_info->cv > UNIFORMITY_CV_THRESHOLD`. Při `-f auto` na neuniformním gridu
to znamenalo:

1. spuštění **6 zkušebních filtfiltů** (`run_filtfilt_trial` pro každého kandidáta),
2. výpis `# Auto cutoff: noise sigma estimate = ...` na stdout,
3. **až poté** tvrdé odmítnutí `ERROR: ... non-uniform grid` na stderr.

Tedy plýtvání výpočtem a matoucí pořadí výstupu (info-hlavička na stdout
před chybou na stderr, kterou uživatel ve výsledku stejně zahodí).

**Oprava:** blok kontroly CV (reject + warnings) přesunut v `butterworth_filtfilt`
**před** sekci auto-cutoff. Neuniformní grid je nyní odmítnut okamžitě, ještě
před jakýmkoli filtrováním. Pro platné (uniformní) vstupy je chování beze změny.

### [RESOLVED v5.11.35] 14. Zavádějící komentáře u přesnosti derivace

Komentáře u `compute_derivatives_5pt()` a v hlavní funkci tvrdily, že použití
`h = h_avg` na neuniformním gridu „adds only O(CV·h²) additional error, well
within filter assumptions". To bylo příliš optimistické: aplikace ekvidistantních
stencilů na neekvidistantní uzly není konzistentní a přesnost degraduje
v nejhorším případě k **prvnímu řádu (O(CV))**, ne O(h⁴).

Vyhlazení samotné je korektní (frekvenční metoda předpokládá uniformní vzorkování
a CV je shora omezeno 0.15). Slabým místem je **derivace** blízko horní hranice CV.

**Oprava:** komentáře přepsány tak, aby pravdivě říkaly, že O(h⁴) platí jen na
uniformním gridu a u CV→0.15 je derivace pouze přibližná. Žádná změna kódu/chování —
jen oprava dokumentační lži. (Implementace neuniformních stencilů by byla větší
zásah a pro frekvenční smoother, který už uniformitu předpokládá, není nutná.)

---

## Nadále otevřené problémy

### [RESOLVED v5.11.40] 3. Kaskádovaná iniciační podmínka se spoléhá na skrytou invariantu

`butterworth.c:225-238` (funkce `apply_cascade`):
```c
double first_val = buf[0];
for (int s = 0; s < NUM_BIQUADS; s++) {
    zi[0] = zi_base[s][0] * first_val;
    zi[1] = zi_base[s][1] * first_val;
    apply_biquad(&sections[s], buf, buf, padded_len, zi);
    first_val = buf[0];  // výstup předchozího biquadu
}
```

Funkční komentář nad funkcí (ř. 220-224) nyní alespoň zmiňuje, že se
spoléhá na unity DC gain. Stále ale chybí vysvětlení, proč je to korektní
i pro první sample, kdy filtr ještě není ustálený.

**Explicitní důkaz (doplněno v re-analýze v5.11.35):** invariant `buf[0]`
plyne z definice TDF-II prvního vzorku a krokového IC. Pro vstup `x` se
škálovaným IC `zi[0] = zi_base[0]·x[0]` platí:

```
y[0] = b0·x[0] + zi[0] = (b0 + zi_base[0])·x[0]
```

Přitom `zi_base[0] = (B0+B1)/det`, `B0 = b1 − a1·b0`, `B1 = b2 − a2·b0`,
`det = 1 + a1 + a2`, takže

```
b0 + zi_base[0] = [b0(1+a1+a2) + b1+b2 − (a1+a2)b0] / (1+a1+a2)
               = (b0+b1+b2) / (1+a1+a2) = Σb/Σa = DC zisk = 1.
```

Tedy `y[0] = x[0]` **přesně**, nezávisle na zbytku signálu a na ustálenosti
filtru. První vzorek je invariantní napříč sekcemi, takže `first_val = buf[0]`
po každé sekci je stále `padded[0]`. Škálování IC per-sekci je proto
matematicky identické se škálováním všech sekcí prvním vzorkem — což je přesně
to, co dělá `scipy.signal.sosfilt_zi` (jednotkový DC zisk každé sekce → `scale`
zůstává 1.0). Implementace je tedy korektní a shoduje se se scipy referencí.

**Stav:** RESOLVED v5.11.40 — důkaz invariantu `buf[0]` (DC zisk = 1, tedy
`y[0] = x[0]` přesně) přenesen do code-komentáře nad `apply_cascade`.

### [RESOLVED v5.11.40] 8. `compute_biquad_ic` tiše nuluje při degenerate case

`butterworth.c:186-193` — komentář říká „should not happen with valid fc",
ale pokud `|det| ≤ 1e-10` nastane, uživatel dostane transient-zatížený výstup
bez varování.

**Doporučení:** Přidat `fprintf(stderr, ...)` nebo vrátit chybu nahoru.
Varování je obzvlášť důležité ve světle nového Morozov-selektoru, který může
iterovat přes extrémní kandidáty.

**Stav:** RESOLVED v5.11.40 — degenerate větev nyní vypisuje
`Warning: biquad initial-condition system is singular (det=...)` na stderr.

### 12. Padding length `3*(order+1)-1 = 14` — POHLCENO nálezem #15

`butterworth.c:54` — blízko scipy defaultu (15 pro monolitický 4. řád),
ale pro biquad by bylo konzistentnější `3 * max(len(a), len(b)) = 9`.
Rozdíl v praxi zanedbatelný.

**Aktualizace 2026-06-10:** otázka 14 vs. 9 je bezpředmětná — pro malá `fc`
jsou obě hodnoty řádově krátké (viz #15). Správný směr je opačný: padding má
růst s klesajícím `fc`, ne se zmenšovat.

### 15. Pevný padding nepokrývá přechodový jev filtru pro malá `fc`

`calculate_pad_length` (`butterworth.c:66-73`) vrací konstantu 14 nezávisle
na `fc`. Délka doznění filtru je ale ~`1/(1-r)`, kde `r` je poloměr
nejpomalejšího pólu (biquad s `theta = pi/8`), a pro malá `fc` roste bez mezí:

| `fc` | pólový poloměr `r` | doznění ~`1/(1-r)` | pad |
|------|--------------------|--------------------|-----|
| 0.05 | 0.930 | ~14 vzorků | 14 |
| 0.02 | 0.971 | ~35 vzorků | 14 |
| 0.01 | 0.986 | ~69 vzorků | 14 |
| 1e-4 (povolené minimum) | 0.99986 | ~6900 vzorků | 14 |

Kroková IC absorbuje přesně DC složku (důkaz u #3), takže buzení pochází jen
z neshody sklonu/křivosti signálu na spoji s odd-reflection paddingem — to
ale stále doznívá přes celou délku `1/(1-r)` a pro `fc < ~0.05` prosakuje za
padding do výsledku jako okrajový artefakt. `FC_MIN_PRACTICAL = 1e-4` hlídá
jen podmíněnost, ne délku transientu. scipy má stejně krátký default, ale
umožňuje `padlen` zvýšit; zde to nejde a nic nevaruje.

**Doporučení:** škálovat `pad_len ~ C/fc` s ořezem na `n-1` (a případně
varovat, když ořez nastane), nebo alespoň vypsat `# WARNING`, když
`1/(1-r_max)` výrazně přesahuje `pad_len`.

**Priorita:** střední (jediný otevřený nález s praktickým dopadem — uživatelé
s `fc < 0.05` dostávají nedokumentované okrajové artefakty).

### [RESOLVED v5.11.40] 16. Morozov fallback jde špatným směrem

`estimate_cutoff_frequency` (`butterworth.c:452-464`): když žádný kandidát
nesplní discrepancy (residual > 1.1·σ̂ i při `fc = 0.5`, tj. signál má
širokopásmový obsah, který i nejslabší filtr poškozuje), použije se fallback
`AUTO_CUTOFF_FALLBACK = 0.2` — tedy **agresivnější vyhlazení** než nejslabší
vyzkoušený kandidát. Logicky by se mělo padnout na největší kandidát (0.5),
který signál poškozuje nejméně, a varovat.

**Stav:** RESOLVED v5.11.40 — fallback nyní padá na největšího kandidáta
(`fc_candidates[N_AUTO_CANDIDATES-1] = 0.5`) a vypisuje `# WARNING` s návrhem
nastavit `fc` ručně. `AUTO_CUTOFF_FALLBACK = 0.2` zůstává pro selhání odhadu
šumu / alokace, kde žádné filtrování neproběhlo.

### [RESOLVED v5.11.40] 17. Konstanta `NOISE_MAD_NORMALIZATION` nesedí na vlastní komentář

`butterworth.c:42` — kód má `1.6528553` s komentářem `sqrt(6) * 0.6745`,
ale `sqrt(6) * 0.6745 = 1.6521808` (rel. chyba 4.1e-4; vypadá to na překlep
`sqrt(6) ~ 2.4505` místo 2.44949). Prakticky neškodné — discrepancy tolerance
je stejně 1.1 — ale konstanta lže o svém původu.

**Stav:** RESOLVED v5.11.40 — hodnota opravena na `1.6521808` (= sqrt(6)*0.6745).

### [RESOLVED v5.11.40] 18. Mrtvý veřejný typ `ButterworthCoeffs`

`butterworth.h:30-32` — typedef definovaný v hlavičce, nikde v projektu
nepoužitý (kód pracuje s polem `BiquadSection[NUM_BIQUADS]` přímo).

**Stav:** RESOLVED v5.11.40 — typedef odstraněn z `butterworth.h`.

### [RESOLVED v5.11.40] 19. Drobné komentářové neshody

- `calculate_pad_length` komentář říká „3 * filter_order" (= 12), formule
  dává `3*(BUTTERWORTH_ORDER+1)-1` = 14 (`butterworth.c:65-68`).
- `residual_std` komentář říká „sample standard deviation", ale dělí `n`
  (populační, `butterworth.c:384-396`). Pro Morozov bez praktického dopadu.
- Mřížka auto-cutoff kandidátů má skoky až 2.5x bez bisekce mezi sousedy —
  vědomé zjednodušení, zmíněno pro úplnost.

**Stav:** RESOLVED v5.11.40 — komentář `calculate_pad_length` opraven na
`3*(order+1)-1 = 14` a `residual_std` na „population standard deviation
(divides by n)". Mřížka kandidátů ponechána (vědomé zjednodušení).

---

## Aktualizované shrnutí priorit

| # | Problém | Priorita | Stav |
|---|---------|----------|------|
| 2 | `-f auto` je stub | **vysoká** | **RESOLVED v5.11.3** (Morozov) |
| 4 | Žádná kontrola stability pólů | střední | **RESOLVED v5.11.2** |
| 5 | Rozpor `fc` rozsah v CLAUDE.md | nízká | **RESOLVED** |
| 6 | Chybí explicitní spodní limit `fc` | nízká | **RESOLVED v5.11.4** |
| 7 | `CUTOFF_FREQ_STABILITY_WARN` jen warning | nízká | **RESOLVED v5.11.5** (rename + docs) |
| 9 | Nepoužívaný parametr `x` | velmi nízká | **RESOLVED v5.11.6** |
| 10 | Matoucí label `sample_rate` | velmi nízká | **RESOLVED v5.11.6** |
| 11 | GB formát pro menší datasety | velmi nízká | **RESOLVED v5.11.6** |
| 1 | Chybí `y_deriv` | **vysoká** | **RESOLVED v5.11.7** (5-bodové stencily O(h⁴)) |
| 13 | Auto-cutoff běží před CV checkem | nízká | **RESOLVED v5.11.35** (přesun CV checku) |
| 14 | Zavádějící komentář přesnosti derivace | nízká | **RESOLVED v5.11.35** (oprava komentářů) |
| 3 | Cascade IC spoléhá na skrytou invariantu | nízká | **RESOLVED v5.11.40** (důkaz v kód-komentáři) |
| 8 | Tichý fallback v `compute_biquad_ic` | nízká | **RESOLVED v5.11.40** (varování na stderr) |
| 12 | Padding length 14 vs. 9 | velmi nízká | pohlceno #15 |
| 15 | Pevný padding vs. transient pro malá `fc` | **střední** | otevřeno |
| 16 | Morozov fallback agresivnější než nejslabší kandidát | nízká | **RESOLVED v5.11.40** (fallback na největší kandidát) |
| 17 | `NOISE_MAD_NORMALIZATION` nesedí na komentář | kosmetická | **RESOLVED v5.11.40** (1.6521808) |
| 18 | Mrtvý typ `ButterworthCoeffs` | kosmetická | **RESOLVED v5.11.40** (odstraněn) |
| 19 | Komentářové neshody (pad délka, residual_std) | kosmetická | **RESOLVED v5.11.40** (komentáře opraveny) |

---

## Závěr

Od původního auditu (v5.11.1) bylo postupně vyřešeno 8 problémů:

1. **v5.11.2** — `check_pole_stability()` jako runtime pojistka proti numerické
   nestabilitě (warn při radius > 0.99, error při radius ≥ 1.0).

2. **v5.11.3** — skutečná implementace `-f auto` přes Morozovovu discrepancy
   principle (MAD odhad šumu + tolerance 1.1·σ̂). CLI slib „automatic cutoff
   selection" nyní odpovídá realitě.

3. **v5.11.4** — explicitní spodní limit `FC_MIN_PRACTICAL = 1e-4` pro `fc`.
   Uživatel s nesmyslně malým `fc` dostane jasnou chybu před designem filtru
   namísto pozdějšího pole-warning na ztrátu přesnosti.

4. **v5.11.5** — oprava pojmenování `CUTOFF_FREQ_STABILITY_WARN` → `_INEFFECTIVE_WARN`,
   zpřesnění textu varování a dokumentace typického užitečného rozsahu `fc`.
   Vyjasněno, že blízkost Nyquistu není numerický, nýbrž UX problém.

5. **v5.11.6** — drobné úklidy (#9, #10, #11): odstranění nepoužívaného
   parametru `x` z `estimate_cutoff_frequency`, vyjasnění labelu „Effective
   sample rate" s komentářem `= 1/h_avg`, adaptivní MB/GB formát v `stderr`.

6. **v5.11.7** — doplněna podpora derivací (#1): `ButterworthResult.y_deriv`
   počítaný 5-bodovým stencilem O(h⁴) z vyhlazeného výstupu. `-d` flag nyní
   funguje pro Butterworth stejně jako pro ostatní metody. 106 testů (3 nové
   derivation unit tests).

7. **v5.11.35** — kritická re-analýza (Opus 4.8): #13 přesun CV checku před
   auto-cutoff (žádné zbytečné zkušební filtfilty ani stdout šum při odmítnutí
   neuniformního gridu), #14 oprava zavádějících komentářů u přesnosti derivace
   (O(h⁴) jen na uniformním gridu, jinak O(CV)). #3 doplněn explicitní důkaz
   invariantu `buf[0]` a shody se `scipy.sosfilt_zi`. Žádná změna chování na
   uniformním gridu; 111 testů.

8. **v5.11.39 (audit 2026-06-10, Fable 5)** — jádro znovu ověřeno jako korektní:
   bilineární koeficienty přepočteny ručně (A0/a1/a2 i čitatel `Wc²(1,2,1)`,
   DC zisk přesně 1), IC shodná se scipy `lfilter_zi`, in-place aliasing
   v `apply_biquad` bezpečný, všechny 5-bodové stencily sedí na standardní
   vzorce, padding indexace v mezích, Morozovův výběr správně orientovaný.
   Nové nálezy #15-#19 (zatím neopravené).

**Po opravách v5.11.40 zbývá jediný otevřený nález: #15** — pevný
padding 14 vzorků nepokrývá přechodový jev filtru pro `fc < ~0.05` (při
povoleném minimu `fc = 1e-4` je doznění ~6900 vzorků), takže malá `fc` dávají
nedokumentované okrajové artefakty. Je střední priority a vyžaduje větší zásah
(škálování `pad_len ~ C/fc` s ořezem). Nízké/kosmetické nálezy #3, #8, #16-#19
byly vyřešeny v v5.11.40.

Matematická stránka (bilineární transformace s prewarpingem, TDF-II biquad,
odd-reflection padding, IC přes Cramerovo pravidlo, Morozov auto-fc,
5-point stencil derivatives) je implementovaná **korektně**. Re-analýza v5.11.35
ověřila návrh filtru, IC i filtfilt řádek po řádku proti scipy referenci a
kolo 2026-06-10 to nezávislým přepočtem potvrdilo — jádro je matematicky
správné. Otevřený dluh: chování na okrajích při malém `fc` (#15) a přesnost
derivace na neuniformním gridu (pravdivě zdokumentovaná od v5.11.35).

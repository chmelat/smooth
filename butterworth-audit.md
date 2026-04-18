# Analýza implementace Butterworthova filtru

**Datum auditu:** 2026-04-18 (aktualizováno)
**Verze projektu:** smooth v5.11.4
**Auditované soubory:** `butterworth.c` (V1.4/2025-12-07), `butterworth.h`, volající část `smooth.c`

**Historie změn dokumentu:**
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

## Nadále otevřené problémy

### 1. Chybějící podpora derivací (narušení API smlouvy)

`ButterworthResult` v `butterworth.h:35-41` stále nemá pole `y_deriv`,
ačkoli CLAUDE.md definuje jako sdílený vzor (`y_smooth`, `y_deriv`, `n`) pro
všechny metody. V `smooth.c` se to řeší ad-hoc varováním.

**Doporučení:** Buď dopočítat derivace centrální diferencí z `y_smooth`
(nulové fázové zpoždění díky filtfilt to umožňuje), nebo to explicitně
dokumentovat v headeru jako architektonické rozhodnutí.

**Priorita:** vysoká (po vyřešení #2 je to jediný zbývající rozpor s architekturou).

### 3. Kaskádovaná iniciační podmínka se spoléhá na skrytou invariantu

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
i pro první sample, kdy filtr ještě není ustálený — v kombinaci s odd-reflection
paddingem (`padded[0] = 2·y[0] - y[pad_len]`) to v praxi funguje, protože
padding transient zatlumí, ale logika by si zasloužila explicitní důkaz
či odkaz na referenci (scipy filtfilt).

**Priorita:** nízká (dokumentační, kód je správně).

### 8. `compute_biquad_ic` tiše nuluje při degenerate case

`butterworth.c:186-193` — komentář říká „should not happen with valid fc",
ale pokud `|det| ≤ 1e-10` nastane, uživatel dostane transient-zatížený výstup
bez varování.

**Doporučení:** Přidat `fprintf(stderr, ...)` nebo vrátit chybu nahoru.
Varování je obzvlášť důležité ve světle nového Morozov-selektoru, který může
iterovat přes extrémní kandidáty.

**Priorita:** nízká.

### 12. Padding length `3*(order+1)-1 = 14`

`butterworth.c:54` — blízko scipy defaultu (15 pro monolitický 4. řád),
ale pro biquad by bylo konzistentnější `3 * max(len(a), len(b)) = 9`.
Rozdíl v praxi zanedbatelný.

**Priorita:** velmi nízká.

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
| 1 | Chybí `y_deriv` | **vysoká** | otevřeno |
| 3 | Cascade IC spoléhá na skrytou invariantu | nízká | otevřeno (dokumentace) |
| 8 | Tichý fallback v `compute_biquad_ic` | nízká | otevřeno |
| 12 | Padding length 14 vs. 9 | velmi nízká | otevřeno |

---

## Závěr

Od původního auditu (v5.11.1) došlo ke čtyřem zlepšením obranné vrstvy:

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

**Jediný zbývající problém s vysokou prioritou** je absence `y_deriv` v
`ButterworthResult` (issue #1) — nekonzistence s API smlouvou sdílenou mezi
metodami (polyfit, savgol, tikhonov). Protože filtfilt je nulově-fázový,
derivace by se dala spolehlivě dopočítat centrální diferencí z `y_smooth`
bez zavlečení dalšího fázového zkreslení.

Matematická stránka (bilineární transformace s prewarpingem, TDF-II biquad,
odd-reflection padding, IC přes Cramerovo pravidlo, Morozov auto-fc) je
implementovaná **korektně**. Zbývající problémy (#3, #8–#12) jsou převážně
dokumentační/kosmetické a neovlivňují typické použití na uniformním gridu.

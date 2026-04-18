# Analýza implementace Butterworthova filtru

**Datum auditu:** 2026-04-18
**Verze projektu:** smooth v5.11.1
**Auditované soubory:** `butterworth.c` (V1.4/2025-12-07), `butterworth.h`, volající část `smooth.c`

---

## Zásadní designové problémy

### 1. Chybějící podpora derivací (narušení API smlouvy)

`ButterworthResult` nemá pole `y_deriv` (`butterworth.h:34-41`), ačkoli CLAUDE.md
definuje jako sdílený vzor, že všechny metody vrací `y_smooth`, `y_deriv`, `n`.
V `smooth.c:564-567` se na to řeší ad-hoc varováním do stderr.

**Doporučení:** Buď dopočítat derivace (např. centrální diferencí z `y_smooth`),
nebo to explicitně dokumentovat v headeru jako architektonické rozhodnutí.

### 2. `estimate_cutoff_frequency` je pouze stub

`butterworth.c:206-210` — vrací napevno 0.2 nezávisle na datech, ale CLI i header
dokumentaci inzerují „automatic cutoff selection" (`-f auto`). Uživatel zadávající
`-f auto` dostane tichý konstantní výsledek — to je funkční lež.

**Doporučení:** Buď implementovat (např. přes odhad SNR / power spectrum), nebo
`-f auto` odmítnout s jasnou chybou, dokud implementace neexistuje.

### 3. Kaskádovaná iniciační podmínka se spoléhá na skrytou invariantu

`butterworth.c:330-343`:
```c
first_val = y_work[0];
for (int s = 0; s < NUM_BIQUADS; s++) {
    zi[0] = zi_base[s][0] * first_val;
    ...
    first_val = y_work[0];  // výstup i-tého biquadu
}
```

Tento postup funguje jen díky tomu, že každý biquad má DC gain = 1, takže výstup
po ustálení rovná vstupu. Jenže na prvním samplu filtr ještě ustálený není —
`y_work[0]` po prvním biquadu není přesně vstupní step. Pro odd-reflection padding,
kde `padded[0] = 2·y[0] - y[pad_len]`, to navíc znamená, že "step amplitude" pro IC
není `y[0]`, ale reflektovaná hodnota. V praxi padding transient zatlumí, ale
logika si zaslouží komentář vysvětlující, proč chain funguje.

### 4. Žádná kontrola stability navrženého filtru

`design_biquad_sections` spočítá koeficienty, ale nikde se neověří, že póly leží
uvnitř jednotkového kruhu (tedy `|a2| < 1` a `|a1| < 1 + a2`). Pro platné
`fc ∈ (0,1)` a numericky zdravý `tan()` by to mělo vždy sedět, ale defenzivní
check za ~4 řádky by zachytil krajní případy (např. `fc` velmi blízko 0 nebo 1
s akumulovanou chybou).

---

## Středně závažné problémy

### 5. Rozpor v konvenci `fc` s CLAUDE.md

CLAUDE.md říká `0 < fc < 0.5`, kód i CLI říkají `0 < fc < 1`
(`butterworth.c:23-24`, `smooth.c:130, 685-686`). CLAUDE.md je zastaralý — buď
sjednotit, nebo dokumentaci opravit.

### 6. Žádný dolní limit na `fc`

Pro `fc → 0` platí `Wc = tan(π·fc/2) → 0` a numerátor biquadu `≈ Wc²/A0`,
kde `A0 ≈ 4`. To dává velmi malé `b` koeficienty a `a2 → 1`, což znamená póly
blízko jednotkového kruhu → ztráta přesnosti.

**Doporučení:** Tvrdý spodní limit (např. `fc > 1e-4`).

### 7. `CUTOFF_FREQ_STABILITY_WARN = 0.95` je jen warning

`butterworth.c:267-269` — pro `fc` v rozsahu (0.95, 1.0) vzniká `Wc` velmi velké,
opět numericky rizikové.

**Doporučení:** Zpřísnit na hard error, nebo alespoň tento rozsah dokumentovat
jako „experimentální".

### 8. `compute_biquad_ic` tiše nuluje při degenerate case

`butterworth.c:122-126` — komentář říká „should not happen with valid fc", ale
pokud k tomu dojde, uživatel dostane zkreslený výstup bez jakéhokoliv varování.

**Doporučení:** Přidat `fprintf(stderr, ...)` nebo vrátit chybu nahoru.

---

## Drobnosti

### 9. Parametr `x` je v hlavní cestě nepoužitý

Vstupuje jen do stub `estimate_cutoff_frequency`. Grid spacing se bere
z `grid_info->h_avg`. Zavádějící signatura.

### 10. `sample_rate` je pouze metadata

Filtr pracuje normalizovaně, `sample_rate = 1/h_avg` slouží jen k tisku hlavičky.
Pro neuniformní grid (CV blízko 0.15) to může být matoucí, protože skutečný
„sample rate" neexistuje. Pro méně matoucí výstup by stálo za zvážení uvést
i `h_avg` s varováním.

### 11. Memory estimate tiskne `GB` i pro menší datasety

`butterworth.c:249-250` — podmínka `n > 50M` to omezuje, ale formát `%.1f GB` je
matoucí — 50M bodů dává ~0.8 GB.

### 12. Padding length `3*(order+1)-1 = 14`

Blízko scipy defaultu (15 pro monolitický 4. řád), ale pro biquad by bylo
konzistentnější `3 * max(len(a), len(b)) = 9`. Rozdíl v praxi zanedbatelný.

---

## Shrnutí priorit

| # | Problém | Priorita |
|---|---------|----------|
| 2 | `-f auto` je stub | **vysoká** (uživatelsky viditelná lež) |
| 1 | Chybí `y_deriv` | **vysoká** (nekonzistence API) |
| 4 | Žádná kontrola stability pólů | střední |
| 3 | Cascade IC spoléhá na skrytou invariantu | střední (dokumentační) |
| 6, 7 | Chybí limity pro extrémní `fc` | střední |
| 5 | Rozpor `fc` rozsah v CLAUDE.md | nízká (jen dokumentace) |
| 8 | Tichý fallback v `compute_biquad_ic` | nízká |
| 9–12 | Kosmetika / drobná UX | nízká |

---

## Závěr

Matematická stránka (bilineární transformace s prewarpingem, TDF-II biquad,
odd-reflection padding, IC přes Cramerovo pravidlo) je implementovaná **korektně**.

Hlavní slabiny jsou v:

- **API konzistenci** — chybí derivace, které ostatní metody poskytují
- **Dodržování slibů CLI** — `-f auto` reálně nefunguje
- **Defenzivních kontrolách** — žádné meze na extrémní `fc`, žádná ověření
  stability, tiché fallbacky

Žádný problém není kritický pro typické použití (uniformní grid, `fc ∈ [0.05, 0.5]`),
ale všechny zvyšují riziko tichého selhání mimo tento „happy path".

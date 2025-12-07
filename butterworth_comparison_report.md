# Porovnání implementací Butterworthova filtru

**Verze:** 1.0  
**Datum:** 2025-12-07  
**Autor:** Technická analýza

---

## 1. Úvod

Tento report porovnává tři implementace 4. řádu Butterworthova dolnopropustného filtru s filtfilt (zero-phase forward-backward filtering) z hlediska použitelnosti na reálná zašuměná data.

### Porovnávané verze

| Verze | Struktura filtru | Počáteční podmínky | Předzpracování | Závislosti |
|-------|------------------|-------------------|----------------|------------|
| V1.3 | Monolitický 4. řád (5 koeficientů) | LAPACK dgesv_ | Žádné | LAPACK |
| V4.0 | Biquad kaskáda (2×2. řád) | Nulové | Lineární detrending | Žádné |
| V1.4 | Biquad kaskáda (2×2. řád) | Analytické pro každý biquad | Žádné | Žádné |

---

## 2. Teoretické pozadí

### 2.1 Biquad kaskáda vs. monolitický filtr

Filtr 4. řádu lze implementovat dvěma způsoby:

**Monolitický přístup (V1.3):**
- Přenosová funkce jako jeden polynom 4. stupně
- 5 koeficientů čitatele (b₀...b₄) a 5 koeficientů jmenovatele (a₀...a₄)
- Potenciální akumulace zaokrouhlovacích chyb

**Biquad kaskáda (V4.0, V1.4):**
- Dva sériově zapojené filtry 2. řádu
- Každá sekce má pouze 3+3 koeficienty
- Lepší numerická stabilita, zejména pro vyšší fc

### 2.2 Počáteční podmínky (IC)

Pro filtfilt je kritické správné nastavení počátečních podmínek filtru, aby se minimalizoval přechodový děj (transient) na okrajích signálu.

**Princip:** Pro konstantní vstup x[n] = c by měl být výstup okamžitě y[n] = c (steady-state).

Řešení: `zi = (I - A)⁻¹ · B`, kde A je companion matrix filtru.

### 2.3 Problém detrendingu V4.0

V4.0 používá lineární detrending založený na dvou bodech (y[0] a y[n-1]):

```
slope = (y[n-1] - y[0]) / (x[n-1] - x[0])
```

**Problém:** U zašuměných dat jsou y[0] a y[n-1] zatíženy šumem, což vede k falešnému odhadu trendu.

---

## 3. Metodika testování

### 3.1 Parametry testů

| Parametr | Hodnota |
|----------|---------|
| Počet bodů (n) | 200 |
| Normalizovaná cutoff frekvence (fc) | 0.15 |
| Standardní odchylka šumu (σ) | 0.5 (pokud není uvedeno jinak) |
| Seed generátoru | 42 |

### 3.2 Metriky

- **RMSE celkové:** Root Mean Square Error přes celý signál
- **RMSE okraje:** RMSE pouze pro prvních a posledních 10 % bodů

### 3.3 Testovací scénáře

1. Konstantní signál s offsetem (y = 100)
2. Lineární trend (y = 50 + 100x)
3. Sinusoida bez offsetu (y = sin(4πx))
4. Sinusoida s velkým offsetem (y = 1000 + sin(4πx))
5. Sinusoida + lineární trend (y = 100x + sin(4πx))
6. Outlier na začátku (y = 50 + outlier na y[0])
7. Vysoký šum (σ = 5)
8. Kvadratický trend (y = 100x²)
9. Step funkce
10. Exponenciální trend (y = e^(2x))

---

## 4. Výsledky

### 4.1 Celkové RMSE

| Test | V1.3 | V4.0 | V1.4 | Nejlepší |
|------|------|------|------|----------|
| Konstanta (y=100) | 0.19175 | **0.16703** | 0.16958 | V4.0 |
| Lineární trend | 0.21883 | **0.21525** | 0.21564 | V4.0 |
| Sinusoida (střed=0) | 0.20907 | 0.12985 | **0.12880** | V1.4 |
| Sinusoida (offset=1000) | 0.16866 | **0.13768** | 0.13812 | V4.0 |
| Sinusoida + trend | 0.23619 | 0.19037 | **0.16887** | V1.4 |
| Outlier na y[0] | **0.96105** | 1.29285 | 1.21452 | V1.3 |
| Vysoký šum (σ=5) | 2.10665 | 1.96016 | **1.86043** | V1.4 |
| Kvadratický trend | **0.15631** | 0.19307 | 0.21925 | V1.3 |
| Step funkce | **5.45837** | 7.73094 | 7.73055 | V1.3 |
| Exponenciální trend | 0.19194 | **0.12315** | 0.12686 | V4.0 |

**Celkový výsledek:** V1.3: 3 vítězství, V4.0: 4 vítězství, V1.4: 3 vítězství

### 4.2 RMSE na okrajích (10 % z každé strany)

| Test | V1.3 | V4.0 | V1.4 | Nejlepší |
|------|------|------|------|----------|
| Konstanta (y=100) | **0.17586** | 0.22201 | 0.23267 | V1.3 |
| Lineární trend | **0.32103** | 0.41445 | 0.41319 | V1.3 |
| Sinusoida (střed=0) | 0.20347 | 0.16073 | **0.15737** | V1.4 |
| Sinusoida (offset=1000) | 0.21286 | **0.18848** | 0.19158 | V4.0 |
| Sinusoida + trend | 0.32493 | 0.33976 | **0.28102** | V1.4 |
| Outlier na y[0] | **2.10199** | 2.87023 | 2.69425 | V1.3 |
| Vysoký šum (σ=5) | **3.19452** | 3.64377 | 3.38602 | V1.3 |
| Kvadratický trend | **0.14056** | 0.34547 | 0.41232 | V1.3 |
| Step funkce | 0.21132 | 0.29136 | **0.20973** | V1.4 |
| Exponenciální trend | 0.26200 | **0.17317** | 0.18410 | V4.0 |

**Celkový výsledek okrajů:** V1.3: 6 vítězství, V4.0: 2 vítězství, V1.4: 2 vítězství

---

## 5. Verifikace počátečních podmínek

### 5.1 Test steady-state odezvy

Pro biquad sekci s fc = 0.15 byl proveden test konstantního vstupu:

**Bez IC (z = [0, 0]):**
```
y[0]=1.2, y[1]=5.4, y[2]=12.5, ..., y[14]=89.7
```

**S IC (z = zi_base × vstup):**
```
y[0]=100.0, y[1]=100.0, y[2]=100.0, ..., y[14]=100.0
```

Výsledek potvrzuje správnost analytického výpočtu IC pro V1.4.

### 5.2 Stabilita pólů

Pro všechny testované fc (0.05 až 0.90) jsou póly biquad sekcí uvnitř jednotkové kružnice:

| fc | Sekce 0 |z| | Sekce 1 |z| | Stav |
|----|----------|----------|------|
| 0.05 | 0.9299 | 0.9704 | Stabilní |
| 0.15 | 0.8007 | 0.9132 | Stabilní |
| 0.30 | 0.6220 | 0.8309 | Stabilní |
| 0.50 | 0.3873 | 0.7288 | Stabilní |
| 0.70 | 0.1991 | 0.6682 | Stabilní |
| 0.90 | 0.5528 | 0.7993 | Stabilní |

---

## 6. Analýza citlivosti detrendingu V4.0

### 6.1 Test s konstantním signálem (y = 50) a šumem σ = 1

20 nezávislých realizací šumu ukázalo:

| Metrika | Hodnota |
|---------|---------|
| Případy s max. chybou trendu > 1.0 | 6 z 20 (30 %) |
| Maximální zaznamenaná chyba trendu | 1.40 |
| Průměrná absolutní chyba slope | 0.015 |

### 6.2 Porovnání metod odhadu trendu

Pro jednu realizaci šumu:

| Metoda | Odhadnutý slope | Chyba |
|--------|-----------------|-------|
| Skutečný | 0.0000 | — |
| V4.0 (2 body) | 0.0027 | 0.0027 |
| Robustní (k=10) | 0.0107 | 0.0107 |
| LSQ (všechna data) | 0.0084 | 0.0084 |

Paradoxně, 2-bodový odhad může být v některých případech nejlepší, ale má vysokou varianci.

---

## 7. Diskuse

### 7.1 Silné stránky jednotlivých verzí

**V1.3:**
- Nejlepší okrajové chování díky správným IC pro celý filtr
- Robustní vůči outlierům a nelineárním trendům
- Nevýhoda: vyžaduje LAPACK

**V4.0:**
- Numericky stabilnější biquad struktura
- Žádné externí závislosti
- Nevýhoda: detrending citlivý na šum v okrajových bodech

**V1.4:**
- Kombinuje výhody obou: biquad + správné IC
- Žádné externí závislosti
- Nejlepší pro sinusoidální signály a vysoký šum
- Nevýhoda: IC se škálují postupně v kaskádě, což může být suboptimální pro některé signály

### 7.2 Doporučení pro výběr verze

| Scénář | Doporučená verze |
|--------|------------------|
| Obecné použití na zašuměná data | V1.4 |
| Signály s nelineárním trendem | V1.3 |
| Prostředí bez LAPACK | V1.4 |
| Lineární trendy | V4.0 nebo V1.4 |
| Sinusoidální signály | V1.4 |
| Kritické okrajové chování | V1.3 |

---

## 8. Závěr

Nová verze V1.4 představuje hybridní řešení kombinující:

1. **Biquad kaskádu** z V4.0 pro lepší numerickou stabilitu
2. **Správné počáteční podmínky** z V1.3 (analyticky počítané bez LAPACK)
3. **Odstranění problematického detrendingu** z V4.0

V1.4 dosahuje nejlepších výsledků pro:
- Sinusoidální signály s trendem (RMSE 0.169 vs. 0.190/0.236)
- Vysoký šum (RMSE 1.860 vs. 1.960/2.107)
- Čisté sinusoidy (RMSE 0.129 vs. 0.130/0.209)

V1.3 zůstává preferovanou volbou pro:
- Nelineární trendy (kvadratické, exponenciální)
- Signály s outliers na okrajích
- Aplikace kde je kritické okrajové chování

Pro praktické použití na reálných zašuměných datech s fc v rozmezí 0.1–0.5 je **V1.4 doporučenou volbou** díky kombinaci robustnosti a absence externích závislostí.

---

## Přílohy

### A. Klíčové vzorce

**Biquad koeficienty (bilineární transformace):**
```
Wc = tan(π·fc/2)
θ₁ = π/8, θ₂ = 3π/8

Pro každou sekci i:
  α = 2·Wc·cos(θᵢ)
  A₀ = 4 + 2α + Wc²
  
  a₁ = (2Wc² - 8) / A₀
  a₂ = (4 - 2α + Wc²) / A₀
  
  gain = Wc² / A₀
  b = [gain, 2·gain, gain]
```

**Počáteční podmínky pro biquad:**
```
det = 1 + a₁ + a₂
B₀ = b₁ - a₁·b₀
B₁ = b₂ - a₂·b₀

zi[0] = (B₀ + B₁) / det
zi[1] = (-a₂·B₀ + (1+a₁)·B₁) / det
```

### B. Testovací prostředí

- Kompilátor: GCC s -O2 optimalizací
- Platforma: Linux x86_64
- Precision: double (IEEE 754)

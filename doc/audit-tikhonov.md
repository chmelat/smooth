# Audit modulu tikhonov.c

Datum: 2026-06-10. Cílený audit modulu `tikhonov.c` (+ `tikhonov.h`) ve verzi
v5.11.38; nálezy TK1 a TK2 opraveny v **v5.11.39** v rámci tohoto auditu.
Čísla řádků odpovídají stavu po opravě (v5.11.39).

Kontext: poslední celokódový audit je `code-audit-v5.11.38.md` (otevřené nálezy
A2/B2/B3 se tohoto modulu netýkají). Výpočetní jádro modulu — Gramova matice
`sum_k w_k d_k^T d_k`, pásové řešení `dpbsv`, správa paměti přes `goto error` —
prošlo čistě: koeficienty 3-bodového stencilu i redukce na
`[1,-4,6,-4,1]*lambda/h^3` na uniformní mřížce ověřeny ručně, error cesty bez
leaků, `free_tikhonov_result` NULL-safe. Slabiny jsou soustředěné v automatické
volbě lambda (GCV / L-curve).

## Opravené nálezy (v5.11.39)

### TK1. ~~Trace aproximace pro n > 5000 měla špatný asymptotický exponent~~ — **FIXED v5.11.39**

`compute_gcv_score_robust` používala pro n > 5000 zkratku
`trace_H = n/(1+sqrt(lambda/h^3))`, tj. ~n*s^(-1/2) pro s = lambda/h^3.
Správná asymptotika z vlastních čísel penalizační matice je ale
~0.35*n*s^(-1/4): trace = sum 1/(1+16*s*sin^4(theta/2)) a integrál pro s >> 1
škáluje jako s^(-1/4). Při s = 1e4 podhodnocení ~3.6x, dále roste jako s^(1/4).
Podhodnocený trace tlačí GCV jmenovatel k 1, GCV degraduje na RSS/n a výběr je
systematicky vychýlen k **podvyhlazení** pro všechna n > 5000.

**Fix:** větev odstraněna; analytická vlastní suma se používá pro všechna n.
Suma je O(n), tedy stejného řádu jako samotný pásový solver — aproximace nic
podstatného nešetřila. Ověřeno na n = 6000: tr(H) při lambda = 1 vychází 68.6,
stará formule by dala ~6 (11x podhodnocení).

### TK2. ~~Tiché přiskřípnutí lambda na hranici pevného rozsahu~~ — **FIXED v5.11.39**

Lambda je rozměrová veličina (škáluje s h^3 a amplitudou y^2), takže pevný
rozsah hledání [1e-8, 1] nemůže sedět každé škále dat. Optimum mimo rozsah se
mlčky přiskříplo na hranici bez jakékoli diagnostiky; refinement navíc
podmínkou `lambda_test > lambda_min && < lambda_max` nesměl za hranici vyjet.

**Fix:** meze `lambda_min`/`lambda_max` vytaženy na začátek
`find_optimal_lambda_gcv` (sdílené oběma větvemi); když vybrané optimum leží na
kraji rozsahu, vypíše se `# WARNING` s doporučením zadat lambda ručně přes
`-l`. README doplněno o odstavec o rozměrovosti lambda. Ověřeno funkčně na
datech s h = 10 (optimum na horní mezi → varování se vypíše; vnitřní optimum →
nevypíše).

## Otevřené nálezy

### TK3. L-curve detekce rohu je trojnásobně křehká — `tikhonov.c:382,400,384-411` — **OPEN** (střední)

Aktivní jen pro n > 20000. Tři nezávislé problémy:

1. Křivost přes `fabs(dx*ddy - dy*ddx)` (`:400`) nerozlišuje orientaci
   zakřivení — může vybrat bod na konvexní bouli místo rohu L-křivky.
   Standardně se bere znaménková křivost správné orientace.
2. Konečné diference jsou parametrizovány indexem a předpokládají rovnoměrný
   krok, ale lambda seznam pro n > 20000 (`:452`) není log-uniformní
   (1e-2 → 5e-2 → 1e-1 → 2e-1 → 5e-1 → 1e0) — křivost na hustším konci se
   uměle nafukuje.
3. Když se nenajde žádná platná křivost, mlčky se vrátí
   `best_idx = n_lambda/2` (`:382`) — arbitrární prostřední lambda bez
   diagnostiky.

**Fix:** znaménková křivost; buď log-uniformní lambda mřížka, nebo diference
vůči log10(lambda) místo indexu; warning při fallbacku.

### TK4. GCV trace je aproximace i na uniformní mřížce — `tikhonov.c:295-308` — **OPEN** (nízká)

Pentadiagonální Gram s přirozenými okrajovými podmínkami (krajní řádky
`[1,-2,1]*lambda/h^3`) není čtverec Dirichletova Laplaciánu, takže vzorec
`16 sin^4(pi k / 2n) / h^3` se odchyluje u nejmenších nenulových vlastních
čísel — a ta dominují trace právě pro velká lambda. Komentář v kódu už
formulačně sedí ("eigenvalue model", v5.11.39), matematická podstata trvá.

**Fix (pokud by byl potřeba):** přesný trace — stochastický Hutchinsonův odhad,
nebo diagonála inverze přes pásovou faktorizaci. Větší zásah; relevantní jen
pokud záleží na přesnosti auto-lambda.

### TK5. Výstupní derivace jen 1. řádu na non-uniformní mřížce — `tikhonov.c:103-108` — **OPEN** (nízká)

`compute_derivatives` používá centrální diferenci
`(y[i+1]-y[i-1])/(x[i+1]-x[i-1])`, která je 2. řádu jen pro h_l = h_r. Modul si
jinak na non-uniformní korektnosti zakládá (celá v5.11 oprava penalizace), ale
`y_deriv` tu péči nemá.

**Fix:** 3-bodový vážený vzorec (2. řád) za stejnou cenu; krajní body zůstávají
1. řádu jednostranné.

### TK6. Dokumentační drift v `tikhonov.h` — `tikhonov.h:36,75-77,91-123` — **OPEN** (nízká)

- Příklad v hlavičce volá `tikhonov_smooth(x, y, n, 0.1)` a
  `find_optimal_lambda_gcv(x, y, n)` bez `grid_info` — podle vlastního příkladu
  se nezkompiluje (`:101,104-105`).
- "Search range: 1e-6 to 1e0" (`:75`) vs. skutečných 1e-8 až 1e0.
- "Warns if regularization term dominates excessively" (`:77`) v kódu
  neexistuje.
- "used for discretization method selection" (`:36`) je relikt z doby před
  sjednocením diskretizace (v5.11.34) — žádná selekce už neexistuje.

**Fix:** triviální úprava komentářů; hlavička aktivně lže, doporučená první.

### TK7. Výkonové redundance v auto-lambda cestě — `tikhonov.c:283-288,460-468` — **OPEN** (kosmetika)

- Pro n > 20000 GCV vyřeší 12 soustav a `find_lambda_lcurve` pak **stejných
  12 lambda řeší znovu** (24 O(n) solverů místo 12) — šlo by předat už
  spočtené `data_term`/`regularization_term`.
- `compute_gcv_score_robust` přepočítává RSS smyčkou (`:286-288`), ač je už
  v `result->data_term`.
- Každý GCV trial počítá a zahazuje `y_deriv` a revaliduje monotonii x
  (~21 zbytečných O(n) průchodů na jedno auto-lambda).

Vše O(n), prakticky neznatelné mimo největší datasety.

### TK8. Drobné nekonzistence — **OPEN** (kosmetika)

- Chybová hláška u `dpbsv` říká "I + lambda*(D2)^T D2" bez W
  (`tikhonov.c:239-241`).
- Selhání mallocu v `find_lambda_lcurve` mlčky vrátí default 0.01 bez warningu
  (`:343-351`).
- Refinement gate `if (n <= 5000)` (`:488`) ztratil po TK1 opodstatnění —
  původně kryl hranici, za níž byl trace jen hrubá aproximace; teď je trace
  stejně přesný pro všechna n, ale pro 5000 < n <= 20000 se refinement pořád
  nedělá. Jednořádkové sjednocení.

## Stav nálezů

| ID  | Místo                          | Závažnost | Stav           |
|-----|--------------------------------|-----------|----------------|
| TK1 | GCV trace (n > 5000)           | střední   | FIXED v5.11.39 |
| TK2 | hranice rozsahu lambda         | střední   | FIXED v5.11.39 |
| TK3 | L-curve roh (n > 20000)        | střední   | OPEN           |
| TK4 | trace vs. natural BC           | nízká     | OPEN           |
| TK5 | derivace na non-uniform mřížce | nízká     | OPEN           |
| TK6 | doc drift v tikhonov.h         | nízká     | OPEN           |
| TK7 | redundantní solvery/průchody   | kosmetika | OPEN           |
| TK8 | hláška, fallback, refine gate  | kosmetika | OPEN           |

Doporučené pořadí oprav: TK6 (triviální), TK8-refinement (přímý důsledek TK1),
TK5 (malý lokální zásah, měřitelný přínos pro `-d` výstup). TK3 a TK4 jen
pokud se reálně pracuje s n > 20000, resp. pokud záleží na přesnosti
auto-lambda — jsou to větší zásahy.

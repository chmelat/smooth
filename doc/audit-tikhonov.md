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

## Opravené nálezy (později)

### TK5. ~~Výstupní derivace jen 1. řádu na non-uniformní mřížce~~ — **FIXED v5.11.49**

`compute_derivatives` používala centrální diferenci
`(y[i+1]-y[i-1])/(x[i+1]-x[i-1])`, 2. řádu jen pro h_l = h_r.

**Fix:** 3-bodový vážený vzorec z neurčitých koeficientů (`tikhonov.c:87-125`),
2. řádu na libovolném rozestupu; krajní body dostaly jednostranný 3-bodový
vzorec, také 2. řádu (audit je měl zůstat 1. řádu — plán 004 je zpřísnil).

Doplněk z plánu 004: **CV je pro tento nález špatný prediktor.** Vedoucí chyba
je `(h_r-h_l)/2 * u''`, tedy úměrná lokální *asymetrii*, ne celkovému rozptylu
rozestupů. Změřeno: alternující mřížka 0.1/0.4 (CV 0.608) se opravou zlepší
**22x**, geometricky odstupňovaná r=1.10 (CV 0.622) **1x**. Jakékoli budoucí
rozhodování typu "je tahle mřížka špatná pro derivace" má sáhnout po
`|h_r-h_l|/h`, ne po `GridAnalysis.cv`.

### TK6. ~~Dokumentační drift v `tikhonov.h`~~ — `d57a778` — **FIXED**

Všechny čtyři body opraveny (rozsah 1e-8, smyšlené varování, relikt o selekci
diskretizace, příklad bez `grid_info`). Plán 003 navíc našel tři další defekty
téže třídy, které audit neměl: dvojitá deklarace `TikhonovResult *result` v
jednom scope, kompilační řádek bez `grid_analysis.c`, a — mimo `tikhonov.h`
úplně — `savgol.h:73,79` volající `analyze_grid(x, 6, 0)` se třemi argumenty
proti dvouargumentové signatuře. Ověřeno extrakcí a kompilací obou příkladů.

## Zaniklé nálezy

### TK3. ~~L-curve detekce rohu je trojnásobně křehká~~ — **DEAD v5.11.44**

Všechny tři problémy (neznaménková křivost, indexové diference nad
non-log-uniformní lambda mřížkou, tichý `best_idx = n_lambda/2` fallback) žily
výhradně ve `find_lambda_lcurve()`, dosažitelné jen z větve `n > 20000`.
Ponytail audit v5.11.43 (nález 1) obojí smazal — funkce ani větev už
neexistují, GCV sweep běží pro všechna n. Nález zaniká bez opravy.

## Otevřené nálezy

### TK4. GCV trace je aproximace i na uniformní mřížce — `tikhonov.c:316-329` — **OPEN, zamítnuto** (nízká)

Pentadiagonální Gram s přirozenými okrajovými podmínkami (krajní řádky
`[1,-2,1]*lambda/h^3`) není čtverec Dirichletova Laplaciánu, takže vzorec
`16 sin^4(pi k / 2n) / h^3` se odchyluje u nejmenších nenulových vlastních
čísel — a ta dominují trace právě pro velká lambda. Komentář v kódu už
formulačně sedí ("eigenvalue model", v5.11.39), matematická podstata trvá.

**Zamítnuto 2026-07-27:** nic na přesnosti auto-lambda nezávisí, takže níže
popsaný zásah nemá co ospravedlnit. Otevřít znovu jen když se to změní.

**Fix (pokud by byl potřeba):** přesný trace — stochastický Hutchinsonův odhad,
nebo diagonála inverze přes pásovou faktorizaci. Větší zásah; relevantní jen
pokud záleží na přesnosti auto-lambda.

### TK7. Výkonové redundance v auto-lambda cestě — `tikhonov.c:298-311,394-411` — **OPEN, zamítnuto** (kosmetika)

**Zamítnuto 2026-07-27:** změřeno na 100k bodech — `-m 2 -l auto` 0.465 s vs.
0.210 s pro `-l 0.1`, tedy ~0.25 s redundantní práce v úloze jinak dominované
I/O. Nestojí to za zásah do auto-lambda cesty.


- ~~Pro n > 20000 GCV vyřeší 12 soustav a `find_lambda_lcurve` pak stejných
  12 lambda řeší znovu~~ — zaniklo s TK3 (v5.11.44).
- `compute_gcv_score_robust` přepočítává RSS smyčkou (`:298-311`), ač je už
  v `result->data_term`.
- Každý GCV trial počítá a zahazuje `y_deriv` a revaliduje monotonii x
  (~21 zbytečných O(n) průchodů na jedno auto-lambda).

Vše O(n), prakticky neznatelné mimo největší datasety.

### TK8. Drobné nekonzistence — **OPEN** (kosmetika)

- Chybová hláška u `dpbsv` říká "I + lambda*(D2)^T D2" bez W
  (`tikhonov.c:263`).
- ~~Selhání mallocu v `find_lambda_lcurve` mlčky vrátí default 0.01~~ — zaniklo
  s TK3 (v5.11.44).
- Refinement gate `if (n <= 5000)` (`:404`) ztratil po TK1 opodstatnění —
  původně kryl hranici, za níž byl trace jen hrubá aproximace; teď je trace
  stejně přesný pro všechna n, ale pro 5000 < n <= 20000 se refinement pořád
  nedělá. Jednořádkové sjednocení.

## Stav nálezů

Aktualizováno 2026-07-29 proti v5.11.53.

| ID  | Místo                          | Závažnost | Stav           |
|-----|--------------------------------|-----------|----------------|
| TK1 | GCV trace (n > 5000)           | střední   | FIXED v5.11.39 |
| TK2 | hranice rozsahu lambda         | střední   | FIXED v5.11.39 |
| TK3 | L-curve roh (n > 20000)        | střední   | DEAD v5.11.44  |
| TK4 | trace vs. natural BC           | nízká     | OPEN (zamítnuto) |
| TK5 | derivace na non-uniform mřížce | nízká     | FIXED v5.11.49 |
| TK6 | doc drift v tikhonov.h         | nízká     | FIXED `d57a778` |
| TK7 | redundantní solvery/průchody   | kosmetika | OPEN (zamítnuto) |
| TK8 | hláška, refine gate            | kosmetika | OPEN           |

Zbývá jen kosmetika. TK4 a TK7 byly při plánování (2026-07-27) posouzeny a
**zamítnuty**, ne odloženy: TK4 (přesný trace přes Hutchinsona nebo pásovou
faktorizaci) je velký zásah opodstatněný jen tehdy, kdyby na přesnosti
auto-lambda reálně něco záviselo — nezávisí; TK7 změřeno na 100k bodech,
`-m 2 -l auto` 0.465 s vs. 0.210 s pro `-l 0.1`, tedy ~0.25 s redundantní práce
v úloze jinak dominované I/O. Z TK8 zbývají dvě jednořádkovky: chybějící W
v chybové hlášce a refinement gate `n <= 5000`, který po TK1 nemá důvod.

Průřezové poučení, které přežilo jednotlivé nálezy: **žádná politika se nemá
věšet na `GridAnalysis.cv`** — TK5 ukázal, že pro derivace rozhoduje lokální
asymetrie, ne CV, a dvě mřížky se shodným CV se v dopadu liší 22x.

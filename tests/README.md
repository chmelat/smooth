# Unit testy pro SMOOTH projekt

Tento adresář obsahuje unit testy pro projekt `smooth` implementované pomocí [Unity testing frameworku](https://github.com/ThrowTheSwitch/Unity).

## Struktura

```
tests/
├── README.md                   # Tento soubor
├── unity.c                     # Unity framework (staženo z GitHub)
├── unity.h                     # Unity header
├── test_main.c                 # Test runner - spouští všechny testy
└── test_grid_analysis.c        # Testy pro grid_analysis.c modul
```

## Jak spustit testy

### Základní spuštění

```bash
# Z hlavního adresáře projektu:
make test
```

To provede:
1. Kompilaci testovacích souborů
2. Linkování s testovanými moduly
3. Spuštění všech testů
4. Vypsání výsledků

### Kontrola memory leaks (Valgrind)

```bash
make test-valgrind
```

To spustí testy s Valgrindem který kontroluje:
- Memory leaky (nevolané `free()`)
- Invalid memory access
- Použití neinicializované paměti

### Vyčištění testovacích souborů

```bash
make test-clean
```

## Co je testováno

### `test_grid_analysis.c` - Testy pro grid_analysis modul

| Test | Co testuje | Účel |
|------|-----------|------|
| `test_grid_perfectly_uniform` | Perfektně uniformní grid | Základní funkcionalita |
| `test_grid_nonuniform` | Výrazně neuniformní grid | Detekce variace |
| `test_grid_minimum_points` | n=2 (minimum bodů) | Edge case |
| `test_grid_nearly_uniform` | Téměř uniformní (malá variace) | Realistický případ |
| `test_grid_null_pointer` | NULL pointer vstup | Error handling |
| `test_grid_large_dataset` | n=1000 bodů | Performance |
| `test_grid_with_outlier` | Grid s jedním outlierem | Realistický problém |

## Příklad výstupu

### Všechny testy prošly:

```
========================================
  SMOOTH PROJECT - UNIT TESTS
========================================
Running tests for grid_analysis module
========================================

--- Basic functionality tests ---
test_grid_analysis.c:44:test_grid_perfectly_uniform:PASS
test_grid_analysis.c:78:test_grid_nonuniform:PASS

--- Edge case tests ---
test_grid_analysis.c:108:test_grid_minimum_points:PASS
test_grid_analysis.c:178:test_grid_null_pointer:PASS

--- Realistic scenario tests ---
test_grid_analysis.c:145:test_grid_nearly_uniform:PASS
test_grid_analysis.c:225:test_grid_with_outlier:PASS

--- Performance tests ---
test_grid_analysis.c:195:test_grid_large_dataset:PASS

-----------------------
7 Tests 0 Failures 0 Ignored
OK
```

### Některý test selhal:

```
test_grid_analysis.c:180:test_grid_null_pointer:FAIL: Expected NULL Was 0x55d2f4a6c2a0
                                                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                    Tento test očekával NULL, ale
                                                    funkce vrátila pointer!

-----------------------
7 Tests 1 Failures 0 Ignored
FAIL
```

## Jak přidat nové testy

### 1. Vytvoř nový testovací soubor (např. `test_savgol.c`)

```c
#include "unity.h"
#include "../savgol.h"

void setUp(void) {
    // Příprava před každým testem
}

void tearDown(void) {
    // Úklid po každém testu
}

void test_savgol_constant(void) {
    // ARRANGE
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[] = {5.0, 5.0, 5.0, 5.0, 5.0};
    int n = 5;

    // ACT
    SavgolResult *result = savgol_smooth(x, y, n, 5, 2, NULL);

    // ASSERT
    TEST_ASSERT_NOT_NULL(result);
    for (int i = 0; i < n; i++) {
        TEST_ASSERT_DOUBLE_WITHIN(0.001, 5.0, result->y_smooth[i]);
    }

    // CLEANUP
    free_savgol_result(result);
}
```

### 2. Přidej nový test do `test_main.c`

```c
// Nahoře přidej deklaraci:
void test_savgol_constant(void);

// V main() přidej:
printf("\n--- Savitzky-Golay tests ---\n");
RUN_TEST(test_savgol_constant);
```

### 3. Aktualizuj Makefile

V `Makefile` přidej nový testovací soubor:

```makefile
TEST_SRC = $(TEST_DIR)/test_main.c \
           $(TEST_DIR)/test_grid_analysis.c \
           $(TEST_DIR)/test_savgol.c \      # <-- nový
           $(TEST_DIR)/unity.c

TEST_MODULES = grid_analysis.o savgol.o  # <-- přidej testované moduly
```

### 4. Spusť testy

```bash
make test
```

## Unity assertions (nejpoužívanější)

```c
// Pointery
TEST_ASSERT_NULL(ptr);
TEST_ASSERT_NOT_NULL(ptr);

// Integers
TEST_ASSERT_EQUAL_INT(expected, actual);

// Doubles/Floats (VŽDY použij WITHIN!)
TEST_ASSERT_DOUBLE_WITHIN(tolerance, expected, actual);

// Porovnání
TEST_ASSERT_LESS_THAN(threshold, actual);
TEST_ASSERT_GREATER_THAN(threshold, actual);

// Arrays
TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, actual, length);

// Boolean
TEST_ASSERT_TRUE(condition);
TEST_ASSERT_FALSE(condition);
```

## Best practices

### ✅ DO (Dělej):

- **Testuj okrajové případy** (n=1, n=2, NULL, prázdné pole)
- **Testuj realistické scénáře** (ne jen ideální případy)
- **Vždy uvolni paměť** na konci testu (zavolej `free_*_result()`)
- **Použij numerickou toleranci** pro double (`WITHIN`)
- **Piš popisné názvy testů** (`test_grid_with_outlier` ne `test1`)
- **Komentuj neobvyklé případy** v testech

### ❌ DON'T (Nedělej):

- **Neporovnávej floaty pomocí `==`** (použij `WITHIN`)
- **Nezapomeň volat free** (memory leak!)
- **Netestuj více věcí najednou** (jeden test = jedna věc)
- **Nepřepisuj setUp/tearDown** pokud to není potřeba
- **Nedělej testy závislé** (každý test musí být samostatný)

## Další kroky

V budoucnu by bylo dobré přidat:

- [ ] Testy pro `savgol.c`
- [ ] Testy pro `tikhonov.c`
- [ ] Testy pro `butterworth.c`
- [ ] Testy pro `polyfit.c`
- [ ] Integrační testy (celý pipeline)
- [ ] Performance benchmarky
- [ ] Coverage report (gcov/lcov)
- [ ] CI/CD integrace (GitHub Actions)

## Reference

- [Unity Framework](https://github.com/ThrowTheSwitch/Unity)
- [Unity Getting Started](http://www.throwtheswitch.org/unity)
- [Test-Driven Development in C](http://www.throwtheswitch.org/build/tdd)

#!/usr/bin/env python3
#
# cte data ze souboru nebo std. vstupu ve formatu (x[min], y[])
# znak # znamena komentar az do konce radky
#
# Použití:
#   ./plot_data.py soubor1.txt              # jeden soubor
#   ./plot_data.py soubor1.txt soubor2.txt  # dva soubory (různé barvy)
#   cat data.txt | ./plot_data.py           # stdin
import sys
import numpy as np
import matplotlib.pyplot as plt

# Určení zdroje dat
if len(sys.argv) == 1:
    # Stdin
    sources = [sys.stdin]
    labels = ['Data']
    colors = ['-b.']
elif len(sys.argv) == 2:
    # Jeden soubor
    sources = [sys.argv[1]]
    labels = [sys.argv[1]]
    colors = ['-b.']
elif len(sys.argv) == 3:
    # Dva soubory
    sources = [sys.argv[1], sys.argv[2]]
    labels = [sys.argv[1], sys.argv[2]]
    colors = ['-b.', '-r.']  # modrá a červená
else:
    print("Použití: ./plot_data.py [soubor1] [soubor2]", file=sys.stderr)
    sys.exit(1)

# Vytvoření grafu (poměr zlatého řezu)
phi = (1 + 5**0.5) / 2  # zlatý řez ≈ 1.618
width = 10
height = width / phi  # ≈ 6.18
plt.figure(figsize=(width, height))

# Načtení a vykreslení dat pro každý zdroj
for source, label, color in zip(sources, labels, colors):
    # Načtení dat (# = komentář)
    data = np.loadtxt(source, comments='#')

    # Rozdělení na x a y
    x = data[:, 0]
    y = data[:, 1]

    # Vykreslení
    # '-b.' = plná modrá čára s body, '-r.' = plná červená čára s body
    plt.plot(x, y, color, linewidth=1, markersize=3, label=label)

plt.xlabel('Čas [min]')
plt.ylabel('Měřená veličina')
# alpha = průhlednost mřížky (0=neviditelné, 1=plná barva, 0.3=lehce průhledné)
plt.grid(True, alpha=0.3)
plt.title('Závislost měřené veličiny na čase')
# Zobrazit legendu pouze pokud jsou dva soubory
if len(sources) > 1:
    plt.legend()
# tight_layout() automaticky upraví rozložení, aby se všechny prvky vešly bez ořezání
plt.tight_layout()

# Zobrazení grafu
plt.show()

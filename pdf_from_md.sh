#!/bin/bash
# Convert README.md to PDF using pandoc with Unicode-compatible fonts
# Uses DejaVu fonts for proper rendering of mathematical symbols (λ, Σ, π, ω, etc.)
# Note: Box-drawing and decorative characters have been replaced with ASCII in README.md

pandoc README.md -o README.pdf \
  --pdf-engine=xelatex \
  --variable geometry:margin=2.5cm \
  --variable fontsize=11pt \
  --variable colorlinks=true \
  --variable linkcolor=blue \
  --variable urlcolor=blue \
  --variable mainfont="DejaVu Serif" \
  --variable monofont="DejaVu Sans Mono"

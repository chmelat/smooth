#!/usr/bin/bash
# Convert README.md to PDF using pandoc with Unicode-compatible fonts
# Uses DejaVu fonts for proper rendering of mathematical symbols (λ, Σ, π, ω, etc.)
# Note: Box-drawing and decorative characters have been replaced with ASCII in README.md

case $# in
  0) echo "Usage: $0 file1.md file2.md ..." 1>&2; exit 2
esac

for i
do

  name="$i"

  if [[ "$name" != *.md ]]; then
    echo "${name} don't like md file!"
    exit 1
  fi

  new_name=$(echo "$name" | sed 's/md$/pdf/')

  echo "Convert ${name} -> ${new_name}"

  pandoc "$name" -o "$new_name" \
    --pdf-engine=xelatex \
    --variable geometry:margin=2.5cm \
    --variable fontsize=11pt \
    --variable colorlinks=true \
    --variable linkcolor=blue \
    --variable urlcolor=blue \
    --variable mainfont="DejaVu Serif" \
    --variable monofont="DejaVu Sans Mono"

done

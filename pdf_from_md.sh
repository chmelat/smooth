#!/bin/bash
pandoc README.md -o README.pdf --pdf-engine=xelatex   --variable geometry:margin=2.5cm  --variable fontsize=11pt   --variable colorlinks=true   --variable linkcolor=blue   --variable urlcolor=blue

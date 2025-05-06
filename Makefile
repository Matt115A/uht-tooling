# Makefile for uht-tooling

# Variables
PYTHON=python
NEXTERA_PRIMERS=scripts/nextera_designer.py
UMI_HUNTER=scripts/umi_hunter.py
DESIGN_SLIM=scripts/design_slim.py

# Rules

.PHONY: all nextera_primers umi_hunter design_slim

all: nextera_primers umi_hunter design_slim

nextera_primers:
	python -m scripts.nextera_designer

umi_hunter:
	python -m scripts.umi_hunter

design_slim:
	python -m scripts.design_slim

clean:
	rm -rf results/*.png
	rm -rf results/*.csv
	echo "Cleaned generated files."

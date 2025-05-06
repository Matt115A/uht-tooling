# Makefile for uht-tooling

# Variables
PYTHON=python
NEXTERA_PRIMERS=scripts/nextera_designer.py
UMI_HUNTER=scripts/umi_hunter.py

# Rules

.PHONY: all nextera_primers umi_hunter

all: nextera_primers umi_hunter

nextera_primers:
	python -m scripts.nextera_designer

umi_hunter:
	python -m scripts.umi_hunter

clean:
	rm -rf results/*.png
	rm -rf results/*.csv
	echo "Cleaned generated files."

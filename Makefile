# Makefile for uht-tooling

# Variables
PYTHON=python
NEXTERA_PRIMERS=scripts/nextera_designer.py

# Rules

.PHONY: all nextera_primers

all: nextera_primers

nextera_primers:
	python -m scripts.nextera_designer

clean:
	rm -rf results/*.png
	rm -rf results/*.csv
	echo "Cleaned generated files."

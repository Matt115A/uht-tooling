# Makefile for uht-tooling

# Variables
PYTHON=python
NEXTERA_PRIMERS=scripts/nextera_designer.py
UMI_HUNTER=scripts/umi_hunter.py
DESIGN_SLIM=scripts/design_slim.py
MUTATION_CALLER=scripts/mutation_caller.py

# Rules

.PHONY: all nextera_primers umi_hunter design_slim mutation_caller

all: nextera_primers umi_hunter design_slim mutation_caller

nextera_primers:
	python -m scripts.nextera_designer

umi_hunter:
	python -m scripts.umi_hunter

design_slim:
	python -m scripts.design_slim

mutation_caller:
	python -m scripts.mutation_caller

ep-library-profile:
	python -m scripts.mut_rate

design_gibson:
	python -m scripts.design_gibson

clean:
	rm -rf results/*.png
	rm -rf results/*.csv
	echo "Cleaned generated files."

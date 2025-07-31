# Makefile for uht-tooling

# Variables
PYTHON=python
NEXTERA_PRIMERS=scripts/nextera_designer.py
UMI_HUNTER=scripts/umi_hunter.py
DESIGN_SLIM=scripts/design_slim.py
MUTATION_CALLER=scripts/mutation_caller.py

# Rules

.PHONY: all nextera_primers umi_hunter design_slim mutation_caller profile_inserts

all: nextera_primers umi_hunter design_slim mutation_caller profile_inserts

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

profile_inserts:
	python -m scripts.profile_inserts

clean:
	rm -rf results/*.png
	rm -rf results/*.csv
	echo "Cleaned generated files."

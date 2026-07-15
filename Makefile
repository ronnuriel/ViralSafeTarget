.PHONY: install test lint demo ui notebook real-hsv2 reproduce-hsv2 clean clean-generated

install:
	python -m pip install -e .[all]

test:
	pytest

lint:
	ruff check src tests scripts

demo:
	bash scripts/run_demo.sh

ui:
	streamlit run app.py

notebook:
	jupyter lab

real-hsv2:
	bash scripts/run_real_hsv2.sh --sample-size 25

reproduce-hsv2:
	vst reproduce hsv2

clean:
	rm -rf .pytest_cache .ruff_cache htmlcov
	find src tests scripts -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -name '.DS_Store' -delete

# Remove disposable workflow outputs while preserving the checked-in public snapshots
# under reports/hsv2_{showcase,genome_wide_exhaustive,evidence_agent,
# virtual_knockout_escape,tool_benchmark}.
clean-generated: clean
	rm -rf reports/demo reports/real_hsv2 reports/hsv2_consensus
	rm -rf reports/hsv2_genome_wide reports/hsv2_gene_function reports/hsv2_project
	rm -rf reports/hsv2_population_heldout reports/hsv2_population_report_balanced
	rm -rf data/processed/*
	touch data/processed/.gitkeep

.PHONY: install test lint demo ui notebook real-hsv2 reproduce-hsv2 clean

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
	rm -rf reports/* data/processed/* .pytest_cache .ruff_cache
	touch reports/.gitkeep data/processed/.gitkeep

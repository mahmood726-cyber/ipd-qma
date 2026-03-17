.PHONY: all data simulate analyze manuscript test clean

all: data simulate analyze manuscript

data:
	python fetch_nhanes.py

simulate:
	python run_simulation_fast.py

analyze:
	python run_analysis.py

manuscript: analyze
	python build_docx.py

test:
	python -m pytest test_ipd_qma.py -v

clean:
	rm -rf output/*.png output/*.csv __pycache__

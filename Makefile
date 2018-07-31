all: week2

week2:
	@echo Running the Week 2 notebook
	cd Week2-GraphTheory && jupyter nbconvert --ExecutePreprocessor.timeout=-1 --inplace --execute GraphTheoryLab.ipynb

getrequirements:
	pip freeze | grep -v conda > requirements.txt

.PHONY: clean docs help
.DEFAULT_GOAL := help

## Remove all build and Python artifacts
clean: clean-build clean-pyc

## Remove build artifacts
clean-build:
	rm -rf build/ dist/ .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

## Remove Python artifacts
clean-pyc:
	find . -type f -name '*.py[co]' -delete
	find . -type d -name '__pycache__' -delete

qa:
	isort --recursive src/sanger_sequencing tests
	black src/sanger_sequencing tests

## Generate Sphinx HTML documentation
docs:
	rm -f docs/sanger-sequencing.rst docs/modules.rst
	sphinx-apidoc -o docs/ src/sanger_sequencing
	$(MAKE) -C docs clean
	$(MAKE) -C docs html

release: clean
	python setup.py sdist bdist_wheel
	twine upload --skip-existing dist/*

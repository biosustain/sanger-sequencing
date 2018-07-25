.PHONY: clean docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import sys
import webbrowser
from os.path import abspath
try:
	from urllib import pathname2url
except ImportError:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

## Remove all build and Python artifacts
clean: clean-build clean-pyc clean-test

## Remove build artifacts
clean-build:
	rm -rf build/ dist/ .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

## Remove Python artifacts
clean-pyc:
	find . -type f -name '*.py[co]' -delete
	find . -type d -name '__pycache__' -delete

## Generate Sphinx HTML documentation
docs:
	rm -f docs/sanger-sequencing.rst docs/modules.rst
	sphinx-apidoc -o docs/ sanger-sequencing
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

## Compile the docs watching for changes
servedocs: docs
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

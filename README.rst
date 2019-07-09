========================
Sanger Sequence Analysis
========================

.. image:: https://img.shields.io/pypi/v/sanger-sequencing.svg
   :target: https://pypi.org/project/sanger-sequencing/
   :alt: Current PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/sanger-sequencing.svg
   :target: https://pypi.org/project/sanger-sequencing/
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/l/sanger-sequencing.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0
   :alt: Apache Software License Version 2.0

.. image:: https://img.shields.io/travis/biosustain/sanger-sequencing/master.svg?label=Travis%20CI
   :target: https://travis-ci.org/biosustain/sanger-sequencing
   :alt: Travis CI

.. image:: https://codecov.io/gh/biosustain/sanger-sequencing/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/biosustain/sanger-sequencing
   :alt: Codecov

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Black

.. summary-start

Semi-automated Sanger sequence analysis for plasmid verification.

This package is the result of an internal hackathon at the Novo Nordisk
Foundation Center for Biosustainability and represents our approach to improving
the workflow of geneticists who need to verify plasmid constructs by Sanger
sequencing.

Getting Started
===============

From a Python environment that has Python 3.7 or later installed you can easily

.. code-block:: console

    $ pip install sanger-sequencing[analysis]

or use ``pip3`` depending on your environment.

When you import the package, two main components are made available to you: a
configuration class that you can instantiate to set some global configuration
values and a high level analysis interface.

.. code-block:: python

    import sanger_sequencing

    config = sanger_sequencing.Configuration()
    print(config.threshold)
    print(config.output)

You can read more about the meaning of those attributes in the configuration
documentation. The main entry point for doing any kind of analysis is the
``sanger_verification`` function. This function requires three arguments: a
template table of what to analyze, a mapping from plasmid identifiers to their
sequence records (typically coming from Genbank files), and a mapping from
sample identifiers to sequence records (``.ab1`` files).

.. summary-end

You can find the complete documentation at: https://sanger-sequencing.readthedocs.io.

Copyright
=========

* Copyright © 2018-2019, Novo Nordisk Foundation Center for Biosustainability,
  Technical University of Denmark
* Free software distributed under the `Apache Software License 2.0
  <https://www.apache.org/licenses/LICENSE-2.0>`_.

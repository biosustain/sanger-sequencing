========================
Sanger Sequence Analysis
========================

.. image:: https://img.shields.io/pypi/v/sanger-sequencing.svg
        :target: https://pypi.python.org/pypi/sanger-sequencing

.. image:: https://img.shields.io/travis/biosustain/sanger-sequencing.svg
        :target: https://travis-ci.org/biosustain/sanger-sequencing

.. image:: https://readthedocs.org/projects/sanger-sequencing/badge/?version=latest
        :target: https://sanger-sequencing.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/biosustain/sanger-sequencing/shield.svg
     :target: https://pyup.io/repos/github/biosustain/sanger-sequencing/
     :alt: Updates

.. summary-start

Semi-automated Sanger sequence analysis for plasmid verification.

This package is the result of an internal hackathon at the Novo Nordisk 
Foundation Center for Biosustainability and represents our approach to 
improving the workflow of geneticists who need to verify plasmid 
constructs by Sanger sequencing.

Getting Started
===============

From a Python environment that has Python 3.6 or later installed you can easily 

.. code-block:: console

    $ pip install sanger-sequencing

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
template table of what to analyze, a mapping from plasmid identifiers to 
their sequence records (typically coming from Genbank files), and a mapping 
from sample identifiers to sequence records (``.ab1`` files).
    
.. summary-end

You can find the complete documentation at: https://sanger-sequencing.readthedocs.io.

Copyright
=========

* Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability, Technical University Denmark licensed
  under the Apache License, Version 2.0

Credits
=======

This package was created using cookiecutter_ and the 
`DD-DeCaF/cookiecutter-decaf-python`_ project template.

.. _cookiecutter: https://github.com/audreyr/cookiecutter
.. _`DD-DeCaF/cookiecutter-decaf-python`: https://github.com/DD-DeCaF/cookiecutter-decaf-python


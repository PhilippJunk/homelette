Welcome to homelette's documentation!
=====================================

homelette is a Python package offering a unified interface to different software for generating and evaluating homology models. This enables users to easily assemble custom homology modelling pipelines. homelette is extensively documented, lightweight and easily extendable.

.. image:: logo.png
   :width: 400


Setting up homelette
====================

This section explains how to set homelette up on your system. homelette is available on GitHub and PyPI. The easiest option to work with homelette is to use a docker container that has all dependencies already installed.

.. toctree::
   :maxdepth: 1

   Installation <installation>
   Docker <docker>


Tutorials
=========

We have prepared a series of 7 tutorials which will teach the interested user everything about using the homelette package. This is a great place to get started with homelette.

For a more interactive experience, all tutorials are available as Jupyter Notebooks through our :ref:`Docker container<access_docker>`.

.. toctree::
   :maxdepth: 1

   Tutorial 1: Basics <Tutorial1_Basics.ipynb>
   Tutorial 2: Model Generation <Tutorial2_Modelling.ipynb>
   Tutorial 3: Model Evaluation <Tutorial3_Evaluation.ipynb>
   Tutorial 4: Extending homelette <Tutorial4_ExtendingHomelette.ipynb>
   Tutorial 5: Parallelization <Tutorial5_Parallelization.ipynb>
   Tutorial 6: Complex Modelling <Tutorial6_ComplexModelling.ipynb>
   Tutorial 7: Assembling Modelling Pipelines <Tutorial7_AssemblingPipelines.ipynb>
   Tutorial 8: Automatic Alignment Generation <Tutorial8_AlignmentGeneration.ipynb>


API Documentation
=================

This is the documentation for all classes, methods and functions in homelette.

.. toctree::
   :maxdepth: 1

   Organization <organization>
   Sequences and Alignments <alignment>
   Routines <routines>
   Evaluation <evaluation>
   PDB File interface <pdb_io>


Extensions
==========

homelette can be extended by new building blocks. This section introduces how extensions work, and where to find them. 

.. toctree::
   :maxdepth: 1

   Extension Overview <extension>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

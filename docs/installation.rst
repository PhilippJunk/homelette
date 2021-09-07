Installation
============

While installing the homelette base package is easy, some of its dependencies are quite complicated to install. If you just want to try out homelette, we would encourage you to start with our :ref:`Docker image<docker>` which has all these dependencies already installed.

homelette
---------

homelette is easily available through our GitHub page (`GitHub homelette <https://github.com/PhilippJunk/homelette/>`_) or through PyPI.

.. code-block:: bash

   python3 -m pip install homelette

Please be aware that homelette requires Python 3.6.12 or newer.


Modelling and Evaluation Software
---------------------------------

homelette doesn't have model generating or model evaluating capabilities on its own. Instead, it provides a unified interface to other software with these capabilities. 

None of the tools and packages listed here are "hard" dependencies in the way that homelette won't work if you have them not installed. Actually, you can still use homelette without any of these packages. However, none of the pre-implemented building blocks would work that way.  It is therefore strongly recommended that, in order to get the most out of homelette, to install as many of these tools and packages.

Again, we want to mention that we have prepared a :ref:`Docker image<docker>` that contains all of these dependencies, and we strongly recommend that you start there if you want to find out if homelette is useful for you. 

MODELLER
^^^^^^^^

Installation instructions for MODELLER can be found here: `Installation MODELLER <https://salilab.org/modeller/download_installation.html>`_. 
Requires a license key (freely available for academic research) which can be requested here: `License MODELLER <https://salilab.org/modeller/registration.html>`_.


altMOD
^^^^^^

altMOD can be installed from here: `GitHub altMOD <https://github.com/pymodproject/altmod>`_. Please make sure that the altMOD directory is in your Python path.


ProMod3
^^^^^^^

ProMod3 has to be compiled from source, instructions can be found here: `Installation ProMod3 <https://openstructure.org/promod3/>`_. Main dependencies are OpenMM (available through conda or from source) and OpenStructure (available here: `Installation OpenStructure <https://openstructure.org/download/>`_).


QMEAN
^^^^^

QMEAN has be compiled from source, instructions can be found here: `GitLab QMEAN <https://git.scicore.unibas.ch/schwede/QMEAN/>`_. Has the same dependencies as ProMod3.


SOAP potential
^^^^^^^^^^^^^^

While the code for evaluation with SOAP is part of MODELLER, some files for SOAP are not included in the standard release and have to be downloaded separately. The files are available here `Download SOAP <https://salilab.org/SOAP/>`_. 

Specifically, you need to have ``soap_protein_od.hdf5`` available in your modlib directory. The modlib directory is placed at ``/usr/lib/modellerXX_XX/modlib/`` if installed with ``dpkg`` or at ``anaconda/envs/yourenv/lib/modellerXX-XX/modlib/`` if installed with ``conda``. These paths might be different on your system.

MolProbity
^^^^^^^^^^

Installation instructures for MolProbity are available here: `Github MolProbity <https://github.com/rlabduke/MolProbity>`_. Pleas make sure that after installation, ``phenix.molprobity`` is in your path.

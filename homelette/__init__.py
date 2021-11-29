'''
homelette
=========

homelette is an interface for various homology modelling tools, enabling the
user to easily assemble custom homology modelling pipelines.

Please check out the documentation and tutorials at
https://homelette.readthedocs.io/.

The docstring examples assume that `homelette` has been imported as `hm`. Code
snipets are indicated by three greater-than signs::

    >>> import homelette as hm

Use the build-in ``help`` function to view the docstring of a function or
class::

    >>> help(hm.Task)

Available subpackages
---------------------

organization
    Classes for organizing workflows and models
alignment
    Classes and functions for handling multiple sequence alignments
routines
    Classes for homology model generation
evaluation
    Classes for homology model evaluation
pdb_io
    Interface for handling and modifying PDB files
extension
    Interface for extending `homelette`
'''

__all__ = ['Task', 'Model', 'Alignment', 'routines', 'evaluation']
__version__ = '1.3'
__author__ = 'Philipp Junk, Christina Kiel'
__email__ = 'philipp.junk@ucdconnect.ie'
__maintainer__ = 'Philipp Junk'
__license__ = 'MIT'

# Standard library imports
import warnings

# Local application imports
from .organization import Task, Model
from .alignment import Alignment
from . import routines
from . import evaluation


# Check third party imports and report missing modules
def _check_imports() -> None:
    '''
    Helper function that checks third-party imports and raises warnings if they
    could not be imported.

    Returns
    -------
    None
    '''
    for (module, imported) in sorted(_IMPORTS.items()):
        if not imported:
            msg = 'Module "{}" could not be imported.'.format(module)
            warnings.warn(msg)

    if not all(_IMPORTS.values()):
        msg = ('Please install the missing modules in order to enjoy the full '
               'functionality of "homology"')
        warnings.warn(msg)


# gather imports from submodules
_IMPORTS = {**evaluation._IMPORTS, **routines._IMPORTS}
_check_imports()

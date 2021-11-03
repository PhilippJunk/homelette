'''
``homelette.organization``
==========================

The :mod:`homelette.organization` submodule contains classes for organizing
workflows.

:class:`Task` is an object orchestrating model generation and evaluation.

:class:`Model` is an object used for storing information about generated
models.

Tutorials
---------

For an introduction to `homelette`'s workflow, :ref:`Tutorial
1</Tutorial1_Basics.ipynb>` is useful.
Assembling custom pipelines is discussed in :ref:`Tutorial
7</Tutorial7_AssemblingPipelines.ipynb>`.


Classes
-------

The following classes are part of this submodule:

    :class:`Task`
    :class:`Model`

------

'''

__all__ = ['Task', 'Model']

# Standard library imports
import concurrent.futures
import contextlib
import glob
import os
import os.path
import shutil
import time
import typing

# Third party imports
import pandas as pd

# Local application imports
from . import pdb_io

# Local imports for type checking
# taken from https://stackoverflow.com/a/39757388/7912251
if typing.TYPE_CHECKING:
    from .alignment import Alignment
    from . import routines
    from . import evaluation


class Task():
    '''
    Class for directing modelling and evaluation.

    It is designed for the modelling of one target sequence from one or
    multiple templates.

    If an already existing folder with models is specified, the Task object
    will load those models in automatically. In this case, it can also be used
    exclusively for evaluation purposes.

    Parameters
    ----------
    task_name : str
        The name of the task
    target : str
        The identifier of the protein to model
    alignment : Alignment
        The alignment object that will be used for modelling
    task_directory : str, optional
        The directory that will be used for this modelling task (default is
        creating a new one based on the task_name)
    overwrite : bool, optional
        Boolean value determining if an already existing task_directory
        should be overwriten. If a directory already exists for a given
        task_name or task_directory, this will determine whether the
        directory and all its contents will be overwritten (True), or
        whether the contained models will be imported (False) (default is
        False)

    Attributes
    ----------
    task_name : str
        The name of the task
    task_directory : str
        The directory that will be used for this modelling task (default is to
        use the task_name)
    target : str
        The identifier of the protein to model
    alignment : Alignment
        The alignment object that will be used for modelling
    models : list
        List of models generated or imported by this task
    routines : list
        List of modelling routines executed by this task

    Returns
    -------
    None
    '''
    def __init__(self, task_name: str, target: str, alignment:
                 typing.Type['Alignment'], task_directory: str = None,
                 overwrite: bool = False) -> None:
        # organization settings
        self.task_name = task_name
        if task_directory is None:
            self.task_directory = os.path.realpath(os.path.expanduser(
                os.path.join(os.getcwd(), task_name)))
        else:
            self.task_directory = os.path.realpath(os.path.expanduser(
                task_directory))

        # modelling settings
        self.target = target
        self.alignment = alignment
        self.models = []
        self.routines = []

        # initialize directory for task
        if os.path.isdir(self.task_directory) and overwrite is False:
            # import models
            for model_file in glob.glob(os.path.join(
                    self.task_directory, '*.pdb')):
                self.models.append(Model(os.path.realpath(os.path.expanduser(
                    model_file)), None, None))
            if len(self.models) != 0:
                print(
                    'Imported {} models from already existing '
                    'task_directory.'.format(len(self.models)))
        elif os.path.isdir(self.task_directory) and overwrite is True:
            shutil.rmtree(self.task_directory)
            os.mkdir(self.task_directory)
        else:
            os.mkdir(self.task_directory)

    @contextlib.contextmanager
    def _cwd_task_folder(self) -> None:
        '''
        Helper function: context manager for executing model generation and
        evaluation inside the task directory
        '''
        # adapted from https://stackoverflow.com/a/37996581/7912251
        curdir = os.getcwd()
        os.chdir(self.task_directory)
        try:
            yield
        finally:
            os.chdir(curdir)

    def execute_routine(self, tag: str, routine:
                        typing.Type['routines.Routine'], templates:
                        typing.Iterable, template_location: str = '.',
                        **kwargs) -> None:
        '''
        Generates homology models using a specified modelling routine

        Parameters
        ----------
        tag : str
            The identifier associated with this combination of routine and
            template(s). Has to be unique between all routines executed by the
            same task object
        routine : Routine
            The routine object used to generate the models
        templates : list
            The iterable containing the identifier(s) of the template(s) used
            for model generation
        template_location : str, optional
            The location of the template PDB files. They should be named
            according to their identifiers in the alignment (i.e. for a
            sequence named "1WXN" to be used as a template, it is expected that
            there will be a PDB file named "1WXN.pdb" in the specified template
            location (default is current working directory)
        **kwargs
            Named parameters passed directly on to the Routine object when the
            modelling is performed. Please check the documentation in order to
            make sure that the parameters passed on are available with the
            Routine object you intend to use

        Returns
        -------
        None
        '''
        template_location = os.path.realpath(
            os.path.expanduser(template_location))
        with self._cwd_task_folder():
            # retrieve template PDB files
            if os.path.realpath(os.getcwd()) == template_location:
                def rm_templates():
                    # helper function that cleans up templates if copied
                    pass
            else:
                for template in templates:
                    shutil.copy(os.path.join(
                        template_location, template + '.pdb'), '.')

                def rm_templates():
                    # helper function that cleans up templates if copied
                    for template in templates:
                        os.remove(template + '.pdb')

            # initialize routine
            try:
                r = routine(self.alignment, self.target, templates, tag,
                            **kwargs)
            except TypeError as err:
                rm_templates()  # rm templates before raising exception
                # extract which routine was called
                routine_class = routine.__name__
                raise TypeError('Routine {0} does not recognize one of the '
                                'keywords you specified. Please review '
                                'help({0}) for information which keywords are '
                                'applicable.\n\nOriginal Error Message:'
                                '\n{1}'.format(routine_class, err))
            except Exception:
                rm_templates()
                raise

            # execute routine
            try:
                r.generate_models()
            except Exception:
                rm_templates()  # rm templates before raising exception
                raise

            # append models
            self.routines.append(r)
            self.models = self.models + r.models

            # clean up templates
            rm_templates()

    def evaluate_models(self, *args: typing.Type['evaluation.Evaluation'],
                        n_threads: int = 1) -> None:
        '''
        Evaluates models using one or multiple evaluation metrics

        Parameters
        ----------
        *args: Evaluation
            Evaluation objects that will be applied to the models
        n_threads : int, optional
            Number of threads used for model evaluation (default is 1, which
            deactivates parallelization)

        Returns
        -------
        None
        '''
        # construct worker functions that runs all evaluation
        # Because of weird interactions of using contextlib in conjunction with
        # threaded applications, it is stronly recommended that all Evaluation
        # objects have redirection of stdout turned off.
        def worker(model):
            for eval_method in args:
                eval_method(model, quiet=False)

        # construct second worker function specifically for parallel execution
        def worker_threaded(model, i):
            # make sure threads start with slight delay in order to avoid any
            # race conditions
            if i <= n_threads:
                time.sleep(i)
            else:
                time.sleep(0.5)
            worker(model)

        # execute evaluations
        with contextlib.redirect_stdout(None), self._cwd_task_folder():
            if n_threads > 1:
                with concurrent.futures.ThreadPoolExecutor(n_threads) as pool:
                    futures = list()
                    for i in range(len(self.models)):
                        futures.append(
                            pool.submit(worker_threaded,
                                        self.models[i], i))
                    pool.shutdown()
                # check for exceptions suppressed by threaded execution
                for future in concurrent.futures.as_completed(futures):
                    if future.exception() is not None:
                        raise future.exception()
            else:
                for model in self.models:
                    worker(model)

    def get_evaluation(self) -> pd.DataFrame:
        '''
        Return evaluation for all models as pandas dataframe.

        Returns
        -------
        pd.DataFrame
            Dataframe containing all model evaluation
        '''
        return pd.DataFrame(
                [m.evaluation for m in self.models])


class Model():
    '''
    Interface used to interact with created protein structure models.

    Parameters
    ----------
    model_file : str
        The file location of the PDB file for this model
    tag : str
        The tag that was used when generating this model (see
        ``Task.execute_routine`` for more details)
    routine : str
        The name of the routine that was used to generate this model

    Attributes
    ----------
    model_file : str
        The file location of the PDB file for this model
    tag : str
        The tag that was used when generating this model (see
        Task.execute_routine for more details)
    routine : str
        The name of the routine that was used to generate this model
    info : dict
        Dictionary that can be used to store metadata about the model (i.e. for
        some evaluation metrics)

    Returns
    -------
    None
    '''
    def __init__(self, model_file: str, tag: str, routine: str) -> None:
        self.model_file = os.path.realpath(os.path.expanduser(model_file))
        self.tag = tag
        self.routine = routine

        self.info = dict()

        # initialize evaluation output
        self.evaluation = {
            'model': os.path.basename(self.model_file),
            'tag': self.tag,
            'routine': self.routine}

    def parse_pdb(self) -> pd.DataFrame:
        '''
        Parses ATOM and HETATM records in PDB file to pandas dataframe
        Useful for giving some evaluations methods access to data from the PDB
        file.

        Returns
        -------
        pd.DataFrame

        Notes
        -----
        Information is extracted according to the PDB file specification
        (version 3.30) and columns are named accordingly. See
        https://www.wwpdb.org/documentation/file-format for more information.
        '''
        return pdb_io.read_pdb(self.model_file).parse_to_pd()

    def get_sequence(self) -> str:
        '''
        Retrieve the 1-letter amino acid sequence of the PDB file associated
        with the Model object.

        Returns
        -------
        str
            Amino acid sequence
        '''
        return pdb_io.read_pdb(self.model_file).get_sequence()

    def rename(self, new_name: str) -> None:
        '''
        Rename the PDB file associated with the Model object.

        Parameters
        ----------
        new_name : str
            New name of PDB file

        Returns
        -------
        None
        '''
        if not new_name[-4:] == '.pdb':
            new_name = new_name + '.pdb'
        new_model_file = os.path.join(
            os.path.dirname(self.model_file), new_name)
        os.rename(self.model_file, new_model_file)
        self.model_file = new_model_file
        self.evaluation['model'] = os.path.basename(self.model_file)

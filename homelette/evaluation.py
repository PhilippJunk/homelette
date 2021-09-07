'''
``homelette.evaluation``
========================

The :mod:`homelette.evaluation` submodule contains different classes for
evaluating homology models.

It is possible to implement custom Evaluation building blocks and use them in
the `homelette` framework.

Tutorials
---------

Working with model evaluations in `homelette` is discussed in detail in
:ref:`Tutorial 3 </Tutorial3_Evaluation.ipynb>`. Implementing custom evaluation
metrics is discussed in
:ref:`Tutorial 4 </Tutorial4_ExtendingHomelette.ipynb>`.
Assembling custom pipelines is discussed in :ref:`Tutorial 7
</Tutorial7_AssemblingPipelines.ipynb>`.

Classes
-------

The following evaluation metrics are implemented:

    :class:`Evaluation_dope`
    :class:`Evaluation_soap_protein`
    :class:`Evaluation_soap_pp`
    :class:`Evaluation_qmean4`
    :class:`Evaluation_qmean6`
    :class:`Evaluation_qmeandisco`
    :class:`Evaluation_mol_probity`

------

'''

__all__ = ['Evaluation_dope', 'Evaluation_soap_protein', 'Evaluation_soap_pp',
           'Evaluation_qmean4', 'Evaluation_qmean6', 'Evaluation_qmeandisco',
           'Evaluation_mol_probity']

# Standard library imports
import contextlib
import re
import subprocess
import typing
# import warnings

# Third party imports
_IMPORTS = dict()
try:
    import modeller
    import modeller.scripts
    import modeller.soap_protein_od
    import modeller.soap_pp
    _IMPORTS['modeller'] = True
except ImportError:
    _IMPORTS['modeller'] = False

try:
    import ost
    _IMPORTS['ost'] = True
except ImportError:
    _IMPORTS['ost'] = False

try:
    import qmean
    _IMPORTS['qmean'] = True
except ImportError:
    _IMPORTS['qmean'] = False

# Local application imports

# Local imports for type checking
# taken from https://stackoverflow.com/a/39757388/7912251
if typing.TYPE_CHECKING:
    from .organization import Model


class Evaluation():
    '''
    Parent class to all evaluation procedures.

    Not supposed to be used by user, used for inheritance. Implements a few
    common attributes shared by all Evaluation objects.
    '''
    # set alibi values
    def __init__(self, model: typing.Type['Model']) -> None:
        self.model = model
        self.output = dict()

    def evaluate(self):
        '''
        Placeholder. Will be implemented by children classes
        '''

    def _evaluate(self, quiet: bool):
        '''
        Performs evaluation. Helper function for children classes.

        Parameters
        ----------
        quiet : bool
            If true, redirect stdout to None
        '''
        if quiet:
            with contextlib.redirect_stdout(None):
                self.evaluate()
        else:
            self.evaluate()
        # update model
        self.model.evaluation.update(self.output)

    @staticmethod
    def _check_dependencies(dependencies: typing.Iterable) -> None:
        '''
        Checks if all dependencies could be loaded. Helper function for
        children classes.

        Parameters
        ----------
        dependencies : Iterable
            Iterable of dependencies

        Returns
        -------
        None

        Raises
        ------
        ImportError
            One or more dependencies were not imported
        '''
        for dependency in dependencies:
            if not _IMPORTS[dependency]:
                raise ImportError(
                    '"{}" is required for this functionality, but could '
                    'not be imported.'.format(dependency))


class Evaluation_dope(Evaluation):
    '''
    Class for evaluating a model with DOPE score.

    Will dump the following entries to the model.evaluation dictionary:

    * dope
    * dope_z_score

    Parameters
    ----------
    model : Model
        The model object to evaluate
    quiet : bool
        If True, will perform evaluation with suppressing stdout (default
        False). Needs to be False for running it asynchronously, as done
        when running Task.evaluate_models with multple cores

    Attributes
    ----------
    model : Model
        The model object to evaluate
    output : dict
        Dictionary that all outputs will be dumped into

    Raises
    ------
    ImportError
        Unable to import dependencies

    Notes
    -----
    DOPE is a staticial potential for the evaluation of homology models [1]_.
    For further information, please check the modeller documentation or the
    associated publication.

    References
    ----------
    .. [1] Shen, M., & Sali, A. (2006). Statistical potential for assessment
       and prediction of protein structures. Protein Science, 15(11),
       2507–2524. https://doi.org/10.1110/ps.062416606
    '''
    def __init__(self, model: typing.Type['Model'], quiet: bool =
                 False) -> None:
        # initialize attributes
        Evaluation.__init__(self, model)
        # check dependencies
        _dependencies = ['modeller']
        self._check_dependencies(_dependencies)
        # execute evaluation
        self._evaluate(quiet)

    def evaluate(self) -> None:
        '''
        Run DOPE evaluation. Automatically called on object initialization

        Returns
        -------
        None
        '''
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        mdl = modeller.scripts.complete_pdb(env, self.model.model_file)
        sele = modeller.selection(mdl.chains[:])
        # set up output
        self.output['dope'] = sele.assess_dope()
        self.output['dope_z_score'] = mdl.assess_normalized_dope()


class Evaluation_soap_protein(Evaluation):
    '''
    Class for evaluating a model with the SOAP protein protential.

    Will dump the following entries to the model.evaluation dictionary:

    * soap_protein

    Parameters
    ----------
    model : Model
        The model object to evaluate
    quiet : bool
        If True, will perform evaluation with suppressing stdout (default
        False). Needs to be False for running it asynchronously, as done when
        running Task.evaluate_models with multple cores

    Attributes
    ----------
    model : Model
        The model object to evaluate
    output : dict
        Dictionary that all outputs will be dumped into

    Raises
    ------
    ImportError
        Unable to import dependencies

    Notes
    -----
    SOAP is a statistical potential for evaluating homology models [2]_. For
    more information, please check the modeller and SOAP documentations or the
    associated publication.

    References
    ----------
    .. [2] Dong, G. Q., Fan, H., Schneidman-Duhovny, D., Webb, B., Sali, A., &
       Tramontano, A. (2013). Optimized atomic statistical potentials:
       Assessment of protein interfaces and loops. Bioinformatics, 29(24),
       3158–3166. https://doi.org/10.1093/bioinformatics/btt560
    '''
    def __init__(self, model: typing.Type['Model'], quiet: bool =
                 False) -> None:
        # initialize attributes
        Evaluation.__init__(self, model)
        # check dependencies
        _dependencies = ['modeller']
        self._check_dependencies(_dependencies)
        # execute evaluation
        self._evaluate(quiet)

    def evaluate(self) -> None:
        '''
        Run SOAP protein evaluation. Automatically called on object
        initialization

        Returns
        -------
        None
        '''
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        mdl = modeller.scripts.complete_pdb(env, self.model.model_file)
        sele = modeller.selection(mdl.chains[:])
        # set up output
        self.output['soap_protein'] = sele.assess(
            modeller.soap_protein_od.Scorer())


class Evaluation_soap_pp(Evaluation):
    '''
    Class for evaluating a model with SOAP interaction potentials. This is used
    for the evaluation of models of protein complexes.

    Will dump the following entries to the model.evaluation dictionary:

    * soap_pp_all
    * soap_pp_atom
    * soap_pp_pair

    Parameters
    ----------
    model : Model
        The model object to evaluate
    quiet : bool
        If True, will perform evaluation with suppressing stdout (default
        False). Needs to be False for running it asynchronously, as done when
        running Task.evaluate_models with multple cores

    Attributes
    ----------
    model : Model
        The model object to evaluate
    output : dict
        Dictionary that all outputs will be dumped into

    Raises
    ------
    ImportError
        Unable to import dependencies

    Notes
    -----
    SOAP is a statistical potential for evaluating homology models [3]_. For
    more information, please check the modeller and SOAP documentations or the
    associated publication.

    References
    ----------
    .. [3] Dong, G. Q., Fan, H., Schneidman-Duhovny, D., Webb, B., Sali, A., &
       Tramontano, A. (2013). Optimized atomic statistical potentials:
       Assessment of protein interfaces and loops. Bioinformatics, 29(24),
       3158–3166. https://doi.org/10.1093/bioinformatics/btt560
    '''
    def __init__(self, model: typing.Type['Model'], quiet: bool =
                 False) -> None:
        # initialize attributes
        Evaluation.__init__(self, model)
        # check dependencies
        _dependencies = ['modeller']
        self._check_dependencies(_dependencies)
        # execute evaluation
        self._evaluate(quiet)

    def evaluate(self) -> None:
        '''
        Run SOAP interaction evaluation. Automatically called on object
        initialization

        Returns
        -------
        None
        '''
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        mdl = modeller.scripts.complete_pdb(env, self.model.model_file)
        sele = modeller.selection(mdl.chains[:])
        # set up output
        self.output['soap_pp_all'] = sele.assess(modeller.soap_pp.Assessor())
        self.output['soap_pp_atom'] = sele.assess(
            modeller.soap_pp.AtomScorer())
        self.output['soap_pp_pair'] = sele.assess(
            modeller.soap_pp.PairScorer())


class Evaluation_qmean4(Evaluation):
    '''
    Class for evaluating a model with the QMEAN4 potential.

    Will dump the following entries to the model.evaluation dictionary:

    * qmean4
    * qmean4_z_score

    Parameters
    ----------
    model : Model
        The model object to evaluate.
    quiet : bool
        If True, will perform evaluation with suppressing stdout (default
        False). Needs to be False for running it asynchronously, as done when
        running Task.evaluate_models with multple cores

    Attributes
    ----------
    model : Model
        The model object to evaluate
    output : dict
        Dictionary that all outputs will be dumped into

    Raises
    ------
    ImportError
        Unable to import dependencies

    See Also
    --------
    Evaluation_qmean6
    Evaluation_qmeandisco

    Notes
    -----
    QMEAN is a statistical potential for evaluating homology models [4]_ [5]_.

    Briefly, QMEAN is a combination of different components. Four compoenents
    (interaction, cbeta, packing and torsion) form the `qmean4` score.

    For more information, please check the QMEAN documentation or the
    associated publications.

    References
    ----------
    .. [4] Benkert, P., Tosatto, S. C. E., & Schomburg, D. (2008). QMEAN: A
       comprehensive scoring function for model quality assessment. Proteins:
       Structure, Function and Genetics, 71(1), 261–277.
       https://doi.org/10.1002/prot.21715

    .. [5] Benkert, P., Biasini, M., & Schwede, T. (2011). Toward the
       estimation of the absolute quality of individual protein structure
       models. Bioinformatics, 27(3), 343–350.
       https://doi.org/10.1093/bioinformatics/btq662
    '''
    def __init__(self, model: typing.Type['Model'], quiet: bool =
                 False) -> None:
        # init attributes
        Evaluation.__init__(self, model)
        # check dependencies
        _dependencies = ['ost', 'qmean']
        self._check_dependencies(_dependencies)
        # execute evaluation
        self._evaluate(quiet)

    def evaluate(self) -> None:
        '''
        Run QMEAN4 protein evaluation. Automatically called on object
        initialization
        Returns
        -------
        None
        '''
        # gather inputs
        model = ost.io.LoadPDB(self.model.model_file)

        # generate score
        qmean_scorer = qmean.QMEANScorer(model)

        # extract scores to model.evaluation dict
        self.output['qmean4'] = qmean_scorer.qmean4_score
        self.output['qmean4_z_score'] = qmean_scorer.qmean4_z_score


class Evaluation_qmean6(Evaluation_qmean4):
    '''
    Class for evaluating a model with the QMEAN6 potential.

    Will dump the following entries to the model.evaluation dictionary:

    * qmean6
    * qmean6_disco

    Requires the following valid entries in the model.info dictionary:

    * accpro_file (.acc file)
    * psipred_file (.horiz file)

    Parameters
    ----------
    model : Model
        The model object to evaluate.
    quiet : bool
        If True, will perform evaluation with suppressing stdout (default
        False). Needs to be False for running it asynchronously, as done when
        running Task.evaluate_models with multple cores

    Attributes
    ----------
    model : Model
        The model object to evaluate
    output : dict
        Dictionary that all outputs will be dumped into

    Raises
    ------
    ImportError
        Unable to import dependencies

    See Also
    --------
    Evaluation_qmean4
    Evaluation_qmeandisco

    Notes
    -----
    QMEAN is a statistical potential for evaluating homology models [6]_ [7]_.

    QMEAN6 is a combination of six different components (interaction, cbeta,
    packing, torsion, ss_agreement, acc_agreement). It is an extension to the
    QMEAN4 score, which additionally evaluates the agreement of the model to
    secondary structur predictions from PSIPRED [8]_ and solvent accessiblity
    predictions from ACCpro [9]_.

    For more information, please check the QMEAN documentation or the
    associated publications.

    References
    ----------
    .. [6] Benkert, P., Tosatto, S. C. E., & Schomburg, D. (2008). QMEAN: A
       comprehensive scoring function for model quality assessment. Proteins:
       Structure, Function and Genetics, 71(1), 261–277.
       https://doi.org/10.1002/prot.21715

    .. [7] Benkert, P., Biasini, M., & Schwede, T. (2011). Toward the
       estimation of the absolute quality of individual protein structure
       models. Bioinformatics, 27(3), 343–350.
       https://doi.org/10.1093/bioinformatics/btq662

    .. [8] Jones, D. T. (1999). Protein secondary structure prediction based on
       position-specific scoring matrices. Journal of Molecular Biology,
       292(2), 195–202. https://doi.org/10.1006/JMBI.1999.3091

    .. [9] Magnan, C. N., & Baldi, P. (2014). SSpro/ACCpro 5: almost perfect
       prediction of protein secondary structure and relative solvent
       accessibility using profiles, machine learning and structural
       similarity. Bioinformatics, 30(18), 2592–2597.
       https://doi.org/10.1093/BIOINFORMATICS/BTU352
    '''
    def evaluate(self) -> None:
        '''
        Run QMEAN6 protein evaluation. Automatically called on object
        initialization

        Returns
        -------
        None
        '''
        # gather inputs
        model = ost.io.LoadPDB(self.model.model_file)
        accpro_handler = self._import_accpro()
        psipred_handler = self._import_psipred()

        # generate score
        qmean_scorer = qmean.QMEANScorer(
                model, accpro=accpro_handler, psipred=psipred_handler)

        # extract scores to model.evaluation dict
        self.output['qmean6'] = qmean_scorer.qmean6_score
        self.output['qmean6_z_score'] = qmean_scorer.qmean6_z_score

    def _import_accpro(self) -> typing.Type['qmean.ACCPROHandler']:
        '''
        Helper function to extract ACCpro annotation from file

        Reads sequence from ``model.file``. Reads ACCpro from
        ``model.info['accpro_file']``.

        If multiple sequences are in the ACCpro file, only the first one will
        be tried out.

        Returns
        -------
        qmean.ACCPROHandler

        Notes
        -----
        Tested with ACCpro version 5.2.
        '''
        try:
            accpro_file = self.model.info['accpro_file']
        except KeyError:
            print("Model.info['accpro_file'] has to be set.")
            raise

        # parse accpro file from accpro_file
        with open(accpro_file, 'r') as f:
            lines = f.readlines()

        # set up data dict
        data = {
                'seq': '',
                'acc': '',
                }

        # extract acc from file
        # get lines that start with >
        new_seqs = [i for i, line in enumerate(lines) if line.startswith('>')]
        # extract everything for the first seq
        if len(new_seqs) == 0:
            raise ValueError('Could not parse file, please check if correct'
                             'file type')
        elif len(new_seqs) == 1:
            accpro_lines = lines[new_seqs[0]+1:]
        else:
            accpro_lines = lines[new_seqs[0]+1:new_seqs[1]]
        # concat and replace '-' with 'b'
        accpro = ''.join([line.strip() for line in accpro_lines])
        accpro = accpro.replace('-', 'b')
        data['acc'] = accpro

        # extract seq from PDB
        data['seq'] = self.model.get_sequence()

        # check if data is not empty
        if (len(data['seq']) == 0 or len(data['acc']) == 0):
            raise ValueError('ACCpro annotation with length 0 not accepted')

        # return accpro handler
        return qmean.ACCPROHandler(data)

    def _import_psipred(self) -> typing.Type['qmean.PSIPREDHandler']:
        '''
        Helper function to extract PSIPRED annotation from file

        Reads PSIPRED predictions, confidence scores and the sequence from the
        ``model.info['psipred_file']``.

        Returns
        -------
        qmean.PSIPREDHandler

        Notes
        -----
        Tested with PSIPRED version 4.01.
        '''
        try:
            psipred_file = self.model.info['psipred_file']
        except KeyError:
            print("Model.info['psipred_file'] has to be set.")
            raise

        with open(psipred_file, 'r') as f:
            lines = f.readlines()

        # extract seq, ss and conf to a dict
        data = {
                'seq': '',
                'ss': '',
                'conf': '',
                }
        for line in lines:
            if line.startswith('  AA'):
                data['seq'] = data['seq'] + line.split(':')[-1].split()[0]
            elif line.startswith('Pred'):
                data['ss'] = data['ss'] + line.split(':')[-1].split()[0]
            elif line.startswith('Conf'):
                data['conf'] = data['conf'] + line.split(':')[-1].split()[0]

        # compare sequence to seq extracted from PDB file
        pdb_seq = self.model.get_sequence()
        if data['seq'] != pdb_seq:
            raise ValueError(
                'Sequence in PDB is not identical to sequence from PSIPRED '
                'file.')

        # check if data is not empty
        if (len(data['seq']) == 0 or len(data['ss']) == 0 or len(data['conf'])
                == 0):
            raise ValueError('PSIPRED annotation with length 0 not accepted')
        # return psipred handler
        return qmean.PSIPREDHandler(data)


class Evaluation_qmeandisco(Evaluation_qmean6):
    '''
    Class for evaluating a model with the QMEAN DisCo potential.

    Will dump the following entries to the model.evaluation dictionary:

    * qmean6
    * qmean6_z_score
    * qmean_local_scores_avg
    * qmean_local_scores_err

    Requires the following valid entries in the model.info dictionary:

    * accpro_file (.acc file)
    * psipred_file (.horiz file)
    * disco_file (generated by ``qmean.DisCoContainer.Save``)

    Parameters
    ----------
    model : Model
        The model object to evaluate.
    quiet : bool
        If True, will perform evaluation with suppressing stdout (default
        False). Needs to be False for running it asynchronously, as done when
        running Task.evaluate_models with multple cores

    Attributes
    ----------
    model : Model
        The model object to evaluate
    output : dict
        Dictionary that all outputs will be dumped into

    Raises
    ------
    ImportError
        Unable to import dependencies

    See Also
    --------
    Evaluation_qmean4
    Evaluation_qmean6

    Notes
    -----
    QMEAN is a statistical potential for evaluating homology models [10]_
    [11]_.

    QMEAN DisCo is an extension of QMEAN by the inclusion of homology derived
    DIStance COnstraints [12]_. These distance contraints do not influence the
    six component of the QMEAN6 score (interaction, cbeta, packing, torsion,
    ss_agreement, acc_agreement), but only the local scores.

    The distance contraints for the target have to be generated before and
    saved to a file.

    For more information, please check the QMEAN documentation or the
    associated publications.

    References
    ----------
    .. [10] Benkert, P., Tosatto, S. C. E., & Schomburg, D. (2008). QMEAN: A
       comprehensive scoring function for model quality assessment. Proteins:
       Structure, Function and Genetics, 71(1), 261–277.
       https://doi.org/10.1002/prot.21715

    .. [11] Benkert, P., Biasini, M., & Schwede, T. (2011). Toward the
       estimation of the absolute quality of individual protein structure
       models. Bioinformatics, 27(3), 343–350.
       https://doi.org/10.1093/bioinformatics/btq662

    .. [12] Studer, G., Rempfer, C., Waterhouse, A. M., Gumienny, R., Haas, J.,
       & Schwede, T. (2020). QMEANDisCo-distance constraints applied on model
       quality estimation. Bioinformatics, 36(6), 1765–1771.
       https://doi.org/10.1093/bioinformatics/btz828

    '''
    def evaluate(self) -> None:
        '''
        Run QMEAN DisCo protein evaluation. Automatically called on object
        initialization

        Returns
        -------
        None
        '''
        # gather inputs
        model = ost.io.LoadPDB(self.model.model_file)
        accpro_handler = self._import_accpro()
        psipred_handler = self._import_psipred()
        dc = self._import_disco()

        # generate score
        qmean_scorer = qmean.QMEANScorer(
                model, accpro=accpro_handler, psipred=psipred_handler,
                dc=dc)

        # extract scores to model.evaluation dict
        self.output['qmean6'] = qmean_scorer.qmean6_score
        self.output['qmean6_z_score'] = qmean_scorer.qmean6_z_score
        self.output['qmean_local_scores_avg'] = (
            qmean_scorer.avg_local_score)
        self.output['qmean_local_scores_err'] = (
            qmean_scorer.avg_local_score_error)

    def _import_disco(self) -> typing.Type['qmean.DisCoContainer']:
        '''
        Helper function to extract distance constraints from Model.

        Requires Model.info['disco_file'] to be set.

        Returns
        -------
        qmean.DisCoContainer
        '''
        try:
            dc = qmean.DisCoContainer.Load(self.model.info['disco_file'])
        except KeyError:
            print("Model.info['disco_file'] has to be set.")
            raise

        # check sequence in disco container against sequence in model PDB
        seq_disco = dc.GetSeqres().string.upper()
        seq_pdb = self.model.get_sequence()

        if seq_disco != seq_pdb:
            raise ValueError(
                'Sequence in PDB is not identical to sequence from DisCo '
                'container.')

        return dc


class Evaluation_mol_probity(Evaluation):
    '''
    Class for evaluating a model with the MolProbity validation service.

    Will dump the following entries to the model.evaluation dictionary:

    * mp_score

    Parameters
    ----------
    model : Model
        The model object to evaluate
    quiet : bool
        If True, will perform evaluation with suppressing stdout (default
        False). Needs to be False for running it asynchronously, as done when
        running Task.evaluate_models with multple cores

    Attributes
    ----------
    model : Model
        The model object to evaluate
    output : dict
        Dictionary that all outputs will be dumped into

    Notes
    -----
    Molprobity is a program that evaluates the quality of 3D structures of
    proteins based on structural features [13]_ [14]_ [15]_. For more
    information, please check the MolProbity webpage or the associated
    publications.

    References
    ----------

    .. [13] Davis, I. W., Leaver-Fay, A., Chen, V. B., Block, J. N., Kapral, G.
       J., Wang, X., Murray, L. W., Arendall, W. B., Snoeyink, J., Richardson,
       J. S., & Richardson, D. C. (2007). MolProbity: all-atom contacts and
       structure validation for proteins and nucleic acids. Nucleic Acids
       Research, 35(suppl_2), W375–W383. https://doi.org/10.1093/NAR/GKM216

    .. [14] Chen, V. B., Arendall, W. B., Headd, J. J., Keedy, D. A.,
       Immormino, R. M., Kapral, G. J., Murray, L. W., Richardson, J. S., &
       Richardson, D. C. (2010). MolProbity: All-atom structure validation for
       macromolecular crystallography. Acta Crystallographica Section D:
       Biological Crystallography, 66(1), 12–21.
       https://doi.org/10.1107/S0907444909042073

    .. [15] Williams, C. J., Headd, J. J., Moriarty, N. W., Prisant, M. G.,
       Videau, L. L., Deis, L. N., Verma, V., Keedy, D. A., Hintze, B. J.,
       Chen, V. B., Jain, S., Lewis, S. M., Arendall, W. B., Snoeyink, J.,
       Adams, P. D., Lovell, S. C., Richardson, J. S., & Richardson, D. C.
       (2018). MolProbity: More and better reference data for improved all-atom
       structure validation. Protein Science, 27(1), 293–315.
       https://doi.org/10.1002/pro.3330
    '''
    def __init__(self, model: typing.Type['Model'], quiet: bool =
                 False) -> None:
        # initialize attributes
        Evaluation.__init__(self, model)
        # execute evaluation
        self._evaluate(quiet)

    def evaluate(self) -> None:
        '''
        Run MolProbity evaluation. Automatically called on object
        initialization

        Returns
        -------
        None
        '''
        # run with subprocess
        command = ['phenix.molprobity', '--coot=False', '--probe_dots=False',
                   '--quiet', self.model.model_file]
        p = subprocess.run(command, stdout=subprocess.PIPE, check=True,
                           universal_newlines=True, shell=False)
        # parse output with regular expressions
        re_mp_score = r'MolProbity score\s*=\s*(\d*[.]*\d*)\n'
        self.output['mp_score'] = float(
            re.search(re_mp_score, p.stdout).group(1))

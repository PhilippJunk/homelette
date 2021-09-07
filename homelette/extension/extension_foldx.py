'''
FoldX extension to homelette
============================

Philipp Junk, 2021

This extension contains evaluation metrics based on FoldX, a force field for
energy calculation and protein design (https://foldxsuite.crg.eu/) [1]_ [2]_.

.. [1] Guerois, R., Nielsen, J. E., & Serrano, L. (2002). Predicting Changes in
   the Stability of Proteins and Protein Complexes: A Study of More Than 1000
   Mutations. Journal of Molecular Biology, 320(2), 369–387.
   https://doi.org/10.1016/S0022-2836(02)00442-4

.. [2] Schymkowitz, J., Borg, J., Stricher, F., Nys, R., Rousseau, F., &
   Serrano, L. (2005). The FoldX web server: an online force field. Nucleic
   Acids Research, 33(Web Server), W382–W388.
   https://doi.org/10.1093/nar/gki387

Usage
-----

.. code-block:: python

   import homelette.extension.extension_foldx as extension_foldx
   help(extension_foldx.Evaluation_foldx_stability)

This extension expects FoldX to be installed and in your path.

Functions and classes
---------------------

Currently contains the following items:
    :class:`Evaluation_foldx_repairmodels`
    :class:`Evaluation_foldx_interaction`
    :class:`Evaluation_foldx_stability`
    :class:`Evaluation_fodlx_alascan_buildmodels`
    :class:`Evaluation_foldx_alascan_interaction`

-----

'''

__all__ = ['Evaluation_foldx_repairmodels', 'Evaluation_foldx_interaction',
           'Evaluation_foldx_stability',
           'Evaluation_foldx_alascan_buildmodels',
           'Evaluation_foldx_alascan_interaction']
__author__ = 'Philipp Junk'

import glob
import os
import os.path
import time
import subprocess

from ..evaluation import Evaluation


class _Evaluation_foldx():
    '''
    Superclass for Foldx-based evaluations. Implements a check for rotabase.
    '''
    def _check_rotabase(self):
        '''
        Checks of rotabase.txt is already present, stops program for 2 seconds
        if only has been created recently (useful for avoiding race conditions
        with FoldX).
        '''
        # check rotabase.txt
        if os.path.exists('rotabase.txt'):
            time_since_modification = (
                    time.time() - os.path.getmtime('rotabase.txt'))
            if time_since_modification < 1:
                time.sleep(2)


class Evaluation_foldx_repairmodels(Evaluation, _Evaluation_foldx):
    '''
    Creates a modified version of the PDB and runs RepairPDB on it

    Will not dump an entry to the model.evaluation dictionary

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

    Notes
    -----
    Most PDBs work fine with FoldX. For a specific use case in which I was
    working with GTP heteroatoms, I had to rename a few atoms to make the PDB
    compliant with FoldX.
    '''
    def __init__(self, model, quiet=False):
        Evaluation.__init__(self, model)
        self.evaluate()

    def evaluate(self):
        '''
        Repairs models with FoldX. Automatically called on object
        initialization

        Returns
        -------
        None
        '''
        # prepare files and directories for foldx input
        model_file = os.path.basename(self.model.model_file)
        model_dir = os.path.realpath(os.path.dirname(self.model.model_file))
        model_name = os.path.splitext(model_file)[0]

        modified_model_file = model_file
        modified_model_dir = os.path.realpath(os.path.join(
            model_dir, 'foldx_modified'))
        if not os.path.exists(modified_model_dir):
            os.makedirs(modified_model_dir, exist_ok=True)

        repair_file = model_name + '_Repair.pdb'
        repair_dir = os.path.realpath(os.path.join(
            model_dir, 'foldx_repair'))
        if not os.path.exists(repair_dir):
            os.makedirs(repair_dir, exist_ok=True)

        # check for rotabase.txt
        _Evaluation_foldx._check_rotabase(self)

        # modify PDB
        # FoldX expects different atom names for some atoms in GTP
        # in order to properly calculate energies, this has be converted
        # ie. "O5'" in PDB has to be modified to "O5*"
        # read PDB
        with open(os.path.join(model_dir, model_file), 'r') as f:
            lines = f.readlines()
        # write modified PDB
        with open(os.path.join(modified_model_dir,
                               modified_model_file), 'w') as f:
            for line in lines:
                if (
                        (line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM') and
                        (line[17:20] == 'GTP' or line[17:20] == 'GDP') and
                        line[15] == "'"):
                    f.write(line[:15] + '*' + line[16:])
                else:
                    f.write(line)

        # RepairPDB
        if not os.path.exists(os.path.join(repair_dir, repair_file)):
            command = [
                'foldx', '--command', 'RepairPDB', '--pdb',
                modified_model_file, '--pdb-dir', modified_model_dir,
                '--output-dir', repair_dir]
            subprocess.run(command, check=True, shell=False,
                           stdout=subprocess.DEVNULL)


class Evaluation_foldx_stability(Evaluation, _Evaluation_foldx):
    '''
    Calculate protein stability with FoldX

    Expects Evaluation_foldx_repairmodels to have been performed beforehand.

    Will dump the following entries to the model.evaluation dictionary:

    * foldx_stability

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
    '''
    def __init__(self, model, quiet=False):
        Evaluation.__init__(self, model)
        self._evaluate(quiet)

    def evaluate(self):
        '''
        Calculates protein stability with FoldX. Automatically called on object
        initialization.

        Returns
        -------
        None
        '''
        model_file = os.path.basename(self.model.model_file)
        model_dir = os.path.realpath(os.path.dirname(self.model.model_file))
        repair_file = os.path.splitext(model_file)[0] + '_Repair.pdb'
        repair_dir = os.path.realpath(os.path.join(model_dir, 'foldx_repair'))
        stability_file = os.path.splitext(repair_file)[0] + '_0_ST.fxout'
        stability_dir = os.path.realpath(os.path.join(
            model_dir, 'foldx_stability'))
        if not os.path.isdir(repair_dir):
            os.makedirs(repair_dir, exist_ok=True)
        if not os.path.isdir(stability_dir):
            os.makedirs(stability_dir, exist_ok=True)

        # Check for rotabase.txt
        _Evaluation_foldx._check_rotabase(self)

        # RepairPDB
        if not os.path.exists(os.path.join(repair_dir, repair_file)):
            Evaluation_foldx_repairmodels(self.model)

        # Stability
        if not os.path.exists(os.path.join(stability_dir, stability_file)):
            command = [
                'foldx', '--command', 'Stability', '--pdb',
                repair_file, '--pdb-dir', repair_dir, '--output-dir',
                stability_dir]
            subprocess.run(command, check=True, shell=False,
                           stdout=subprocess.DEVNULL)

        # collect output
        with open(os.path.join(stability_dir, stability_file), 'r') as f:
            line = f.read()
            self.output['foldx_stability'] = float(
                line.split('\t')[1])


class Evaluation_foldx_interaction(Evaluation, _Evaluation_foldx):
    '''
    Calculates interaction energy with FoldX

    Requires a protein-protein complex. Expects Evaluation_foldx_repairmodels
    to have been performed beforehand.

    Will dump the following entries to the model.evaluation dictionary:

    * foldx_interaction

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
    '''
    def __init__(self, model, quiet=False):
        Evaluation.__init__(self, model)
        self._evaluate(quiet)

    def evaluate(self):
        '''
        Calculates protein interaction energy with FoldX. Automatically called
        on object initialization.

        Returns
        -------
        None
        '''
        # prepare files and directories for foldx input
        model_file = os.path.basename(self.model.model_file)
        model_dir = os.path.realpath(os.path.dirname(self.model.model_file))

        repair_file = os.path.splitext(model_file)[0] + '_Repair.pdb'
        repair_dir = os.path.realpath(os.path.join(model_dir, 'foldx_repair'))

        interaction_file = ('Interaction_' + os.path.splitext(repair_file)[0] +
                            '_AC.fxout')
        interaction_dir = os.path.realpath(os.path.join(
            model_dir, 'foldx_interaction'))
        if not os.path.isdir(interaction_dir):
            os.makedirs(interaction_dir, exist_ok=True)

        # check for rotabase
        _Evaluation_foldx._check_rotabase(self)

        # perform repairmodels if repaired file does not already exists
        if not os.path.exists(os.path.join(repair_dir, repair_file)):
            Evaluation_foldx_repairmodels(self.model)

        # Interaction
        if not os.path.exists(os.path.join(interaction_dir, interaction_file)):
            command = [
                'foldx', '--command', 'AnalyseComplex', '--pdb',
                repair_file, '--pdb-dir', repair_dir, '--output-dir',
                interaction_dir]
            subprocess.run(command, check=True, shell=False,
                           stdout=subprocess.DEVNULL)

        # collect output
        with open(os.path.join(interaction_dir, interaction_file), 'r') as f:
            lines = f.readlines()
            self.output['foldx_interaction'] = float(
                lines[9].split('\t')[5])


class Evaluation_foldx_alascan_buildmodels(Evaluation, _Evaluation_foldx):
    '''
    Generates alanine point mutations for all positions in the given model
    using FoldX. Automatically called on object initialization.

    Expects Evaluation_foldx_repairmodels to have been performed beforehand.

    Will not dump an entry to the model.evaluation dictionary.

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

    See Also
    --------
    Evaluation_foldx_alascan_interaction

    Notes
    -----
    This Evaluation is very RAM intensive, so expect only to run 1 or 2 threads
    ni parallel.
    '''
    def __init__(self, model, quiet=False):
        Evaluation.__init__(self, model)
        self.evaluate()

    def evaluate(self):
        '''
        Generates alanine point mutations  for all positions in the given
        model. Automatically called on object initialization.

        Returns
        -------
        None
        '''
        # prepare files and directories for foldx input
        model_file = os.path.basename(self.model.model_file)
        model_dir = os.path.realpath(os.path.dirname(self.model.model_file))
        model_name = os.path.splitext(model_file)[0]

        repair_file = model_name + '_Repair.pdb'
        repair_dir = os.path.realpath(os.path.join(model_dir, 'foldx_repair'))
        if not os.path.isdir(repair_dir):
            os.makedirs(repair_dir, exist_ok=True)

        mutants_dir = os.path.realpath(os.path.join(
            model_dir, 'foldx_alascan/models'))
        mutant_file = os.path.realpath(os.path.join(
            model_dir, 'foldx_alascan', 'individual_list_' + model_name +
            '.txt'))
        if not os.path.isdir(mutants_dir):
            os.makedirs(mutants_dir, exist_ok=True)

        # check for rotabase
        _Evaluation_foldx._check_rotabase(self)

        # generate input for mutagenesis
        pdb = self.model.parse_pdb()
        # collapse dataframe
        residues = (
            pdb[pdb.record.eq('ATOM')][['chainID', 'resSeq', 'resName']]
            .groupby(['chainID', 'resSeq', 'resName'])
            .count()
            .reset_index()
        )

        # change 3 letter code to 1 letter code
        def _321(aa):
            '''
            Transform 3 letter amino acid code to 1 letter code
            '''
            aa_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                       'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                       'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                       'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                       'SEC': 'U'}
            return aa.map(aa_code)
        residues = (
            residues.assign(resName=_321(residues['resName']))
        )
        # concat
        mutlist = list(residues.resName + residues.chainID +
                       residues.resSeq.astype('str') + 'A')
        with open(mutant_file, 'w') as f:
            for mut in mutlist:
                f.write(mut + ';\n')

        # BuildModel
        command = [
            'foldx', '--command', 'BuildModel', '--pdb', repair_file,
            '--pdb-dir', repair_dir, '--output-dir', mutants_dir,
            '--mutant-file', mutant_file, '--moveNeighbours', 'false'
        ]
        subprocess.run(command, check=True, shell=False,
                       stdout=subprocess.DEVNULL)


class Evaluation_foldx_alascan_interaction(Evaluation, _Evaluation_foldx):
    '''
    Calculates protein interaction energy with FoldX for all alanine point
    mutations generated by Evaluation_foldx_alascan_buildmodels.

    Expects Evaluation_foldx_alascan_buildmodels to have been run before.

    Will dump the following entry to the model.evaluation dictionary:

    * foldx_alascan: Dictionary of all interaction energies for all alanine
                     scan mutations.

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

    See Also
    --------
    Evaluation_foldx_alascan_buildmodels
    '''
    def __init__(self, model, quiet=False):
        Evaluation.__init__(self, model)
        self._evaluate(quiet)

    def evaluate(self):
        '''
        Calculates protein interaction energy with FoldX for all alanine point
        mutations generated by Evaluation_foldx_alascan_buildmodels.

        Returns
        -------
        None
        '''
        # prepare files and directories for foldx input
        model_file = os.path.basename(self.model.model_file)
        model_dir = os.path.realpath(os.path.dirname(self.model.model_file))
        model_name = os.path.splitext(model_file)[0]

        mutants_dir = os.path.realpath(os.path.join(
            model_dir, 'foldx_alascan/models'))
        mutant_file = os.path.realpath(os.path.join(
            model_dir, 'foldx_alascan', 'individual_list_' + model_name +
            '.txt'))

        interaction_dir = os.path.realpath(os.path.join(
            model_dir, 'foldx_alascan/interaction'))
        if not os.path.isdir(interaction_dir):
            os.makedirs(interaction_dir, exist_ok=True)

        # check for rotabase
        _Evaluation_foldx._check_rotabase(self)

        # recreate mutlist
        pdb = self.model.parse_pdb()
        # collapse dataframe
        residues = (
            pdb[pdb.record.eq('ATOM')][['chainID', 'resSeq', 'resName']]
            .groupby(['chainID', 'resSeq', 'resName'])
            .count()
            .reset_index()
        )

        # change 3 letter code to 1 letter code
        def _321(aa):
            '''
            Transform 3 letter amino acid code to 1 letter code
            '''
            aa_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                       'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                       'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                       'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                       'SEC': 'U'}
            return aa.map(aa_code)
        residues = (
            residues.assign(resName=_321(residues['resName']))
        )
        # concat
        mutlist = list(residues.resName + residues.chainID +
                       residues.resSeq.astype('str') + 'A')

        # Interaction
        for i in range(len(mutlist)):
            model_mut_file = model_name + '_Repair_' + str(i + 1) + '.pdb'
            command = [
                'foldx', '--command', 'AnalyseComplex', '--pdb',
                model_mut_file, '--pdb-dir', mutants_dir, '--output-dir',
                interaction_dir
            ]
            subprocess.run(command, check=True, shell=False,
                           stdout=subprocess.DEVNULL)

        # collect output
        output = dict()
        for i in range(len(mutlist)):
            interaction_file = ('Interaction_' + model_name + '_Repair_' +
                                str(i + 1) + '_AC.fxout')
            with open(
                    os.path.join(interaction_dir, interaction_file), 'r') as f:
                lines = f.readlines()
                output[mutlist[i]] = float(
                    lines[9].split('\t')[5])
        self.output['foldx_alascan'] = output

        # delete model files
        for f in glob.glob(
                os.path.join(mutants_dir, model_name + '_Repair_*.pdb')):
            os.remove(f)
        for f in glob.glob(
                os.path.join(mutants_dir, '*_' + model_name +
                             '_Repair.fxout')):
            os.remove(f)
        # delete interaction files
        for f in glob.glob(
                os.path.join(interaction_dir, '*_' + model_name +
                             '_Repair_*.fxout')):
            os.remove(f)
        # delete mutfile
        if os.path.exists(mutant_file):
            os.remove(mutant_file)

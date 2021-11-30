'''
``homelette.alignment``
=======================

The :mod:`homelette.alignment` submodule contains a selection of tools for
handling sequences and alignments, as well as for the automatic generation of
sequences from a target sequence.

Tutorials
---------

Basic handing of alignments with `homelette` is demonstrated in :ref:`Tutorial
1 </Tutorial1_Basics.ipynb>`. The assembling of alignments for complex
modelling is discussed in
:ref:`Tutorial 6 </Tutorial6_ComplexModelling.ipynb>`. The automatic generation
of alignments is shown in :ref:`Tutorial 8
</Tutorial8_AlignmentGeneration.ipynb>`.

Functions and classes
---------------------

Functions and classes present in `homelette.alignment` are listed below:

    :class:`Alignment`
    :class:`Sequence`
    :class:`AlignmentGenerator`
    :class:`AlignmentGenerator_pdb`
    :class:`AlignmentGenerator_hhblits`
    :class:`AlignmentGenerator_from_aln`
    :func:`assemble_complex_aln`

-----

'''

__all__ = ['Alignment', 'Sequence', 'AlignmentGenerator',
           'AlignmentGenerator_pdb', 'AlignmentGenerator_hhblits',
           'AlignmentGenerator_from_aln', 'assemble_complex_aln']

# Standard library imports
import abc
import contextlib
import itertools
import json
import os.path
import re
import shutil
import subprocess
import typing
import urllib.error
import urllib.parse
import urllib.request
import warnings

#  Third party imports
import pandas as pd

# Local application imports
from . import pdb_io
from .organization import Task


class Sequence():
    '''
    Class that contains individual sequences and miscellaneous information
    about them.

    Parameters
    ----------
    name : str
        Identifier of the sequence
    sequence : str
        Sequence in 1 letter amino acid code
    **kwargs :
        Annotations, for acceptable keys see :func:`Sequence.annotate`

    Attributes
    ----------
    name : str
        Identifier of the sequence
    sequence : str
        Sequence in 1 letter amino acid code
    annotation : dict
        Collection of annotation for this sequence

    Notes
    -----
    See :func:`Sequence.annotate` for more information on the annotation of
    sequences.
    '''
    def __init__(self, name: str, sequence: str, **kwargs) -> None:
        self.sequence = sequence
        self.name = name
        self.annotation = {
            'seq_type': str(),
            'pdb_code': str(),
            'begin_res': str(),
            'begin_chain': str(),
            'end_res': str(),
            'end_chain': str(),
            'prot_name': str(),
            'prot_source': str(),
            'resolution': str(),
            'r_factor': str()
        }
        self.annotate(**kwargs)

    def annotate(self, **kwargs: str):
        '''
        Change annotation for sequence object.

        Keywords not specified in the Notes section will be ignored.

        Parameters
        ----------
        kwargs : str
            Annotations. For acceptible values, see Notes.

        Returns
        -------
        None

        Notes
        -----

        Annotations are important for MODELLER in order to properly process
        alignment in PIR format. The following annotations are supported and
        can be modified.

        +---------------+---------------------------------------------------+
        | annotation    | explanation                                       |
        +===============+===================================================+
        | seq_type      | Specification whether sequence should be treated  |
        |               | as a template (set to 'structure') or as a target |
        |               | (set to 'sequence')                               |
        +---------------+---------------------------------------------------+
        | pdb_code      | PDB code corresponding to sequence (if available) |
        +---------------+---------------------------------------------------+
        | begin_res     | Residue number for the first residue of the       |
        |               | sequence in the corresponing PDB file             |
        +---------------+---------------------------------------------------+
        | begin_chain   | Chain identifier for the first residue of the     |
        |               | sequence in the corresponding PDB file            |
        +---------------+---------------------------------------------------+
        | end_res       | Residue number for the last residue of the        |
        |               | sequence in the corresponding PDB file            |
        +---------------+---------------------------------------------------+
        | end_chain     | Chain identifier for the last residue of the      |
        |               | sequence in the corresponding PDB file            |
        +---------------+---------------------------------------------------+
        | prot_name     | Protein name, optional                            |
        +---------------+---------------------------------------------------+
        | prot_source   | Protein source, optional                          |
        +---------------+---------------------------------------------------+
        | resolution    | Resolution of PDB structure, optional             |
        +---------------+---------------------------------------------------+
        | R_factor      | R-factor of PDB structure, optional               |
        +---------------+---------------------------------------------------+

        Different types of annotations are required, depending whether a target
        or a template is annotated. For *targets*, it is sufficient to seq the
        seq_type to 'sequence'. For *templates*, it is required by MODELLER
        that seq_type and pdb_code are annotated. begin_res, begin_chain,
        end_res and end_chain are recommended. The rest can be left unannoted.

        Examples
        --------

        Annotation for a target sequence.

        >>> target = hm.alignment.Sequence(name = 'target', sequence =
        ...     'TARGET')
        >>> target.annotation
        {'seq_type': '', 'pdb_code': '', 'begin_res': '', 'begin_chain': '',
        'end_res': '', 'end_chain': '', 'prot_name': '', 'prot_source': '',
        'resolution': '', 'r_factor': ''}
        >>> target.annotate(seq_type = 'sequence')
        >>> target.annotation
        {'seq_type': 'sequence', 'pdb_code': '', 'begin_res': '',
        'begin_chain': '', 'end_res': '', 'end_chain': '', 'prot_name': '',
        'prot_source': '', 'resolution': '', 'r_factor': ''}

        Annotation for a template structure.

        >>> template = hm.alignment.Sequence(name = 'template', sequence =
        ...     'TEMPLATE')
        >>> template.annotation
        {'seq_type': '', 'pdb_code': '', 'begin_res': '', 'begin_chain': '',
        'end_res': '', 'end_chain': '', 'prot_name': '', 'prot_source': '',
        'resolution': '', 'r_factor': ''}
        >>> template.annotate(seq_type = 'structure', pdb_code = 'TMPL',
        ...     begin_res = '1', begin_chain = 'A', end_res = '8', end_chain =
        ...     'A')
        >>> template.annotation
        {'seq_type': 'structure', 'pdb_code': 'TMPL', 'begin_res': '1',
        'begin_chain': 'A', 'end_res': '8', 'end_chain': 'A', 'prot_name': '',
        'prot_source': '', 'resolution': '', 'r_factor': ''}
        '''
        for key, value in kwargs.items():
            # only include keys that are already present
            if key in self.annotation.keys():
                self.annotation[key] = value

    def get_annotation_pir(self) -> str:
        '''
        Return annotation in the colon-separated format expected from the PIR
        alignment format used by MODELLER.

        Returns
        -------
        str
            Annotation in PIR format

        Examples
        --------
        >>> template = hm.alignment.Sequence(name = 'template', sequence =
        ...     'TEMPLATE', seq_type = 'structure', pdb_code = 'TMPL',
        ...     begin_res = '1', begin_chain = 'A', end_res = '8', end_chain =
        ...     'A')
        >>> template.get_annotation_pir()
        'structure:TMPL:1:A:8:A::::'
        '''
        output = ('{seq_type}:{pdb_code}:{begin_res}:{begin_chain}:{end_res}:'
                  '{end_chain}:{prot_name}:{prot_source}:{resolution}:'
                  '{r_factor}'.format_map(self.annotation)
                  )
        return output

    def get_annotation_print(self) -> None:
        '''
        Print annotation to console

        Returns
        -------
        None

        Examples
        --------
        >>> template = hm.alignment.Sequence(name = 'template', sequence =
        ...     'TEMPLATE', seq_type = 'structure', pdb_code = 'TMPL',
        ...     begin_res = '1', begin_chain = 'A', end_res = '8', end_chain =
        ...     'A')
        >>> template.get_annotation_print()
        Sequence Type	structure
        PDB ID		TMPL
        Start Residue	1
        Start Chain	A
        End Residue	8
        End Chain	A
        Protein Name
        Protein Source
        Resolution
        R-Factor
        '''
        output = ('Sequence Type\t{seq_type}\nPDB ID\t\t{pdb_code}\nStart '
                  'Residue\t{begin_res}\nStart Chain\t{begin_chain}\nEnd '
                  'Residue\t{end_res}\nEnd Chain\t{end_chain}\nProtein '
                  'Name\t{prot_name}\nProtein Source\t{prot_source}\n'
                  'Resolution\t{resolution}\nR-Factor\t{r_factor}'.format_map(
                      self.annotation)
                  )
        print(output)

    def get_gaps(self) -> tuple:
        '''
        Find gap positions in sequence

        Returns
        -------
        tuple
            Positions of gaps in sequence

        Examples
        --------
        >>> seq = hm.alignment.Sequence(name = 'seq', sequence = 'SEQ-UEN--CE')
        >>> seq.get_gaps()
        (3, 7, 8)
        '''
        i = 0
        gaps = []
        while i < len(self.sequence):
            i = self.sequence.find('-', i)
            if i == -1:
                break
            gaps.append(i)
            i += 1
        return tuple(gaps)

    def remove_gaps(self, remove_all: bool = False,
                    positions: typing.Optional[typing.Iterable[int]] =
                    None) -> None:
        '''
        Remove gaps from the sequence.

        Gaps in the alignment are symbolized by '-'. Removal can either happen
        at specific or all positions. Indexing for specific positions is
        zero-based and checked before removal (raises Warning if the attempted
        removal of a non-gap position is detected)

        Parameters
        ----------
        remove_all : bool
            Remove all gaps (default False)
        positions : iterable
            Positions to remove (zero-based indexing)

        Returns
        -------
        None

        Warns
        -----
        UserWarning
            Specified position is not a gap

        Examples
        --------

        Example 1: remove all

        >>> seq = hm.alignment.Sequence(name = 'seq', sequence = 'SEQ-UEN--CE')
        >>> seq.remove_gaps(remove_all = True)
        >>> seq.sequence
        'SEQUENCE'

        Example 2: selective removal

        >>> seq = hm.alignment.Sequence(name = 'seq', sequence = 'SEQ-UEN--CE')
        >>> seq.remove_gaps(positions = (7, 8))
        >>> seq.sequence
        'SEQ-UENCE'
        '''
        if remove_all is True:
            self.sequence = self.sequence.replace('-', '')
        else:
            sequence = list(self.sequence)
            new_sequence = list()
            for i, res in enumerate(sequence):
                if i in positions:
                    # check if position is actually a gap
                    if res != '-':
                        warnings.warn(('Position {} is not a gap and can not'
                                       ' be removed').format(i))
                        new_sequence.append(sequence[i])
                else:
                    new_sequence.append(res)
            self.sequence = ''.join(new_sequence)


class Alignment():
    '''
    Class for managing sequence alignments.

    Parameters
    ----------
    file_name : str, optional
        The file to read the alignment from. If no file name is given, an empty
        alignment object will be created (default None)
    file_format : str, optional
        The format of the alignment file. Can be 'fasta' or 'pir' (default
        'fasta')

    Attributes
    ----------
    sequences : dict
        Collection of sequences. Sequences names are the dictionary keys,
        Sequence objects the values

    Raises
    ------
    ValueError
        File_format specified is not 'fasta' or 'pir'
    '''
    def __init__(self, file_name: str = None, file_format: str =
                 'fasta') -> None:
        self.sequences = dict()
        # if no file is given, initialize empty
        if file_name is None:
            pass
        # else import file
        elif file_format == 'fasta':
            self._read_fasta_aln(file_name)
        elif file_format == 'pir':
            self._read_pir(file_name)
        else:
            raise ValueError(
                    ('Unknown format for the alignment file: {}\nPlease use'
                     '"fasta" or "pir"".').format(
                                  file_format))

    def _read_pir(self, file_name: str) -> None:
        '''
        Parser for PIR format alignment file.

        Parameters
        ----------
        file_name : str
            File name of PIR file to parse

        Returns
        -------
        None
        '''
        with open(file_name, 'r') as file_handler:
            lines = file_handler.readlines()

        sequences = dict()  # dict that will replace self.sequences
        seq, seq_name = str(), str()  # temporary storage for name and sequence
        annotations = dict()  # temporary storage for annotations

        def add_sequence(seq_name, seq, annotations):
            '''
            Helper function for adding sequences to the alignment after
            performing sanity checks

            Checks if sequence length is not 0, checks if last character is PIR
            alignment format stopping character "*" and removes it if present.
            Then transform to sequence and adds to sequences.
            '''
            if len(seq) != 0:  # add if not empty
                # check for last character, if * remove
                if seq[-1] == '*':
                    seq = seq[:-1]
                sequences[seq_name] = Sequence(
                    seq_name, seq, **annotations)

        # iterate over file
        for line in lines:
            if line.startswith('>P1;'):  # new sequence
                add_sequence(seq_name, seq, annotations)
                # reset seq, seq_name and annotations
                seq, seq_name, annotations = str(), str(), dict()
                seq_name = line.replace('>P1;', '').strip()
            elif len(seq) == 0 and len(annotations) == 0:  # annotations
                values = line.split(':') + ['']
                keys = ['seq_type', 'pdb_code', 'begin_res', 'begin_chain',
                        'end_res', 'end_chain', 'prot_name', 'prot_source',
                        'resolution', 'r_factor']
                for key, value in zip(keys, values):
                    annotations[key] = value
            else:  # sequence
                seq += line.strip().upper()
        if len(seq) != 0:  # last sequence
            add_sequence(seq_name, seq, annotations)

        # update self.sequences
        self.sequences = sequences

    def _read_fasta_aln(self, file_name: str) -> None:
        '''
        Parser for FASTA format alignment file.

        Parameters
        ----------
        file_name : str
            File name of FASTA file to parse

        Returns
        -------
        None
        '''
        with open(file_name, 'r') as file_handler:
            lines = file_handler.readlines()

        sequences = dict()  # dict that will replace self.sequences
        seq, seq_name = str(), str()  # temporary storage for name and sequence

        for line in lines:
            if line.startswith('>'):  # new sequence
                if len(seq) != 0:  # add if not empty
                    sequences[seq_name] = Sequence(seq_name, seq)
                    seq, seq_name = str(), str()  # reset seq and seq_name
                seq_name = line.replace('>', '').strip()
            else:
                seq += line.strip().upper()
        if len(seq) != 0:  # last sequence
            sequences[seq_name] = Sequence(seq_name, seq)

        self.sequences = sequences

    def get_sequence(self, sequence_name: str) -> typing.Type['Sequence']:
        '''
        Retrieve sequence object by sequence name.

        Parameters
        ----------
        sequence_name : str
            Name of sequence to retrieve

        Returns
        -------
        Sequence
        '''
        return self.sequences[sequence_name]

    def select_sequences(self, sequence_names: typing.Iterable) -> None:
        '''
        Select sequences to remain in the alignment by sequence name

        Parameters
        ----------
        sequence_names : iterable
            Iterable of sequence names

        Returns
        -------
        None

        Raises
        ------
        KeyError
            Sequence name not found in alignment
        '''
        sequences = dict()
        for seq_name in sequence_names:
            sequences[seq_name] = self.sequences[seq_name]
        self.sequences = sequences

    def remove_sequence(self, sequence_name: str) -> None:
        '''
        Remove a sequence from the alignment by sequence name.

        Parameters
        ----------
        sequence_name : str
            Sequence name to remove from alignment

        Returns
        -------
        None
        '''
        del self.sequences[sequence_name]

    def rename_sequence(self, old_name: str, new_name: str) -> None:
        '''
        Rename sequence in the alignment

        Parameters
        ----------
        old_name : str
            Old name of sequence
        new_name : str
            New name of sequence

        Returns
        -------
        None
        '''
        self.sequences[new_name] = self.sequences.pop(old_name)

    def write_pir(self, file_name: str, line_wrap: int = 50) -> None:
        '''
        Write alignment to file in the PIR file format.

        Parameters
        ----------
        file_name : str
            File name to write to
        line_wrap : int
            Characters per line (default 50)

        Returns
        -------
        None
        '''
        with open(file_name, 'w') as file_handler:
            for sequence_name, sequence in self.sequences.items():
                # write name and annotation
                file_handler.write('>P1;{}\n'.format(sequence_name))
                file_handler.write(
                    '{}\n'.format(sequence.get_annotation_pir()))
                # write sequence with newline every nth character
                file_handler.write(
                    re.sub(''.join(['(.{', str(line_wrap), '})']), '\\1\n',
                           sequence.sequence, 0, re.DOTALL))
                file_handler.write('\n*\n')

    def write_fasta(self, file_name: str, line_wrap: int = 80) -> None:
        '''
        Write alignment to file in the FASTA alignment file format.

        Parameters
        ----------
        file_name : str
            File name to write to
        line_wrap : int
            Characters per line (default 80)

        Returns
        -------
        None
        '''
        with open(file_name, 'w') as file_handler:
            for sequence_name, sequence in self.sequences.items():
                # write name
                file_handler.write('>{}\n'.format(sequence_name))
                # write sequence with newline every nth character
                file_handler.write(
                        re.sub(''.join(['(.{', str(line_wrap), '})']), '\\1\n',
                               sequence.sequence, 0, re.DOTALL))
                file_handler.write('\n')

    def print_clustal(self, line_wrap: int = 80) -> None:
        '''
        Print alignment to console in the clustal file format.

        Parameters
        ----------
        line_wrap : int
            Characters per line (default 80)

        Returns
        -------
        None
        '''
        # get sequence names and sequences and trim sequence_names
        sequences = [sequence.sequence for sequence in self.sequences.values()]
        sequence_names = [(name + '          ')[:10] for name in
                          self.sequences]
        # assemble output
        while len(sequences[0]) != 0:
            new_sequences = list()
            for sequence, sequence_name in zip(sequences, sequence_names):
                new_sequences.append(sequence[line_wrap:])
                print('{}  {}'.format(sequence_name, sequence[:line_wrap]))
            print('\n')
            sequences = new_sequences

    def write_clustal(self, file_name: str, line_wrap: int = 50) -> None:
        '''
        Write alignment to file in the clustal file format.

        Parameters
        ----------
        file_name : str
            File name to write to
        line_wrap : int
            Characters per line (default 50)

        Returns
        -------
        None
        '''
        # redirect output of self.print_clustal to file
        with open(file_name, 'w') as file_handler:
            with contextlib.redirect_stdout(file_handler):
                self.print_clustal(line_wrap=line_wrap)

    def remove_redundant_gaps(self) -> None:
        '''
        Remove gaps in the alignment that are present in every column.

        Returns
        -------
        None
        '''
        all_gaps = [sequence.get_gaps() for sequence in
                    self.sequences.values()]
        # get intersection of all gaps
        # adapted from https://stackoverflow.com/a/10066921/7912251
        intersection_gaps = list(set.intersection(*map(set, all_gaps)))
        # remove gap
        if len(intersection_gaps) != 0:
            for sequence in self.sequences.values():
                sequence.remove_gaps(positions=intersection_gaps)

    def replace_sequence(self, seq_name: str, new_sequence: str) -> None:
        '''
        Targeted replacement of sequence in alignment.

        Parameters
        ----------
        seq_name : str
            The identifier of the sequence that will be replaced.
        new_sequence : str
            The new sequence.

        Notes
        -----
        This replacement is designed to introduce missing residues from
        template structures into the alignment and therefore has very strict
        requirements. The new and old sequence have to be identical, except
        that the new sequence might contain unmodelled residues. These are
        indicated by the letter 'X' in the new sequence, and will result in a
        gap '-' in the alignment after replacement. It is important that all
        unmodelled residues, even at the start or beginning of the template
        sequence are correctly labeled as 'X'.

        Examples
        --------
        >>> aln = hm.Alignment(None)
        >>> aln.sequences = {
        ...     'seq1': hm.alignment.Sequence('seq1', 'AAAACCCCDDDD'),
        ...     'seq2': hm.alignment.Sequence('seq2', 'AAAAEEEEDDDD'),
        ...     'seq3': hm.alignment.Sequence('seq3', 'AAAA----DDDD')
        ...     }
        >>> replacement_seq1 = 'AAAAXXXXXDDD'
        >>> replacement_seq3 = 'AAXXXXDD'
        >>> aln.replace_sequence('seq1', replacement_seq1)
        >>> aln.print_clustal()
        seq1        AAAA-----DDD
        seq2        AAAAEEEEDDDD
        seq3        AAAA----DDDD
        >>> aln.replace_sequence('seq3', replacement_seq3)
        >>> aln.print_clustal()
        seq1        AAAA-----DDD
        seq2        AAAAEEEEDDDD
        seq3        AA--------DD
        '''
        # check if sequences fully match
        if re.fullmatch(
                new_sequence.upper().replace('X', r'\w'),
                self.sequences[seq_name].sequence.upper().replace('-', '')
                ) is None:
            raise ValueError(
                f'{seq_name}: New sequence does not match with sequence in '
                'alignment.')

        new_sequence = list(new_sequence.replace('-', '').upper())
        old_sequence = list(self.sequences[seq_name].sequence.upper())
        replaced_seq = list()
        for old_res in old_sequence:
            if old_res == '-':
                # transfer gaps
                replaced_seq.append(old_res)
            else:
                new_res = new_sequence.pop(0)
                if old_res == new_res:
                    # match
                    replaced_seq.append(new_res)
                elif old_res != new_res and new_res == 'X':
                    # introduce gaps for unmodelled residues
                    replaced_seq.append('-')
                else:
                    # after checking with re.fullmatch, this should never run
                    raise ValueError(
                        '{}: Mismatch detected'.format(seq_name))

        # update sequence
        self.sequences[seq_name].sequence = ''.join(replaced_seq)

    def calc_identity(self, sequence_name_1: str,
                      sequence_name_2: str) -> float:
        '''
        Calculation of sequence identity between two sequences in the
        alignment.

        Parameters
        ----------
        sequence_name_1, sequence_name_2 : str
            Sequence pair to calculate identity for

        Returns
        -------
        identity : float
            Sequence identity between the two sequences

        See Also
        --------
        calc_identity_target
        calc_pairwise_identity_all

        Notes
        -----
        There are mutiple ways of calculating sequence identity, which can be
        useful in different situations. Here implemented is one way which makes
        a lot of sence for evaluating templates for homology modelling. The
        sequence identity is calculated by dividing the number of matches by
        the length of sequence 1 (mismatches and gaps are handled identically,
        no gap compression).

        .. math::

             \\text{seqid} = \\frac{\\text{matches}}
             {\\text{length}(\\text{sequence1})}

        Examples
        --------

        Gaps and mismatches are treated equally.

        >>> aln = hm.Alignment(None)
        >>> aln.sequences = {
        ...     'seq1': hm.alignment.Sequence('seq1', 'AAAACCCCDDDD'),
        ...     'seq2': hm.alignment.Sequence('seq2', 'AAAAEEEEDDDD'),
        ...     'seq3': hm.alignment.Sequence('seq3', 'AAAA----DDDD')
        ...     }
        >>> aln.calc_identity('seq1', 'seq2')
        66.67
        >>> aln.calc_identity('seq1', 'seq3')
        66.67

        Normalization happens for the length of sequence 1, so the order of
        sequences matters.

        >>> aln = hm.Alignment(None)
        >>> aln.sequences = {
        ...     'seq1': hm.alignment.Sequence('seq1', 'AAAACCCCDDDD'),
        ...     'seq2': hm.alignment.Sequence('seq3', 'AAAA----DDDD')
        ...     }
        >>> aln.calc_identity('seq1', 'seq2')
        66.67
        >>> aln.calc_identity('seq2', 'seq1')
        100.0
        '''

        def _calc_identity(sequence_1: list, sequence_2: list) -> float:
            '''
            Helper function for calculating sequence identity

            seq_id = matches / length(sequence_1)

            Parameters
            ----------
            sequence_1, sequence_2: list
                Sequence string transformed to list

            Returns
            -------
            float
            '''
            matches = 0
            length = len([s for s in sequence_1 if s != '-'])
            for res_1, res_2 in zip(sequence_1, sequence_2):
                if res_1 == res_2 and res_1 != '-':
                    matches += 1
            return round(100 * matches / length, 2)

        # extract sequences of interest
        sequence_1 = list(self.sequences[sequence_name_1].sequence)
        sequence_2 = list(self.sequences[sequence_name_2].sequence)
        # calculate identity
        return _calc_identity(sequence_1, sequence_2)

    def calc_pairwise_identity_all(self) -> typing.Type['pd.DataFrame']:
        '''
        Calculate identity between all sequences in the alignment.

        Returns
        -------
        identities : pd.DataFrame
            Dataframe with pairwise sequence identites

        See Also
        --------
        calc_identity
        calc_identity_target

        Notes
        -----
        Calculates sequence identity as descripted for calc_identity:

        .. math::

             \\text{seqid} = \\frac{\\text{matches}}
             {\\text{length}(\\text{sequence1})}
        '''
        output = {'sequence_1': [], 'sequence_2': [], 'identity': []}
        # iterate over all pairs of sequences
        for sequence_name_1, sequence_name_2 in itertools.product(
                self.sequences.keys(), repeat=2):
            if not sequence_name_1 == sequence_name_2:
                output['sequence_1'].append(sequence_name_1)
                output['sequence_2'].append(sequence_name_2)
                output['identity'].append(self.calc_identity(
                    sequence_name_1, sequence_name_2))
        return pd.DataFrame(output)

    def calc_identity_target(
            self, sequence_name: str) -> typing.Type['pd.DataFrame']:
        '''
        Calculate identity of all sequences in the alignment to specified
        target sequence.

        Parameters
        ----------
        sequence_name : str
            Target sequence

        Returns
        -------
        identities : pd.DataFrame
            Dataframe with pairwise sequence identities

        See Also
        --------
        calc_identity
        calc_pairwise_identity_all

        Notes
        -----
        Calculates sequence identity as descripted for calc_identity:

        .. math::

             \\text{seqid} = \\frac{\\text{matches}}
             {\\text{length}(\\text{sequence1})}
        '''
        output = {'sequence_1': [], 'sequence_2': [], 'identity': []}
        for sequence_name_2 in self.sequences:
            if not sequence_name == sequence_name_2:
                output['sequence_1'].append(sequence_name)
                output['sequence_2'].append(sequence_name_2)
                output['identity'].append(self.calc_identity(
                    sequence_name, sequence_name_2))
        return pd.DataFrame(output)

    def calc_coverage(self, sequence_name_1: str,
                      sequence_name_2: str) -> float:
        '''
        Calculation of coverage of sequence 2 to sequence 1 in the alignment.

        Parameters
        ----------
        sequence_name_1, sequence_name_2 : str
            Sequence pair to calculate coverage for

        Returns
        -------
        coverage : float
            Coverage of sequence 2 to sequence 1

        See Also
        --------
        calc_coverage_target
        calc_pairwise_coverage_all

        Notes
        -----
        Coverage in this context means how many of the residues in sequences 1
        are assigned a residue in sequence 2. This is useful for evaluating
        potential templates, because a low sequence identity (as implemented in
        homelette) could be caused either by a lot of residues not being
        aligned at all, or a lot of residues being aligned but not with
        identical residues.

        .. math::

            \\text{coverage} = \\frac{\\text{aligned residues}}
            {\\text{length}(\\text{sequence1})}

        Examples
        --------

        Gaps and mismatches are not treated equally.

        >>> aln = hm.Alignment(None)
        >>> aln.sequences = {
        ...     'seq1': hm.alignment.Sequence('seq1', 'AAAACCCCDDDD'),
        ...     'seq2': hm.alignment.Sequence('seq2', 'AAAAEEEEDDDD'),
        ...     'seq3': hm.alignment.Sequence('seq3', 'AAAA----DDDD')
        ...     }
        >>> aln.calc_coverage('seq1', 'seq2')
        100.0
        >>> aln.calc_coverage('seq1', 'seq3')
        66.67

        Normalization happens for the length of sequence 1, so the order of
        sequences matters.

        >>> aln = hm.Alignment(None)
        >>> aln.sequences = {
        ...     'seq1': hm.alignment.Sequence('seq1', 'AAAACCCCDDDD'),
        ...     'seq2': hm.alignment.Sequence('seq3', 'AAAA----DDDD')
        ...     }
        >>> aln.calc_coverage('seq1', 'seq2')
        66.67
        >>> aln.calc_coverage('seq2', 'seq1')
        100.0
        '''
        def _calc_coverage(sequence_1: typing.Iterable, sequence_2:
                           typing.Iterable) -> float:
            '''
            Helper function for calculating sequence coverage

            coverage = aligned_res / length(sequence_1)

            Parameters
            ----------
            sequence_1, sequence_2 : Iterable
                Sequence string transformed to list

            Returns
            -------
            float
            '''
            aligned_res = 0
            length = len([s for s in sequence_1 if s != '-'])
            for res_1, res_2 in zip(sequence_1, sequence_2):
                if res_1 != '-' and res_2 != '-':
                    aligned_res += 1
            return round(100 * aligned_res / length, 2)

        # extract sequences of interest
        sequence_1 = list(self.sequences[sequence_name_1].sequence)
        sequence_2 = list(self.sequences[sequence_name_2].sequence)
        return _calc_coverage(sequence_1, sequence_2)

    def calc_coverage_target(
            self, sequence_name: str) -> typing.Type['pd.DataFrame']:
        '''
        Calculate coverage of all sequences in the alignment to specified
        target sequence.

        Parameters
        ----------
        sequence_name : str
            Target sequence

        Returns
        -------
        coverages : pd.DataFrame
            Dataframe with pairwise coverage

        See Also
        --------
        calc_coverage
        calc_pairwise_coverage_all

        Notes
        -----
        Calculates coverage as described for calc_coverage:

        .. math::

            \\text{coverage} = \\frac{\\text{aligned residues}}
            {\\text{length}(\\text{sequence1})}
        '''
        output = {
                'sequence_1': [],
                'sequence_2': [],
                'coverage': [],
                }
        for sequence_name_2 in self.sequences:
            if not sequence_name == sequence_name_2:
                output['sequence_1'].append(sequence_name)
                output['sequence_2'].append(sequence_name_2)
                output['coverage'].append(self.calc_coverage(
                    sequence_name, sequence_name_2))
        return pd.DataFrame(output)

    def calc_pairwise_coverage_all(self) -> typing.Type['pd.DataFrame']:
        '''
        Calculate coverage between all sequences in the alignment.

        Returns
        -------
        coverages : pd.DataFrame
            Dataframe with pairwise coverage

        See Also
        --------
        calc_coverage
        calc_coverage_target

        Notes
        -----
        Calculates coverage as described for calc_coverage:

        .. math::

            \\text{coverage} = \\frac{\\text{aligned residues}}
            {\\text{length}(\\text{sequence1})}
        '''
        output = {
                'sequence_1': [],
                'sequence_2': [],
                'coverage': [],
                }
        # iterate over all pairs of sequences
        for sequence_name_1, sequence_name_2 in itertools.product(
                self.sequences.keys(), repeat=2):
            if not sequence_name_1 == sequence_name_2:
                output['sequence_1'].append(sequence_name_1)
                output['sequence_2'].append(sequence_name_2)
                output['coverage'].append(self.calc_coverage(
                    sequence_name_1, sequence_name_2))
        return pd.DataFrame(output)


def assemble_complex_aln(*args: typing.Type['Alignment'], names:
                         dict) -> typing.Type['Alignment']:
    '''
    Assemble complex alignments compatible with MODELLER from individual
    alignments.

    Parameters
    ----------
    *args : Alignment
        The input alignments
    names : dict
        Dictionary instructing how sequences in the different alignment objects
        are supposed to be arranged in the complex alignment. The keys are the
        names of the sequences in the output alignments. The values are
        iterables of the sequence names from the input alignments in the order
        they are supposed to appaer in the output alignment. Any value that can
        not be found in the alignment signals that this position in the complex
        alignment should be filled with gaps.

    Returns
    -------
    Alignment
        Assembled complex alignment

    Examples
    --------
    >>> aln1 = hm.Alignment(None)
    >>> aln1.sequences = {
    ...     'seq1_1': hm.alignment.Sequence('seq1_1', 'HELLO'),
    ...     'seq2_1': hm.alignment.Sequence('seq2_1', 'H---I'),
    ...     'seq3_1': hm.alignment.Sequence('seq3_1', '-HI--')
    ...     }
    >>> aln2 = hm.Alignment(None)
    >>> aln2.sequences = {
    ...     'seq2_2': hm.alignment.Sequence('seq2_2', 'KITTY'),
    ...     'seq1_2': hm.alignment.Sequence('seq1_2', 'WORLD')
    ...     }
    >>> names = {'seq1': ('seq1_1', 'seq1_2'),
    ...          'seq2': ('seq2_1', 'seq2_2'),
    ...          'seq3': ('seq3_1', 'gaps')
    ...     }
    >>> aln_assembled = hm.alignment.assemble_complex_aln(
    ...     aln1, aln2, names=names)
    >>> aln_assembled.print_clustal()
    seq1        HELLO/WORLD
    seq2        H---I/KITTY
    seq3        -HI--/-----
    '''
    output_aln = Alignment(None)

    # every entry in the names dict is processed as follows:
    # {output_name: (input_names)}
    for output_name in names.keys():
        sequence = list()
        input_names = names[output_name]
        # for each output_name, extract sequences from input alignments that
        # correspond to input_names
        for i, aln in enumerate(args):
            try:
                sequence.append(aln.get_sequence(input_names[i]).sequence)
            except (KeyError, IndexError):
                sequence_len = len(next(iter(aln.sequences.values())).sequence)
                sequence.append('-' * sequence_len)
        # assemble sequence from parts and add to output alignment
        sequence = '/'.join(sequence)
        output_aln.sequences[output_name] = Sequence(output_name, sequence)
    return output_aln


class AlignmentGenerator(abc.ABC):
    '''
    Parent class for the auto-generation of alignments and template selection
    based on sequence input.

    Parameters
    ----------
    sequence : str
        Target sequence in 1 letter amino acid code.
    target : str
        The name of the target sequence (default "target").
    template_location : str
        Directory where processed templates will be stored (default
        "./templates/").

    Attributes
    ----------
    alignment : Alignment
        The alignment.
    target_seq : str
        The target sequence.
    target : str
        The name of the target sequence.
    template_location : str
        Directory where processed templates will be stored.
    state
        Dictionary describing the state of the AlignmentGenerator object

    Returns
    -------
    None
    '''
    def __init__(self, sequence: str, target: str = 'target',
                 template_location: str = './templates/') -> None:
        self.alignment = None
        self.target_seq = sequence
        self.target = target
        self.template_location = os.path.abspath(template_location)
        # Simple state machine, together with _check_state
        self.state = {
            'has_alignment': False,
            'is_processed': False,
            }

    @abc.abstractmethod
    def get_suggestion(self):
        '''
        Generate suggestion for templates and alignment
        '''

    @classmethod
    def from_fasta(cls, fasta_file: str, template_location: str =
                   './templates/') -> 'AlignmentGenerator':
        '''
        Generates an instance of the AlignemntGenerator with the first sequence
        in the fasta file.

        Parameters
        ----------
        fasta_file : str
            Fasta file from which the first sequence will be read.
        template_location : str
            Directory where processed templates will be stored (default
            "./templates/").

        Returns
        -------
        AlignmentGenerator

        Raises
        ------
        ValueError
            Fasta file not properly formatted
        '''
        with open(fasta_file, 'r') as file_handler:
            lines = file_handler.readlines()

        # check if proper file
        if not any((line.startswith('>') for line in lines)):
            raise ValueError(
                'Could not identify any line that marks beginning of a '
                'sequence block (staring with ">").')

        target_name, target_seq = None, str()
        for line in lines:
            if line.startswith('>') and target_name is None:
                target_name = line.replace('>', '').strip()
            elif line.startswith('>') and target_name is not None:
                break
            else:
                target_seq += line.replace('-', '').strip().upper()

        if len(target_seq) == 0:
            raise ValueError(
                'No sequence found in the first sequence block.')

        return cls(target_seq, target_name, template_location)

    def _check_state(self, has_alignment: bool = None, is_processed: bool =
                     None) -> None:
        '''
        Check state of the AlignmentGenerator object.

        Parameters
        ----------
        has_alignment : bool
            Assesses whether an alignment has already been generated.
        is_processedb : bool
            Assesses whether, based on the alignment, templates have already
            been downloaded and processed from the PDB.

        Raises
        ------
        RuntimeError
            Current state does not support requested behaviour.

        Notes
        -----
        Some functionality can only be performed with the object being
        in a certain state.
        Giving None as an input prevents the check for that particular state to
        occur.
        '''
        def perform_check(required, state):
            '''
            Helper function for checking the state

            Returns true if required is None
            '''
            if (required is None) or (required is state):
                return True
            return False

        if not (perform_check(has_alignment, self.state['has_alignment']) and
                perform_check(is_processed, self.state['is_processed'])):
            msg = (
                f'Current state does not support requested behaviour.\n'
                f'has_alignment: current: {self.state["has_alignment"]} '
                f'required: {has_alignment}\n'
                f'is_processed: current: {self.state["is_processed"]} '
                f'required: {is_processed}\n'
                )
            # add more detail if has_alignment has mismatch
            if not perform_check(has_alignment, self.state['has_alignment']):
                if not has_alignment and self.state['has_alignment']:
                    msg = msg + (
                        '\nYou have already generated an alignment.'
                        )
                if has_alignment and not self.state['has_alignment']:
                    msg = msg + (
                        '\nRequested action requires an alignment. Please call'
                        ' get_suggestion to generate an alignment.'
                        )
            # add more detail if is_processed as mismatch
            if not perform_check(is_processed, self.state['is_processed']):
                if not is_processed and self.state['is_processed']:
                    msg = msg + (
                        'You have already downloaded the template structures '
                        'from the PDB and updated the alignment.\n'
                        )
                if is_processed and not self.state['is_processed']:
                    msg = msg + (
                        'Requested action requires the templates to be '
                        'downloaded and the alignment to be processed. Please '
                        'call get_pdbs to perform these steps.\n'
                        )
            # raise error with assembled details
            raise RuntimeError(msg)

    def _guess_pdb_format_from_aln(self) -> str:
        '''
        Guess which organization layer of PDB is used in alignment.

        Returns
        -------
        pdb_format : str
            Can be one of `entry`, `polymer_entity`, or
            `polymer_entity_instance`.

        Raises
        ------
        ValueError
            PDB format could not be guessed.

        Notes
        -----
        The following organizational layers exist within the RCSB:

        * entry: The full PDB entry (i.e. 3NY5)
        * polymer_entity: On of the polymer entities under the entry (i.e.
        3NY5_1)
        * polymer_entity_instance: One of the instances of a polymer entity
        under the entry (i.e. 3NY5.A)
        '''
        # check state
        self._check_state(has_alignment=True, is_processed=None)
        # parse template ids from aln
        templates = [t for t in self.alignment.sequences.keys() if t !=
                     self.target]
        # identify pattern
        r_entry = re.compile(r'^[A-Za-z0-9]{4}')
        r_entity = re.compile(r'^[A-Za-z0-9]{4}[\W_][0-9]')
        r_instance = re.compile(r'^[A-Za-z0-9]{4}[\W_][A-Za-z]')
        if all((re.fullmatch(r_entry, template) for template in
                templates)):
            pdb_format = 'entry'
        elif all((re.fullmatch(r_entity, template) for template in
                  templates)):
            pdb_format = 'polymer_entity'
        elif all((re.fullmatch(r_instance, template) for template in
                  templates)):
            pdb_format = 'polymer_entity_instance'
        else:
            raise ValueError(
                'Unable to guess pdb_format from template names. Please '
                'make sure all template names in the alignment follow one of'
                ' the proposed naming schemes.')
        return pdb_format

    def show_suggestion(self, get_metadata: bool = False
                        ) -> typing.Type['pd.DataFrame']:
        '''
        Shows which templates have been suggested by the AlignmentGenerator, as
        well as some useful statistics (sequence identity, coverage).

        Parameters
        ----------
        get_metadata : bool
            Retrieve additional metadata (experimental method, resolution,
            structure title) from the RCSB.

        Returns
        -------
        suggestion : pd.DataFrame
            DataFrame with calculated sequence identity and sequence coverage
            for target

        Raises
        ------
        RuntimeError
            Alignment has not been generated yet

        See Also
        --------
        Alignment.calc_identity
        Alignment.calc_coverage

        Notes
        -----
        The standard output lists the templates in the alignment and shows both
        coverage and sequence identity to the target sequence. The templates
        are ordered by sequence identity.

        In addition, the experimental method (Xray, NMR or Electron
        Microscopy), the resolution (if applicable) and the title of the
        template structure can be retrieved from the RCSB. Retrieving metadata
        from the PDB requires a working internet connecction.
        '''
        self._check_state(has_alignment=True, is_processed=None)

        # calculate coverage and identity
        df_coverage = self.alignment.calc_coverage_target(self.target)
        df_identity = self.alignment.calc_identity_target(self.target)

        if get_metadata:
            # Fetch annotation from RCSB
            templates = list(
                (t[0:4] for t in self.alignment.sequences
                 if t != self.target))
            url = 'https://data.rcsb.org/graphql?'

            # query for structure annotation
            # Documentation of query API:
            # https://data.rcsb.org/index.html
            query = (
                f'query={{'
                f' entries(entry_ids: {templates!r}) {{'
                f'  rcsb_id'
                f'  struct {{'
                f'   title'
                f'  }}'
                f'  exptl {{'
                f'   method'
                f'  }}'
                f'  rcsb_entry_info {{'
                f'   resolution_combined'
                f'  }}'
                f' }}'
                f'}}'
                )
            # format query
            query = query.replace("'", '"')
            # encode URL
            query = urllib.parse.quote(query, safe='=():,')

            # access query
            with urllib.request.urlopen(url + query) as response:
                # check status
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                else:
                    raise urllib.error.URLError(
                        f'Unknown URL status: expected 200, got '
                        f'{response.status}')

            # extract data
            annot_id = []
            for r in response_decoded['data']['entries']:
                annot_id.append([
                    r['rcsb_id'],
                    r['exptl'][0]['method'],
                    r['rcsb_entry_info']['resolution_combined'][0],
                    r['struct']['title'],
                    ])
            df_annotation = pd.DataFrame(annot_id, columns=[
                'pdbid', 'method', 'resolution', 'title'])

        # combine data frames
        output = (
            pd.merge(df_coverage, df_identity, on=('sequence_1', 'sequence_2'))
            # remove column with sequence_1 and rename sequence_2
            .drop('sequence_1', axis=1)
            # sort values
            .sort_values(by=['identity', 'coverage'], ascending=[False, False])
            .rename({'sequence_2': 'template'}, axis=1)
            )
        if get_metadata:
            output = (
                output
                .assign(pdbid=lambda df: df['template'].map(
                    lambda template: template[0:4]))
                )
            output = (
                pd.merge(output, df_annotation, on=('pdbid', 'pdbid'))
                # sort values
                .sort_values(by=['identity', 'coverage'],
                             ascending=[False, False])
                # remove merge column
                .drop('pdbid', axis=1)
                )

        return output

    def select_templates(self, templates: typing.Iterable) -> None:
        '''
        Select templates from suggested templates by identifier.

        Parameters
        ----------
        templates : iterable
            The selected templates as an interable.

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            Alignment has not been generated yet
        '''
        self._check_state(has_alignment=True, is_processed=None)
        selection = [self.target] + list(templates)
        self.alignment.select_sequences(selection)
        self.alignment.remove_redundant_gaps()

    def get_pdbs(self, pdb_format: str = 'auto', verbose: bool = True) -> None:
        '''
        Downloads and processes templates present in alignment.

        Parameters
        ----------
        pdb_format : str
            Format of PDB identifiers in alignment (default auto)
        verbose : bool
            Explain what operations are performed

        Raises
        ------
        RuntimeError
            Alignment has not been generated yet
        ValueError
            PDB format could not be guessed

        Notes
        -----
        pdb_format tells the function how to parse the template identifiers in
        the alignment:

        * auto: Automatic guess for pdb_format
        * entry: Sequences are named only be their PDB identifier (i.e. 4G0N)
        * entity: Sequences are named in the format PDBID_ENTITY (i.e. 4G0N_1)
        * instance: Sequences are named in the format PDBID_CHAIN (i.e. 4G0N_A)

        Please make sure that all templates follow one naming convention, and
        that there are no sequences in the alignment that violate the naming
        convention (except the target sequence).

        During the template processing, all hetatms will be remove from the
        template, as well as all other chains. All chains will be renamed to
        "A" and the residue number will be set to 1 on the first residue. The
        corresponding annotations are automatically made in the alignment
        object.
        '''
        # check state
        self._check_state(has_alignment=True, is_processed=False)

        # Helper functions
        # verbose behaviour
        # adapted from https://stackoverflow.com/a/5980173/7912251
        if verbose:
            def vprint(*args, **kwargs):
                for arg in args:
                    print(arg, **kwargs)
        else:
            def vprint(*args):
                pass

        def get_entities_from_entries(templates, alignment) -> (dict, dict):
            '''
            For an alignment with PDB entry identifiers (i.e. 1LFD), create
            mapping of entry identifiers to the entity identifiers (i.e.
            1LFD_1) and download entity information.
            '''
            # get number of entities for each entry
            url = 'https://data.rcsb.org/graphql?'
            query = (
                f'query={{'
                f' entries(entry_ids: '
                f' {templates!r}) {{'
                f'  rcsb_id'
                f'  rcsb_entry_info {{'
                f'   polymer_entity_count_protein'
                f'  }}'
                f' }}'
                f'}}'
                )
            query = query.replace("'", '"')
            query = urllib.parse.quote(query, safe='=():,')

            with urllib.request.urlopen(url + query) as response:
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                else:
                    raise urllib.error.URLError(
                        f'Unknown URL status: expected 200, got '
                        f'{response.status}')

            entities = list()
            for r in response_decoded['data']['entries']:
                entry = r['rcsb_id']
                entity_count = (
                    r['rcsb_entry_info']['polymer_entity_count_protein'])
                for entity_num in range(1, entity_count+1):
                    entities.append(f'{entry}_{entity_num}')

            # pull sequence information for all entities
            query = (f'query={{'
                     f' polymer_entities (entity_ids: {entities!r}) {{'
                     f'  rcsb_id'
                     f'  entity_poly {{'
                     f'   pdbx_seq_one_letter_code_can'
                     f'   pdbx_strand_id'
                     f'  }}'
                     f' }}'
                     f'}}'
                     )
            query = query.replace("'", '"')
            query = urllib.parse.quote(query, safe='=():,')

            with urllib.request.urlopen(url + query) as response:
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                else:
                    raise urllib.error.URLError(
                            f'Unknown URL status: expected 200, got '
                            f'{response.status}')

            entities = dict()
            for r in response_decoded['data']['polymer_entities']:
                entity = r['rcsb_id']
                seq = r['entity_poly']['pdbx_seq_one_letter_code_can']
                chains = r['entity_poly']['pdbx_strand_id'].split(',')
                entities[entity] = {
                        'sequence': seq,
                        'chains': chains}

            # match entries with entities based on entities
            mapping = dict()
            for template in alignment.sequences:
                for entity in entities:
                    if template != entity[0:4]:
                        continue
                    seq_alignment = (
                        alignment.sequences[template].sequence
                        .replace('-', ''))
                    seq_entity = entities[entity]['sequence']
                    seq_match = re.search(seq_alignment, seq_entity)
                    if seq_match is not None:
                        mapping[entity] = [{'entry': template,
                                            'aln_name': template}]

            # filter entities so that only matched entities remain
            entities = {
                entity: seq for entity, seq in entities.items() if entity in
                mapping.keys()}

            return (mapping, entities)

        def get_entities_from_instances(templates, alignment) -> (dict, dict):
            '''
            For a list of templates with PDB polymer instance identifiers (i.e.
            1LFD_A), create mapping of instance identifiers to the entity
            identifiers (i.e. 1LFD_1) and download entity information.
            '''
            # make sure instances are correctly formatted
            templates_formatted = {f'{t[0:4]}.{t[5]}': t for t in templates}
            # pull entity information for all instances
            url = 'https://data.rcsb.org/graphql?'
            query = (
                f'query={{'
                f' polymer_entity_instances(instance_ids: '
                f' {list(templates_formatted.keys())!r}) {{'
                f'  rcsb_id'
                f'  rcsb_polymer_entity_instance_container_identifiers {{'
                f'   entity_id'
                f'  }}'
                f' }}'
                f'}}'
                )

            query = query.replace("'", '"')
            query = urllib.parse.quote(query, safe='=():,')

            with urllib.request.urlopen(url + query) as response:
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                else:
                    raise urllib.error.URLError(
                        f'Unknown URL status: expected 200, got '
                        f'{response.status}')

            # parse entity information
            mapping = dict()
            for r in response_decoded['data']['polymer_entity_instances']:
                instance = r['rcsb_id']
                entity = (
                    instance[0:4] + '_' +
                    r['rcsb_polymer_entity_instance_container_identifiers']
                    ['entity_id'])
                if entity in mapping.keys():
                    mapping[entity].append(
                        {'instance': instance,
                         'aln_name': templates_formatted[instance]})
                else:
                    mapping[entity] = [
                            {'instance': instance,
                             'aln_name': templates_formatted[instance]}]

            # pull sequences for entites from RCSB
            query = (f'query={{'
                     f' polymer_entities (entity_ids:'
                     f' {list(mapping.keys())!r}) {{'
                     f'  rcsb_id'
                     f'  entity_poly {{'
                     f'   pdbx_seq_one_letter_code_can'
                     f'   pdbx_strand_id'
                     f'  }}'
                     f' }}'
                     f'}}'
                     )
            query = query.replace("'", '"')
            query = urllib.parse.quote(query, safe='=():,')

            with urllib.request.urlopen(url + query) as response:
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                else:
                    raise urllib.error.URLError(
                            f'Unknown URL status: expected 200, got '
                            f'{response.status}')

            entities = dict()
            for r in response_decoded['data']['polymer_entities']:
                entity = r['rcsb_id']
                seq = r['entity_poly']['pdbx_seq_one_letter_code_can']
                chains = r['entity_poly']['pdbx_strand_id'].split(',')
                entities[entity] = {
                        'sequence': seq,
                        'chains': chains}

            # check that, if multiple instances are assigned a common identity,
            # their alignment is identical
            for entity, instances in dict(mapping).items():
                if len(instances) > 1:
                    seqs_instances = [
                        alignment.sequences[instance['aln_name']].sequence for
                        instance in instances]
                    if len(set(seqs_instances)) != 1:
                        # if multiple instance with different alignments:
                        # choose alignment with highest sequence identity
                        # and discard alternative alignments
                        temp_aln = Alignment(None)
                        temp_aln.sequences[self.target] = (
                                alignment.sequences[self.target])
                        for instance in instances:
                            aln_name = instance['aln_name']
                            sequence = alignment.get_sequence(aln_name)
                            temp_aln.sequences[aln_name] = sequence
                        # calculate identities
                        best_instance = (
                            temp_aln.calc_identity_target(self.target)
                            .sort_values('identity', ascending=False)
                            ['sequence_2'][0]
                            )
                        mapping[entity] = [
                                instance for instance in instances
                                if instance['aln_name'] == best_instance]
                        warnings.warn(
                            f'Entity: {entity}\n'
                            f'Found instances of the same PDB entity, but with'
                            f' different alignments. Since all instances of an'
                            f' entity share a sequence, this should not '
                            f'happen.\n'
                            f'Automatically selected instance with highest '
                            f'sequence identity: {best_instance}\n'
                            f'To avoid this behaviour, please select '
                            f'one of the instances before running '
                            f'get_pdbs.')

            return (mapping, entities)

        def get_entities_from_entities(templates, alignment) -> (dict, dict):
            '''
            Process templates and alignment to mapping of entities to entries
            in the alignment and download entity information.
            '''
            # make sure entity names are properly formatted
            mapping = {f'{t[0:4]}_{t[5]}': [{
                'entity': f'{t[0:4]}_{t[5]}',
                'aln_name': t}] for t in templates}

            # pull sequences for entites from RCSB
            url = 'https://data.rcsb.org/graphql?'
            query = (f'query={{'
                     f' polymer_entities (entity_ids:'
                     f' {list(mapping.keys())!r}) {{'
                     f'  rcsb_id'
                     f'  entity_poly {{'
                     f'   pdbx_seq_one_letter_code_can'
                     f'   pdbx_strand_id'
                     f'  }}'
                     f' }}'
                     f'}}'
                     )
            query = query.replace("'", '"')
            query = urllib.parse.quote(query, safe='=():,')

            with urllib.request.urlopen(url + query) as response:
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                else:
                    raise urllib.error.URLError(
                            f'Unknown URL status: expected 200, got '
                            f'{response.status}')

            entities = dict()
            for r in response_decoded['data']['polymer_entities']:
                entity = r['rcsb_id']
                seq = r['entity_poly']['pdbx_seq_one_letter_code_can']
                chains = r['entity_poly']['pdbx_strand_id'].split(',')
                entities[entity] = {
                        'sequence': seq,
                        'chains': chains}

            return (mapping, entities)

        def adjust_template_seq(seq_alignment: str,
                                seq_template: str):
            '''
            Adjusts the sequence from the template to the sequence present in
            the alignment

            Check match
            Pad right and left with missing residues
            '''
            # check if template sequence is in alignment sequence
            seq_match = re.search(seq_template.replace('X', r'\w'),
                                  seq_alignment)

            if seq_match is None:
                raise RuntimeError(
                    'Could not match template sequence with aligned sequence.')

            # pad template sequence if necessary
            seq_template_padded = (
                    seq_template
                    .rjust(len(seq_template) + seq_match.start(), 'X')
                    .ljust(len(seq_alignment), 'X'))

            return seq_template_padded

        def add_seq_to_aln(aln, seq_name, old_sequence, new_sequence):
            '''
            Adds sequence to alignment, then performs sequence replacement and
            annotates the sequence.
            '''
            annotations = {
                'seq_type': 'structure',
                'pdb_code': seq_name,
                'begin_res': '1',
                'begin_chain': 'A',
                }
            aln.sequences.update({
                seq_name: Sequence(seq_name, old_sequence, **annotations)})
            aln.replace_sequence(seq_name, new_sequence)

            return aln

        # list of 3 letter amino acid code
        amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
                       'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                       'THR', 'TRP', 'TYR', 'VAL', 'SEC']

        # check pdb_format
        templates = [template for template in self.alignment.sequences.keys()
                     if template != self.target]
        if pdb_format == 'auto':
            vprint('Guessing template naming format...')
            pdb_format = self._guess_pdb_format_from_aln()
            vprint(f'Template naming format guessed: {pdb_format}!\n')
        elif pdb_format not in ['entry', 'polymer_entity',
                                'polymer_entity_instance']:
            raise ValueError(
                f'Invalid value to pdb_format: {pdb_format}\nHas to be one of '
                f'"auto", "entry", "polymer_entity", '
                f'"polymer_entity_instance".')

        # Initialize template dir
        vprint('Checking template dir...')
        if not os.path.exists(self.template_location):
            vprint('Template dir not found...')
            os.makedirs(self.template_location, exist_ok=True)
            vprint(f'New template dir created at\n'
                   f'"{self.template_location}"!\n')
        else:
            vprint('Template dir found!\n')

        # depending on the format, extract mapping of templates to pdb entities
        # as well as sequences and chains for these entities from PDB
        vprint('Processing templates:\n')
        if pdb_format == 'entry':
            mapping, entities = get_entities_from_entries(
                    templates, self.alignment)
        elif pdb_format == 'polymer_entity':
            mapping, entities = get_entities_from_entities(
                    templates, self.alignment)
        elif pdb_format == 'polymer_entity_instance':
            mapping, entities = get_entities_from_instances(
                    templates, self.alignment)

        # initialize changed alignment
        new_aln = Alignment(None)
        new_aln.sequences = {
                self.target: self.alignment.sequences[self.target]}
        new_aln.get_sequence(self.target).annotate(seq_type='sequence')

        # compare sequences from the alignment and the entities to find out if
        # template structures need to be clipped
        # (i.e. hhblits does not necessary include the full template sequence
        # in its alignment, so to properly process it, we need to find out
        # which parts are included and remove the rest from the template)
        for entity in mapping:
            # get sequences
            seq_alignment = (
                self.alignment.sequences[mapping[entity][0]['aln_name']]
                .sequence.replace('-', ''))
            seq_entity = entities[entity]['sequence']
            # find first set of indexes: where to clip the alignment
            try:
                index_first_res_clipped, index_last_res_clipped = re.search(
                    seq_alignment, seq_entity).span()
            except AttributeError:
                # Could not find the sequence from the alignment in in entity
                # sequence
                msg = (
                    f'{mapping[entity][0]["aln_name"]}\n'
                    f'Could not find sequence from alignment in sequence of '
                    f'assigned entity. This can happen if the chains were '
                    f'reassigned in the PDB after the pdb70 database was '
                    f'created.\n'
                    f'This sequence will be skipped. In order to use it, '
                    f' please save your alignment and change the name of the'
                    f' sequence')
                warnings.warn(msg)
                continue
            # download template
            vprint(f'{entity[0:4]} downloading from PDB...')
            try:
                pdb = pdb_io.download_pdb(entity[0:4])
            except urllib.error.URLError:
                vprint(f'{entity[0:4]} could not be downloaded...\n')
                continue
            vprint(f'{entity[0:4]} downloaded!')
            # iterate over chains
            for chain in entities[entity]['chains']:
                pdb_chain = pdb.transform_extract_chain(chain)
                vprint(f'{entity[0:4]}_{chain}: Chain extracted!')
                seq_template = pdb_chain.get_sequence(ignore_missing=False)
                try:
                    seq_template_padded = adjust_template_seq(
                            seq_entity, seq_template)
                except RuntimeError:
                    msg = (
                        f'{entity[0:4]}_{chain}: Could not match sequence in '
                        f'alignment to sequence from structure.')
                    # raise RuntimeError(msg) from exc
                    warnings.warn(msg)
                    continue
                # get first and last residue from template pdb
                pdb_chain_onlyprotein = pdb_chain.transform_filter_res_name(
                        amino_acids, mode='in')
                template_res_span = (
                    int(pdb_chain_onlyprotein.lines[0][22:26]),
                    int(pdb_chain_onlyprotein.lines[-1][22:26]))
                # get index of first and last residue in alignment
                for i, res in enumerate(seq_template_padded):
                    if res != 'X':
                        index_first_res_template = i
                        break
                for i, res in enumerate(seq_template_padded[::-1]):
                    if res != 'X':
                        index_last_res_template = len(seq_template_padded) - i
                        break
                # calculate residue range to clip from indices
                lower = (template_res_span[0] + index_first_res_clipped -
                         index_first_res_template)
                upper = (template_res_span[1] + index_last_res_clipped -
                         index_last_res_template)
                # adjust template structure
                pdb_chain = pdb_chain.transform_filter_res_seq(lower, upper)
                # update alignment
                new_aln = add_seq_to_aln(
                    new_aln, f'{entity[0:4]}_{chain}',
                    (self.alignment.sequences[mapping[entity][0]['aln_name']]
                     .sequence),
                    (seq_template_padded
                     [index_first_res_clipped:index_last_res_clipped]))
                vprint(f'{entity[0:4]}_{chain}: Alignment updated!')

                # process template
                pdb_chain = (
                    pdb_chain
                    .transform_remove_hetatm()
                    .transform_renumber_residues(starting_res=1)
                    .transform_change_chain_id(new_chain_id='A'))
                pdb_chain.write_pdb(os.path.join(
                    self.template_location, f'{entity[0:4]}_{chain}.pdb'))
                vprint(f'{entity[0:4]}_{chain}: PDB processed!')

        # check alignment for empty sequences (might happen during processing)
        for template, sequence in dict(new_aln.sequences).items():
            sequence = sequence.sequence
            if sequence == len(sequence) * '-':
                vprint(f'{template}: Adjusting the sequence to the template '
                       f'resulted in an empty sequence.\nRemoving sequence '
                       f'from the alignment.')
                del new_aln.sequences[template]
                if os.path.exists(os.path.join(
                        self.template_location, template + '.pdb')):
                    os.remove(os.path.join(
                        self.template_location, template + '.pdb'))

        # update alignment
        self.alignment = new_aln
        # update state
        self.state['is_processed'] = True
        vprint(f'\nFinishing... All templates successfully\ndownloaded and '
               f'processed!\nTemplates can be found in\n'
               f'"{self.template_location}".')

    def initialize_task(self, task_name: str = None, overwrite: bool = False,
                        task_class: Task = Task) -> Task:
        '''
        Initialize a homelette Task object for model generation and evaluation.

        Parameters
        ----------
        task_name : str
            The name of the task to initialize. If None, initialize as
            models_{target}.
        overwrite : bool
            Whether to overwrite the task directory if a directory of the same
            name already exists (default False).
        task_class : Task
            The class to initialize the Task with. This makes it possible to
            define custom child classes of Task and construct them from this
            function (default Task)

        Returns
        -------
        Task

        Raises
        ------
        RuntimeError
            Alignment has not been generated or templates have not been
            downloaded and processed.
        '''
        # check state
        self._check_state(has_alignment=True, is_processed=True)
        if task_name is None:
            task_name = f'models_{self.target}'
        return task_class(
                task_name=task_name,
                overwrite=overwrite,
                target=self.target,
                alignment=self.alignment)


class AlignmentGenerator_pdb(AlignmentGenerator):
    '''
    Identification of templates using the RCSB search API, generation of
    alignment using Clustal Omega and download and processing of template
    structures.

    Parameters
    ----------
    sequence : str
        Target sequence in 1 letter amino acid code.
    target : str
        The name of the target sequence (default "target").
    template_location : str
        Directory where processed templates will be stored (default
        "./templates/").

    Attributes
    ----------
    alignment : Alignment
        The alignment.
    target_seq : str
        The target sequence.
    target : str
        The name of the target sequence.
    template_location : str
        Directory where processed templates will be stored.
    state
        Dictionary describing the state of the AlignmentGenerator object

    Returns
    ----------
    None

    Notes
    -----
    The AlignmentGenerator uses the RCSB Search API [1]_ to identify potential
    template structures given the target sequence using MMseq2 [2]_. The
    sequences of the potentially downloaded and locally aligned using Clustal
    Omega [3]_ [4]_.

    References
    ----------
    .. [1] Rose, Y., Duarte, J. M., Lowe, R., Segura, J., Bi, C., Bhikadiya,
       C., Chen, L., Rose, A. S., Bittrich, S., Burley, S. K., & Westbrook, J.
       D. (2021). RCSB Protein Data Bank: Architectural Advances Towards
       Integrated Searching and Efficient Access to Macromolecular Structure
       Data from the PDB Archive. Journal of Molecular Biology, 433(11),
       166704. https://doi.org/10.1016/J.JMB.2020.11.003

    .. [2] Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive
       protein sequence searching for the analysis of massive data sets.
       Nature Biotechnology 2017 35:11, 35(11), 1026–1028.
       https://doi.org/10.1038/nbt.3988

    .. [3] Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K.,
       Li, W., Lopez, R., McWilliam, H., Remmert, M., Söding, J., Thompson, J.
       D., & Higgins, D. G. (2011). Fast, scalable generation of high-quality
       protein multiple sequence alignments using Clustal Omega. Molecular
       Systems Biology, 7(1), 539. https://doi.org/10.1038/MSB.2011.75

    .. [4] Sievers, F., & Higgins, D. G. (2018). Clustal Omega for making
       accurate alignments of many protein sequences. Protein Science, 27(1),
       135–145. https://doi.org/10.1002/PRO.3290
    '''
    def get_suggestion(self, seq_id_cutoff: float = 0.5, min_length: int = 30,
                       max_results: int = 50, xray_only: bool = True,
                       verbose: bool = True) -> None:
        '''
        Identifies potential templates, retrieves their sequences and aligns
        them locally using Clustal Omega.

        Parameters
        ----------
        seq_id_cutoff : float
            The sequence identity cutoff for the identification of template
            structures. Templates below this threshold will be ignored (default
            0.5).
        min_length : int
            The minimum length of template sequence to be included in the
            results (default 30 amino acids).
        max_results : int
            The number of results returned (default 50).
        xray_only : bool
            Only consider templates structures generated with X-ray
            crystallography (default True).
        verbose : bool
            Explain what is done (default True).

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            Alignment already generated.
        '''
        # check state
        self._check_state(has_alignment=False, is_processed=False)

        # Helper functions
        # verbose behaviour
        # adapted from https://stackoverflow.com/a/5980173/7912251
        if verbose:
            def vprint(*args, **kwargs):
                for arg in args:
                    print(arg, **kwargs)
        else:
            def vprint(*args):
                pass

        def query_pdb(sequence: str, seq_id_cutoff: float = 0.5,
                      min_length: int = 30, max_results: int = 50) -> list:
            '''
            Queries sequence with sequence identity cutoff against the PDB data
            base.

            Parameters
            ----------
            sequence : str
                Sequence as string
            seq_id_cutoff : float
                Sequence identity cutoff, between 0 and 1 (default 0.5)
            min_length : int
                Minimal length for potential templates (default 30)
            max_results : int
                Maximum number of results to return (default 50)

            Returns
            -------
            templates : list
            '''
            # assembly query
            url = 'https://search.rcsb.org/rcsbsearch/v1/query?json='
            # documentation for query structures can be found at
            # https://search.rcsb.org/index.html
            if xray_only:
                query_experimentalmethod = '''
                  {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                      "attribute": "exptl.method",
                      "operator": "exact_match",
                      "negation": false,
                      "value": "X-RAY DIFFRACTION"
                    }
                  },
                '''
            else:
                query_experimentalmethod = ''

            query = f'''
            {{
              "query": {{
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                  {{
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {{
                      "evalue_cutoff": 1,
                      "identity_cutoff": {seq_id_cutoff},
                      "target": "pdb_protein_sequence",
                      "value": "{sequence}"
                    }}
                  }},
                  {query_experimentalmethod}
                  {{
                    "type": "terminal",
                    "service": "text",
                    "parameters": {{
                      "attribute": "entity_poly.rcsb_sample_sequence_length",
                      "operator": "greater_or_equal",
                      "negation": false,
                      "value": {min_length}
                    }}
                }}
                ]
              }},
              "request_options": {{
                "pager": {{
                  "start": 0,
                  "rows": {max_results}
                }},
                "scoring_strategy": "sequence"
              }},
              "return_type": "polymer_entity"
            }}
            '''
            # format query
            # remove whitespaces
            query = re.sub(r'\s\s+', '', query)
            # encode URL
            query = urllib.parse.quote_plus(query)

            # access query
            with urllib.request.urlopen(url + query) as response:
                # check status
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                elif response.status == 204:
                    return []
                else:
                    raise urllib.error.URLError(
                            f'Unknown URL status: expected 200 or 204, got '
                            f'{response.status}')

            # extract templates from response
            templates = list()
            for hit in response_decoded['result_set']:
                templates.append(hit['identifier'])

            return templates

        def generate_alignment(templates):
            '''
            Generate alignment based on list of PDB entities.
            '''
            vprint('Retrieving sequences...')
            # download fastas for entities
            url = 'https://data.rcsb.org/graphql?'
            query = (f'query={{'
                     f' polymer_entities (entity_ids: {templates!r}) {{'
                     f'  rcsb_id'
                     f'  entity_poly {{'
                     f'   pdbx_seq_one_letter_code_can'
                     f'  }}'
                     f' }}'
                     f'}}'
                     )
            # format query
            query = query.replace("'", '"')
            # encode URL
            query = urllib.parse.quote(query, safe='=():,')

            # access query
            with urllib.request.urlopen(url + query) as response:
                # check status
                if response.status == 200:
                    response_decoded = json.loads(response.read().decode(
                        'utf-8'))
                else:
                    raise urllib.error.URLError(
                            f'Unknown URL status: expected 200, got '
                            f'{response.status}')

            # parse response to file
            fasta_file_name = os.path.realpath(
                    f'.sequences_{self.target}.fa')
            with open(fasta_file_name, 'w') as file_handle:
                file_handle.write(f'>{self.target}\n')
                file_handle.write(f'{self.target_seq}\n')
                for r in response_decoded['data']['polymer_entities']:
                    template = r['rcsb_id']
                    sequence = r['entity_poly']['pdbx_seq_one_letter_code_can']
                    file_handle.write(f'>{template}\n')
                    file_handle.write(f'{sequence}\n')
            vprint('Sequences succefully retrieved!\n')

            # perform alignment with clustalo
            vprint('Generating alignment...')
            aln_file_name = os.path.realpath(
                    f'aln_{self.target}.fasta_aln')
            command = [
                'clustalo', '-i', fasta_file_name, '-o', aln_file_name,
                '--force']
            subprocess.run(command, stdout=None, check=True, shell=False)

            aln = Alignment(aln_file_name)
            # remove files
            for file_name in [fasta_file_name, aln_file_name]:
                if os.path.exists(file_name):
                    os.remove(file_name)

            return aln

        # send query
        vprint('Querying PDB...')
        templates = query_pdb(self.target_seq, seq_id_cutoff, min_length,
                              max_results)
        if len(templates) == 0:
            print(f'Query found no potential templates with current '
                  f'parameters.\nseq_id_cutoff:\t{seq_id_cutoff}\n'
                  f'min_length:\t{min_length}')
            return None
        vprint(f'Query successful, {len(templates)} found!\n')
        # generate alignment
        self.alignment = generate_alignment(templates)
        vprint('Alignment generated!\n')
        # update state
        self.state['has_alignment'] = True
        # call show_suggestion
        vprint(f'Query successful.\n{len(self.alignment.sequences.keys())} '
               f'sequences have been found.')


# NOTES
# github https://github.com/soedinglab/hh-suite
# databases http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs
class AlignmentGenerator_hhblits(AlignmentGenerator):
    '''
    Identification of templates using hhblits to search a local PDB database,
    generation of alignment by combining pairwise alignments of target and
    template together.

    Parameters
    ----------
    sequence : str
        Target sequence in 1 letter amino acid code.
    target : str
        The name of the target sequence (default "target").
    template_location : str
        Directory where processed templates will be stored (default
        "./templates/").

    Attributes
    ----------
    alignment : Alignment
        The alignment.
    target_seq : str
        The target sequence.
    target : str
        The name of the target sequence.
    template_location : str
        Directory where processed templates will be stored.
    state
        Dictionary describing the state of the AlignmentGenerator object.

    Returns
    -------
    None

    Notes
    -----
    HHblits from the HHsuite [5]_ is used to query the databases. The resulting
    pairwise sequence alignments of template to target are combined using the
    target sequence as the master sequence. The resulting alignment is
    therefore, strictly speaking, not a proper multiple sequence alignment.
    However, all information from the pairwise alignments is preserved, and for
    homology modelling, the alignments of templates among each others do not
    have any influence.

    References
    ----------

    .. [5] Sievers, F., & Higgins, D. G. (2018). Clustal Omega for making
       accurate alignments of many protein sequences. Protein Science, 27(1),
       135–145. https://doi.org/10.1002/PRO.3290

    '''
    def get_suggestion(
            self, database_dir: str = './databases/',
            use_uniref: bool = False, evalue_cutoff: float = 0.001,
            iterations: int = 2, n_threads: int = 2, neffmax: float = 10.0,
            verbose: bool = True) -> None:
        '''
        Use HHblits to identify template structures and create a multiple
        sequence alignment by combination of pairwise alignments on target
        sequence.

        Parameters
        ----------
        database_dir : str
            The directory where the pdb70 (and the UniRef30) database are
            stored (default ./databases/).
        use_uniref : bool
            Use UniRef30 to create a MSA before querying the pdb70 database
            (default False). This leads to better results, however it takes
            longer and requires the UniRef30 database on your system.
        evalue_cutoff : float
            E-value cutoff for inclusion in result alignment (default 0.001)
        iterations : int
            Number of iterations when querying the pdb70 database.
        n_threads : int
            Number of threads when querying the pdb70 (or UniRef30) database
            (default 2).
        neffmax : float
            The neffmax value used when querying the pdb70 database
            (default 10.0).
        verbose : bool
            Explain which operations are performed (default True).

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            Alignment has already been generated.

        Notes
        -----
        This function expects "hhblits" to be installed and in the path. In
        addition, the pdb70 database needs to be downloaded and extracted in
        the database_dir. The files need to be called "pdb70_*" for hhblits to
        correctly find the database. If UniRef30 is used to create a
        pre-alignment for better results, the UniRef30 database needs to be
        downloaded and extracted in the database_dir. The files need to be
        called "UniRef30_*".

        For more information on neffmax, please check the hhblits
        documentation.

        If UniRef30 is used to generate a prealignment, then hhblits will be
        called for one iteration with standard parameters.
        '''
        # check state
        self._check_state(has_alignment=False, is_processed=False)

        # init verbose
        if verbose:
            def vprint(*args, **kwargs):
                for arg in args:
                    print(arg, **kwargs)
        else:
            def vprint(*args):
                pass

        # Helper functions and classes for parsing HHR files
        class Assembler():
            '''
            Helper class for assembling a MSA-like data structure from multiple
            pairwise alignments.

            Performs assembly based on multiple buffers and a master sequence.

            All pairwise columns between query (master) and templates (slaves)
            are conserved.
            '''
            def __init__(self, master_seq, master_name, aln_buffers):
                self.master_seq = master_seq
                self.master_name = master_name
                self.aln_buffers = aln_buffers

                self.index = 0

                # initialize data structure
                self.alignment = {master_name: list()}
                for buffer in self.aln_buffers:
                    template_name = buffer.slave_name
                    self.alignment[template_name] = list()

                # initialize indexes in buffers
                for buffer in self.aln_buffers:
                    buffer_seq = buffer.master_seq.replace('-', '')
                    buffer.set_index(re.search(buffer_seq, master_seq).start())

            def next_col(self):
                # select buffers which are in index
                selected_buffers = [
                    buffer for buffer in self.aln_buffers if (
                        self.index >= buffer.master_index)
                    and not buffer.is_empty()]

                # get next positions
                next_position = self.master_seq[self.index]
                next_positions_buffers = [
                    buffer.present_next_master() for buffer in
                    selected_buffers]

                # if all agree, insert
                if all([next_position_buffer == next_position for
                        next_position_buffer in next_positions_buffers]):
                    processed = [self.master_name]
                    # insert position in master
                    self.alignment[self.master_name].append(next_position)
                    # insert position in templates
                    for buffer in selected_buffers:
                        template_name = buffer.slave_name
                        self.alignment[template_name].append(buffer.pop_next())
                        processed.append(template_name)
                    # insert gap everywhere else
                    non_processed = [name for name in self.alignment if name
                                     not in processed]
                    for name in non_processed:
                        self.alignment[name].append('-')
                    # increment index
                    self.index += 1

                # if not, get non agreers and pop them
                elif len(selected_buffers) > 0:
                    selected_buffers = [
                        buffer for buffer in selected_buffers
                        if buffer.present_next_master() != next_position]
                    processed = list()
                    # insert in selected
                    for buffer in selected_buffers:
                        template_name = buffer.slave_name
                        self.alignment[template_name].append(buffer.pop_next())
                        processed.append(template_name)
                    # insert gap everywhere else
                    non_processed = [name for name in self.alignment if name
                                     not in processed]
                    for name in non_processed:
                        self.alignment[name].append('-')

                # if all buffers are finished, but there is still query
                # sequence to process
                else:
                    self.alignment[self.master_name].append(next_position)
                    self.index += 1
                    non_processed = [name for name in self.alignment
                                     if name != self.master_name]
                    for name in non_processed:
                        self.alignment[name].append('-')

            def construct_aln(self):
                while not (
                        all([buffer.is_empty() for buffer in self.aln_buffers])
                        and self.index + 1 >= len(self.master_seq)):
                    self.next_col()

                out = Alignment(None)
                out.sequences = {name: Sequence(name, ''.join(sequence))
                                 for name, sequence in self.alignment.items()}
                out.remove_redundant_gaps()
                return out

        class AlnBuffer():
            '''
            Helper class for assembling MSA-like objects from multiple pairwise
            sequence alignments.

            Flexible data structure for removing elements from the first
            position of a list.
            '''
            def __init__(self, master_name, slave_name, alignment):
                self.master_name = master_name
                self.slave_name = slave_name

                self.master_seq = alignment.sequences[master_name].sequence
                self.slave_seq = alignment.sequences[slave_name].sequence

                # index at which to start considering this sequence
                self.master_index = 0

            def present_next_master(self):
                return self.master_seq[0]

            def pop_next(self):
                next_slave = self.slave_seq[0]

                self.master_seq = self.master_seq[1:]
                self.slave_seq = self.slave_seq[1:]

                return next_slave

            def is_empty(self):
                return len(self.master_seq) == 0

            def set_index(self, master_index):
                self.master_index = master_index

        def parse_hhr(hhr_file, evalue_cutoff=0.001):
            '''
            Extract information from HHR file into a list of AlnBuffer Objects.
            '''
            target = self.target

            with open(hhr_file, 'r') as file_handler:
                content = file_handler.read()

            buffers = list()
            for block in content.split('\n>')[1:]:
                evalue = float(re.search(r'(?<=E-value=)\S+', block).group())

                if evalue < evalue_cutoff:
                    template_name = re.search(r'\w+', block).group()
                    query_seq = ''.join([
                        re.split(r'\s+', line)[-1] for line in
                        re.findall(fr'Q\s+{target}\s+\d+\s+[\w-]+', block)])
                    template_seq = ''.join([
                        re.split(r'\s+', line)[-1] for line in
                        re.findall(fr'T\s+{template_name}\s+\d+\s+[\w-]+',
                                   block)])
                    # make sure template name is only 6 characters long!
                    template_name = template_name[0:6]

                    aln = Alignment(None)
                    aln.sequences = {
                        target: Sequence(target, query_seq),
                        template_name: Sequence(template_name, template_seq)
                    }

                    buffers.append(AlnBuffer(target, template_name, aln))

            if verbose:
                n_hits = len(content.split('\n>')[1:])
                vprint(
                    f'Identified {n_hits} potential templates.\n'
                    f'Applying E-value cutoff of {evalue_cutoff}...\n'
                    f'{len(buffers)} potential templates remaining\n')

            return buffers

        # create query file
        query_file = os.path.realpath(f'{self.target}.fa')
        with open(query_file, 'w') as file_handler:
            file_handler.write(f'>{self.target}\n')
            file_handler.write(f'{self.target_seq}')

        # Uniref enrichment
        if use_uniref:
            vprint('UniRef prealignment...')
            vprint('This might take some time...')
            database_uniref30 = os.path.join(database_dir, 'UniRef30')
            a3m_file = os.path.realpath(f'{self.target}_query.a3m')
            command = [
                'hhblits', '-i', query_file, '-d', database_uniref30, '-oa3m',
                a3m_file, '-n', '1', '-cpu', str(n_threads), '-v', '1']
            subprocess.run(command, stdout=None, check=True, shell=False)
            # replace query file with output alignment
            shutil.move(a3m_file, query_file)
            vprint('UniRef prealignment completed!\n')

        # run hhblits
        database_pdb70 = os.path.join(database_dir, 'pdb70')
        hhr_file = os.path.realpath(f'{self.target}.hhr')

        vprint('Performing PDB database search...')
        command = [
            'hhblits', '-i', query_file, '-d', database_pdb70, '-o', hhr_file,
            '-n', str(iterations), '-cpu', str(n_threads), '-neffmax',
            str(neffmax), '-v', '1']
        subprocess.run(command, stdout=None, check=True, shell=False)
        vprint('PDB database search completed!')

        # parse pairwise alignments from hhr file
        vprint('Parse results...')
        buffers = parse_hhr(hhr_file, evalue_cutoff)

        # assemble multiple alignment from pairwise alignments
        vprint('Assemble combined alignment...')
        assembler = Assembler(self.target_seq, self.target, buffers)
        self.alignment = assembler.construct_aln()
        vprint('Alignment assembled!')

        # update state
        self.state['has_alignment'] = True

        # remove temporary files
        for file_name in [query_file, hhr_file]:
            if os.path.exists(file_name):
                os.remove(file_name)


class AlignmentGenerator_from_aln(AlignmentGenerator):
    '''
    Reads an alignment from file into the AlignmentGenerator workflow.

    Parameters
    ----------
    alignment_file : str
        The file to read the alignment from.
    target : str
        The name of the target sequence in the alignment.
    template_location : str
        Directory where processed templates will be stored (default
        './templates/').
    file_format : str, optional
        The format of the alignment file. Can be 'fasta' or 'pir' (default
        'fasta').

    Attributes
    ----------
    alignment : Alignment
        The alignment.
    target_seq : str
        The target sequence.
    target : str
        The name of the target sequence.
    template_location : str
        Directory where processed templates will be stored.
    state : dict
        Dictionary describing the state of the AlignmentGenerator object.

    Returns
    -------
    None

    Notes
    -----
    Useful for making use of the PDB download and processing functions that
    come with the AlignmentGenerator classes.
    '''
    def __init__(self, alignment_file: str, target: str,
                 template_location: str = './templates/',
                 file_format: str = 'fasta') -> None:
        self.alignment = Alignment(alignment_file, file_format)
        self.target = target
        self.target_seq = (self.alignment.sequences[target].sequence
                           .replace('-', ''))
        self.template_location = template_location
        self.state = {
            'has_alignment': True,
            'is_processed': False,
            }

    def get_suggestion(self):
        '''
        Not implemented, since alignment is read from file on initialization.

        Raises
        ------
        NotImplementedError
        '''
        msg = ''
        raise NotImplementedError(msg)

    def from_fasta(self, *args, **kwargs):
        '''
        Not implemented, since alignment is read from file on initialization.


        Raises
        ------
        NotImplementedError
        '''
        msg = ''
        raise NotImplementedError(msg)

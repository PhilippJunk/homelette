'''
``homelette.alignment``
=======================

The :mod:`homelette.alignment` submodule contains a selection of tools for
handling sequences and alignments.

Tutorials
---------

Basic handing of alignments with `homelette` is demonstrated in :ref:`Tutorial
1 </Tutorial1_Basics.ipynb>`. The assembling of alignments for complex
modelling is discussed in
:ref:`Tutorial 6 </Tutorial6_ComplexModelling.ipynb>`.

Functions and classes
---------------------

Functions and classes present in `homelette.alignment` are listed below:

    :class:`Alignment`
    :class:`Sequence`
    :func:`assemble_complex_aln`

-----

'''

__all__ = ['Alignment', 'Sequence', 'assemble_complex_aln']

# Standard library imports
import abc
import contextlib
import itertools
import json
import os.path
import re
import subprocess
import time
import typing
import urllib.error
import urllib.parse
import urllib.request
import warnings

#  Third party imports
import pandas as pd

# Local application imports
from . import pdb_io


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
        with open(file_name, 'r') as f:
            lines = f.readlines()

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
        with open(file_name, 'r') as f:
            lines = f.readlines()

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
        with open(file_name, 'w') as f:
            for sequence_name, sequence in self.sequences.items():
                # write name and annotation
                f.write('>P1;{}\n'.format(sequence_name))
                f.write('{}\n'.format(sequence.get_annotation_pir()))
                # write sequence with newline every nth character
                f.write(re.sub(''.join(['(.{', str(line_wrap), '})']), '\\1\n',
                        sequence.sequence, 0, re.DOTALL))
                f.write('\n*\n')

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
        with open(file_name, 'w') as f:
            for sequence_name, sequence in self.sequences.items():
                # write name
                f.write('>{}\n'.format(sequence_name))
                # write sequence with newline every nth character
                f.write(re.sub(''.join(['(.{', str(line_wrap), '})']), '\\1\n',
                        sequence.sequence, 0, re.DOTALL))
                f.write('\n')

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
                          self.sequences.keys()]
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
        with open(file_name, 'w') as f:
            with contextlib.redirect_stdout(f):
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

    # TODO think about name:
    # maybe something like 'remove_missing_res'?
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
        '''  # TODO examples
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
        for sequence_name_2 in self.sequences.keys():
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
        for sequence_name_2 in self.sequences.keys():
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
    based on sequence input

    Parameters
    ----------
    sequence : str
        Target sequence in 1 letter amino acid code
    '''
    def __init__(self, sequence: str, target: str = 'target',
                 template_location: str = './templates'):
        self.alignment = None
        self.target_seq = sequence
        self.target = 'target'
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
        pass

    def _check_state(self, has_alignment: bool = None, is_processed: bool =
                     None) -> None:
        '''
        Check state of the AlignmentGenerator object.

        Parameters
        ----------
        has_alignment: bool
            Assesses whether an alignment has already been generated.
        is_processed: bool
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
            if required is None:
                return True
            if required is state:
                return True
            if required is not state:
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

    def show_suggestion(self) -> typing.Type['pd.DataFrame']:
        '''
        Shows which templates have been suggested by the AlignmentGenerator, as
        well as some useful statistics (sequence identity, coverage).

        Returns
        -------
        suggestion : pd.DataFrame
            DataFrame with calculates sequence identity and sequence coverage
            for target

        Raises
        ------
        RuntimeError
            Alignment has not been generated yet
        '''
        self._check_state(has_alignment=True, is_processed=None)

        df_coverage = self.alignment.calc_coverage_target(self.target)
        df_identity = self.alignment.calc_identity_target(self.target)
        # TODO maybe get method for PDB structure? and resolution?
        # but that would require multiple web requests, so not really ideal..
        # TODO maybe propose ranking? (borda) or just sort by seq_id?

        output = (
            # combine data frame
            pd.merge(df_coverage, df_identity, on=('sequence_1', 'sequence_2'))
            # sort values
            .sort_values(by=['identity', 'coverage'], ascending=[False, False])
            # remove column with sequence_1 and rename sequence_2
            .drop('sequence_1', axis=1)
            .rename({'sequence_2': 'template'}, axis=1)
            )
        return output

    def select_templates(self, templates: typing.Iterable) -> None:
        '''
        Select templates from suggested templates by identifier.

        Raises
        ------
        RuntimeError
            Alignment has not been generated yet
        '''
        self._check_state(has_alignment=True, is_processed=False)
        selection = ['target'] + list(templates)
        self.alignment.select_sequences(selection)

    # TODO name of the function? download_templates?
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


        Notes
        -----
        pdb_format tells the function how to parse the template identifiers in
        the alignment:

        * auto: Automatic guess for pdb_format
        * chain: Sequences are named in the format PDBID_CHAIN (i.e. 4G0N_A)
        * entity: Sequences are named in the format PDBID_ENTITY (i.e. 4G0N_1)
        * identifier: Sequences are named only be their PDB identifier (i.e.
        4G0N)
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

        def guess_pdb_format(templates):
            '''
            Guess PDB identifier format from template names
            '''
            r_identifier = re.compile(r'^[A-Za-z0-9]{4}')
            r_chain = re.compile(r'^[A-Za-z0-9]{4}[\W_][A-Za-z]')
            r_entity = re.compile(r'^[A-Za-z0-9]{4}[\W_][0-9]')
            if all([re.fullmatch(r_identifier, template) for template in
                   templates]):
                pdb_format = 'identifier'
            elif all([re.fullmatch(r_chain, template) for template in
                     templates]):
                pdb_format = 'chain'
            elif all([re.fullmatch(r_entity, template) for template in
                     templates]):
                pdb_format = 'entity'
            else:
                raise ValueError(
                    'Unable to guess pdb_format from template names. Please '
                    'set pdb_format manually and make sure all template names '
                    'in the alignment follow one of the proposed naming '
                    'schemes.')
            return pdb_format

        def parse_alignment(alignment: Alignment, templates: typing.Iterable,
                            pdb_format) -> dict:
            '''
            '''  # TODO
            out = {}
            for template in templates:
                # get PDBID
                pdbid = template[0:4]

                # include in out dict
                if pdbid not in out.keys():
                    out[pdbid] = list()

                # collect information about alignment entry
                aln_entry = {
                        'template': template,
                        'sequence': alignment.sequences[template],
                        }
                if pdb_format == 'chain':
                    aln_entry['chain'] = template[5]

                out[pdbid].append(aln_entry)

            return out

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

        def process_template_pdb(pdb, pdb_file_name, template_location):
            '''
            Remove HOH, renumber residues, rename chain ID, save to file.
            '''
            pdb = (
                pdb
                .transform_filter_res_name(['HOH'])
                .transform_renumber_residues(starting_res=1)
                .transform_change_chain_id(new_chain_id='A'))
            pdb.write_pdb(os.path.join(
                template_location, pdb_file_name))

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

        # check pdb_format
        templates = [template for template in self.alignment.sequences.keys()
                     if template != self.target]
        if pdb_format == 'auto':
            vprint('Guessing template naming format...')
            pdb_format = guess_pdb_format(templates)
            vprint(f'Template naming format guessed: {pdb_format}!\n')
        elif pdb_format not in ['chain', 'identifier', 'entity']:
            raise ValueError(
                f'Invalid value to pdb_format: {pdb_format}\nHas to be one of '
                f'"auto", "chain", "entity", "identifier".')

        # Initialize template dir
        vprint('Checking template dir...')
        if not os.path.exists(self.template_location):
            vprint('Template dir not found...')
            os.makedirs(self.template_location, exist_ok=True)
            vprint(f'New template dir created at\n'
                   f'"{self.template_location}"!\n')
        else:
            vprint('Template dir found!\n')

        # parse templates
        templates_parsed = parse_alignment(self.alignment, templates,
                                           pdb_format)
        vprint('Processing templates:\n')

        # initialize changed alignment
        new_aln = Alignment(None)
        new_aln.sequences = {
                self.target: self.alignment.sequences[self.target]}

        # iterate over PDB identifier
        for pdbid in templates_parsed:
            # download template
            vprint(f'{pdbid} downloading from PDB...')
            pdb = pdb_io.download_pdb(pdbid)
            vprint(f'{pdbid} downloaded!')

            # different behaviour whether we know which chains to extract
            # (pdb_format == 'chain'), or need to check which chains fit the
            # entity (pdb_format != 'chain')
            if pdb_format == 'chain':
                for chain_info in templates_parsed[pdbid]:
                    template_name = chain_info['template']
                    chain = chain_info['chain']
                    seq_alignment = chain_info['sequence'].sequence.replace(
                            '-', '')
                    pdb_chain = pdb.transform_extract_chain(chain)
                    vprint(f'{pdbid}_{chain}: Chain extracted!')
                    seq_template = pdb_chain.get_sequence(ignore_missing=False)
                    try:
                        seq_template_padded = adjust_template_seq(
                                seq_alignment, seq_template)
                    except RuntimeError as e:
                        raise RuntimeError(f'Template: {pdbid}_{chain}') from e

                    new_aln = add_seq_to_aln(
                        new_aln, f'{pdbid}_{chain}',
                        self.alignment.sequences[template_name].sequence,
                        seq_template_padded)
                    vprint(f'{pdbid}_{chain}: Alignment updated!')

                    process_template_pdb(pdb_chain, f'{pdbid}_{chain}.pdb',
                                         self.template_location)
                    vprint(f'{pdbid}_{chain}: PDB processed!')

            # if template were not given with chain identifiers, iterate over
            # all combinations of entities and chains and try to match them
            else:
                for entity_info, chain in itertools.product(
                        templates_parsed[pdbid], pdb.get_chains()):
                    template_name = entity_info['template']
                    seq_alignment = entity_info['sequence'].sequence.replace(
                            '-', '')
                    pdb_chain = pdb.transform_extract_chain(chain)
                    seq_template = pdb_chain.get_sequence(
                            ignore_missing=False)
                    try:
                        seq_template_padded = adjust_template_seq(
                                seq_alignment, seq_template)
                    except RuntimeError:
                        # if not matching, dont continue loop
                        continue

                    vprint(f'{pdbid}_{chain}: Chain extracted!')

                    # update alignment
                    new_aln = add_seq_to_aln(
                        new_aln, f'{pdbid}_{chain}',
                        self.alignment.sequences[template_name].sequence,
                        seq_template_padded)
                    vprint(f'{pdbid}_{chain}: Alignment updated!')

                    # process template pdb
                    process_template_pdb(pdb_chain, f'{pdbid}_{chain}.pdb',
                                         self.template_location)
                    vprint(f'{pdbid}_{chain}: PDB processed!')

        # update alignment
        self.alignment = new_aln
        # update state
        self.state['is_processed'] = True
        vprint(f'\nFinishing... All templates successfully\ndownloaded and '
               f'processed!\nTemplates can be found in\n'
               f'"{self.template_location}".')


# TODO remove after testing
class TestAlignmentGenerator(AlignmentGenerator):
    '''
    TEST ONLY
    '''
    def get_suggestion(self):
        # check state
        self._check_state(has_alignment=False, is_processed=False)
        # For testing purpose, just retrieve sequence alignment from
        # examples/data/single/aln_2.fasta_aln
        self.alignment = Alignment('/home/junkpp/work/programs/homelette/'
                                   'test_alngen.fasta_aln')
        self.alignment.rename_sequence('ARAF', 'target')
        self.target_seq = self.alignment.sequences['target'].sequence

        # update state
        self.state['has_alignment'] = True


class AlignmentGenerator_pdb(AlignmentGenerator):
    '''
    Auto-generation of alignments based on a pdbblast search for a target
    sequence

    Parameters
    ----------
    sequence : str
        Target sequence
    '''
    def get_suggestion(self, seq_id_cutoff: float = 0.5, min_length: int = 30,
                       verbose=True) -> None:
        '''
        '''  # TODO
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
            # TODO consider exp. method? maybe xray_only?
            # documentation for query structures can be found at
            # https://search.rcsb.org/index.html
            print(f'{sequence}, {seq_id_cutoff}, {min_length}, {max_results}')
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
                  {{
                    "type": "terminal",
                    "service": "text",
                    "parameters": {{
                      "attribute": "exptl.method",
                      "operator": "exact_match",
                      "negation": false,
                      "value": "X-RAY DIFFRACTION"
                    }}
                  }},
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
                "return_all_hits": true,
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
            # download fastas for entities
            # TODO potentially outsource to sequence_collection class?
            sequences = dict()
            vprint('Retrieving sequences...')
            # set up status updates during downloads
            for template in templates:
                url = (
                    f'https://rcsb.org/fasta/entity/{template}/'
                    f'download')
                with urllib.request.urlopen(url) as response:
                    fasta_file = response.read().decode('utf-8')
                sequences[template] = fasta_file.split('\n')[1]
                # give short delay between requests
                time.sleep(0.5)
            vprint('Sequences succefully retrieved!\n')

            # write to output file
            fasta_file_name = os.path.realpath(
                    f'.sequences_{self.target}.fa')
            with open(fasta_file_name, 'w') as file_handle:
                file_handle.write(f'>{self.target}\n')
                file_handle.write(f'{self.target_seq}\n')
                for template, seq in sequences.items():
                    file_handle.write(f'>{template}\n')
                    file_handle.write(f'{seq}\n')

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
            for file_name in [fasta_file_name]:
                if os.path.exists(file_name):
                    os.remove(file_name)

            return aln

        # send query
        vprint('Querying PDB...')
        templates = query_pdb(self.target_seq, seq_id_cutoff, min_length)
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
        vprint('Query successful.\nThe following sequences have been found')
        vprint(self.show_suggestion())


# NOTES
# github https://github.com/soedinglab/hh-suite
# databases http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs
class AlignmentGenerator_hhblits(AlignmentGenerator):
    '''
    '''
    def get_suggestion(self):
        '''
        '''  # TODO
        # check state
        self._check_state(has_alignment=False, is_processed=False)
        # TODO
        # update state
        self.state['has_alignment'] = True
        pass


class AlignmentGenerator_hmmer(AlignmentGenerator):
    '''
    '''
    # NOTES
    # github https://github.com/EddyRivasLab/hmmer
    # documentation https://github.com/EddyRivasLab/hmmer
    # not sure if I can just used a downloaded PDB sequence database?
    # if so, these would contain ID_CHAIN nomenclature
    # is there a HMMER docker container? try out?

    def get_suggestion(self):
        '''
        '''  # TODO
        # check state
        self._check_state(has_alignment=False, is_processed=False)
        # TODO
        # update state
        self.state['has_alignment'] = True
        pass
# TODO maybe implement dependency check for non-python packages?

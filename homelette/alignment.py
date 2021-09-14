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
import contextlib
import itertools
import re
import typing
import warnings

#  Third party imports
import pandas as pd


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
        sequences matters

        >>> aln = hm.Alignment(None)
        >>> aln.sequences = {
        ...     'seq1': hm.alignment.Sequence('seq1', 'AAAACCCCDDDD'),
        ...     'seq3': hm.alignment.Sequence('seq3', 'AAAA----DDDD')
        ...     }
        >>> aln.calc_identity('seq1', 'seq2')
        66.67
        >>> aln.calc_identity('seq3', 'seq1')
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
            output['sequence_1'].append(sequence_name)
            output['sequence_2'].append(sequence_name_2)
            output['identity'].append(self.calc_identity(
                sequence_name, sequence_name_2))
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

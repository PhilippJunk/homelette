'''
``homelette.pdb_io``
====================

The :mod:`homelette.pdb_io` submodule contains an object for parsing and
manipulating PDB files. There are several constructor function that can read
PDB files or download them from the internet.


Functions and classes
---------------------

Functions and classes present in `homelette.pdb_io` are listed below:

    :class:`PdbObject`
    :func:`read_pdb`
    :func:`download_pdb`

-----

'''

__all__ = ['read_pdb', 'download_pdb', 'PdbObject']

# Standard library imports
import gzip
import itertools
import typing
import urllib.request
import urllib.error
import warnings

# Third party imports
import pandas as pd


class PdbObject:
    '''
    Object encapsulating functionality regarding the processing of PDB files

    Parameters
    ----------
    lines : Iterable
        The lines of the PDB

    Attributes
    ----------
    lines
        The lines of the PDB, filtered for ATOM and HETATM records

    Returns
    -------
    None

    See Also
    --------
    read_pdb
    download_pdb

    Notes
    -----
    Please contruct instances of PdbObject using the constructor functions.

    If a PDB file with multiple MODELs is read, only the first model will be
    conserved.
    '''
    def __init__(self, lines: typing.Iterable) -> None:
        # filter for ATOM, HETATM, MODEL and TER records
        lines = [line for line in lines if
                 line.startswith('ATOM') or
                 line.startswith('HETATM') or
                 line.startswith('TER') or
                 line.startswith('MODEL')]
        # if multiple models are present, choose the first one
        models = 0
        new_lines = list()
        for line in lines:
            if line.startswith('MODEL'):
                models += 1
                if models > 1:
                    break
            else:
                new_lines.append(line)

        self.lines = new_lines

    def write_pdb(self, file_name) -> None:
        '''
        Write PDB to file.

        Parameters
        ----------
        file_name : str
            The name of the file to write the PDB to.

        Returns
        -------
        None
        '''
        with open(file_name, 'w') as file_handler:
            file_handler.writelines(self.lines)

    def parse_to_pd(self) -> pd.DataFrame:
        '''
        Parses PDB to pandas dataframe.

        Returns
        -------
        pd.DataFrame

        Notes
        -----
        Information is extracted according to the PDB file specification
        (version 3.30) and columns are named accordingly. See
        https://www.wwpdb.org/documentation/file-format for more information.
        '''
        # parse
        out = []
        for line in [line for line in self.lines if not
                     line.startswith('TER')]:
            out.append({
                'record': line[0:6].strip(),
                'serial': int(line[6:11]),
                'name': line[12:16].strip(),
                'altLoc': line[16].strip(),
                'resName': line[17:20].strip(),
                'chainID': line[21].strip(),
                'resSeq': int(line[22:26]),
                'iCode': line[26].strip(),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54]),
                'occupancy': float(line[54:60]),
                'tempFactor': float(line[60:66]),
                'element': line[76:78].strip(),
                'charge': line[78:80].strip(),
            })
        # concat to pd.DataFrame
        return pd.DataFrame(out)

    def get_sequence(self, ignore_missing: bool = True) -> str:
        '''
        Retrieve the 1-letter amino acid sequence of the PDB, grouped by
        chain.

        Parameters
        ----------
        ignore_missing : bool
            Changes behaviour with regards to unmodelled residues. If True,
            they will be ignored for generating the sequence (default). If
            False, they will be represented in the sequence with the character
            X.

        Returns
        -------
        str
            Amino acid sequence
        '''
        def _321(aminoacid):
            '''
            Transform 3 letter amino acid code to 1 letter code
            '''
            aa_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                       'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                       'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                       'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                       'SEC': 'U'}
            return aminoacid.map(aa_code)

        pdb_df = self.parse_to_pd()
        # extract residues from pdb df
        residues = (
            pdb_df[pdb_df.record.eq('ATOM')][['chainID', 'resSeq', 'resName']]
            .groupby(['chainID', 'resSeq', 'resName'])
            .count()
            .reset_index()
        )
        # transform residues from 3 letter to 1 letter code
        residues = (
            residues.assign(resName=_321(residues['resName']))
        )
        # check for multiple chains
        if len(residues.chainID.unique()) > 1:
            warnings.warn(
                'Multiple chains present. Chains will be separated by "/".')
        # handle missing residues
        if not ignore_missing:
            for chainID in residues.chainID.unique():
                # retrieve min and max residues
                min_resSeq = int(
                        residues[residues['chainID'] == chainID]
                        ['resSeq']
                        .min())
                max_resSeq = int(
                        residues[residues['chainID'] == chainID]
                        ['resSeq']
                        .max())
                # create new dataframe with complete records from min to max
                new_df = pd.DataFrame({
                    'chainID': chainID,
                    'resSeq': list(range(min_resSeq, max_resSeq + 1))
                    })
                # merge new with old data frame, replace missing values with X
                # adapted from https://stackoverflow.com/a/54906079/7912251
                residues = (
                    residues.merge(new_df, on=('chainID', 'resSeq'),
                                   how='outer', sort=True)
                    .fillna('X'))
        # extract chains on their own, then concat with / as separator
        sequences = (
            [''.join(residues[residues['chainID'] == chainID]['resName']
                     .tolist())
             for chainID in residues.chainID.unique()])
        return '/'.join(sequences).upper()

    def get_chains(self) -> list:
        '''
        Extract all chains present in the PDB.

        Returns
        -------
        list
        '''
        return sorted(list({line[21] for line in self.lines}))

    def transform_extract_chain(self, chain) -> 'PdbObject':
        '''
        Extract chain from PDB.

        Parameters
        ----------
        chain : str
            The chain ID to be extracted.

        Returns
        -------
        PdbObject
        '''
        filtered_lines = [line for line in self.lines if line[21] == chain]
        return PdbObject(filtered_lines)

    def transform_renumber_residues(
            self, starting_res: int = 1) -> 'PdbObject':
        '''
        Renumber residues in PDB.

        Parameters
        ----------
        starting_res : int
            Residue number to start renumbering at (default 1)

        Returns
        -------
        PdbObject

        Notes
        -----
        Missing residues in the PDB (i.e. unmodelled) will not be considered in
        the renumbering. If multiple chains are present in the PDB, numbering
        will be continued from one chain to the next one.
        '''
        #
        curr_resSeq, prev_resSeq = str(), str()
        index = int(starting_res) - 1
        output_lines = []
        for line in self.lines:
            curr_resSeq = line[22:26]
            if curr_resSeq != prev_resSeq:
                index += 1
                if len(str(index)) > 4:
                    raise RuntimeError(
                        'Digits of residues exceed available columns (4) in '
                        'the PDB file format. Please choose a smaller number.')
            output_lines.append(line[:22] + str(index).rjust(4) + line[26:])
            prev_resSeq = curr_resSeq

        return PdbObject(output_lines)

    def transform_change_chain_id(self, new_chain_id) -> 'PdbObject':
        '''
        Replace chain ID for every entry in PDB.

        Parameters
        ----------
        new_chain_id : str
            New chain ID.

        Returns
        -------
        PdbObject
        '''
        changed_lines = [(line[:21] + new_chain_id + line[22:])
                         for line in self.lines]
        return PdbObject(changed_lines)

    def transform_remove_hetatm(self) -> 'PdbObject':
        '''
        Remove all HETATM entries from PDB.

        Returns
        -------
        PdbObject
        '''
        filtered_lines = [line for line in self.lines if not
                          line.startswith('HETATM')]
        return PdbObject(filtered_lines)

    def transform_filter_res_name(
            self, selection: typing.Iterable,
            mode: str = 'out') -> 'PdbObject':
        '''
        Filter PDB by residue name.

        Parameters
        ----------
        selection : Iterable
            For which residue names to filter
        mode : str
            Filtering mode. If mode = "out", the selection will be filtered out
            (default). If mode = "in", everything except the selection will be
            filtered out.

        Returns
        -------
        PdbObject
        '''
        # implement different behaviour based on mode
        if mode == 'out':
            def check_line(line):
                return line[17:20].strip() not in selection
        elif mode == 'in':
            def check_line(line):
                return line[17:20].strip() in selection
        else:
            raise ValueError('Mode has to be "out" or "in". "' + mode +
                             '" is not an accepted value.')
        filtered_lines = [line for line in self.lines if check_line(line)]
        return PdbObject(filtered_lines)

    def transform_filter_res_seq(self, lower: int, upper: int) -> 'PdbObject':
        '''
        Filter PDB by residue number.

        Parameters
        ----------
        lower : int
            Lower bound of range to filter with.
        upper : int
            Upper bound of range to filter with, inclusive.

        Returns
        -------
        PdbObject
        '''
        filtered_lines = [line for line in self.lines
                          if lower <= int(line[22:26]) <= upper]
        return PdbObject(filtered_lines)

    def transform_concat(self, *others: 'PdbObject') -> 'PdbObject':
        '''
        Concat PDB with other PDBs.

        Parameters
        ----------
        *others : 'PdbObject
            Any number of PDBs.

        Returns
        -------
        PdbObject
        '''
        list_of_lines = [self.lines] + [pdb.lines for pdb in others]
        return PdbObject(itertools.chain(*list_of_lines))


# Constructor functions for PdbObject
def read_pdb(file_name: str) -> PdbObject:
    '''
    Reads PDB from file.

    Parameters
    ----------
    file_name : str
        PDB file name

    Returns
    -------
    PdbObject

    Notes
    -----
    If a PDB file with multiple MODELs is read, only the first model will be
    conserved.
    '''
    with open(file_name, 'r') as file_handle:
        return PdbObject(file_handle.read().splitlines(keepends=True))


def download_pdb(pdbid: str) -> PdbObject:
    '''
    Download PDB from the RCSB.

    Parameters
    ----------
    pdbid : str
        PDB identifier

    Returns
    -------
    PdbObject

    Notes
    -----
    If a PDB file with multiple MODELs is read, only the first model will be
    conserved.
    '''
    # adapted from https://stackoverflow.com/a/7244263/7912251
    url = 'https://files.rcsb.org/download/' + pdbid + '.pdb.gz'
    with urllib.request.urlopen(url) as response:
        if response.status == 200:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                pdb = uncompressed.read().decode('utf-8')
        else:
            raise urllib.error.URLError()
    return PdbObject(pdb.splitlines(keepends=True))

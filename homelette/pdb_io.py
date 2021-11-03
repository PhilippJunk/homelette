'''
``homelette.pdb_io``
====================

The :mod:`homelette.pdb_io` submodule contains an object for parsing and
manipulating PDB files. There are several constructor function that can read
PDB files or download them from the internet.


Tutorials
---------
TODO?

Functions and classes
---------------------

Functions and classes present in `homelette.pdb_io` are listed below:

    :class:`PdbObject`
    :func:`read_pdb`
    :func:`download_pdb`

-----

'''  # TODO

__all__ = ['read_pdb', 'download_pdb', 'PdbObject']

# Standard library imports
import gzip
import typing
import urllib.request

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

    See Also
    --------
    read_pdb
    download_pdb

    Notes
    -----
    Please contruct instances of PdbObject using the constructor functions.

    Information is extracted according to the PDB file specification (version
    3.30) and columns are named accordingly. See
    https://www.wwpdb.org/documentation/file-format for more information.
    '''  # TODO
    def __init__(self, lines: typing.Iterable) -> None:
        # filter for ATOM and HETATM records
        self.lines = [line for line in lines if line.startswith('ATOM') or
                      line.startswith('HETATM')]

    def write_pdb(self, file_name) -> None:
        '''
        Write PDB to file.

        Parameters
        ----------

        file_name : str
            The name of the file to write the PDB to.
        '''
        with open(file_name, 'w') as file_handler:
            file_handler.writelines(self.lines)

    def parse_to_pd(self) -> pd.DataFrame:
        '''
        Parses PDB to pandas dataframe.

        Returns
        -------
        pd.DataFrame
        '''
        # parse
        out = []
        for line in self.lines:
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

    def get_sequence(self) -> str:
        '''
        Retrieve the 1-letter amino acid sequence of the PDB, grouped by
        chain.

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

        # TODO deal with multiple chains
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
        return ''.join(residues['resName'].tolist()).upper()

    def transform_extract_chain(self, chain) -> 'PdbObject':
        '''
        '''  # TODO
        pass

    def transform_renumber_residues(self, starting_res) -> 'PdbObject':
        '''
        '''  # TODO
        pass

    def transform_change_chain_id(self, new_chain_id) -> 'PdbObject':
        '''
        '''  # TODO
        pass


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
    '''
    with open(file_name, 'r') as file_handle:
        return PdbObject(file_handle.readlines())


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
    '''
    # adapted from https://stackoverflow.com/a/7244263/7912251
    url = 'https://files.rcsb.org/download/' + pdbid + '.pdb.gz'
    with urllib.request.urlopen(url) as response:
        with gzip.GzipFile(fileobj=response) as uncompressed:
            pdb = uncompressed.read().decode('utf-8')
    return PdbObject(pdb.splitlines())

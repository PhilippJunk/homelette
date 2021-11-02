'''
'''
import gzip
import urllib.request

import pandas as pd

class PdbObject:
    # TODO change name to camel case
    '''
    Object encapsulating functionality regarding PDB files and transformations

    Parameters
    ----------

    Attributes
    ----------
    '''  # TODO
    def __init__(self, lines):
        # TODO maybe filter for ATOM and HETATM records only?
        self.lines = [line for line in lines if line.startswith('ATOM') or
                line.startswith('HETATM')]

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
        Retrieve the 1-letter amino acid sequence of the PDB, grouped by chains.

        Returns
        -------
        str
            Amino acid sequence
        '''
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

        # TODO deal with multiple chains?
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
        pass

    def transform_renumber_residues(self, starting_res) -> 'PdbObject':
        pass

    def transform_change_chain_id(self, new_chain_id) -> 'PdbObject':
        pass


# Constructor functions for PdbObject
def read_pdb(file_name) -> PdbObject:
    '''
    Reads PDB from file.
    '''  # TODO
    with open(file_name, 'r') as file_handle:
        return PdbObject(file_handle.readlines())

def download_pdb(pdbid) -> PdbObject:
    '''
    Download PDB from the RCSB.
    '''
    url = 'https://files.rcsb.org/download/' + pdbid + '.pdb.gz'
    with urllib.request.urlopen(url) as response:
                with gzip.GzipFile(fileobj=response) as uncompressed:
                    pdb = uncompressed.read().decode('utf-8')
    return PdbObject(pdb.splitlines())

from Bio.PDB.MMCIFParser import MMCIFParser
import sys
from Bio.PDB import Selection, NeighborSearch, PDBParser, PDBIO, Select
import os
import collections
import json

dict_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
           'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

class SelectChains(Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

def write_sel_chain_pdb(pdbfile,
                        chain_id,
                        overwrite=False):
    writer = PDBIO()
    struct = str_object(pdbfile)
    writer.set_structure(struct)
    outfle = (os.getcwd() 
                + '/' 
                + pdbfile.rsplit('.')[0]
                + '_'
                + ''.join(chain_id) 
                + '.pdb')
    if (not overwrite) and (os.path.isfile(outfle)):
        print("Chain %s of '%s' already extracted to '%s'." %
                (", ".join(chain_id), pdbfile, outfle))
        return outfle
    print("Extracting chain %s from %s..." % (", ".join(chain_id), pdbfile))
    writer.save(outfle, select=SelectChains(chain_id))

def str_object(pdbfile):
    '''
    Given a pdb file creates a structure object
    '''
    if pdbfile.rsplit('.')[-1] == 'pdb':
        parser = PDBParser(QUIET=True)
    elif pdbfile.rsplit('.')[-1] == 'ent':
        parser = PDBParser(QUIET=True)
    elif pdbfile.rsplit('.')[0] == 'cif':
        parser = MMCIFParser()
    else:
        print 'Not a valid PDB format'
        sys.exit()
    structure = parser.get_structure('id', pdbfile)
    return structure

def get_num_chains(structure):
    '''
    Returns a list of chain objects in a structure
    '''
    ch_list = Selection.unfold_entities(structure, 'C')
    return ch_list

def get_CAatoms(chain_object):
    '''
    for a given chain object, returns a list of CA atoms
    '''
    ca_list = []
    atm_list = Selection.unfold_entities(chain_object, 'A')
    for a in atm_list:
        if a.get_id() == 'CA':
            ca_list.append(a)
    return ca_list
    
def interface_residues(pdbfile,
                    len_chains,
                    numres_cut,
                    dist_cutoff):
    '''
    Given a PDB file calculates interface residues
    '''
#     numres_cut = int(numres_cut)
    structure = str_object(pdbfile)
    ch_list = get_num_chains(structure)
    int_reslist_j = []
    int_reslist_i = []
    dict_intf = {}
    dict_residue_contacts = {}
    for i in xrange(0,len(ch_list)-1):
        for j in xrange(i+1,len(ch_list)):
            int_reslist_j = []
            int_reslist_i = []
            intf_ch = ch_list[j].get_id() + '_' + ch_list[i].get_id()
            ca_i = get_CAatoms(ch_list[i])
            ca_j = get_CAatoms(ch_list[j])
            # print ca_i, intf_ch
            if len(ca_i) <len_chains or len(ca_j) <len_chains: 
                # dict_intf = {} 
                # dict_residue_contacts = {}
                continue
            ns_i = NeighborSearch(ca_i)
            for atms in ca_j:
                n = ns_i.search(atms.get_coord(), dist_cutoff,'R')
                if len(n) >0 :
                    res = (str(atms.get_parent().get_resname()) 
                             + str(atms.get_parent().get_id()[1]))
                    int_reslist_j.append(res)
                    t = []
                    for at in n:
                        res1 = (str(at.get_resname()) 
                             + str(at.get_id()[1]))
                        int_reslist_i.append(res1)
                        t.append(res1)
                    if intf_ch in dict_residue_contacts.keys():  
                        dict_residue_contacts[intf_ch].update({res:t})
                    else:
                        dict_residue_contacts[intf_ch] = {res:t}
            #Number of interface residues from each chain is >=10
            if len(set(int_reslist_j)) >= numres_cut and len(set(int_reslist_i)) >=numres_cut:
                dict_intf[intf_ch] = [list(set(int_reslist_j)),list(set(int_reslist_i))]   
    intf_outfle_prefix = pdbfile.split('.')[0]
    intf_outfle_name = intf_outfle_prefix + '_interface_residues_dict.json'
    
    if dict_intf: 
        with open(intf_outfle_name,'w') as outfle:
            json.dump(dict_intf,outfle)
    return dict_intf, dict_residue_contacts
    
    
def polar_residues(reslist):
    assert len(reslist) >0 
    polar_ref = ['SER', 
                 'THR',
                 'ASN',
                 'GLN',
                 'HIS',
                 'TYR']
    polar_reslist = []
    for r in reslist:
        if r in polar_ref:
            polar_reslist.append(r)
    frac_polar_residue = len(polar_reslist)/float(len(reslist))
    return frac_polar_residue

def hydrophobic_residues(reslist):
    assert len(reslist) >0 
    hydro_ref = ['ALA',
                 'LEU',
                 'ILE',
                 'VAL',
                 'PHE',
                 'TRP',
                 'CYS',
                 'MET'
                ]
    hydro_reslist = []
    for r in reslist:
        if r in hydro_ref:
            hydro_reslist.append(r)
    frac_hydro_residue = len(hydro_reslist)/float(len(reslist))
    return frac_hydro_residue
            
def charged_residues(reslist):
    assert len(reslist) >0 
    charged_ref = ['LYS',
                   'ARG',
                   'ASP',
                   'GLU']
    charged_reslist = []
    for r in reslist:
        if r in charged_ref:
            charged_reslist.append(r)
    frac_charged_residue = len(charged_reslist)/float(len(reslist))
    return frac_charged_residue

# def conserved_residues(reslist):

def run_freesasa(pdbfile):
    '''
    Runs FreeSASA and writes a RSA file
    '''
    cmd = ('freesasa '
            + pdbfile
            +' --format=rsa >'
            + pdbfile.rsplit('.')[0].rsa)
    os.system(cmd)
    
def parse_rsafile(rsafile):
    '''
    Reads the rsafile of freesasa and returns a dictionary
    where keys are the chain IDs and value is the list of acceesible residues
    '''
    dict_solvent_accessible = {}
    with open(rsafile) as rsa:
        for lne in rsa:
            if lne.startswith('RES'):
                fl = lne.split()
                if fl[7] >= 5:
                    if dict_solvent_accessible[fl[2]]:
                        dict_solvent_accessible[fl[2]].append(fl[1] + fl[3])
                    else:
                        dict_solvent_accessible[fl[2]] = fl[1] + fl[3]
    return dict_solvent_accessible

def accessible_interface_residues(reslist,reflist):
    '''
    reslist is list of all accessible residues calcualted using freesasa
    reflist is list of interface residues
    this function checks for is interface residue is solvent accessible    
    '''
    accessible_reslist = []
    for r in reslist:
        if r in reflist:
            accessible_reslist.append(r)
    return accessible_reslist
                    
def get_nr_lst_interface_intra_pdb():
    '''
    takes the text file as input generated after running ialign and getting the non-redundant
    list of interfaces within a PDB
    '''
    nr_intf_fle = ''
    to_check = os.listdir(os.getcwd())
    nr_intf_lst = []

    for fl in to_check:
        if 'nr_interface_lst.txt' in fl:
            nr_intf_fle = fl
            break

    if nr_intf_fle:
        with open(nr_intf_fle) as nr:
            for nri in nr:
                nri = nri.strip()
                if nri not in nr_intf_lst:
                    nr_intf_lst.append(nri[-2] + '_' + nri[-1])
    else:
        for fl in to_check:
            if 'interface_residues_dict.json' in fl:
                with open(fl) as json_file:  
                    dict_intf = json.load(json_file)
                    nr_intf_lst = dict_intf.keys()
    
    return nr_intf_lst

def measure_contact_counts(dict_contact):
    '''
    Takes the output from the interface_residues and calculate contact counts
    calculates Cij
    now the function checks for the non redundant interfaces list if present after intra pdb clusters 
    '''
    nr_intf_lst = get_nr_lst_interface_intra_pdb()
    contact_pair = collections.OrderedDict()
    for ch in dict_contact.keys():
        if nr_intf_lst:
            if ch in nr_intf_lst:
                for res in dict_contact[ch]:
                    if res[0:3] in dict_aa.keys():
                        for i in dict_contact[ch][res]:
                            pair_name = res[0:3]+ '_' + i[0:3]
                            pair_name1 = i[0:3]+ '_' + res[0:3]
                            if pair_name in contact_pair.keys() or pair_name1 in contact_pair.keys():
                                try:
                                    contact_pair[pair_name] += 1
                                except:
                                    contact_pair[pair_name1] += 1
                            else:
                                contact_pair[pair_name] = 1
        else:
            for res in dict_contact[ch]:
                    if res[0:3] in dict_aa.keys():
                        for i in dict_contact[ch][res]:
                            pair_name = res[0:3]+ '_' + i[0:3]
                            pair_name1 = i[0:3]+ '_' + res[0:3]
                            if pair_name in contact_pair.keys() or pair_name1 in contact_pair.keys():
                                try:
                                    contact_pair[pair_name] += 1
                                except:
                                    contact_pair[pair_name1] += 1
                            else:
                                contact_pair[pair_name] = 1
    #     OrderedDict([('PRO_PHE', 3), ('PRO_VAL', 1), ('PRO_TYR', 1), ('ASN_PHE', 2)])
    return contact_pair
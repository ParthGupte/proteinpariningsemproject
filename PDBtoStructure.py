from PDB import PDB
from proteinhelper import *

def slice_cord_line(line):    
    return (int(line[6:11]),line[12:16].lstrip().rstrip(),line[16:21].lstrip().rstrip(),int(line[24:30].lstrip().rstrip()),float(line[30:38]),float(line[38:46]),float(line[46:54]))

# amino_acids = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER']

def PDBtoStructure(pdb: PDB):
    i = 0
    res_D = {1:{}}
    for line in pdb.ATOM_lines:
        data = slice_cord_line(line)
        amino_num = data[3]
        if amino_num == i:
            continue
        else:
            atom = data[1]
            if {'N','CA','C','O'} == set(res_D[i+1].keys()):
                i += 1
                res_D[i+1] = {} 
                continue

            if atom == 'N':
                N_coord = data[4:]
                res_D[i+1]['N'] = N_coord
            elif atom == 'CA':
                CA_coord = data[4:]
                res_D[i+1]['CA'] = CA_coord
            elif atom == 'C':
                C_coord = data[4:]
                res_D[i+1]['C'] = C_coord
            elif atom == 'O':
                O_coord = data[4:]
                res_D[i+1]['O'] = O_coord
    
    points = []
    vectors = []
    for j in range(1,i):
        res_info = res_D[j]
        points.append(list(res_info['CA']))
        v1 = np.array(res_info['N'])- np.array(res_info['CA'])
        v2 = np.array(res_info['C']) - np.array(res_info['CA'])
        v3 = np.array(res_info['O']) - np.array(res_info['C'])
        vectors.append([v1,v2,v3])
    return Structure(np.array(points,dtype=np.float64),np.array(vectors,dtype=np.float64))

# anishpdb = PDB("PDBfiles/anish.pdb")
# print(anishpdb.file_name)

# anishstruc = PDBtoStructure(anishpdb)
# anishstruc.display('r')


        
        



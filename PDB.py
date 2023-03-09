# from .ss_id_module import Identify_ss
class PDB:
    '''Creates a new object type pdb that parses and process a pdbfile for easy use later
    by Mukundan S

    # objects
    file_name: file name
    all_lines : all lines
    ATOM_lines : all lines starting with ATOM
    HETATM_lines : All lines starting with HETATM
    ATOM_HETATM_lines : All lines starting with ATOM and HETATM
    res_pdb_line_dict : Dictionary of format {res_id:{atom:cordinate_line}}  res_id = line[21:26] atom = line[12:16].strip()
    res_pdb_cord_dict : Dictionary of format {res_id:{atom:tuple of x,y,z cordinates}} res_id = line[21:26] atom = line[12:16].strip()
    res_id_res_name_di : Dictionary of format {res_id:res_name}
    atom_line_di: Dictionary of format {res_idatom:cordinate_line}
    atom_cord_di: Dictionary of format {res_idatom:tuple of x,y,z cordinates}

    # functions
    get_ss() : calculates secondary structures
        returns : (helix,sheet)
        helix = {'helix id':[residue ids]}
        sheet = {sheet id_strand id:[residue id]}

    get_distance(res_id_1,atm_1,res_id_2,atm_2):
        Gives distance between two atoms sepcified by res_id and atom name

    eucl_dist(cord1,cord2):
        gives euclidian distance between two points
    
    get_res_id(ATOM_line):
        gives residue id of a given line

    write_file(file_name)
        writes file with all atom lines
    '''
    def __init__(self,file_name:str,lines = False):
        # save file name
        try:
            self.file_name = file_name
            with open(file_name) as inf:
                self.all_lines = inf.readlines()
        except:
            self.all_lines = lines
        # Only atom lines
        try:
            self._process_lines()
        except Exception as e:
            import sys
            if file_name:
                print('ERROR processing file to make PDB object! File name: '+file_name)
                print(e)
                sys.exit()
            else:
                print('ERROR while making PDB object!')
                print(e)
                sys.exit()

    def get_ss(self):        
        iss = Identify_ss(self.all_lines)
        iss.find_helix()
        iss.find_sheet()
        return (iss.helix,iss.sheet)

    def get_distance(self,res_id_1:str,atm_1:str,res_id_2:str,atm_2:str):
        '''
        Gives euclidian distance between res_id_1 atm_1 and res_id_2 atm_2
        '''
        return self.eucl_dist(self.res_pdb_cord_dict[res_id_1][atm_1],self.res_pdb_cord_dict[res_id_2][atm_2])
    
    def eucl_dist(self,cord_1,cord_2):
        return (sum([(cord_1[x]-cord_2[x])**2 for x in range(len(cord_1))]))**0.5
    
    def get_res_id(self,ATOM_line):
        return ATOM_line[21:26]
    
    def write_file(self,file_name):
        # writes file as file name from all lines
        with open(file_name,'w') as inf:
            inf.writelines(self.all_lines)
        self.file_name = file_name
    
    def write_assembled_file(self,file_name,atom_no_start = 1):
        '''
        Assembles inner values res_pdb_line_dict and writes lines.
        WARNING: Renumbers atoms
        '''
        lines = []
        for res_id in self.res_pdb_line_dict:
            lines+=list(self.res_pdb_line_dict[res_id].values())
        with open(file_name,'w') as inf:
            inf.writelines(self.rewrite_atom_id(atom_no_start,lines))
        self.file_name = file_name

    def rewrite_atom_id(self,start_atom_id,PDB_lines): # rewrites atom and res id based on starting vals
        new_lines = []        
        for line in PDB_lines:
#             atom_id,res_id =  
            new_line = line[0:6]+'{:5d} '+line[12:]            
            new_lines.append(new_line.format(start_atom_id))
            start_atom_id+=1
            # print(start_atom_id)
        return new_lines

    def rewrite_atom_id(self,start_atom_id,PDB_lines): # rewrites atom and res id based on starting vals
        new_lines = []        
        for line in PDB_lines:
#             atom_id,res_id =  
            new_line = line[0:6]+'{:5d} '+line[12:]            
            new_lines.append(new_line.format(start_atom_id))
            start_atom_id+=1
            # print(start_atom_id)
        return new_lines

    def _process_lines(self):
        self.ATOM_lines = []
        self.ATOM_HETATM_lines = []
        self.HETATM_lines = []
        # add more objects here to add more kind of lines
        
        # other_objects
        self.res_pdb_line_dict = {}
        self.res_pdb_cord_dict = {}
        self.res_id_res_name_di = {}
        self.atom_line_di = {}
        self.atom_cord_di = {}

        for line in self.all_lines:
            # print(line)
            if line[0:6] == 'ATOM  ':
                self.ATOM_lines.append(line)
                self.ATOM_HETATM_lines.append(line)
            
            if line[0:6] == 'HETATM':
                self.ATOM_HETATM_lines.append(line)
                self.HETATM_lines.append(line)
            
            # Processing dict
            if line in self.ATOM_HETATM_lines:
                res = line[21:26]
                atom = line[12:16].strip()
                self.atom_line_di[res+atom] = line
                try:
                    self.atom_cord_di[res+atom] = (float(line[30:38]),float(line[38:46]),float(line[46:54]))
                except Exception as e:
                    print('Skipping line:',line,'Due to ERROE:',e)
                if res in self.res_pdb_line_dict.keys():
                    self.res_pdb_line_dict[res][atom] = line
                    self.res_pdb_cord_dict[res][atom] = (float(line[30:38]),float(line[38:46]),float(line[46:54]))
                else:
                    self.res_pdb_line_dict[res] = {}
                    self.res_pdb_line_dict[res][atom] = line
                    self.res_pdb_cord_dict[res] = {}
                    self.res_pdb_cord_dict[res][atom] = (float(line[30:38]),float(line[38:46]),float(line[46:54]))
                # res name
                self.res_id_res_name_di[res] =  line[17:20]       

class Identify_ss():
    '''
    Module to identify residues belonging to helices and sheets 
    Mukundan S
    '''
    def __init__(self,pdb):
        '''        
        Identifies the residues crropsoinding to alpha helix and beta sheets from pdb records
        
        Steps:
        1) read HELIX and SHEET lines
        2) For each helix in helix line identify residues
        3) For each sheet identify residues for each strand
        
        Helices can be found as a dictionary {'helix id':[residue ids]} as self.helix
        and sheets as {sheet id_strand id:[residue id]} as self.sheet
        '''        
        self.pdb = pdb
        self.helix_record = [x for x in pdb if x.startswith('HELIX')]
        self.sheet_record = [x for x in pdb if x.startswith('SHEET')]        
        
        self.find_helix()
        self.find_sheet()
        
    def find_sheet(self):
        self.sheet = {} # dictionary of the form {sheet id_strand id:[residue id]}
        for x in self.sheet_record:
            sheet_id = x[11:14].strip()
            strand_id = x[7:10].strip()
            # making key for dict
            key = sheet_id+'_'+strand_id
            
            init_res_chain =x[21]
            init_res_no = int(x[22:26])
            term_res_no = int(x[33:37])
            res_id_list = [self.get_re_id(init_res_chain,x) for x in range(init_res_no,term_res_no+1)]
            self.sheet[key] = res_id_list
    
    def find_helix(self):
        self.helix = {} # dictionary of the form {'helix id':[residue ids]}
        for x in self.helix_record:
            helix_id = x[11:14].strip()
            init_res_chain= x[19]
            init_res_no = int(x[21:25])            
            term_res_no = int(x[33:37])
            
#             print(init_res_no,term_res_no)
            res_id_list = [self.get_re_id(init_res_chain,x) for x in range(init_res_no,term_res_no+1)]
            self.helix[helix_id] = res_id_list
            
    def get_re_id(self,chain,res):
        return '{:1s}{:4d}'.format(chain,res) 
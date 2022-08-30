from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np
import extract_result_ansys as ERA
import util
from scipy.optimize import linear_sum_assignment


SOLID_COMP = "SOLID_NODES"
FLUID_COMP = "FLUID_NODES"

@dataclass
class Result_propreties:
    """
    Notes : The result propreties class is a way to store the result 
    information that we want.

    Attributes:
        modes_to_extract( type:list ): 
            The index of the list of frequencies to extract

        n_ddl( type:int ):
            The number of degree of freedom that the analysis contains

        solid_comp( type:str ):
            Indicates the name of the solid component indicated in the
            .mac file
        
        fluid_comp( type:str ):
            Indicates the name of the fluid component indicated in the
            .mac file

        disp_ddls( type:list ):
            Indicates which degree of freedoms are the displacement one
        
        press_ddls( type:list ):
            Indicates which degree of freedoms are the pressure one
        
        extract_matrix( type:bool ):
            Indicates if we want to extract the matrix from the .full file
        

    """
    modes_to_extract:list
    n_ddl:int
    solid_comp:str = SOLID_COMP
    fluid_comp:str = FLUID_COMP
    disp_ddls:list = field(default_factory=lambda: [0,1,2])
    press_ddls:list = field(default_factory=list)
    extract_matrix:bool = False
    
    def __str__(self):
        return str(vars(self))

# class Result_propreties:
#     def __init__(self, modes_to_extract, n_ddl, solid_comp = SOLID_COMP, fluid_comp = FLUID_COMP, 
#                 disp_ddls = [0,1,2], press_ddls = None, extract_matrix = False):
#         self.modes_to_extract = modes_to_extract
#         self.n_ddl = n_ddl
#         self.solid_comp = solid_comp
#         self.fluid_comp = fluid_comp
#         self.disp_ddls = disp_ddls
#         self.press_ddls = press_ddls
#         self.extract_matrix = extract_matrix
    
#     def __str__(self):
#         return str(vars(self))

class Result_class:
    def __init__(self, nodal_position:np.ndarray, nodal_results:np.ndarray, 
                 id_to_index:dict, components:dict = None, name:str = "results"):
        """
        Notes : The Result_class is used to represent a general result from Ansys.

        Attributes:
            nodes ( :obj:np.array ): 
                Cartesian positions of the nodes of the analysis, 
                nodes[index_node,:] = [index_position_x, index_position_y, index_position_z]

            nodal_results( :obj:np.array ): 
                Raw results of the analysis classed like: 
                nodal_results[index_node, index_ddl, other_dimension_index] = ddl_value
            
            components( type:dict ):
                dictionnary of the components of the model. For example, if the components
                are seperated into a solid component "SOLID_NODES" and a fluid component 
                "FLUID_NODES" :
                components = {"SOLID_NODES":solid_indexs_nodes, "FLUID_NODES":fluid_indexs_nodes}

            id_to_index( type:dic ): 
                For a node id given from the APDL model, id_to_index[id] gives the index of 
                all the array which corresponds to the nodes identified by the id

        """
        self.nodes = nodal_position
        self.nodal_results = nodal_results
        self.id_to_index = id_to_index
        if components != None:
            self.components = {}
            for comp, array in components.items():
                self.components[comp] = ERA.list_id_to_index(array, id_to_index)
        self.name = name


class Modal_result(Result_class):
    """
    Notes : The Modal_result class is used to represent modal results from an 
    APDL simulation.

    Attributes:
        frequencies( type:list ): 
            List of the frequencies of the results

        nodes ( :obj:np.array ): 
            Cartesian positions of the nodes of the analysis, 
            nodes[index_node,:] = [index_position_x, index_position_y, index_position_z]

        nodal_results( :obj:np.array ): 
            Raw results of the analysis classed like: 
            nodal_results[index_node, index_ddl, index_freq] = ddl_value
        
        components( type:dict ):
            dictionnary of the components of the model. For example, if the components
            are seperated into a solid component "SOLID_NODES" and a fluid component 
            "FLUID_NODES" :
            components = {"SOLID_NODES":solid_indexs_nodes, "FLUID_NODES":fluid_indexs_nodes}

        id_to_index( type:dic ): 
            For a node id given from the APDL model, id_to_index[id] gives the index of 
            all the array which corresponds to the nodes identified by the id
        
        freq_index( type:list ): 
            For an index of the freqs results available, freq_index[index] = index of the 
            list of all the frequency available for the APDL results

        freq_id_to_index( type:list ): 
            The inverse of freq index 
        
        displacement( :obj: np.array)
            Gives all the displacement of the nodes :
            displacement[node_index,:, freq_index] = [x,y,z]
        
        pressure( :obj: np.array)
            Gives all the displacement of the nodes :
            displacement[node_index,:, freq_index] = pressure

    """
    def __init__(self, frequencies:list, nodal_position:np.ndarray, components:dict,
                nodal_results:np.ndarray, id_to_index:dict, mode_index_to_freq_index:list, 
                disp_ddls:list = None, press_ddls:list = None, name:str = "Modal results"):
        super().__init__(nodal_position, nodal_results, id_to_index, components, name)
        self.frequencies = frequencies
        self.nodes = nodal_position
        self.nodal_results = nodal_results #deplacements et pression
        # self.components = {}
        # for comp, array in components.items():
        #     self.components[comp] = ERA.list_id_to_index(array, id_to_index)
        # self.id_to_index = id_to_index
        self.freq_index = mode_index_to_freq_index
        if disp_ddls != None:
            self.set_disp(disp_ddls)
        if press_ddls != None:
            self.set_press(press_ddls)
        self.index_id_freq()
        
    
    def index_id_freq(self):
        """ the method index_id_freq

        Note: 
            assign the attributes freq_id_to_index
        """
        L = [self.freq_index[i] for i in range(self.n_freq())]
        self.freq_id_to_index = L

    def disp_pol(self, id_freq):
        """the method disp_pol

        Note: 
            Convert the displacements to polar displacement
        
        """
        index_freq = self.freq_id_to_index[id_freq]
        displacement = self.displacement[:,:,index_freq]
        nodes_pos = self.nodes
        total_disp = nodes_pos + displacement
        XY = total_disp[:,0:2]
        Z = total_disp[:,2]
        pol = util.cart2pol_array(XY)
        pol = np.insert(pol, 2, Z, axis=1)
        pol_pos = self.pol_position()
        return pol - pol_pos

    def n_freq(self):
        """
        return the number of modes in the results
        """
        return len(self.freq_index)

    def n_nodes(self):
        """
        return the number of nodes in the results
        """
        n_ligns,  = np.shape(self.nodes[:,0])
        return n_ligns

    def node_indexs(self):
        """
        return the list of nodes indexs
        """
        indexs = [i for i in range(self.n_nodes())]
        return indexs

    def set_disp(self, disp_ddls:list):
        """
        set the displacement from the modal results
        """
        displacement = self.nodal_results[:,disp_ddls,:]
        zeros = np.zeros((self.n_nodes(),1))
        i = len(disp_ddls)
        while i < 3:
            displacement = np.insert(displacement, i, zeros, axis=1)
            i += 1
        self.displacement = displacement
    
    def set_press(self, press_ddls:list):
        """
        set the pressure from the modal result
        """
        self.pressure = self.nodal_results[:,press_ddls,:]

    def pol_position(self):
        """
        return the position of nodes in polar coordinates
        """
        positions = self.nodes
        positions_xy = positions[:,0:2]
        positions_z = positions[:,2]
        pol = util.cart2pol_array(positions_xy)
        pol = np.insert(pol, 2, positions_z, axis=1)
        return pol
    
    def select(self, by:str = "position", type:str = "equal",
                coord:str = "cart", X = None, Y = None, 
                Z = None, tol:float = 1e-8):
        """
        It's the selection method to select methodically the nodes

        Args:
            by( str ): 
                Indicates where the selecton focus :
                    "position"(default) : the selection focus on the position of the nodes

            type( str ): 
                Indicates how to perform the selection :
                    "equal"(default) : selection the nodes where coordinates are equal (with tolerance),
                                    to the value entered.
                    "between" : selection the nodes where coordinates are between 2 different
                                values entered.

            coord( str ):
                To choose whether we select in cartesian or polar coordinates :
                    "cart"(default) : we select in cartesian coordinates

                    "pol": we select in polar coordinates, X = R, Y = Theta, Z = Z

            X, Y, Z(float of list):
                The values for the selection
            
            tol( float ):
                The tolerance of the selection
            
        Returns:
            indexs( list ): list of the indexs of the selected nodes
            

        """
        indexs = self.node_indexs()
        if by == "position":
            if coord == "cart":
                positions = self.nodes
            if coord == "pol":
                positions = self.pol_position()
            if X != None:
                if type == "equal":
                    _, indexs_x = ERA.select_equal(positions, X, 
                                                0, tol)
                if type == "between":
                    _, indexs_x = ERA.select_between(positions, X[0], X[1], 
                                                0, tol)
            else:
                indexs_x = indexs
    
            if Y != None:
                if type == "equal":
                    _, indexs_y = ERA.select_equal(positions, Y, 
                                                1, tol)
                if type == "between":
                    _, indexs_y = ERA.select_between(positions, Y[0], Y[1], 
                                                1, tol)
            else:
                indexs_y = indexs

            if Z != None:
                if type == "equal":
                    _, indexs_z = ERA.select_equal(positions, Z, 
                                                2, tol)
                if type == "between":
                    _, indexs_z = ERA.select_between(positions, Z[0], Z[1], 
                                                2, tol)
            else:
                indexs_z = indexs
            
            indexs = list(set(list(indexs_x)).intersection(set(list(indexs_y))).intersection(set(list(indexs_z))))
            indexs.sort()
            return indexs
    
    def unique_frequencies(self, decimal:int = None, ids_ = False, neg = False):
        """
        The methode get_frequencies
        Note : 
            return the frequency of the analysis
        
        Args :
            decimal( int ):
                it's the number of decimal we want for the separation of the unicity
                selection of the frequencies.
            
            ids_( bool ):
                Default = False. Indicates whether to return the indexs of self.frequencies
                chosen to be unique.
            
            neg( bool ):
                Default = False. Indicates if we want negative frequency or not.
            
        Returns:
            frequencies( np.array ):
                The list of frequency selected
            
            ids( list ):
                The list of indexs of the frequencies selected
            
        """
        frequencies = self.frequencies
        if decimal == None:
            id_pos = np.where(frequencies >= 0)
            return frequencies[id_pos], id_pos
        else:
            frequencies_, ids = np.unique(frequencies.round(decimals=decimal), return_index=True)
            if neg == False:
                id_pos = np.where(frequencies_ >= 0)
                ids = ids[id_pos]
            if ids_:
                return frequencies[ids], ids
            else:
                return frequencies[ids]
    
    def associate_frequencies(self, approx_freqs, decimal = None, 
    score = False, family_id = False):
        freqs, ids = self.unique_frequencies(decimal=decimal)
        freqs_assigned = np.zeros(np.size(approx_freqs))
        cost_matrix = util.diff_cost(approx_freqs, freqs)
        row_inds, col_inds = linear_sum_assignment(cost_matrix)
        for row_ind, col_ind in zip(row_inds, col_inds):
            freqs_assigned[row_ind] = freqs[col_ind]
        if score:
            score = np.sum(cost_matrix[row_inds, col_inds])
            if family_id:
                return freqs_assigned, row_inds, score
            else:
                return freqs_assigned, score
        return freqs_assigned

            

    
    def get_frequency(self, approx_freq:float, id = False, decimal = None, 
                    inferior = False):
        """
        The methode get_frequency
        Note : 
            To select the frequency with an approximation of the frequency
        
        Args :
            approx_freq( float ):
                The approximation of the frequency
            
            id( bool ):
                Default = False. Indicates whether we want the index of self.frequencies
                wich corresponds to the frequency selected
            
        Returns:
            frequencies( float ):
                The frequency selected
            
            ids( list ):
                The index of the frequency selected
            
        """
        frequencies, idxs = self.unique_frequencies(decimal, True)
        if inferior:
            id_0, = np.where((approx_freq - frequencies) >= 0)
            if len(id_0) == 0:
                condition = abs(frequencies-approx_freq) == min(abs(frequencies-approx_freq))
                id_ = np.where(condition)
            else:
                frequencies = frequencies[id_0]
                id_, = np.where( abs(frequencies[id_0]-approx_freq) == min(abs(frequencies[id_0]-approx_freq)))
                id_ = id_0[id_]
        else:
            condition = abs(frequencies-approx_freq) == min(abs(frequencies-approx_freq))
            id_ = np.where(condition)
        frequency = frequencies[id_]
        if id:
            return frequency, idxs[0][id_]
        else : 
            return frequency
    
    




class Modal_result_w_matrixs(Modal_result):
    """
    Notes : The Modal_result class is used to represent modal results from an 
    APDL simulation.

    Attributes:
        frequencies( type:list ): 
            List of the frequencies of the results

        nodes ( :obj:np.array ): 
            Cartesian positions of the nodes of the analysis, 
            nodes[index_node,:] = [index_position_x, index_position_y, index_position_z]

        nodal_results( :obj:np.array ): 
            Raw results of the analysis classed like: 
            nodal_results[index_node, index_ddl, index_freq] = ddl_value
        
        components( type:dict ):
            dictionnary of the components of the model. For example, if the components
            are seperated into a solid component "SOLID_NODES" and a fluid component 
            "FLUID_NODES" :
            components = {"SOLID_NODES":solid_indexs_nodes, "FLUID_NODES":fluid_indexs_nodes}

        id_to_index( type:dic ): 
            For a node id given from the APDL model, id_to_index[id] gives the index of 
            all the array which corresponds to the nodes identified by the id
        
        freq_index( type:list ): 
            For an index of the freqs results available, freq_index[index] = index of the 
            list of all the frequency available for the APDL results

        freq_id_to_index( type:list ): 
            The inverse of freq index 
        
        displacement( :obj: np.array)
            Gives all the displacement of the nodes :
            displacement[node_index,:, freq_index] = [x,y,z]
        
        pressure( :obj: np.array)
            Gives all the displacement of the nodes :
            displacement[node_index,:, freq_index] = pressure

    """
    def __init__(self, frequencies:list, nodal_position:np.ndarray, components:dict,
                nodal_results:np.ndarray, id_to_index:dict, mode_index_to_freq_index:list, 
                disp_ddls:list = None, press_ddls:list = None, K:np.array = None, M:np.array = None,
                eigen_vectors = None, ev = None, name:str = "Modal results"):
        
        super().__init__(frequencies, nodal_position, components, nodal_results, id_to_index,
        mode_index_to_freq_index, disp_ddls, press_ddls)
        self.K = K
        self.M = M
        self.eigen_vectors = eigen_vectors
        self.ev = ev
        self.modal_mass_dic = {}
        self.modal_stiff_dic = {}
    
    def get_eigen_vec(self, idx_mode):
        return self.eigen_vectors[:, idx_mode]

    def eigen_val(self, idx_mode):
        ev = self.ev[idx_mode]
        ev = (2 * np.pi * ev) ** 2
        return ev

    def modal_mass(self, idx_mode):
        modal_mass_dic = self.modal_mass_dic
        if idx_mode not in modal_mass_dic:
            eigen_vec = self.get_eigen_vec(idx_mode)
            modal_mass_coeff = util.compute_modal_coeff(self.M, eigen_vec)
            modal_mass_dic[idx_mode] = modal_mass_coeff
        else:
            modal_mass_coeff = modal_mass_dic[idx_mode]
        return modal_mass_coeff

    def modal_stiff(self, idx_mode):
        modal_stiff_dic = self.modal_stiff_dic
        if idx_mode not in modal_stiff_dic:
            eigen_vec = self.get_eigen_vec(idx_mode)
            modal_stiff_coeff = util.compute_modal_coeff(self.K, eigen_vec)
            modal_stiff_dic[idx_mode] = modal_stiff_coeff
        else:
            modal_stiff_coeff = modal_stiff_dic[idx_mode]
        return modal_stiff_coeff
    
    def compute_added_mass(self, acoustic_freq, idx_mode):
        ev_acc = util.freq_to_ev(acoustic_freq)
        modal_mass_coeff = self.modal_mass(idx_mode)
        modal_stiff_coeff = self.modal_stiff(idx_mode)
        added_mass = 1 / ev_acc * modal_stiff_coeff - modal_mass_coeff
        return added_mass

    def supp_matrixs(self):
        delattr(self, "K")
        delattr(self, "M")
        

    def get_frequencies(self, decimal:int = None, ids_ = False, neg = False):
        """
        The methode get_frequencies
        Note : 
            return the frequency of the analysis
        
        Args :
            decimal( int ):
                it's the number of decimal we want for the separation of the unicity
                selection of the frequencies.
            
            ids_( bool ):
                Default = False. Indicates whether to return the indexs of self.frequencies
                chosen to be unique.
            
            neg( bool ):
                Default = False. Indicates if we want negative frequency or not.
            
        Returns:
            frequencies( np.array ):
                The list of frequency selected
            
            ids( list ):
                The list of indexs of the frequencies selected
            
        """
        if decimal == None:
            return self.ev
        else:
            frequencies = self.ev
            frequencies_, ids = np.unique(frequencies.round(decimals=decimal), return_index=True)
            if neg == False:
                id_pos = np.where(frequencies_ >= 0)
                ids = ids[id_pos]
            if ids_:
                return frequencies[ids], ids
            else:
                return frequencies[ids]

                
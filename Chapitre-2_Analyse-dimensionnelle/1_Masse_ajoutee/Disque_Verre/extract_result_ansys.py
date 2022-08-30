import matplotlib.pyplot as plt
from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.core import launch_mapdl
import ansys.dpf.core as dpf
import numpy as np
import pathlib
import pickle as pkl
import analyse_result_ansys as ARA
import util
import pprint as pp

def extract_save_modal(APDL_folder:pathlib.Path, result_file_name:pathlib.Path,
                data_file_name:pathlib.Path, modes_to_extract:list = None, n_ddl:int = None, 
                saving_type:str = "result-class", disp_ddls:list = [0,1,2], press_ddls:list = [3], 
                mapdl = None, save:bool = True, save_file:pathlib.Path = None, 
                result_propreties:ARA.Result_propreties = None, 
                extract_matrix_:bool = False, file_name_full:pathlib.Path = None):
        
    """ the function extract_save

        Note: 
            This function combines all the belows functions to extract the information
            from a modal analysis. Needs to be extended for other type of analysis.

        Args:
            APDL_folder( pathlib.Path ): 
                The location of the APDL working folder.

            result_filde_name( pathlib.Path ): 
                The relative location of the .rst file into the working folder

            data_file_name( pathlib.Path ):
                The absolute location of the data file .dat

            modes_to_extract( list ):
                The index of the modes to extract. Indexs corresponds to the list indexs
                of the frequency list.
            
            n_ddl( int ):
                The exact number of ddl taken into account in the APDL analysis. The most difficult
                information to have but must be exact.

            saving_type( str ):
                Gives the type object to save the result of the extraction :
                    "result-class"( default ): all the information is stored into
                        a result object.
                    "dic": all the information is stored into a dictionnary
            
            disp_ddls( list ):
                the indexs of the ddls which correspond to the displacement ddls of the analysis
            
            press_ddls( list ):
                the indexs of the ddl which correspond to the pressure ddl of the analysis
            
            mapdl( mapdl ):
                The mapdl session if already charged by launch_mapdl
            
            save( bool ):
                If the results need to be saved or not

            save_file( pathlib.Path ):
                The location of the file to save the results
            
            result_propreties( ARA.Result_propreties ):
                An object where all the propreties are already stored
            
            extract_matrix_( bool ):
                If the function needs to extract the matrixs of the problem or not.
                Default = False
            
            file_name_full( pathlib.Path ):
                The location of the .full file where are the matrixs informations.
            

        Returns:
            results( object or dict ):
                returns all the result under a dict or an object form. 
    """

    if result_propreties != None:
        modes_to_extract = result_propreties.modes_to_extract
        n_ddl = result_propreties.n_ddl
        disp_ddls = result_propreties.disp_ddls
        press_ddls = result_propreties.press_ddls
        extract_matrix_ = result_propreties.extract_matrix

    #extraction of the frequencies :
    frequencies = extract_modal_frequencies(result_file_name, mapdl = mapdl)
    

    #extraction of the nodal position :
    nodal_position = extract_nodal_position(data_file_name)

    #extraction of the components file
    components = extract_components(data_file_name)

    #extraction of the nodal results
    nodal_results, id_to_index, sol_index_to_freq_index = extract_modal_result(APDL_folder, result_file_name,
                                                modes_to_extract, n_ddl, mapdl = mapdl)
    frequencies = frequencies[sol_index_to_freq_index]
    if extract_matrix_:
        M_matrix, K_matrix, ev, eigen_vectors = extract_solve(APDL_folder, file_name_full, 
        mapdl=mapdl, n_ev = len(frequencies), array = True)
        dictionnary = {"freq": frequencies, "comp": components, 
                        "nodal-position":nodal_position, 
                        "nodal-result": nodal_results, 
                        "id-to-index": id_to_index,
                        "freq-index": sol_index_to_freq_index,
                        "K":K_matrix, "M":M_matrix}
    else:
        dictionnary = {"freq": frequencies, "comp": components, 
                        "nodal-position":nodal_position, 
                        "nodal-result": nodal_results, 
                        "id-to-index": id_to_index,
                        "freq-index": sol_index_to_freq_index}
    if saving_type == "dic":
        if save:
            save_object(save_file, dictionnary, "object")
        return dictionnary
    if saving_type == "result-class":
        if extract_matrix_:
            result = ARA.Modal_result_w_matrixs(frequencies, nodal_position, components,
                                nodal_results, id_to_index, sol_index_to_freq_index,
                                disp_ddls=disp_ddls, press_ddls=press_ddls, K = K_matrix,
                                M = M_matrix, ev=ev, eigen_vectors=eigen_vectors)
        else:
            result = ARA.Modal_result(frequencies, nodal_position, components,
                                nodal_results, id_to_index, sol_index_to_freq_index,
                                disp_ddls=disp_ddls, press_ddls=press_ddls)
        if save:
            save_object(save_file, result, "object")
        return result

def save_object(file_name:pathlib.Path, object, array_type = False):
    """ the function save_object 

        Note: 
                Allows to save arrays or object

        Args:
            filename( pathlib.Path ): 
                Gives the path of the file where to save the object
            object: 
                can be any type of object
            array_type( boolean ) : 
                indicates if it's an array or not


    """
    file_name.parents[1].mkdir(parents = True, exist_ok = True)
    with open(file_name, "wb") as f:
        if array_type:
            np.save(f, object)
        else:
            pkl.dump(object, f)

def load_object(file_name:pathlib.Path, array_type = False):
    """ the function load_object

        Note: 
            Allows to load object from file where save_object has been used before

        Args:
            filename( pathlib.Path ): Gives the path of the file where to save the object
            type_( boolean ) : indicates if it's an array or not


    """
    with open(file_name, "rb") as f:
        if array_type:
            object = np.load(f, allow_pickle = True)
        else:
            object = pkl.load(f)
        return object

def select_equal(U:np.array, u:float, ddl:int, tol:float):
    """ the function select_equal

        Note: 
            Selects the nodes where his ddl is equal or close to u

        Args:
            U( np.array ): 
                It's a vector where every ligns refer to a node and every column
                refers to a ddl

            u( float ): 
                It's the value that we want to select

            ddl( int ):
                It's the integer that refers to the ddl in question 

            tol( float ):
                It's the tolerance which give the range of the selections :
                abs(U - u) <= tol are selected
        
        Returns:
            U_selected( np.array ): 
                The ligns selected from U 
            ligns( int ):
                The index of the ligns selected
    """
    (ligns,) = np.where(np.abs(U[:,ddl] - u) <= tol )
    return U[ligns, :], ligns

def select_between(U:np.array, u_0:float, u_1:float, 
                    ddl:int):
    """ the function select_between

        Note: 
            Selects the nodes where the ddl is between 2 limits

        Args:
            U( np.array ): 
                It's a vector where every ligns refer to a node and every column
                refers to a ddl

            u_0( float ): 
                It's where U[:,ddl] >= u_0

            u_1( float ): 
                It's where U[:,ddl] <= u_1

            ddl( int ):
                It's the integer that refers to the ddl in question
        
        Returns:
            U_selected( np.array ): 
                The ligns selected from U 
            ligns( int ):
                The index of the ligns selected
    """
    ligns = np.where(u_0 <= U[:,ddl] and U[:,ddl] <= u_1)
    return U[ligns, :], ligns

def list_id_to_index(ids:list, id_to_index:dict):
    """ the function select_between

        Note: 
            From a list of APDL ids, gives a list of 
            python indexs

        Args:
            ids( list ): 
                List of the ids

            id_to_index( dictionnary ):
                A dictionnary where the keys are the ids and
                the associated items are the ids
        
        Returns:
            indexs( list of int ):
                A list of indexs.
    """
    indexs = []
    for id in ids:
        indexs.append(id_to_index[id])
    return indexs

def sup_nan_data(array:np.array, tresh:float = 1e20):
    """ method sup_nan_data

        Note: 
            supress the erroned datas represented by high number by APDL by
            putting it to 0.

        Args:
            array( np.array ): 
                The array to treat the nan_data

            tresh( float ):
                The treshold where the data above should be fixed to 0
        
        Returns:
            array ( np.array ):
                The array processed
    """
    array = np.copy(array)
    ligns, columns = np.where(array > tresh)
    for lign, column in zip(ligns, columns):
            array[lign, column] = 0
    return array

def order_array(array:np.array, id_to_index:dict, index_to_id:list):
    """ method order_array

        Note: 
            Order an array with the ids order.

        Args:
            array( np.array ): 
                The array to be ordered

            id_to_index( dictionnary ):
                A dictionnary where the keys are the ids and
                the associated items are the ids
            
            index_to_id( list ):
                A list that gives the id which corresponds to the index
                of the list
        
        Returns:
            ordered_array ( np.array ):
                The array ordered
    """
    array = np.copy(array)
    n_lign, n_col = np.shape(array)
    ordered_array = np.zeros((n_lign, n_col))
    for index, id in enumerate(index_to_id):
        ordered_array[id_to_index[id],:] = array[index, :]
    return ordered_array

def extract_modal_frequencies(result_file:pathlib.Path, save:bool = False, 
                            saving_file:pathlib.Path = None, mapdl = None):
    """ function extract_modal_frequencies

        Note: 
            From a given result file extract the modal frequencies
            using the ansys module DPF 

        Args:
            result_file( pathlib.Path ): 
                It's the .rst file after the simulation

            save( bool ):
                If the needs to save to a file this information. Default is False.
            
            saving_file( pathlib.Path ):
                The location of the file where saving this information
            
            mapdl( mapdl ):
                The mapdl session if already charged by launch_mapdl
        
        Returns:
            frequencies( np.array ):
                list of the different frequencies
    """
    #extractions des fréquences :
    result_file = str(result_file)
    model = dpf.Model(result_file)
    metadatas = model.metadata

    if mapdl == None:
        tf = metadatas.time_freq_support
        frequencies = tf.time_frequencies.data
    else:
        xpl = mapdl.xpl
        pure_path = pathlib.PurePath(result_file)
        folder_result_file = pure_path.parts[-1]
        xpl.open(folder_result_file)
        freqs = xpl.read('TIM').asarray()
        frequencies = freqs

    if save:
        save_object(saving_file, frequencies)
    return frequencies

def extract_nodal_position(model_file:pathlib.Path, save:bool = False, 
                    saving_file:pathlib.Path = None):
    """ function extract_modal_position

        Note: 
            From a given model file extract the position of the nodes

        Args:
            model_file( pathlib.Path ): 
                It's the .dat when the model is created

            save( bool ):
                If the needs to save to a file this information. Default is False.
            
            saving_file( pathlib.Path ):
                The location of the file where saving this information
        
        Returns:
            nodes( np.array ):
                nodes[index,:] = [x, y, z]
    """
    model_file = str(model_file)
    archive = pymapdl_reader.Archive(model_file)
    if save:
        save_object(saving_file, archive.nodes)
    return archive.nodes

def extract_components(model_file:pathlib.Path, save = False, 
                    saving_file:pathlib.Path = None):
    """ function extract_components

        Note: 
            From a given model file extract the components of the model.
            For example, the component could be FLUID_NODES and SOLID_NODES

        Args:
            model_file( pathlib.Path ): 
                It's the .dat when the model is created

            save( bool ):
                If the needs to save to a file this information. Default is False.
            
            saving_file( pathlib.Path ):
                The location of the file where saving this information
        
        Returns:
            nodes_components( dict ):
                dictionnary where keys are the components name and items are the
                ids of the node
    """
    model_file = str(model_file)
    archive = pymapdl_reader.Archive(model_file)
    if save:
        save_object(saving_file, archive.node_components, "object")
    return archive.node_components

def extract_modal_result(folder:pathlib.Path, result_file:pathlib.Path, 
                    mode_indexs:list, n_ddl:int, save:bool = False, saving_file:pathlib.Path = None, 
                    mapdl = None):
    """ function extract_components

        Note: 
            From a given working APDL folder location and the result file location, extract the 
            informations from the modal analysis

        Args:
            folder( pathlib.Path ): 
                It's the working folder of the analysis
            
            result_file( pathlib.Path ): 
                It's the location of the .rst file into the working folder

            mode_indexs( list ): 
                The modes indexs (the index are the index of the frequency list 
                extracted with extract_nodal_frequency)
            
            n_ddl( int ):
                The exact number of ddl taken into account into the analysis. (It's may be the most 
                difficult part to know and needs to be exact)
            
            save( bool ):
                If the needs to save to a file this information. Default is False.
            
            saving_file( pathlib.Path ):
                The location of the file where saving this information

            mapdl( mapdl ):
                The mapdl session if already charged by launch_mapdl
        
        Returns:
            U( np.array ):
                U has 3 different components, U[ligns, columns, depths], where the ligns
                refer to the node index, the columns refer to the ddl and the depths, to the
                index of the modes 
                WARNING : the index of the mode or different from the index from the frequency list

            id_to_index( dict ):
                From a given id, gives the correspondant index

            sol_index_to_freq_index( list ):
                For a given mode index of the solution U, gives the correspondant
                index from the frequency list.  
            
    """
    #extractions of index to id :
    result_file = str(result_file)
    model = dpf.Model(result_file)
    metadatas = model.metadata
    nodes = metadatas.meshed_region.nodes
    id_to_index = nodes.mapping_id_to_index
    
    #launch of mapdl
    if mapdl == None:
        folder = str(folder)
        mapdl = launch_mapdl(run_location = folder)
    xpl = mapdl.xpl
    pure_path = pathlib.PurePath(result_file)
    folder_result_file = pure_path.parts[-1]
    xpl.open(folder_result_file)
    index_to_id = xpl.read("NOD").asarray()
    (n_lign,) = np.shape(index_to_id)
    n_modes = len(mode_indexs)
    U = np.zeros((n_lign, n_ddl, n_modes))
    sol_index_to_freq_index = []
    nb_modes_rst_file = util.nb_modes(xpl)
    if len(mode_indexs) > nb_modes_rst_file:
        mapdl.exit()
        string = "Not enough modes in the RST file, only : " + str(nb_modes_rst_file) 
        raise ValueError(string)

    for i, mode_index in enumerate(mode_indexs):
        sol_index_to_freq_index.append(mode_index)
        xpl.goto("DSI::SET" + str(mode_index + 1))
        u = xpl.read("NSL")
        u = u.asarray()
        u = np.reshape(u, (-1, n_ddl))
        u = order_array(u, id_to_index, index_to_id)
        u = sup_nan_data(u)
        U[:,:,i] = u

    xpl.close()
    if save:
        dic = {"disp": U, "id_to_index": id_to_index,
                "mode-index":sol_index_to_freq_index}
        save_object(saving_file, dic, "object")
    if mapdl == None:
        mapdl.exit()
    return U, id_to_index, sol_index_to_freq_index

def extract_matrix(folder, file_name, mapdl = None, array = False, 
                    damping = False):
    if mapdl == None:
        folder = str(folder)
        mapdl = launch_mapdl(run_location = folder)
    mm = mapdl.math
    mapdl.finish()
    mm.free()
    K = mm.stiff(fname=file_name)
    M = mm.mass(fname=file_name)
    if damping:
        D = mm.damp(fname=file_name)
    if array:
        K = K.asarray()
        K = util.from_unsym_to_sym(K)
        M = M.asarray()
        M = util.from_unsym_to_sym(M)
        if damping:
            D = D.asarray()
            D = util.from_unsym_to_sym(D)
    if damping:
        return M, K, D
    else:
        return M, K

def eigs(M, K, mapdl, n_ev, array = True):
    mm = mapdl.math
    eigen_vectors = mm.zeros(K.nrow, n_ev)
    ev = mm.eigs(n_ev, K, M, phi=eigen_vectors, fmin = 1)
    if array:
        eigen_vectors = eigen_vectors.asarray()
        ev = ev.asarray()

    return ev, eigen_vectors

def extract_solve(folder, file_name, mapdl, n_ev, array = False):
    M, K = extract_matrix(folder, file_name, mapdl)
    ev, eigen_vectors = eigs(M, K, mapdl, n_ev, array)
    if array:
        K = K.asarray()
        K = util.from_unsym_to_sym(K)
        M = M.asarray()
        M = util.from_unsym_to_sym(M)

    return M, K, ev, eigen_vectors






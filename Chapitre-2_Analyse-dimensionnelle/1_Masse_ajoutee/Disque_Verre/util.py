import numpy as np
import scipy as spy
import scipy.interpolate as interp
import scipy.integrate as integ
import scipy.optimize as opt

#------------------ Exceptions ------------------
#----------------------------------------------------------




#------------------ Integrations ------------------
#----------------------------------------------------------

def scalar_prod_samples_func_1D(x_samples, y_samples, product_function, bounds, integration_product = None):
    inter_func = interp.interp1d(x_samples, y_samples, kind='cubic')
    if integration_product == None:
        product = lambda x : inter_func * product_function
        a, b = bounds
        result, error_est = integ.quad(product, a, b)
    else:
        result = integration_product(inter_func, product_function)
    return result

def phase_computation(func_1, func_2, scalar_product):
    phase = 0
    def func_to_minimize(phase):
        delta_func_2 = lambda x : func_2(x + delta_func_2)
        return - scalar_product(func_1, func_2)
    final_phase = opt.minimize(func_to_minimize, phase)
    return final_phase    
    

#------------------ MAPDL - XPL ------------------
#----------------------------------------------------------

def nb_modes(xpl):
    dic = xpl.json()
    list_informations = dic['children']
    dic = list_informations[3]
    list_of_modes = dic['children']
    return len(list_of_modes)


#------------------ Matrixs treatment ------------------
#----------------------------------------------------------

def from_unsym_to_sym(M):
    diag = spy.sparse.diags(M.diagonal())
    return M + M.transpose() - diag

def compute_modal_coeff(matrix, eigen_vector):
    coeff = eigen_vector.T.dot(matrix.dot(eigen_vector)) \
        / eigen_vector.T.dot(eigen_vector)
    return coeff

#------------------Eigen value following ------------------
#----------------------------------------------------------
def diff_cost(freqs_1, freqs_2):
    n_1, n_2 = np.size(freqs_1), np.size(freqs_2)
    cost_matrix = np.zeros((n_1, n_2))
    for i, freq_1 in enumerate(freqs_1):
        for j, freq_2 in enumerate(freqs_2):
            cost_matrix[i,j] = abs(freq_1 - freq_2)
    return cost_matrix

#------------------Added mass computation------------------
#----------------------------------------------------------

def freq_to_ev(freq):
    return (2 * np.pi * freq) ** 2











#---------------------Polar coordinates---------------------
#-----------------------------------------------------------
def cart2pol(x, y, positive_only=True):
    """Converts 2D cartesian coordinates to polar coordinates

    Parameters
    ----------
    x (float): x coordinate of point
    y (float): y coordinate of point
    positive_only (boolean): to have angles from 0 to 2*pi rather than -pi to pi

    Returns
    -------
    radius (float): calculated radius of point
    angle (float): calculated angle of point

    """
    radius = np.linalg.norm([x,y])
    angle = np.arctan2(y, x)
    if positive_only:
        if angle < 0:
            angle = 2*np.pi + angle

    return (radius, angle)

def cart2pol_array(U, positive_only = True):
    """
    Applies cart2pol to an array.
    """
    X = U[:,0]
    Y = U[:,1]
    (n_ligns,) = np.shape(X)
    U_pol = np.zeros((n_ligns, 2))
    for i in range(n_ligns):
        U_pol[i,0], U_pol[i,1] = cart2pol(X[i], Y[i], positive_only)
    return U_pol
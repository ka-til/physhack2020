import numpy as np
import matplotlib.pylab as plt
import scipy.constants as spc
import pandas as pd
from scipy import interpolate

class spectrum_line(object):
    '''
    class spectrum line contains useful functions to calculate the transmittance of a give wavenumber
    '''

    def __init__(self):
        '''
        As soon as the class is called, pandas is used to read the excel file to get some useful parameters
        '''
        print('Reading the excel file.')
        self.excel_file = pd.read_excel("Venus Atmosphere.xlsx")

    def getParams(self, text_file):
        '''
        Input: text_file, the text_file is dowloaded from HITRAN, contains information about a single gas
        Output: molecule_id: according to the HITRAN database
                v_center: center of the spectral line shape
                S_0: Spectral line intensity at 296K
                delta_air: Pressure shift at T = 296K and p = 1atm
                alphaLorentz_a: air -broadened half width half maximum at T = 296K and p = 1atm
                alphaDoppler_s: self-broadened half width half maximum at T = 296K and p = 1atm
                gamma: The coefficient of the temperature dependence of the air-broadened half width
                E_l: The lower-state energy of the transition (cm-1)
        '''
        f = np.loadtxt(text_file, delimiter=',')
        molecule_id = f[:, 0]
        print('I am molecule ID is '+str(molecule_id[0])+ ' taken from HITRAN database')

        isotopologue_id = f[:, 1]
        print('My isotopologue ID is '+str(isotopologue_id[0])+ ' taken from HITRAN database')

        v_center = f[:, 2]
        S_0 = f[:, 3]
        delta_air = f[:, 4]
        alphaLorentz_a = f[:, 5]
        alphaDoppler_s = f[:, 6]
        gamma = f[:, 7]
        E_l = f[:, 8]

        return molecule_id, v_center, S_0, delta_air, alphaLorentz_a, alphaDoppler_s, gamma, E_l

    def setParams(self, molecule_id, v_center, S_0, delta_air, alphaLorentz_a, alphaDoppler_s, gamma, E_l):

        '''
        setParams sets common parameters to be further used in the class
        The inputs are same as the outputs from getParams
        '''

        self.v_center = v_center
        self.S_0 = S_0
        self.delta_air = delta_air
        self.alphaLorentz_a = alphaLorentz_a
        self.alphaDoppler_s = alphaDoppler_s
        self.gamma = gamma
        self.E_l = E_l

        if molecule_id == 5.0:
            molecule = self.excel_file.CO
        if molecule_id == 22.0:
            molecule = self.excel_file.N2
        if molecule_id == 9.0:
            molecule = self.excel_file.SO2
        if molecule_id == 1.0:
            molecule = self.excel_file.H2O
        if molecule_id == 2.0:
            molecule = self.excel_file.CO2
        if molecule_id == 7.0:
            molecule = self.excel_file.O2

        self.q = molecule[0]
        self.m = molecule[1]
        self.coeff1_L = molecule[2]
        self.coeff2_L = molecule[3]
        self.coeff3_L = molecule[4]
        self.coeff4_L = molecule[5]
        self.coeff1_M = molecule[6]
        self.coeff2_M = molecule[7]
        self.coeff3_M = molecule[8]
        self.coeff4_M = molecule[9]

    def setLayersParams(self, pressure, temperature):
        '''
        Pressure and temepratures vary according to altitude
        This function is written for convinence to set pressures and temperatures when needed
        These pressures and temperatures will be used for the below functions
        '''
        self.pressure = pressure
        self.temperature = temperature
        self.path_length = 5 #[km] this is fixed parameter

    def doppler_broadening(self):
        '''
        Half width half maximum is defined and returned
        This is only applicable for pressures over 0.1 atm
        '''
        T = self.temperature
        v_center = self.v_center
        m = self.m
        sqrt = np.sqrt(2*spc.k*T/(m*spc.c**2))
        alpha_doppler = v_center*sqrt
        return alpha_doppler

    def lorentz_broadening(self):
        '''
        Half width half maximum is defined and returned
        This is only applicable for pressures below or equal to 0.1 atm
        '''
        P, T = self.pressure, self.temperature
        P_0, T_0 = 1, 296                         #[atm, K]
        alphaLorentz_a = self.alphaLorentz_a
        alphaLorentz_s = self.alphaDoppler_s
        q = self.q
        gamma = self.gamma
        pressureRatio = P/P_0
        tempRatio = T_0/T
        alpha_lorentz = ((1-q)*alphaLorentz_a+q*alphaLorentz_s)*pressureRatio*tempRatio**gamma
        return alpha_lorentz

    def voigt_profile_substitute(self,v):
        '''
        This function defines the shape of the spectral line
        For pressures over 0.1 the shape is a gaussian
        For pressures below or equal to the shape is a lorentzian
        The return parameter f is the shape of the spectral line
        '''
        v_c = self.v_center
        if self.pressure > 0.1:
            alpha = self.lorentz_broadening()
            a = ((v - v_c)/alpha)**2
            f = 1/(np.pi*alpha*(1+a))
            return f
        else:
            alpha = self.doppler_broadening()*1e17
            sigma = alpha*2/2.355
            a = ((v - v_c)/4*sigma)**2
            f = np.exp(-a)/(sigma*np.sqrt(2*np.pi))
            return f

    def line_shape(self, v):
        '''
        Making additional adjustment to the lineshape for the far winged effect
        The return parameter g is the final line shape
        '''
        T = self.temperature
        v_center = self.v_center
        f = self.voigt_profile_substitute(v)
        tanh = np.tanh(spc.h*spc.c*v/(2*spc.k*T))
        tanh_center = np.tanh(spc.h*spc.c * v_center/(2*spc.k*T))
        g = v * tanh * f/(v_center*tanh_center)
        return g

    def total_internal_partition(self, temperature):
        '''
        Total internal partition sum(TIPS) are taken from the paper -
        https://www.sciencedirect.com/science/article/abs/pii/S0022286099002665#:~:text=Total%20internal%20partition%20sums%20
        The paper defines TIPS as a polynominal
        The data from the papaer is stored in the excel file
        Note: Valid only for temeratures between 70-1500K(Low and medium temerature coefficients are only stored on the excel sheet)
        Returns the partition sum.
        '''
        T = temperature

        if T >= 70 and T <= 500:
            a, b, c, d = self.coeff1_L, self.coeff2_L, self.coeff3_L, self.coeff4_L
            Q = a + b*T + c*T**2 + d*T**3
            return Q

        if T > 500 and T <= 1500:
            a, b, c, d = self.coeff1_M, self.coeff2_M, self.coeff3_M, self.coeff4_M
            Q = a + b*T + c*T**2 + d*T**3
            return Q


    def line_intensity(self):
        '''
        Line intensity uses the total partition sum
        Line intensity is returned here
        '''
        T = self.temperature
        v_center = self.v_center
        S_0, T_0 = self.S_0, 296
        Q_T = self.total_internal_partition(T)
        Q_T0 = self.total_internal_partition(296)
        E_l = self.E_l
        exp_lower = np.exp(spc.h*spc.c*E_l/(spc.k*T))
        exp_lower0 = np.exp(spc.h*spc.c*E_l/(spc.k*T_0))
        exp_center = np.exp(spc.h*spc.c*v_center/(spc.k*T))
        exp_center0 = np.exp(spc.h*spc.c*v_center/(spc.k*T))
        S = S_0*Q_T0*exp_lower*exp_center/(Q_T*exp_lower0*exp_center0)
        return S

    def transmittance(self, v):
        '''
        Using the functions and parameters defined above the transmittance is finally calculated
        Transmittance formula used here is given by Beer-Lambert Law
        The function returns transmittance value for a given wavenumber
        '''
        P = self.pressure
        T = self.temperature
        x = self.path_length
        q = self.q
        g = self.line_shape(v)
        S = self.line_intensity()
        numerator = -q*P*x*S*g
        denominator = spc.Boltzmann * T
        tau = np.exp(numerator/denominator)
        return tau

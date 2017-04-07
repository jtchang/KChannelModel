""" @package MultiCompartmentNeuronHelpers
    This package contains helper functions for handling data output from NEURON simulations as found in Chang & Higley (2016).
    

"""
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import ipywidgets as widgets
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pandas
from scipy.integrate import simps
import scipy

def import_GABARun_files(path, run_id):
    """ Loads the raw data for uninhibited bAP, inhibited bAP, and inhibition alone into dataframes stored in a dictionary.
    
        Args: 
					path (str):File Path
					run_id (str):Run ID #
       Returns: 
					dictionary (dict): dictionary of the raw data in pandas dataframes.
    """
    
    dictionary={}
    dictionary["uninhibitedbAP"]=[]
    dictionary["inhibitedbAP"]=[]
    dictionary["inhibited"]=[]
    
    pathLoad= path % run_id
        
        
    filename=pathLoad+"controlbapShaft.dat"
    dictionary["uninhibitedbAP"]=pandas.read_csv(filename, " ", header=None, index_col=0, names=file_column_headers())
        
    filename=pathLoad+"inhibbapShaft.dat"
    dictionary["inhibitedbAP"]=pandas.read_csv(filename, " ", header=None, index_col=0, names=file_column_headers())    
       
    filename=pathLoad+"inhibShaft.dat"
    dictionary["inhibited"]=pandas.read_csv(filename, " ", header=None, index_col=0, names=file_column_headers())    
            
            
    return dictionary

def batch_load_files(path, start, stop):
    """ Loads raw data for Run #'s start to stop. and stores these in an array.
    
        Args: 
					path (str): File Path 
					start (int): Start Run Id # 
					stop (int): Stop Run Id #
        Returns: 
					dictionary (array): Returns the array of the raw data dictionarys for the associated runs.
    """
    length=(stop-start)+1
    dictionary=[None]* length
   
    
    for i in range(0, length):
        dictionary[i]=import_GABARun_files(path, i+start)
    
    
    return dictionary

def file_column_headers():
    """ Defines the column headers for .dat files from NEURON simulations 
    
        Args: 
        Returns: 
					(array): Returns the column headers for the raw data found in the .dat files.
    """
    return ["time",
            "soma_v",
            "dendrite_v",
            "spine_v",
            "dendrite_ca",
            "spine_ca",
            "dendrite_na",
            "dendrite_kv",
            "dendrite_kad",
            "dendrite_pas",
            "dendrite_syn",
            "dendrite_ica"]

def calculate_calcium_inhibition(Uninhibited, Inhibited, time, start, stop, blstart, blstop):
    """ Calculates Calcium Inhibition, Calcium Peak Uninhibited, Calcium Peak Inhibited, and baseline Calcium.
    
        Args: 
					Uninhibited (array): Raw AP Ca 
					Inhibited (array): Raw Inhib-AP Ca
					time (array): time data
					start (double): start time
					stop (double): stop time
					blstart (double): baseline start time
					blstop (double): baseline stop time
        Returns: 
					(array): Array of double values for Calcium Inhibition, Calcium Peak Uninhibited, Calcium Peak Inhibited, and baseline Calcium
    """
    start_index = time.tolist().index(start)
    stop_index = time.tolist().index(stop)
    blstart_index= time.tolist().index(blstart)
    blstop_index= time.tolist().index(blstop)
    
    baseline_calcium= np.mean(Uninhibited.as_matrix()[blstart_index:blstop_index])
    U= Uninhibited.as_matrix()[start_index:stop_index] - baseline_calcium
    I=Inhibited.as_matrix()[start_index:stop_index] - baseline_calcium
    U_area=simps(U)
    I_area=simps(I)
    percent_cainh= (I_area-U_area)/U_area * 100
    U_capk= U.max()
    I_capk= I.max()
    
    return [percent_cainh, U_capk, I_capk, baseline_calcium]

def calculate_voltage_inhibition(Uninhibited, Inhibited, time, start, stop, blstart, blstop):
    """ Calculates Voltage Inhibition, Voltage Peak Uninhibited, Voltage Peak Inhibited, and baseline Vm.
    
         Args: 
					Uninhibited (array): Raw AP V
					Inhibited (array): Raw Inhib-AP V
					time (array): time data
					start (double): start time
					stop (double): stop time
					blstart (double): baseline start time
					blstop (double): baseline stop time
        Returns: 
					array: Voltage Inhibition, Voltage Peak Uninhibited, Voltage Peak Inhibited, and baseline Vm.
    """
    start_index = time.tolist().index(start)
    stop_index = time.tolist().index(stop)
    blstart_index= time.tolist().index(blstart)
    blstop_index= time.tolist().index(blstop)
    
    baseline_v= np.mean(Uninhibited.as_matrix()[blstart_index:blstop_index])
    U= Uninhibited.as_matrix()[start_index:stop_index] - baseline_v
    I=Inhibited.as_matrix()[start_index:stop_index] - baseline_v
    U_vpk= Uninhibited.as_matrix()[start_index:stop_index].max()
    I_vpk= Inhibited.as_matrix()[start_index:stop_index].max()
    
    percent_vinh= (I_vpk-U_vpk)/U_vpk * 100
    
    return [percent_vinh, U_vpk, I_vpk, baseline_v]


def batch_cainh(run_dict, start, stop, blstart, blstop):
    """ Calculates Calcium Inhibition, Calcium Peak Uninhibited, Calcium Peak Inhibited, and baseline Calcium for multiple runs. Calculates Voltage Inhibition, Voltage Peak Uninhibited, Voltage Peak Inhibited, and baseline Vm for multiple runs.
    
        Args: 
					run_dict (dict): Dictionary that include rawdata
					start (double): start time
					stop (double): stop time
					blstart (double): baseline start time
					blstop (double): baseline stop time      
        Returns:
					null
					Adds dictionary entries to run_dict.
    """
    run_array=run_dict["RawData"]
    cainh_array=[None]*len(run_array)
    caUpk_array=[None]*len(run_array)
    caIpk_array=[None]*len(run_array)
    cabl_array=[None]*len(run_array)
    vinh_array=[None]*len(run_array)
    vUpk_array=[None]*len(run_array)
    vIpk_array=[None]*len(run_array)
    vbl_array=[None]*len(run_array)
    
    for i in range(0, len(run_array)):
        U=run_array[i]["uninhibitedbAP"]["dendrite_ca"]
        I=run_array[i]["inhibitedbAP"]["dendrite_ca"]
        T=run_array[i]["uninhibitedbAP"]["time"]
        [cainh_array[i],caUpk_array[i],caIpk_array[i],cabl_array[i] ]= calculate_calcium_inhibition(U, I, T, start, stop, blstart, blstop)
        U=run_array[i]["uninhibitedbAP"]["dendrite_v"]
        I=run_array[i]["inhibitedbAP"]["dendrite_v"]
        T=run_array[i]["uninhibitedbAP"]["time"]
        [vinh_array[i],vUpk_array[i],vIpk_array[i],vbl_array[i] ]= calculate_voltage_inhibition(U, I, T, start, stop, blstart, blstop)
                  
                       
    run_dict["CaInh"] =cainh_array
    run_dict["CaPkU"]=caUpk_array
    run_dict["CaPkI"]=caIpk_array
    run_dict["CaBl"]=cabl_array
    run_dict["VInh"]=vinh_array
    run_dict["VPkU"]=vUpk_array
    run_dict["VPkI"]=vIpk_array
    run_dict["VBl"]=vbl_array
    
def sign(number):
    """ Compares number to zero, if greater returns 1 if less returns -1
    
        Args: 
				number (num): takes a number
        Returns: 
					int: +1 or -1
    """
    if number > 0:
        return 1
    else:
        return -1
    

def AP_FWHM(X,Y):
    """ Finds the FWHM of Y series data in X series data. Linear interpolation is used to calculate the FWHM
    
        Args: 
					X (array): X series data
					Y (array): Y series data
        Returns: 
					FWHM (int): FWHM of X series data
    """
    
    half_max = max(Y) / 2.

    d=[None]*len(Y)
    for i in range(1, len(Y)):
        d[i]= half_max - Y[i]

    j=0
    for i in range(1, len(Y)-1):
        if sign(d[i]) != sign(d[i+1]):
            if j==0:
                j=1
                
                y_t=(d[i+1]-d[i])/(X[i+1]-X[i])
                x0= d[i]/y_t
                left_x= x0+X[i]
            else:
                                
                y_t=(d[i+1]-d[i])/(X[i+1]-X[i])
                x0= d[i]/y_t
                right_x= x0+X[i]
        else:
            d[i]=0
    return right_x-left_x #return the difference (full width)

def batch_AP_width(run):
    """ Batch computes action-potential FWHM for dendritic compartment
    
        Args: 
					run (dict): Run dictionary of raw data
        Output: 
					null
					Updates run dictionary with AP width values (AP_widthU, AP_widthI)
    """
    run_array=run["RawData"]
    
    run["AP_widthU"]=np.array([None]*len(run["RawData"]))
    run["AP_widthI"]=np.array([None]*len(run["RawData"]))
    run["AP_widthDelta"]=np.array([None]*len(run["RawData"]))
    for i in range(0, len(run["RawData"])):
        run["AP_widthU"][i]=AP_FWHM(run_array[i]["uninhibitedbAP"].time.as_matrix(),run_array[i]["uninhibitedbAP"].dendrite_v.as_matrix()-run["VBl"][i])
        run["AP_widthI"][i]=AP_FWHM(run_array[i]["inhibitedbAP"].time.as_matrix(),run_array[i]["inhibitedbAP"].dendrite_v.as_matrix()-run["VBl"][i])
        run["AP_widthDelta"][i]=  run["AP_widthI"][i]- run["AP_widthU"][i] 
 
    
def ica_area(time,ICA, start, stop):
    """ Computes calcium current measures for dendritic compartment over period start to stop
        
        Args: 
					time (array): time data
					ICA (array): calcium data
					start (double): start time
					stop (double): stop time
        Returns: 
				 array: Array of flux (area), peak current (peak), and fwhm (width) 
    """
    start_index=time.tolist().index(start)
    stop_index= time.tolist().index(stop)
    
    area=simps(ICA[start_index:stop_index], x=time[start_index:stop_index])
    peak=ICA[start_index:stop_index].min()
    width=AP_FWHM(time,-ICA)
    return [area, peak, width]
    
def batch_ica_area(run, start, stop):
    """ Batch computes total calcium current flux, peak, and fwhm over window start to stop for inhibited and noninhibited case
        
        Args: 
					run (dict): Run dictionary of raw data
					start (double): start time
					stop (double): stop time
        Returns: 
					null
					Updates run dictionary with calcium flux (ICa_AreahU, AP_AreaI), calcium peak (ICa_PeakU, ICa_PeakI), and fwhm (ICa_WidthU, ICa_WidthI).
    
    """
    run_array=run["RawData"]
    
    run["ICa_PeakU"]=np.array([None]*len(run["RawData"]))
    run["ICa_PeakI"]=np.array([None]*len(run["RawData"]))
    run["ICa_AreaU"]=np.array([None]*len(run["RawData"]))
    run["ICa_AreaI"]=np.array([None]*len(run["RawData"]))
    run["ICa_AreaDiff"]=np.array([None]*len(run["RawData"]))
    run["ICa_WidthU"]=np.array([None]*len(run["RawData"]))
    run["ICa_WidthI"]=np.array([None]*len(run["RawData"]))
    for i in range(0, len(run["RawData"])):
        [run["ICa_AreaU"][i],run["ICa_PeakU"][i], run["ICa_WidthU"][i]]=ica_area(run_array[i]["uninhibitedbAP"].time.as_matrix(),
                                  run_array[i]["uninhibitedbAP"].dendrite_ica.as_matrix(), 
                                  start, 
                                  stop)
        [run["ICa_AreaI"][i],run["ICa_PeakI"][i], run["ICa_WidthI"][i]]=ica_area(run_array[i]["inhibitedbAP"].time.as_matrix(),
                                  run_array[i]["inhibitedbAP"].dendrite_ica.as_matrix(), 
                                  start, 
                                  stop)     
    
    

def load_ena_data():
    """ Loads ENA data from NEURON modeling run found in 'data' directory
    
        Args: 
					None
        Returns:
					ENA_test (dict): Dictionary of ENA data
    """
    
    path_super="data/GABAg_2em3/ena/ena_%s/shift_%s/"
    ena_start=50
    ena_stop=75
    ENA_test={}
    ENA_test["car_1"]={}

    path=path_super % ('%d', 'car_1')
    ENA_test["car_1"]["RawData"]=batch_load_files(path, ena_start, ena_stop)

    batch_cainh(ENA_test["car_1"], 500, 600, 450, 499)
    ENA_test["car_1"]["ena"]=np.arange(ena_start, ena_stop+1)
    batch_AP_width(ENA_test["car_1"])
    batch_ica_area(ENA_test["car_1"], 500, 520)

    return ENA_test

def load_gk_data():
    """ Loads GK data from NEURON modeling run found in 'data' directory
    
        Args: 
					None
        Returns: 
					gk_test (dict): Dictionary of GK data
    """
    path_super="data/GABAg_5em4/gk/gk_%s/shift_%s/"
    gk_start=1

    gk_stop=10
    gk_test={}
    gk_test["car_1"]={}

    path=path_super % ('%d', 'car_1')
    gk_test["car_1"]["RawData"]=batch_load_files(path, gk_start, gk_stop)
    batch_cainh(gk_test["car_1"], 500, 600, 450, 499)
    batch_AP_width(gk_test["car_1"])
    gk_test["car_1"]["gk"]=np.arange(gk_start, gk_stop+1)*24/99 
    batch_ica_area(gk_test["car_1"], 500, 520)
    
    return gk_test



def make_figure(ena, gk):
    """ Generates Figures 6 found in Chang & Higley (2017) and saves pdfs to '/Figures' folder.
        
        Args: 
					ena (dict): ena data
					gk (dict): gk data
        Returns: 
					fig1: Fig6 handle
    """
    # Create Our Figure Layout
    fig1 = plt.figure(1, figsize=[9,9])
    ax1=plt.subplot2grid((3,12), (0,0), colspan=4)
    ax2=plt.subplot2grid((3,12), (0,4), colspan=4)
    ax3=plt.subplot2grid((3,12), (0,8), colspan=4)

    ax4=plt.subplot2grid((3,12), (1,0), colspan=3)
    ax5=plt.subplot2grid((3,12), (1,3), colspan=3)

    ax6=plt.subplot2grid((3,12), (1,6), colspan=3)
    ax7=plt.subplot2grid((3,12), (1,9), colspan=3)
    
    title="ENA= %d mV \n gK=%f mS/cm2 \n Ca Inh= %f (%%)"
    gk_index=0
    ax1.plot(gk["car_1"]["RawData"][gk_index]["uninhibitedbAP"].time,gk["car_1"]["RawData"][gk_index]["uninhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][gk_index]["uninhibitedbAP"].dendrite_ica.min())
    ax1.plot(gk["car_1"]["RawData"][gk_index]["inhibitedbAP"].time,gk["car_1"]["RawData"][gk_index]["inhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][gk_index]["uninhibitedbAP"].dendrite_ica.min()
    )
    ax1.set_title(title % (ena["car_1"]["ena"][10],gk["car_1"]["gk"][gk_index],gk["car_1"]["CaInh"][gk_index] ))
    ax1.set_xlabel("time (ms)")
    ax1.set_ylabel("Normalized $I_{Ca^{2+}}$ (a.u.)")
    ax1.set_xlim(498,510)
    ax1.set_ylim(0,1)
    format_axis(ax1)
    
    ax2.plot(gk["car_1"]["RawData"][40]["uninhibitedbAP"].time,gk["car_1"]["RawData"][40]["uninhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][40]["uninhibitedbAP"].dendrite_ica.min())
    ax2.plot(gk["car_1"]["RawData"][40]["inhibitedbAP"].time,gk["car_1"]["RawData"][40]["inhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][40]["uninhibitedbAP"].dendrite_ica.min())
    ax2.set_title(title % (ena["car_1"]["ena"][10],gk["car_1"]["gk"][40],gk["car_1"]["CaInh"][40] ))
    ax2.set_xlabel("time (ms)")
    ax2.set_ylabel("Normalized $I_{Ca^{2+}}$ (a.u.)")
    ax2.set_xlim(498,510)
    ax2.set_ylim(0,1)
    format_axis(ax2)
    
    
    ena_index=25
    
    ax3.plot(ena["car_1"]["RawData"][ena_index]["uninhibitedbAP"].time,ena["car_1"]["RawData"][ena_index]["uninhibitedbAP"].dendrite_ica/ena["car_1"]["RawData"][ena_index]["uninhibitedbAP"].dendrite_ica.min())
    ax3.plot(ena["car_1"]["RawData"][ena_index]["inhibitedbAP"].time,ena["car_1"]["RawData"][ena_index]["inhibitedbAP"].dendrite_ica/ena["car_1"]["RawData"][ena_index]["uninhibitedbAP"].dendrite_ica.min())
    ax3.set_title(title  % (ena["car_1"]["ena"][ena_index],gk["car_1"]["gk"][40],ena["car_1"]["CaInh"][ena_index] ))
    ax3.set_xlabel("time (ms)")
    ax3.set_ylabel("Normalized $I_{Ca^{2+}}$ (a.u.)")
    ax3.set_xlim(498,510)
    ax3.set_ylim(0,1)
    format_axis(ax3)
    
    ax4.plot(ena["car_1"]["ena"],ena["car_1"]["VPkU"], 'k')
    ax4.plot(ena["car_1"]["ena"],ena["car_1"]["VPkI"], 'r')
    ax4.plot(ena["car_1"]["ena"],ena["car_1"]["VBl"], 'g')
    ax4.set_xlabel("$E_{Na}$ (mV)")
    ax4.set_ylabel("Peak Voltage (mV)")
    ax4.set_xlim(49,76)
    ax4.set_ylim(-70,10)
    format_axis(ax4)
    
    
    ax5.plot(ena["car_1"]["ena"], ena["car_1"]["CaInh"], 'ko', ms=5)
    ax5.set_xlabel("$E_{Na}$ (mV)")
    ax5.set_ylabel("Calcium Inhibition (%)")
    ax5.set_xlim(49,76)
    ax5.set_ylim(-30, -15)
    format_axis(ax5)
    
    
    ax6.plot(gk["car_1"]["gk"],gk["car_1"]["VPkU"], 'k')
    ax6.plot(gk["car_1"]["gk"],gk["car_1"]["VPkI"], 'r')
    ax6.plot(gk["car_1"]["gk"],gk["car_1"]["VBl"], 'g')
    ax6.set_xlabel("$\\bar{g}_{KAD}$ ($mS/cm^2$)")
    ax6.set_ylabel("Peak Voltage (mV)")
    ax6.set_xlim(2,25)
    ax6.set_ylim(-70,10)
    format_axis(ax6)
    
    ax7.plot(gk["car_1"]["gk"], gk["car_1"]["CaInh"], 'wo', ms=5)
    ax7.set_xlabel("$\\bar{g}_{KAD}$ ($mS/cm^2$)")
    ax7.set_ylabel("Calcium Inhibition(%)")
    ax7.set_xlim(2,25)
    ax7.set_ylim(-30,-15)
    format_axis(ax7)
    
   
    ax1a=plt.subplot2grid((3,12), (2,1), colspan=4)
    ax2a=plt.subplot2grid((3,12), (2,7), colspan=4)

    ax1a.plot(gk["car_1"]["VPkU"],gk["car_1"]["CaInh"], 'wo', ms=5)
    ax1a.plot(ena["car_1"]["VPkU"],ena["car_1"]["CaInh"], 'ko', ms=5)
    ax1a.set_xlabel("Uninhibited AP Peak Voltage (mV)")
    ax1a.set_ylabel("Calcium Inhibition(%)")
    ax1a.set_ylim(-30, -15) 
    format_axis(ax1a)
    
    ax2a.plot(gk["car_1"]["AP_widthU"],gk["car_1"]["CaInh"], 'wo', ms=5, label="$\\bar{g}_{KAD}$")
    ax2a.plot(ena["car_1"]["AP_widthU"],ena["car_1"]["CaInh"], 'ko', ms=5, label="$E_{Na}$")
    ax2a.set_xlabel("AP FWHM (ms)")
    ax2a.set_ylabel("Calcium Inhibition(%)")
    ax2a.set_ylim(-30, -15) 
    format_axis(ax2a)
    
    ax2a.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints=1, frameon=False)
    
    fig1.tight_layout()
    fig1.savefig('Figures/Figure6.pdf', format ='pdf')
    #fig2.savefig('Figures/Figure6-supp1.pdf', format ='pdf')
    return fig1

def format_axis(ax):
    """ Simple Helper Function that formats axis to only display left and bottom
			Args:
				ax (matplotlib.pyplot.axis): axis handle
			Returns:
				null
		
		"""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
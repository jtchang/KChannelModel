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
import MultiCompartmentNeuronHelpers as mcnh

def genFig(ena, gk, figname):
    # Create Our Figure Layout
  fig1 = plt.figure(1, figsize=[9,4])
  ax1=plt.subplot2grid((1,12), (0,0), colspan=4)
  ax2=plt.subplot2grid((1,12), (0,4), colspan=4)
  ax3=plt.subplot2grid((1,12), (0,8), colspan=4)

  # ax4=plt.subplot2grid((3,12), (1,0), colspan=6)
  #ax5=plt.subplot2grid((3,12), (2,0), colspan=6)

  #ax6=plt.subplot2grid((3,12), (1,6), colspan=6)
  #ax7=plt.subplot2grid((3,12), (2,6), colspan=6)

  title="ENA= %d mV \n gK=%f mS/cm2 \n Ca Inh= %f (%%)"
  gk_index=0
  ax1.plot(gk["car_1"]["RawData"][gk_index]["uninhibitedbAP"].time,gk["car_1"]["RawData"][gk_index]["uninhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][gk_index]["uninhibitedbAP"].dendrite_ica.min())
  ax1.plot(gk["car_1"]["RawData"][gk_index]["inhibitedbAP"].time,gk["car_1"]["RawData"][gk_index]["inhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][gk_index]["inhibitedbAP"].dendrite_ica.min())
  ax1.set_title(title % (ena["car_1"]["ena"][10],gk["car_1"]["gk"][gk_index],gk["car_1"]["CaInh"][gk_index] ))
  ax1.set_xlabel("time (ms)")
  ax1.set_ylabel("NORM $I_{Ca^{2+}}$ (a.u.)")
  ax1.set_xlim(498,600)
  ax1.set_ylim(0,1)
  mcnh.format_axis(ax1)

  ax2.plot(gk["car_1"]["RawData"][30]["uninhibitedbAP"].time,gk["car_1"]["RawData"][30]["uninhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][30]["uninhibitedbAP"].dendrite_ica.min())
  ax2.plot(gk["car_1"]["RawData"][30]["inhibitedbAP"].time,gk["car_1"]["RawData"][30]["inhibitedbAP"].dendrite_ica/gk["car_1"]["RawData"][30]["inhibitedbAP"].dendrite_ica.min())
  ax2.set_title(title % (ena["car_1"]["ena"][10],gk["car_1"]["gk"][30],gk["car_1"]["CaInh"][30] ))
  ax2.set_xlabel("time (ms)")
  ax2.set_ylabel("NORM $I_{Ca^{2+}}$ (a.u.)")
  ax2.set_xlim(498,600)
  ax2.set_ylim(0,1)
  mcnh.format_axis(ax2)


  ena_index=25

  ax3.plot(ena["car_1"]["RawData"][ena_index]["uninhibitedbAP"].time,ena["car_1"]["RawData"][ena_index]["uninhibitedbAP"].dendrite_ica/ena["car_1"]["RawData"][ena_index]["uninhibitedbAP"].dendrite_ica.min())
  ax3.plot(ena["car_1"]["RawData"][ena_index]["inhibitedbAP"].time,ena["car_1"]["RawData"][ena_index]["inhibitedbAP"].dendrite_ica/ena["car_1"]["RawData"][ena_index]["inhibitedbAP"].dendrite_ica.min())
  ax3.set_title(title  % (ena["car_1"]["ena"][ena_index],gk["car_1"]["gk"][30],ena["car_1"]["CaInh"][ena_index] ))
  ax3.set_xlabel("time (ms)")
  ax3.set_ylabel("NORM $I_{Ca^{2+}}$ (a.u.)")
  ax3.set_xlim(498,600)
  ax3.set_ylim(0,1)
  mcnh.format_axis(ax3)


  fig1.tight_layout()
  fig1.savefig('Figures/'+figname+'.pdf', format ='pdf')

def genArbFig(ena, gk,x, y, figname):
    fig1 = plt.figure(1, figsize=[4,4])
    ax1=plt.subplot2grid((1,1), (0,0), colspan=4)
    #ax2=plt.subplot2grid((1,12), (0,4), colspan=4)
    mec=['k', 'r']

    for i in range(0,len(x)):
        #   ax1.plot(gk["car_1"][x[i]],gk["car_1"][y[i]], 'ko', markeredgecolor=mec[i])
        ax1.plot(ena["car_1"][x[i]],ena["car_1"][y[i]], 'wo',markeredgecolor=mec[i])
        ax1.set_xlabel(x[i])
        ax1.set_ylabel(y[i])
        #ax2.set_xlabel(x[i])
        #ax2.set_ylabel(y[i])
    fig1.tight_layout()
    fig1.savefig('Figures/'+figname+'.pdf', format ='pdf')

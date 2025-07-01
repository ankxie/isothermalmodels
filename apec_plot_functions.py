"""
This code contains functions for reading in .dat files with isothermal models from xspec, and plotting the data. 
These functions are modified from the basic_plot.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp

tnr="Times New Roman"
tnrb="Times New Roman Bold"

xmax = 30

def read_dat_file_spacemodel(filename):
    """
    Specialized function for reading .dat files and outputting a 3-tiered
    list (described below).

    Specially for reading in a data in model space.
    """
    datfile = open(filename, "r")
    line = datfile.readline().strip()
    
    #skipping header
    while line[0] not in "0123456789":
        line = datfile.readline().strip()
    
    """
    a 2-tiered list:
    tier 1: the data set (MEG or HEG)
    tier 2: 0   central wavelength in angstroms
            1   half-width of bin
            2   model, in photons/sec/cm^2/angstroms
    """
    
    data = []
    eg = [[],[],[]] #to hold MEG, then HEG, as we go through appending
    while line != "": #until end of file, keep reading

#if we find a divider (the NO NO NO NO NO line), then append MEG data to data list and reset eg
#to hold HEG data
        if line[0] == "N":
            data.append(eg)
            eg = [[],[],[]]
        
        else:
            line = line.split(" ") #otherwise, split line into list of strings
            eg[0].append(float(line[0])) #append float of wavelength to wavelength list
            eg[1].append(float(line[1])) #append half bin width
            eg[2].append(float(line[2])) #append model

        
        #grab the next line
        line = datfile.readline().strip()
    
    #at end of file, append HEG data to data list
    data.append(eg)

    if data[0][0][0] > data[0][0][1]: #if MEG data in order of decreasing wavelength
        for i in range(3):
            data[0][i].reverse() #reverse order
    if data[1][0][0] > data[1][0][1]: #if HEG data in order of decreasing wavelength
        for i in range(3):
            data[1][i].reverse() #reverse order
    
    return data

def read_dat_file_spacedata(filename):
    """
    Specialized function for reading .dat files and outputting a 3-tiered
    list (described below).

    Specially for reading in a model in data space.
    """
    datfile = open(filename, "r")
    line = datfile.readline().strip()
    
    #skipping header
    while line[0] not in "0123456789":
        line = datfile.readline().strip()
    
    """
    a 2-tiered list:
    tier 1: the data set (MEG or HEG)
    tier 2: 0   central wavelength in angstroms
            1   half-width of bin
            2   data value, in counts/sec/angstrom
            3   formal uncertainty, same units
            4   model 
    """
    
    data = []
    eg = [[],[],[],[],[]] #to hold MEG, then HEG, as we go through appending
    while line != "": #until end of file, keep reading

#if we find a divider (the NO NO NO NO NO line), then append MEG data to data list and reset eg
#to hold HEG data
        if line[0] == "N":
            data.append(eg)
            eg = [[],[],[],[],[]]
        
        else:
            line = line.split(" ") #otherwise, split line into list of strings
            eg[0].append(float(line[0])) #append float of wavelength to wavelength list
            eg[1].append(float(line[1])) #append half bin width
            eg[2].append(float(line[2])) #append data value
            eg[3].append(float(line[3])) #append error bar value
            eg[4].append(float(line[4])) #append model value

        
        #grab the next line
        line = datfile.readline().strip()
    
    #at end of file, append HEG data to data list
    data.append(eg)

    if data[0][0][0] > data[0][0][1]: #if MEG data in order of decreasing wavelength
        for i in range(5):
            data[0][i].reverse() #reverse order
    if data[1][0][0] > data[1][0][1]: #if HEG data in order of decreasing wavelength
        for i in range(5):
            data[1][i].reverse() #reverse order
    
    return data

    
def histogram_data(data, set, space):
    """
    Takes a data set and returns a data set that will produce a histogram.
    This involves, essentially, doubling the size of the model. Instead
    of plotting each value once at the center of the bin, we want to plot
    twice on either side of the bin.
    
    set = 0 if both datasets come from same grating and can be added together, 
    1 to look at MEG data, 2 to look at HEG data

    space = 0 if working in data space
    space = 1 if working in model space
    """
    try:
        if space == 0:
            #each of these has two parts: one for the first data set,
            #one for the second
            histwvs = [[],[]]
            histmods = [[],[]]
            histvals = [[],[]]
            
            for i in range(len(data)):
                for j in range(len(data[i][1])):
                    histwvs[i].append(data[i][0][j]-data[i][1][j]) #append left wavelength of bin
                    histwvs[i].append(data[i][0][j]+data[i][1][j]) #append right wavelength of bin
                    histvals[i].append(data[i][2][j]) #append data value twice. This actually isn't necessary for any of these
                    #plots, so you can delete it if you want.
                    histvals[i].append(data[i][2][j])
                    histmods[i].append(data[i][4][j]) #append model value twice
                    histmods[i].append(data[i][4][j])
            
            if set == 0: #co-adds counts together 
                tot_histwvs = histwvs[0]
                tot_histvals = np.add(histvals[0], histvals[1])
                #tot_histvals = [1, 1]
                tot_histmods = np.add(histmods[0], histmods[1])
                #tot_histmods = [1, 1]
                return tot_histwvs,tot_histvals,tot_histmods
            if set == 1:
                histwvs_MEG = histwvs[0]
                tot_histvals_MEG = histvals[0]
                tot_histmods_MEG = histmods[0]
                return histwvs_MEG,tot_histvals_MEG,tot_histmods_MEG
            if set == 2:
                histwvs_HEG = histwvs[1]
                tot_histvals_HEG = histvals[1]
                tot_histmods_HEG = histmods[1]
                return histwvs_HEG, tot_histvals_HEG, tot_histmods_HEG
            else:
                print("please input 0, 1, or 2 for the set parameter")
        
        if space == 1:
            #each of these has two parts: one for the first data set,
            #one for the second
            histwvs = [[],[]]
            histmods = [[],[]]
            
            for i in range(len(data)):
                for j in range(len(data[i][1])):
                    histwvs[i].append(data[i][0][j]-data[i][1][j]) #append left wavelength of bin
                    histwvs[i].append(data[i][0][j]+data[i][1][j]) #append right wavelength of bin
                    histmods[i].append(data[i][2][j]) #append model value twice
                    histmods[i].append(data[i][2][j])
            
            if set == 0: #adds HEG and MEG counts together 
                tot_histwvs = histwvs[0]
                tot_histmods = np.add(histmods[0], histmods[1])
                #tot_histvals = [1, 1]
                return tot_histwvs,tot_histmods
            if set == 1:
                histwvs_MEG = histwvs[0]
                tot_histmods_MEG = histmods[0]
                return histwvs_MEG,tot_histmods_MEG,
            if set == 2:
                histwvs_HEG = histwvs[1]
                tot_histmods_HEG = histmods[1]
                return histwvs_HEG, tot_histmods_HEG
            else:
                print("please input 0, 1, or 2 for the set parameter")
    
    except:
        print("error: did you use the right read_dat_file?")

def plot_model(data, set, c, space = 0, temp = False, log = False):
    
    """set = 0 if both datasets come from same grating and can be added together, 
    1 to look at MEG data, 2 to look at HEG data

    space = 0 if working in data space
    space = 1 if working in model space
    """

    #plot
    plt.rcParams['figure.figsize'] = (15.0, 8.0)
    plt.figure()
    plt.xlabel("Wavelength (Å)",fontsize=32,fontname=tnr)
    plt.xticks(fontname=tnr)
    plt.yticks(fontname=tnr)
    plt.tick_params(labelsize=24)

    if space == 0:
        histwvs,histvals,histmods = histogram_data(data, set, space)
        plt.ylabel("Counts s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)
    
    if space == 1:
        histwvs,histmods = histogram_data(data, set, space)
        plt.ylabel("Photons cm$^\mathregular{-2}$ s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)


    plt.plot(histwvs, histmods, color = c, label = temp)
    #plt.xlim(np.min(histwvs[j]),np.max(histwvs[j]))
    
    plt.xlim(4.5, xmax)

    if log == True:
        plt.yscale('log')

    if temp != False:
        plt.legend()


def isothermal_model_loop(file_format, n, colors, set, space = 0, temp = False, log = False, save = False):
    """this function inputs three arrays of the same length: 
        an array of integers to retrieve datfiles 
        an array of strings to choose colors
        an array of strings to indicate temperatures
    if save = True, the plots will be saved
    if save = False, the plots will not be saved
    The txt files must have a consistant format

    Outputs plots of isothermal models in model space
    """
    for i in range(len(n)):
        if space == 0:
            data = read_dat_file_spacedata(file_format + str(n[i]) + ".qdp")
            space_type = "spacedata"
        if space == 1:
            data = read_dat_file_spacemodel(file_format + str(n[i]) + ".qdp")
            space_type = "spacemodel"

        if temp == False:
            plot_model(data, set, colors[i], space, log)
        else:
            plot_model(data, set, colors[i], space, temp[i], log)
            plt.legend()
        
        if log == True:
            scale_type = "log"
        else:
            scale_type = "linear"

        if save == True:
            plt.savefig("isothermalmodel_" + str(n[i]) + "_" + space_type + "_" + scale_type + ".png")

def isothermal_model_subplots(file_format, n, colors, set, space = 0, temp = False, log = False, save = False):
    
    #set up figure 
    plt.clf() 
    fig1 = plt.figure(1)
    fig1.set_size_inches(8, 15)
    plt.box(False)

    fig1.supxlabel("Wavelength (Å)",fontsize=32,fontname=tnr, y = 0.05)
    #plt.xticks(fontname=tnr)
    #plt.yticks(fontname=tnr)
    plt.tick_params(labelsize=24, bottom = False, left = False, labelbottom = False, labelleft = False)

    if space == 0:
        fig1.supylabel("Counts s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)
    
    if space == 1:
        fig1.supylabel("Photons cm$^\mathregular{-2}$ s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)

    ax1 = fig1.add_subplot(611)
    ax2 = fig1.add_subplot(612)
    ax3 = fig1.add_subplot(613)
    ax4 = fig1.add_subplot(614)
    ax5 = fig1.add_subplot(615)
    ax6 = fig1.add_subplot(616)

    ax = [ax1, ax2, ax3, ax4, ax5, ax6]

    max = 0
    min = 10e-8 

    for i in range(len(n)):
        if space == 0:
            data = read_dat_file_spacedata(file_format + str(n[i]) + ".qdp")
            histwvs,histvals,histmods = histogram_data(data, set, space)
            space_type = "spacedata"

        if space == 1:
            data = read_dat_file_spacemodel(file_format + str(n[i]) + ".qdp")
            histwvs, histmods = histogram_data(data, set, space)
            space_type = "spacemodel"
        
        if temp == False:
            ax[i].plot(histwvs, histmods, color = colors[i]) 
        else:
            ax[i].plot(histwvs, histmods, color = colors[i], label = temp[i])
            ax[i].legend()

        nmax = np.max(histmods)
        nmin = np.min(histmods)
        if nmax > max:
            max = nmax + (max * 0.5)

    
    #for the weighted isothermal plots, in order to actually see the higher temp spectra
    ax1.set_xlim(4.5, xmax)
    ax2.set_xlim(4.5, xmax)
    ax3.set_xlim(4.5, xmax)
    ax4.set_xlim(4.5, xmax)
    ax5.set_xlim(4.5, xmax)
    ax6.set_xlim(4.5, xmax)

    if log == True:
        #min = 0.1
        #print(min)
        ax1.set_yscale('log')
        ax2.set_yscale('log')
        ax3.set_yscale('log')
        ax4.set_yscale('log')
        ax5.set_yscale('log')
        ax6.set_yscale('log')
        ax1.set_ylim(min, max)
        ax2.set_ylim(min, max)
        ax3.set_ylim(min, max)
        ax4.set_ylim(min, max)
        ax5.set_ylim(min, max)
        ax6.set_ylim(min, max)
        scale_type = "log"
    else:
        #print(max)
        ax1.set_ylim(0, max)
        ax2.set_ylim(0, max)
        ax3.set_ylim(0, max)
        ax4.set_ylim(0, max)
        ax5.set_ylim(0, max)
        ax6.set_ylim(0, max)
        scale_type = "linear"

    if save == True:
        plt.savefig("isothermal_model_subplots" + space_type + "_" + scale_type + ".png")


def isothermal_model_composite(file_format, n, colors, set, space = 0, temp = False, log = False, save = False):
    """this function outputs a composite plot of isothermal models
    """

    plt.rcParams['figure.figsize'] = (15.0, 8.0)
    plt.figure()
    plt.xlabel("Wavelength (Å)",fontsize=32,fontname=tnr)
    plt.xticks(fontname=tnr)
    plt.yticks(fontname=tnr)
    plt.tick_params(labelsize=24)

    for i in range(len(n)):

        if space == 0:
            data = read_dat_file_spacedata(file_format + str(n[i]) + ".qdp")
            histwvs, histvals, histmods = histogram_data(data, set, space)
            plt.ylabel("Counts s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)
            space_type = "spacedata"

        if space == 1:
            data = read_dat_file_spacemodel(file_format + str(n[i]) + ".qdp")
            histwvs,histmods = histogram_data(data, set, space)
            plt.ylabel("Photons cm$^\mathregular{-2}$ s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)
            space_type = "spacemodel"

        if temp == False:
            plt.plot(histwvs, histmods, color = colors[i])
        else:
            plt.plot(histwvs, histmods, color = colors[i], label = temp[i])
            plt.legend()


    plt.xlim(4.5, xmax)
    
    if log == True:
        plt.yscale('log')
        plt.ylim(bottom = 10e-8)
        scale_type = "log"
    else:
        scale_type = "linear"

    if save == True:
        plt.savefig("isothermal_model_composite" + space_type + "_" + scale_type + ".png")

def isothermal_model_combined(file_format, n, set, space = 0, log = False, save = False):
    '''this function takes in different isothermal models and then combines their fluxes to create one spectrum'''

    plt.rcParams['figure.figsize'] = (15.0, 8.0)
    plt.figure()
    plt.xlabel("Wavelength (Å)",fontsize=32,fontname=tnr)
    plt.xticks(fontname=tnr)
    plt.yticks(fontname=tnr)
    plt.tick_params(labelsize=24)

    histmods_list = []

    for i in range(len(n)):

        if space == 0:
            data = read_dat_file_spacedata(file_format + str(n[i]) + ".qdp")
            histwvs, histvals, histmods = histogram_data(data, set, space)
            histmods_list.append(histmods)
            plt.ylabel("Counts s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)
            space_type = "spacedata"

        if space == 1:
            data = read_dat_file_spacemodel(file_format + str(n[i]) + ".qdp")
            histwvs,histmods = histogram_data(data, set, space)
            histmods_list.append(histmods)
            plt.ylabel("Photons cm$^\mathregular{-2}$ s$^\mathregular{-1}$ Å$^\mathregular{-1}$",fontsize=32,fontname=tnr)
            space_type = "spacemodel"

    histmods_combined = np.zeros_like(histmods)
    for i in range(len(histmods_list[0])):
            for j in range(len(histmods_list)):
                histmods_combined[i] += histmods_list[j][i]

    plt.plot(histwvs, histmods_combined, color = 'k')

    plt.xlim(4.5, xmax)
    
    if log == True:
        plt.yscale('log')
        plt.ylim(bottom = 10e-8)
        scale_type = "log"
    else:
        scale_type = "linear"

    if save == True:
        plt.savefig("isothermal_model_composite" + space_type + "_" + scale_type + ".png")

print("importing fin! 1")
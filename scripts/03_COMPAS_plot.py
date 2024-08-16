####### PARENT SCRIPT: tertius_plotting_081224.py
#------> this version applies orbital period masks instead of semimaj axis

# DESCRIPTION:
# ---> Visualizes the processed files extraxcted from the sÄ±mulation output by the previous script.
# ---> Creates occurence or percentage heatmaps for black hole - star bÄ±naries as a function of orbital period and semimajor axis.

import os, sys
import numpy as np               # for handling arrays
import seaborn as sns
import numpy.ma as ma
import pandas as pd
import time                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import matplotlib

sns.set_theme(style="ticks")

c = datetime.now()
# Displays Time
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')
print("current time :", current_time)

if not os.path.exists("/data/a.saricaoglu/Plots/" +  str(c.strftime("%m.%d"))): 
    os.makedirs("/data/a.saricaoglu/Plots/"  +  str(c.strftime("%m.%d"))) 

# Choose the IMF to visualise outputs.
IMFs = ["Kroupa"]
mode = "SingleBH"

# Plug in the date of the output data.
date = "08.15"

# Choose plotting option: occurance rates if occrate = True, percentages if occrate = False.
occrate = True

# Choose the grÄ±d sÄ±ze. large = False : 5 x 5 , large = True : 15 x 15
large = False

if large is False:
    binsp = 5
    bindcop = [np.log10([0.3, 0.8, 1.9, 4.8, 11.9, 30.0]), np.log10([5.0,9.1, 16.6, 30.2, 54.9, 100.0])]
    bindcsemax = 5
    fsize = 40
    lsize = 40 
    fgsize = (20,16)
    sz = 30
else:
    bin = 15
    binsp = 15
    bindcsemax = 15
    bindcop = 15
    fsize = 50
    lsize = 50 
    fgsize = (50,36)
    sz = 30

def mylog(x,y):
    return np.log10(x[(x > 0) & (y > 0)]), np.log10(y[(x > 0) & (y > 0)])

if occrate:
    def percentage(size, vals):
        return vals * (1/size)
    rate = " Occurence Rate"
else:
    def percentage(size, vals):
        return vals * (100/size)
    rate = " Percentage"
    


# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
# Choose an output hdf5 file to work with
for IMF in IMFs:

    pathToData = '/data/a.saricaoglu/Files/' + IMF + date #change the date accordingly the date of the files created via Sec

    # Keeps corresponding numerical values 
    SPs = []
    DCs = []

    STELLARTYPEZAMSSP1   =  []
    STELLARTYPEZAMSSP2   =  []
    STELLARTYPESP1   =  []
    STELLARTYPESP2   =  []       
    STELLARTYPEDC1   =  []
    STELLARTYPEDC2   =  [] 

    MASSZAMSSP1 = []
    MASSZAMSSP2 = []
    MASSDC1 = []
    MASSDC2 = []

    SEMIMAJORAXISSP = []
    SEMIMAJORAXISDC = []

    ORBITALPERIODSP = []
    ORBITALPERIODDC = []

    # Boolean values for masking
    MASKSPBH1 = np.array([], dtype='bool')
    MASKSPBH2 = np.array([], dtype='bool')
    MASKDCBH1 = np.array([], dtype='bool')
    MASKDCBH2 = np.array([], dtype='bool')

    MASKSPunb = np.array([], dtype='bool')
    MASKSPdco = np.array([], dtype='bool')
    MASKSPmrgr = np.array([], dtype='bool')

    MASKDCinSP = np.array([], dtype='bool')
    MASKDCBHNS1 = np.array([], dtype='bool')
    MASKDCBHNS2 = np.array([], dtype='bool')

    MASKDCnonBH =  np.array([], dtype='bool')
    MASKDCSEMAJ = np.array([], dtype='bool')
    MASKDCORBPER = np.array([], dtype='bool')
    MASKDCCOSI = np.array([], dtype='bool')
    MASKDCnegCE = np.array([], dtype='bool')
    MASKDCposCE = np.array([], dtype='bool')

   
    # Dataframes to use with sns plotting package
    data_SP = {}
    SPdf = pd.DataFrame(data_SP)
    data_DC = {}
    DCdf = pd.DataFrame(data_DC)
    
    # Reads files produced by ... and contructs the dataframes. SPs and DCs are separate since they have different sizes.
    arrays_to_save = [SPs,DCs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEDC1,STELLARTYPEDC2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSDC1,MASSDC2,SEMIMAJORAXISSP,SEMIMAJORAXISDC,MASKSPBH1,MASKSPBH2,MASKDCBH1,MASKDCBH2,
                    MASKSPunb,MASKSPdco,MASKSPmrgr,MASKDCinSP,MASKDCBHNS1,MASKDCBHNS2, MASKDCnonBH, MASKDCSEMAJ,MASKDCnegCE, MASKDCposCE, 
                    ORBITALPERIODSP, ORBITALPERIODDC, MASKDCORBPER, MASKDCCOSI]
    filenames = ["SPs","DCs","STELLARTYPEZAMSSP1","STELLARTYPEZAMSSP2","STELLARTYPESP1","STELLARTYPESP2","STELLARTYPEDC1","STELLARTYPEDC2",
                "MASSZAMSSP1","MASSZAMSSP2","MASSDC1","MASSDC2","SEMIMAJORAXISSP","SEMIMAJORAXISDC","MASKSPBH1","MASKSPBH2","MASKDCBH1","MASKDCBH2",
                "MASKSPunb","MASKSPdco","MASKSPmrgr","MASKDCinSP","MASKDCBHNS1","MASKDCBHNS2", "MASKDCnonBH", "MASKDCSEMAJ","MASKDCnegCE","MASKDCposCE",
                "ORBITALPERIODSP", "ORBITALPERIODDC", "MASKDCORBPER", "MASKDCCOSI"]
    i = 0
    runs= [x[2] for x in os.walk(pathToData)][0]
    i=0
    lenSP = 0
    lenDC = 0
    for run in runs:
        print("run :", run)
        for var in filenames:
            print("var :", var)
            print("if check :" , arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]])
            if var in run and len(arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]]) == 0 :
                index = np.where(np.asarray(filenames) == var)[0][0]
                arrays_to_save[np.where(np.asarray(filenames) == var)[0][0]] = 1

        print(pathToData  + run)
        data = np.loadtxt(pathToData  + run)
        print("index: ", index)
        print("File name: ", run, " Variable name: ", filenames[index])
        if len(data) != 0:
            print("data reached.")
        print("len: ", len(data))
        print("df col name:", filenames[index])     
        if index == 0:
            lenSP = len(data)
        if index == 1:
            lenDC = len(data)

        if len(data) == lenSP:

            SPdf[filenames[index]] = data
            print(SPdf[filenames[index]] .isnull().sum(), SPdf[filenames[index]] .max())
            print("SUCCESS (SP)")
            i = i+1

        if len(data) == lenDC:

            DCdf[filenames[index]] = data
            print(DCdf[filenames[index]].isnull().sum(), DCdf[filenames[index]].max())
            print("SUCCESS (DC)")
            i = i+1


    print(SPdf.keys())
    print(DCdf.keys())

    systemSize = lenSP
    print("System size: ", lenSP, " 100/systemsize: ", 100/lenSP)
    DCdf_singleBH = pd.DataFrame()
    DCdf_tot = pd.DataFrame()
    #Systems which primary is BH
    DCdf_singleBH["BH1"] = DCdf["MASSDC1"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]
    DCdf_singleBH["CP2"] = DCdf["MASSDC2"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]
    DCdf_singleBH["SA1"] = DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]
    DCdf_singleBH["OP1"] = DCdf["ORBITALPERIODDC"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]


    #Systems which secondary is BH
    DCdf_singleBH["BH2"] = DCdf["MASSDC2"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]
    DCdf_singleBH["CP1"] = DCdf["MASSDC1"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]
    DCdf_singleBH["SA2"] = DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]
    DCdf_singleBH["OP2"] = DCdf["ORBITALPERIODDC"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]


########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    DCdf_tot["BlackHole"] = np.concatenate([DCdf_singleBH["BH1"], DCdf_singleBH["BH2"]])
    DCdf_tot["Companion"] = np.concatenate([DCdf_singleBH["CP1"], DCdf_singleBH["CP2"]])
    DCdf_tot["Semax"] = np.concatenate([DCdf_singleBH["SA1"] , DCdf_singleBH["SA2"] ])
    DCdf_tot["OrbitalPeriod"] = np.concatenate([DCdf_singleBH["OP1"] , DCdf_singleBH["OP2"] ])

##########
    DCdf_singleBH["BH1_negCE"]= DCdf["MASSDC1"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]
    DCdf_singleBH["CP2_negCE"]= DCdf["MASSDC2"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]
    DCdf_singleBH["SA1_negCE"]= DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]
    DCdf_singleBH["OP1_negCE"] = DCdf["ORBITALPERIODDC"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]

    #Systems which secondary is BH
    DCdf_singleBH["BH2_negCE"]= DCdf["MASSDC2"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]
    DCdf_singleBH["CP1_negCE"]= DCdf["MASSDC1"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]
    DCdf_singleBH["SA2_negCE"]= DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]
    DCdf_singleBH["OP2_negCE"] = DCdf["ORBITALPERIODDC"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    DCdf_tot["BlackHole_negCE"]= np.concatenate([DCdf_singleBH["BH1_negCE"], DCdf_singleBH["BH2_negCE"]])
    DCdf_tot["Companion_negCE"]= np.concatenate([DCdf_singleBH["CP1_negCE"], DCdf_singleBH["CP2_negCE"]])
    DCdf_tot["Semax_negCE"]= np.concatenate([DCdf_singleBH["SA1_negCE"], DCdf_singleBH["SA2_negCE"]])
    DCdf_tot["OrbitalPeriod_negCE"] = np.concatenate([DCdf_singleBH["OP1_negCE"] , DCdf_singleBH["OP2_negCE"] ])

    DCdf_singleBH["BH1_posCE"]= DCdf["MASSDC1"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]
    DCdf_singleBH["CP2_posCE"]= DCdf["MASSDC2"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]
    DCdf_singleBH["SA1_posCE"]= DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]
    DCdf_singleBH["OP1_posCE"] = DCdf["ORBITALPERIODDC"]*DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]

    #Systems which secondary is BH
    DCdf_singleBH["BH2_posCE"]= DCdf["MASSDC2"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]
    DCdf_singleBH["CP1_posCE"]= DCdf["MASSDC1"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]
    DCdf_singleBH["SA2_posCE"]= DCdf["SEMIMAJORAXISDC"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]
    DCdf_singleBH["OP2_posCE"] = DCdf["ORBITALPERIODDC"]*DCdf["MASKDCBH2"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    DCdf_tot["BlackHole_posCE"]= np.concatenate([DCdf_singleBH["BH1_posCE"], DCdf_singleBH["BH2_posCE"]])
    DCdf_tot["Companion_posCE"]= np.concatenate([DCdf_singleBH["CP1_posCE"], DCdf_singleBH["CP2_posCE"]])
    DCdf_tot["Semax_posCE"]= np.concatenate([DCdf_singleBH["SA1_posCE"], DCdf_singleBH["SA2_posCE"]])
    DCdf_tot["OrbitalPeriod_posCE"] = np.concatenate([DCdf_singleBH["OP1_posCE"] , DCdf_singleBH["OP2_posCE"] ])

##########

    SPdf_singleBH = pd.DataFrame()
    SPdf_tot = pd.DataFrame()
     #Systems which primary is BH
    SPdf_singleBH["BH1"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH1"]
    SPdf_singleBH["CP2"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH1"]
    SPdf_singleBH["SA1"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH1"]
    SPdf_singleBH["OP1"] = SPdf["ORBITALPERIODSP"]*SPdf["MASKSPBH1"]


    #Systems which secondary is BH
    SPdf_singleBH["BH2"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH2"]
    SPdf_singleBH["CP1"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH2"]
    SPdf_singleBH["SA2"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH2"]
    SPdf_singleBH["OP2"] = SPdf["ORBITALPERIODSP"]*SPdf["MASKSPBH2"]


########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    SPdf_tot["BlackHole"] = np.concatenate([SPdf_singleBH["BH1"], SPdf_singleBH["BH2"]])
    SPdf_tot["Companion"] = np.concatenate([SPdf_singleBH["CP1"], SPdf_singleBH["CP2"]])
    SPdf_tot["Semax"] = np.concatenate([SPdf_singleBH["SA1"] , SPdf_singleBH["SA2"] ])
    SPdf_tot["OrbitalPeriod"] = np.concatenate([SPdf_singleBH["OP1"] , SPdf_singleBH["OP2"] ])

##########
    
     #Systems which primary is BH (DCs in SP)
    SPdf_singleBH["BH1_DC"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH1"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["CP2_DC"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH1"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["SA1_DC"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH1"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["OP1_DC"] = SPdf["ORBITALPERIODSP"]*SPdf["MASKSPBH1"]*SPdf["MASKDCinSP"]


    #Systems which secondary is BH (DCs in SP)
    SPdf_singleBH["BH2_DC"] = SPdf["MASSZAMSSP2"]*SPdf["MASKSPBH2"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["CP1_DC"] = SPdf["MASSZAMSSP1"]*SPdf["MASKSPBH2"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["SA2_DC"] = SPdf["SEMIMAJORAXISSP"]*SPdf["MASKSPBH2"]*SPdf["MASKDCinSP"]
    SPdf_singleBH["OP2_DC"] = SPdf["ORBITALPERIODSP"]*SPdf["MASKSPBH2"]*SPdf["MASKDCinSP"]


########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary (DCs in SP)
    SPdf_tot["BlackHole_DC"] = np.concatenate([SPdf_singleBH["BH1_DC"], SPdf_singleBH["BH2_DC"]])
    SPdf_tot["Companion_DC"] = np.concatenate([SPdf_singleBH["CP1_DC"], SPdf_singleBH["CP2_DC"]])
    SPdf_tot["Semax_DC"] = np.concatenate([SPdf_singleBH["SA1_DC"] , SPdf_singleBH["SA2_DC"] ])
    SPdf_tot["OrbitalPeriod_DC"] = np.concatenate([SPdf_singleBH["OP1_DC"] , SPdf_singleBH["OP2_DC"] ])

##########
    print("\n Number of systems with orbital period < 30 days AND with orbital inclination (cosi) <  2.228e-6 AND undergoing CE (DC): ", str(np.sum(DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCposCE"])),
          "\n Number of systems with orbital period < 30 days AND with orbital inclination (cosi) <  2.228e-6 AND NOT undergoing CE (DC): ", str(np.sum(DCdf["MASKDCBH1"]*DCdf["MASKDCORBPER"]*DCdf["MASKDCCOSI"]*DCdf["MASKDCnegCE"])),)
    print(max(DCdf_tot["Semax_negCE"]))
    print(max(DCdf_tot["BlackHole_negCE"]))
    print(max(DCdf_tot["OrbitalPeriod_negCE"]))
 
    print(min(DCdf_tot["Semax_negCE"]))
    print(min(DCdf_tot["BlackHole_negCE"]))
    print(min(DCdf_tot["OrbitalPeriod_negCE"]))

    print(max(DCdf_tot["Semax_posCE"]))
    print(max(DCdf_tot["BlackHole_posCE"]))
    print(max(DCdf_tot["OrbitalPeriod_posCE"]))
 
    print(min(DCdf_tot["Semax_posCE"]))
    print(min(DCdf_tot["BlackHole_posCE"]))
    print(min(DCdf_tot["OrbitalPeriod_posCE"]))
    
    print(max(SPdf_tot["Semax_DC"]))
    print(max(SPdf_tot["BlackHole_DC"]))
    print(max(SPdf_tot["OrbitalPeriod_DC"]))
 
    print(min(SPdf_tot["Semax_DC"]))
    print(min(SPdf_tot["BlackHole_DC"]))
    print(min(SPdf_tot["OrbitalPeriod_DC"]))

    # print(len(*mylog(DCdf_tot["Semax"],DCdf_tot["BlackHole"])[0]))
    # print(len(*mylog(DCdf_tot["OrbitalPeriod_CE"],DCdf_tot["BlackHole_CE"])[0]))
    
######## UPDATE 07.11.24
    

    # Produces heatmaps. Only black holes, companions are commented since not needed.
    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax"],SPdf_tot["BlackHole"]),bins=binsp)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Semimajor Axis inital [AU]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)"+ rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_semax_" +  IMF + "_SP_Mp_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax_DC"],SPdf_tot["BlackHole_DC"]), bins=binsp)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_DC": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax_DC": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax_DC", columns="BlackHole_DC", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.set_xticks(hmapplot.get_xticks())
    ax.set_xticklabels(hmapplot.get_xticklabels())
    ax.set_yticks(hmapplot.get_yticks())
    ax.set_yticklabels(hmapplot.get_yticklabels())
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Semimajor Axis inital [AU]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$) - Star ($M_s$) Binary Systems (DCs in SP)"+ rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_semax_" +  IMF + "_SP_Mp_DC_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax"],SPdf_tot["Companion"]), bins=bin)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="Semax", values="frequency")
    # f, ax = plt.subplots(figsize=fgsize)
    # xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=lsize)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=lsize)
    # ax.set_xlabel("Semimajor Axis inital [AU]"+ rate,   fontsize=fsize, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)",   fontsize=fsize, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_semax_" +  IMF + "_SP_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["OrbitalPeriod_DC"],SPdf_tot["Companion_DC"]), bins=bin)
    # df = pd.DataFrame({
    # "Companion_DC": np.repeat(xbins[:-1], len(ybins)-1),
    # "OrbitalPeriod_DC": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion_DC", columns="OrbitalPeriod_DC", values="frequency")
    # f, ax = plt.subplots(figsize=fgsize)
    # xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=lsize)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=lsize)
    # ax.set_xlabel("Semimajor Axis inital [AU]"+ rate,   fontsize=fsize, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs in SP)",   fontsize=fsize, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_semax_" +  IMF + "_SP_Ms_DC_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["Semax_negCE"],DCdf_tot["BlackHole_negCE"]), bins=bindcsemax)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_negCE": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax_negCE": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax_negCE", columns="BlackHole_negCE", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Semimajor Axis [AU]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs without CE)"+ rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_semax_" +  IMF + "_DC_Mp_nonCE" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["Semax_posCE"],DCdf_tot["BlackHole_posCE"]), bins=bindcsemax)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_posCE": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax_posCE": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax_posCE", columns="BlackHole_posCE", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Semimajor Axis [AU]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs with CE)"+ rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_semax_" +  IMF + "_DC_Mp_CE_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    #### ORBITAL PERIOD ON X AXIS

    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["OrbitalPeriod"],SPdf_tot["BlackHole"]),bins=binsp)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "OrbitalPeriod": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="OrbitalPeriod", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Orbital Period inital [days]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)"+ rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_orbperiod_" +  IMF + "_SP_Mp_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["OrbitalPeriod_DC"],SPdf_tot["BlackHole_DC"]), bins=binsp)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_DC": np.repeat(xbins[:-1], len(ybins)-1),
    "OrbitalPeriod_DC": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="OrbitalPeriod_DC", columns="BlackHole_DC", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Orbital Period inital [days]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$) - Star ($M_s$) Binary Systems (DCs in SP)"+ rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_orbperiod_" +  IMF + "_SP_Mp_DC_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["OrbitalPeriod"],SPdf_tot["Companion"]), bins=bin)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "OrbitalPeriod": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="OrbitalPeriod", values="frequency")
    # f, ax = plt.subplots(figsize=fgsize)
    # xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=lsize)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=lsize)
    # ax.set_xlabel("Orbital Period inital [days]"+ rate,   fontsize=fsize, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)",   fontsize=fsize, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_orbperiod_" +  IMF + "_SP_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["OrbitalPeriod_DC"],SPdf_tot["Companion_DC"]), bins=bin)
    # df = pd.DataFrame({
    # "Companion_DC": np.repeat(xbins[:-1], len(ybins)-1),
    # "OrbitalPeriod_DC": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion_DC", columns="OrbitalPeriod_DC", values="frequency")
    # f, ax = plt.subplots(figsize=fgsize)
    # xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=lsize)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=lsize)
    # ax.set_xlabel("Orbital Period inital [days]"+ rate,   fontsize=fsize, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs in SP)",   fontsize=fsize, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_orbperiod_" +  IMF + "_SP_Ms_DC_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["OrbitalPeriod_posCE"],DCdf_tot["BlackHole_posCE"]), bins=bindcop)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_posCE": np.repeat(xbins[:-1], len(ybins)-1),
    "OrbitalPeriod_posCE": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="OrbitalPeriod_posCE", columns="BlackHole_posCE", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Orbital Period [days]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs with CE)" + rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_orbperiod_" +  IMF + "_DC_Mp_CE" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["OrbitalPeriod_negCE"],DCdf_tot["BlackHole_negCE"]), bins=bindcop)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_negCE": np.repeat(xbins[:-1], len(ybins)-1),
    "OrbitalPeriod_negCE": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="OrbitalPeriod_negCE", columns="BlackHole_negCE", values="frequency")
    f, ax = plt.subplots(figsize=fgsize)
    xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=lsize)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=lsize, left=True, labelleft=True)
    ax.set_xlabel("Orbital Period [days]",   fontsize=fsize, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DCs without CE)"+ rate,   fontsize=fsize, pad=30)
    f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_orbperiod_" +  IMF + "_DC_Mp_nonCE_" + current_time + ".png",   bbox_inches='tight')
    plt.close()



    # values, xbins, ybins = np.histogram2d(*mylog(DCdf_tot["OrbitalPeriod_DC"],DCdf_tot["Companion"]), bins=bin)
    # values = percentage(systemSize, values)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "OrbitalPeriod": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="OrbitalPeriod", values="frequency")
    # f, ax = plt.subplots(figsize=fgsize)
    # xlabl = ["{:.1f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl, annot_kws={"size":sz})
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=lsize)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=lsize)
    # ax.set_xlabel("Orbital Period inital [days]"+ rate,   fontsize=fsize, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",   fontsize=fsize, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (DC)",   fontsize=fsize, pad=30)
    # f.savefig("/data/a.saricaoglu/Plots/" + str(c.strftime("%m.%d")) + "/" + mode  + "_orbperiod_" +  IMF + "_DC_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()


    #####

    

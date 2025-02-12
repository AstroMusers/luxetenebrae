####### PARENT SCRIPT: tertius_plotting_071724.py
#------> this version plots for contrained semimajor range < 0.5 AU for MTs
#%%
import os, sys
import numpy as np               # for handling arrays
import seaborn as sns
import numpy.ma as ma
import pandas as pd
import time                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import matplotlib
import seaborn as sns
from astropy.io import fits

def mylog(x, y):
    # Compute the logarithm of x and y
    log_x = np.log10(x[y > 0])
    log_y = np.log10(y[x > 0])
    
    # Filter out NaN values
    mask = ~np.isnan(log_x) & ~np.isnan(log_y)
    log_x_filtered = log_x[mask]
    log_y_filtered = log_y[mask]
    
    return log_x_filtered, log_y_filtered
def percentage(size, vals):
    return vals * (100/size)
sns.set_theme(style="ticks")

c = datetime.now()
# Displays Time
current_time = c.strftime('%H%M')
print("current time :", current_time)

mode = ["Default_WD_Enabled/","Limited_WD_Enabled/"] #"Default/","Limited/", 

path = '/data/a.saricaoglu/repo/COMPAS'
# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
# Choose an output hdf5 file to work with
for mod in mode:
    matplotlib.rcParams['figure.figsize'] = (21,14)
    matplotlib.rcParams['lines.markersize'] = 1
    matplotlib.rcParams['font.size'] = 17
    matplotlib.rcParams['legend.loc'] = "upper right"

    pathToData = path + '/Files/' + mod + "02.12/" #change the date accordingly the date of the files created via Sec

    # Keeps corresponding numerical values 
    SPs = []
    MTs = []
    CEs = []

    STELLARTYPEZAMSSP1   =  []
    STELLARTYPEZAMSSP2   =  []
    STELLARTYPESP1   =  []
    STELLARTYPESP2   =  []       
    STELLARTYPEMT1   =  []
    STELLARTYPEMT2   =  [] 

    MASSZAMSSP1 = []
    MASSZAMSSP2 = []
    MASSMT1 = []
    MASSMT2 = []

    SEMIMAJORAXISZAMSSP = []
    SEMIMAJORAXISMT = []

    # Boolean values for masking
    MASKSPBH1 = np.array([], dtype='bool')
    MASKSPBH2 = np.array([], dtype='bool')
    MASKMTBH1 = np.array([], dtype='bool')
    MASKMTBH2 = np.array([], dtype='bool')

    MASKSPunb = np.array([], dtype='bool')
    MASKSPMTo = np.array([], dtype='bool')
    MASKSPmrgr = np.array([], dtype='bool')

    MASKMTinSP = np.array([], dtype='bool')
    MASKMTBHNS1 = np.array([], dtype='bool')
    MASKMTBHNS2 = np.array([], dtype='bool')

    MASKMTnonBH =  np.array([], dtype='bool')
    MASKMTSEMAJ = np.array([], dtype='bool')
   
    # Dataframes to use with sns plotting package
    data_SP = {}
    SPdf = pd.DataFrame(data_SP)
    data_MT = {}
    MTdf = pd.DataFrame(data_MT)
    
    # Reads files produced by ... and contructs the dataframes. SPs and MTs are separate since they have different sizes.
    arrays_to_save = [SPs,MTs,CEs,STELLARTYPEZAMSSP1,STELLARTYPEZAMSSP2,STELLARTYPESP1,STELLARTYPESP2,STELLARTYPEMT1,STELLARTYPEMT2,
                    MASSZAMSSP1,MASSZAMSSP2,MASSMT1,MASSMT2,SEMIMAJORAXISZAMSSP,SEMIMAJORAXISMT,MASKSPBH1,MASKSPBH2,MASKMTBH1,MASKMTBH2,
                    MASKSPunb,MASKSPMTo,MASKSPmrgr,MASKMTinSP,MASKMTBHNS1,MASKMTBHNS2, MASKMTnonBH, MASKMTSEMAJ]
    filenames = ["SPs","MTs","CEs","STELLARTYPEZAMSSP1","STELLARTYPEZAMSSP2","STELLARTYPESP1","STELLARTYPESP2","STELLARTYPEMT1","STELLARTYPEMT2",
                "MASSZAMSSP1","MASSZAMSSP2","MASSMT1","MASSMT2","SEMIMAJORAXISZAMSSP","SEMIMAJORAXISMT","MASKSPBH1","MASKSPBH2","MASKMTBH1","MASKMTBH2",
                "MASKSPunb","MASKSPMTo","MASKSPmrgr","MASKMTinSP","MASKMTBHNS1","MASKMTBHNS2", "MASKMTnonBH", "MASKMTSEMAJ"]

    fits_data = pathToData + 'secundus.fits'
    with fits.open(fits_data) as hdul:
        hdul.info()
        SP = hdul[1].data
        MT = hdul[2].data
        CE = hdul[3].data
    
    fits_data = pathToData + 'tertius.fits'
    with fits.open(fits_data) as hdul:
        hdul.info()
        SP_mask = hdul[1].data
        MT_mask = hdul[2].data


    events = MT['EVENTSMT']
    print(min(events), max(events), np.mean(events))


    systemSize = len(SP['SPs'])
    print(f"System size: {systemSize}")
    MTdf_pre_singleBH = pd.DataFrame()
    MTdf_pre_tot = pd.DataFrame()
    #Systems which primary is BH
    MTdf_pre_singleBH["BH1"] = MT["MASSPREMT1"]*MT_mask["MASKPREMTBH1"]*MT_mask["MASKPREMTORBPER"]
    MTdf_pre_singleBH["CP2"] = MT["MASSPREMT2"]*MT_mask["MASKPREMTBH1"]*MT_mask["MASKPREMTORBPER"]
    MTdf_pre_singleBH["SA1"] = MT["SEMIMAJORAXISPREMT"]*MT_mask["MASKPREMTBH1"]*MT_mask["MASKPREMTORBPER"]

    #Systems which secondary is BH
    MTdf_pre_singleBH["BH2"] = MT["MASSPREMT2"]*MT_mask["MASKPREMTBH2"]*MT_mask["MASKPREMTORBPER"]
    MTdf_pre_singleBH["CP1"] = MT["MASSPREMT1"]*MT_mask["MASKPREMTBH2"]*MT_mask["MASKPREMTORBPER"]
    MTdf_pre_singleBH["SA2"] = MT["SEMIMAJORAXISPREMT"]*MT_mask["MASKPREMTBH2"]*MT_mask["MASKPREMTORBPER"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    MTdf_pre_tot["BlackHole"] = np.concatenate([MTdf_pre_singleBH["BH1"], MTdf_pre_singleBH["BH2"]])
    MTdf_pre_tot["Companion"] = np.concatenate([MTdf_pre_singleBH["CP1"], MTdf_pre_singleBH["CP2"]])
    MTdf_pre_tot["Semax"] = np.concatenate([MTdf_pre_singleBH["SA1"] , MTdf_pre_singleBH["SA2"] ])


    MTdf_pst_singleBH = pd.DataFrame()
    MTdf_pst_tot = pd.DataFrame()
    #Systems which primary is BH
    MTdf_pst_singleBH["BH1"] = MT["MASSPSTMT1"]*MT_mask["MASKPSTMTBH1"]*MT_mask["MASKPSTMTORBPER"]
    MTdf_pst_singleBH["CP2"] = MT["MASSPSTMT2"]*MT_mask["MASKPSTMTBH1"]*MT_mask["MASKPSTMTORBPER"]
    MTdf_pst_singleBH["SA1"] = MT["SEMIMAJORAXISPSTMT"]*MT_mask["MASKPSTMTBH1"]*MT_mask["MASKPSTMTORBPER"]

    #Systems which secondary is BH
    MTdf_pst_singleBH["BH2"] = MT["MASSPSTMT2"]*MT_mask["MASKPSTMTBH2"]*MT_mask["MASKPSTMTORBPER"]
    MTdf_pst_singleBH["CP1"] = MT["MASSPSTMT1"]*MT_mask["MASKPSTMTBH2"]*MT_mask["MASKPSTMTORBPER"]
    MTdf_pst_singleBH["SA2"] = MT["SEMIMAJORAXISPSTMT"]*MT_mask["MASKPSTMTBH2"]*MT_mask["MASKPSTMTORBPER"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    MTdf_pst_tot["BlackHole"] = np.concatenate([MTdf_pst_singleBH["BH1"], MTdf_pst_singleBH["BH2"]])
    MTdf_pst_tot["Companion"] = np.concatenate([MTdf_pst_singleBH["CP1"], MTdf_pst_singleBH["CP2"]])
    MTdf_pst_tot["Semax"] = np.concatenate([MTdf_pst_singleBH["SA1"] , MTdf_pst_singleBH["SA2"] ])

    change = (MT_mask['MASKTYPECHANGE1'] == 1) | (MT_mask['MASKTYPECHANGE2'] == 1)
    MTdf_pst_tot["BlackHole_changemsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*change, MTdf_pst_singleBH["BH2"]*change])
    MTdf_pst_tot["Companion_changemsk"] = np.concatenate([MTdf_pst_singleBH["CP1"]*change, MTdf_pst_singleBH["CP2"]*change])
    MTdf_pst_tot["Semax_changemsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*change , MTdf_pst_singleBH["SA2"]*change ])  

    first = MT_mask['MASKFIRSTMT']
    MTdf_pre_tot["BlackHole_firstmsk"] = np.concatenate([MTdf_pre_singleBH["BH1"]*first, MTdf_pre_singleBH["BH2"]*first])
    MTdf_pre_tot["Companion_firstmsk"] = np.concatenate([MTdf_pre_singleBH["CP1"]*first, MTdf_pre_singleBH["CP2"]*first])
    MTdf_pre_tot["Semax_firstmsk"] = np.concatenate([MTdf_pre_singleBH["SA1"]*first , MTdf_pre_singleBH["SA2"]*first ]) 

    last = MT_mask['MASKLASTMT']
    MTdf_pst_tot["BlackHole_lastmsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*last, MTdf_pst_singleBH["BH2"]*last])
    MTdf_pst_tot["Companion_lastmsk"] = np.concatenate([MTdf_pst_singleBH["CP1"]*last, MTdf_pst_singleBH["CP2"]*last])
    MTdf_pst_tot["Semax_lastmsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*last , MTdf_pst_singleBH["SA2"]*last ])

    print("Checking for NaN values in MTdf_pst_tot['Semax_lastmsk']:", np.isnan(*mylog(MTdf_pst_tot["Semax_lastmsk"],MTdf_pst_tot["BlackHole_lastmsk"])).any())


    mainsequence_pre = (MT_mask['MASKPREMTBHMS1'] | MT_mask['MASKPREMTBHMS2']) & first
    MTdf_pre_tot["BlackHole_mainpremsk"] = np.concatenate([MTdf_pre_singleBH["BH1"]*mainsequence_pre, MTdf_pre_singleBH["BH2"]*mainsequence_pre])
    MTdf_pre_tot["Companion_mainpremsk"] = np.concatenate([MTdf_pre_singleBH["CP1"]*mainsequence_pre, MTdf_pre_singleBH["CP2"]*mainsequence_pre])
    MTdf_pre_tot["Semax_mainpremsk"] = np.concatenate([MTdf_pre_singleBH["SA1"]*mainsequence_pre , MTdf_pre_singleBH["SA2"]*mainsequence_pre ])

    mainsequence_pst = (MT_mask['MASKPSTMTBHMS1'] | MT_mask['MASKPSTMTBHMS2']) & last
    MTdf_pst_tot["BlackHole_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst, MTdf_pst_singleBH["BH2"]*mainsequence_pst])
    MTdf_pst_tot["Companion_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["CP1"]*mainsequence_pst, MTdf_pst_singleBH["CP2"]*mainsequence_pst])
    MTdf_pst_tot["Semax_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst , MTdf_pst_singleBH["SA2"]*mainsequence_pst ])




##########

    SPdf_singleBH = pd.DataFrame()
    SPdf_tot = pd.DataFrame()
     #Systems which primary is BH
    SPdf_singleBH["BH1"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH1"]
    SPdf_singleBH["CP2"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH1"]
    SPdf_singleBH["SA1"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH1"]

    #Systems which secondary is BH
    SPdf_singleBH["BH2"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH2"]
    SPdf_singleBH["CP1"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH2"]
    SPdf_singleBH["SA2"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH2"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    SPdf_tot["BlackHole"] = np.concatenate([SPdf_singleBH["BH1"], SPdf_singleBH["BH2"]])
    SPdf_tot["Companion"] = np.concatenate([SPdf_singleBH["CP1"], SPdf_singleBH["CP2"]])
    SPdf_tot["Semax"] = np.concatenate([SPdf_singleBH["SA1"] , SPdf_singleBH["SA2"] ])
##########
    
     #Systems which primary is BH (MTs in SP)
    SPdf_singleBH["BH1_MT"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH1"]*SP_mask["MASKMTinSP"]
    SPdf_singleBH["CP2_MT"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH1"]*SP_mask["MASKMTinSP"]
    SPdf_singleBH["SA1_MT"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH1"]*SP_mask["MASKMTinSP"]

    #Systems which secondary is BH (MTs in SP)
    SPdf_singleBH["BH2_MT"] = SP["MASSZAMSSP2"]*SP_mask["MASKSPBH2"]*SP_mask["MASKMTinSP"]
    SPdf_singleBH["CP1_MT"] = SP["MASSZAMSSP1"]*SP_mask["MASKSPBH2"]*SP_mask["MASKMTinSP"]
    SPdf_singleBH["SA2_MT"] = SP["SEMIMAJORAXISZAMSSP"]*SP_mask["MASKSPBH2"]*SP_mask["MASKMTinSP"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary (MTs in SP)
    SPdf_tot["BlackHole_MT"] = np.concatenate([SPdf_singleBH["BH1_MT"], SPdf_singleBH["BH2_MT"]])
    SPdf_tot["Companion_MT"] = np.concatenate([SPdf_singleBH["CP1_MT"], SPdf_singleBH["CP2_MT"]])
    SPdf_tot["Semax_MT"] = np.concatenate([SPdf_singleBH["SA1_MT"] , SPdf_singleBH["SA2_MT"] ])
##########

    print(max(MTdf_pre_tot["Semax"]))
    print(max(MTdf_pre_tot["BlackHole"]))
    print(min(MTdf_pre_tot["Semax"]))
    print(min(MTdf_pre_tot["BlackHole"]))
    

    
######## UPDATE 07.11.24
    

 
    if not os.path.exists(path+ "/Plots/" + mod +  str(c.strftime("%m.%d")+ "/" + current_time + "/") ): 
        os.makedirs(path + "/Plots/"  + mod +  str(c.strftime("%m.%d")+ "/" + current_time + "/") ) 
    directoryp = path + "/Plots/"  + mod +  str(c.strftime("%m.%d") + "/" + current_time + "/")  
    # Produces heatmaps. Only black holes, companions are commented since not needed.
    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax"], SPdf_tot["BlackHole"]), bins=20)
    values = percentage(systemSize, values)

    # Check the shape of values to ensure it matches the expected dimensions
    print("Shape of values:", values.shape)
    print("Length of xbins:", len(xbins))
    print("Length of ybins:", len(ybins))

    df = pd.DataFrame({
        "BlackHole": np.repeat(xbins[:-1], len(ybins) - 1),
        "Semax": np.tile(ybins[:-1], len(xbins) - 1),
        'frequency': values.T.flatten()
    })

    # Check the first few rows of the DataFrame to ensure the data is correct
    print(df.head())

    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")

    # Check the shape of the pivot table to ensure it matches the expected dimensions
    print("Shape of pivot table:", vals.shape)

    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]

    # Check the labels to ensure they are correct
    print("X-axis labels:", xlabl)
    print("Y-axis labels:", ylabl)

    hmapplot = sns.heatmap(vals, annot=True, fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]", fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]", fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)", fontsize=40, pad=30)
    f.savefig(directoryp + "SP_Mp_" + current_time + ".png", bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax_MT"],SPdf_tot["BlackHole_MT"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole_MT": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax_MT": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax_MT", columns="BlackHole_MT", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (MTs in SP)",  fontsize=40, pad=30)
    f.savefig(directoryp + "SP_Mp_MT_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax"],SPdf_tot["Companion"]), bins=20)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="Semax", values="frequency")
    # f, ax = plt.subplots(figsize=(36, 32))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (SP)",  fontsize=40, pad=30)
    # f.savefig(directoryp + "SP_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax_MT"],SPdf_tot["Companion_MT"]), bins=20)
    # df = pd.DataFrame({
    # "Companion_MT": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax_MT": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion_MT", columns="Semax_MT", values="frequency")
    # f, ax = plt.subplots(figsize=(36, 32))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (MTs in SP)",  fontsize=40, pad=30)
    # f.savefig(directoryp + "SP_Ms_MT_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax"],MTdf_pre_tot["BlackHole"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title(f"Black Hole ($M_p$)  - Star ($M_s$) Binary Systems(MTO) N:{len(df['BlackHole'])}",  fontsize=40, pad=30)
    f.savefig(directoryp + "MT_All_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    #%%
    
    len(change)
    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax_firstmsk"],MTdf_pre_tot["BlackHole_firstmsk"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title(f"Black Hole($M_p$)-Star($M_s$) Binaries before first MT N:{len(df['BlackHole'])}",  fontsize=40, pad=30)
    f.savefig(directoryp + "MT_preFirstMT_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["Semax_lastmsk"],MTdf_pst_tot["BlackHole_lastmsk"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title(f"Black Hole($M_p$)-Star($M_s$) Binaries after last MT N:{len(df['BlackHole'])}",  fontsize=40, pad=30)
    f.savefig(directoryp + "MT_postLastMT_" + current_time + ".png",   bbox_inches='tight')
    plt.close()


    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax_mainpremsk"],MTdf_pre_tot["BlackHole_mainpremsk"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title(f"Black Hole($M_p$)-Main Sequence Star($M_s$) Binaries before first MT N:{len(df['BlackHole'])}",  fontsize=40, pad=30)
    f.savefig(directoryp + "MT_MainSequence_preFirstMT_" + current_time + ".png",   bbox_inches='tight')    
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["Semax_mainpstmsk"],MTdf_pst_tot["BlackHole_mainpstmsk"]), bins=20)
    values = percentage(systemSize, values)
    df = pd.DataFrame({
    "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    "Semax": np.tile(ybins[:-1], len(xbins)-1),
    'frequency': values.T.flatten()})
    vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    f, ax = plt.subplots(figsize=(36, 32))
    xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    ax.set_title(f"Black Hole($M_p$)-Main Sequence Star($M_s$) Binaries after last MT N:{len(df['BlackHole'])}",  fontsize=40, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_" + current_time + ".png",   bbox_inches='tight')
    plt.close()


    # values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax_MT"],MTdf_pre_tot["Companion"]), bins=20)
    # values = percentage(systemSize, values)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="Semax", values="frequency")
    # f, ax = plt.subplots(figsize=(36, 32))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major', labelsize=20)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=40, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=40, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (MT)",  fontsize=40, pad=30)
    # f.savefig(directoryp + "_MT_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()


    #####

    
# %%

####### PARENT SCRIPT: tertius_plotting_020525
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
import calculations as calc

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

mode = ["Default_WD_Enabled/"] #"Default/","Limited/", ,"Limited_WD_Enabled/"

path = '/data/a.saricaoglu/repo/COMPAS'
# # Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print(sys.path)
# Choose an output hdf5 file to work with
for mod in mode:
    matplotlib.rcParams['figure.figsize'] = (21,14)
    matplotlib.rcParams['lines.markersize'] = 1
    matplotlib.rcParams['font.size'] = 20
    matplotlib.rcParams['legend.loc'] = "upper right"

    pathToData = path + '/Files/' + mod + "02.19/" #change the date accordingly the date of the files created via Sec

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
    MTdf_pre_singleBH["BH1"] = MT["MASSPREMT1"]*MT_mask["MASKPREMTBH1"]
    MTdf_pre_singleBH["CP2"] = MT["MASSPREMT2"]*MT_mask["MASKPREMTBH1"]
    MTdf_pre_singleBH["SA1"] = MT["SEMIMAJORAXISPREMT"]*MT_mask["MASKPREMTBH1"]

    #Systems which secondary is BH
    MTdf_pre_singleBH["BH2"] = MT["MASSPREMT2"]*MT_mask["MASKPREMTBH2"]
    MTdf_pre_singleBH["CP1"] = MT["MASSPREMT1"]*MT_mask["MASKPREMTBH2"]
    MTdf_pre_singleBH["SA2"] = MT["SEMIMAJORAXISPREMT"]*MT_mask["MASKPREMTBH2"]

########## UPDATE 7.11.24
    #disregarding initially primary or secondary stars, BH objects are considered as primary and the star as secondary.
    MTdf_pre_tot["BlackHole"] = np.concatenate([MTdf_pre_singleBH["BH1"], MTdf_pre_singleBH["BH2"]])
    MTdf_pre_tot["Companion"] = np.concatenate([MTdf_pre_singleBH["CP1"], MTdf_pre_singleBH["CP2"]])
    MTdf_pre_tot["Semax"] = np.concatenate([MTdf_pre_singleBH["SA1"] , MTdf_pre_singleBH["SA2"] ])


    MTdf_pst_singleBH = pd.DataFrame()
    MTdf_pst_tot = pd.DataFrame()
    #Systems which primary is BH
    MTdf_pst_singleBH["BH1"] = MT["MASSPSTMT1"]*MT_mask["MASKPSTMTBH1"]
    MTdf_pst_singleBH["CP2"] = MT["MASSPSTMT2"]*MT_mask["MASKPSTMTBH1"]
    MTdf_pst_singleBH["SA1"] = MT["SEMIMAJORAXISPSTMT"]*MT_mask["MASKPSTMTBH1"]
    MTdf_pst_singleBH["OP1"] = MT_mask["ORBITALPERIODPSTMT"]*MT_mask["MASKPSTMTBH1"]

    #Systems which secondary is BH
    MTdf_pst_singleBH["BH2"] = MT["MASSPSTMT2"]*MT_mask["MASKPSTMTBH2"]
    MTdf_pst_singleBH["CP1"] = MT["MASSPSTMT1"]*MT_mask["MASKPSTMTBH2"]
    MTdf_pst_singleBH["SA2"] = MT["SEMIMAJORAXISPSTMT"]*MT_mask["MASKPSTMTBH2"]
    MTdf_pst_singleBH["OP2"] = MT_mask["ORBITALPERIODPSTMT"]*MT_mask["MASKPSTMTBH2"]
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
    MTdf_pst_tot["OP_mainpstmsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst , MTdf_pst_singleBH["OP2"]*mainsequence_pst ])

    mainsequence_pre_Tlim = mainsequence_pre & MT_mask['MASKPREMTORBPER']
    MTdf_pre_tot["BlackHole_mainpreTlimmsk"] = np.concatenate([MTdf_pre_singleBH["BH1"]*mainsequence_pre_Tlim, MTdf_pre_singleBH["BH2"]*mainsequence_pre_Tlim])
    MTdf_pre_tot["Companion_mainpreTlimmsk"] = np.concatenate([MTdf_pre_singleBH["CP1"]*mainsequence_pre_Tlim, MTdf_pre_singleBH["CP2"]*mainsequence_pre_Tlim])
    MTdf_pre_tot["Semax_mainpreTlimmsk"] = np.concatenate([MTdf_pre_singleBH["SA1"]*mainsequence_pre_Tlim , MTdf_pre_singleBH["SA2"]*mainsequence_pre_Tlim ])

    mainsequence_pst_Tlim = mainsequence_pst & MT_mask['MASKPSTMTORBPER']
    MTdf_pst_tot["BlackHole_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst_Tlim, MTdf_pst_singleBH["BH2"]*mainsequence_pst_Tlim])
    MTdf_pst_tot["Companion_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["CP1"]*mainsequence_pst_Tlim, MTdf_pst_singleBH["CP2"]*mainsequence_pst_Tlim])
    MTdf_pst_tot["Semax_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst_Tlim , MTdf_pst_singleBH["SA2"]*mainsequence_pst_Tlim ])
    MTdf_pst_tot["OP_mainpstTlimmsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst_Tlim , MTdf_pst_singleBH["OP2"]*mainsequence_pst_Tlim ])

    searchability = calc.searchability(MTdf_pst_singleBH["SA1"])
    N = np.sum(searchability)
    searchability_index = N/len(searchability)
    # searchability_index = np.sum(MT_mask['MASKPSTMTSEARCHABILITY'])/len(MT_mask['MASKPSTMTSEARCHABILITY'])
    print('Searchability index:', searchability_index)
    mainsequence_pst_alim_masuda = mainsequence_pst & MT_mask['MASKPSTMTSEMAJ_masuda'] 
    MTdf_pst_tot["BlackHole_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst_alim_masuda, MTdf_pst_singleBH["BH2"]*mainsequence_pst_alim_masuda])
    MTdf_pst_tot["Companion_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["CP1"]*mainsequence_pst_alim_masuda, MTdf_pst_singleBH["CP2"]*mainsequence_pst_alim_masuda])
    MTdf_pst_tot["Semax_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst_alim_masuda , MTdf_pst_singleBH["SA2"]*mainsequence_pst_alim_masuda ])
    MTdf_pst_tot["OP_mainpstalimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst_alim_masuda , MTdf_pst_singleBH["OP2"]*mainsequence_pst_alim_masuda ])

    # print('Number of systems past last MT event with Orbital period less than 30 days:', np.sum(mainsequence_pst_alim_masuda))
    print('Number of searchable systems past last MT event with semaj less than 0.4 au:', searchability_index*np.sum(mainsequence_pst_alim_masuda))
    print('Number of systems past last with semaj less than 0.4 au:', np.sum(mainsequence_pst_alim_masuda))

    mainsequence_pst_Tlim_masuda = mainsequence_pst & MT_mask['MASKPSTMTORBPER_masuda']
    MTdf_pst_tot["BlackHole_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["BH1"]*mainsequence_pst_Tlim_masuda, MTdf_pst_singleBH["BH2"]*mainsequence_pst_Tlim_masuda])
    MTdf_pst_tot["Companion_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["CP1"]*mainsequence_pst_Tlim_masuda, MTdf_pst_singleBH["CP2"]*mainsequence_pst_Tlim_masuda])
    MTdf_pst_tot["Semax_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["SA1"]*mainsequence_pst_Tlim_masuda , MTdf_pst_singleBH["SA2"]*mainsequence_pst_Tlim_masuda ])
    MTdf_pst_tot["OP_mainpstTlimmasudamsk"] = np.concatenate([MTdf_pst_singleBH["OP1"]*mainsequence_pst_Tlim_masuda , MTdf_pst_singleBH["OP2"]*mainsequence_pst_Tlim_masuda ])

    print('Number of searchable systems past last MT event with Orbital period less than 30 days:', searchability_index* np.sum(mainsequence_pst_Tlim_masuda))
    print('Number of systems past last MT event with Orbital period less than 30 days:', np.sum(mainsequence_pst_Tlim_masuda))
    print('Number of systems past last MT with Orbital period less than 30 days:', np.sum((MT_mask['MASKPSTMTORBPER_masuda']) * last))
    print('Number of systems past last MT with orbital period less than 70 days:', np.sum((MT_mask['MASKPSTMTORBPER']) * last))
    print('max Orbital period in days masuda:', np.max(MTdf_pst_tot["OP_mainpstTlimmasudamsk"]))
    print('mean Orbital period in days masuda:', np.mean(MTdf_pst_tot["OP_mainpstTlimmasudamsk"]))
    print('max semimajor axis in au masuda:', np.max(MTdf_pst_tot["Semax_mainpstTlimmasudamsk"]))
    print('max orbital period in days:', np.max(MTdf_pst_tot["OP_mainpstTlimmsk"]))
    print('mean orbital period in days:', np.mean(MTdf_pst_tot["OP_mainpstTlimmsk"]))
    print('max semimajor axis in au:', np.max(MTdf_pst_tot["Semax_mainpstTlimmsk"]))
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

    # print(max(MTdf_pre_tot["Semax"]))
    # print(max(MTdf_pre_tot["BlackHole"]))
    # print(min(MTdf_pre_tot["Semax"]))
    # print(min(MTdf_pre_tot["BlackHole"]))
    

    
######## UPDATE 07.11.24
    

 
    if not os.path.exists(path+ "/Plots/" + mod +  str(c.strftime("%m.%d")+ "/" + current_time + "/") ): 
        os.makedirs(path + "/Plots/"  + mod +  str(c.strftime("%m.%d")+ "/" + current_time + "/") ) 
    directoryp = path + "/Plots/"  + mod +  str(c.strftime("%m.%d") + "/" + current_time + "/")  
    # Produces heatmaps. Only black holes, companions are commented since not needed.
    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax"], SPdf_tot["BlackHole"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # Check the shape of values to ensure it matches the expected dimensions
    print("Shape of values:", values.shape)
    print("Length of xbins:", len(xbins))
    print("Length of ybins:", len(ybins))
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    print('values:', values)
    print('vals:', vals)

    # Check the shape of the pivot table to ensure it matches the expected dimensions
    print("Shape of pivot table:", vals.shape)

    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True, fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=30)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    print(ax.get_xticks())
    print(ax.get_yticks()) 
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)", fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)", fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Star($M_s$) Binary Syst at ZAMS (SP) N:{len(SP['MASSZAMSSP1'])}", fontsize=30, pad=30)
    f.savefig(directoryp + "SP_atZAMS_inMT" + current_time + ".png", bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(SPdf_tot["Semax_MT"],SPdf_tot["BlackHole_MT"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    print('values:', values)
    print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=30)    

    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Star($M_s$) Binaries at ZAMS undergoing MT (MTs in SP) N:{np.sum(SP_mask['MASKMTinSP'])})",  fontsize=30, pad=30)
    f.savefig(directoryp + "SP_atZAMS_inMT_" + current_time + ".png",   bbox_inches='tight')
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
    # ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=45, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=45, labelpad=25)
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
    # ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=45, labelpad=25)
    # ax.set_ylabel("$M_s$  [$M_{\odot}$]",  fontsize=45, labelpad=25)
    # ax.set_title("Black Hole ($M_p$)  - Star ($M_s$) Binary Systems (MTs in SP)",  fontsize=40, pad=30)
    # f.savefig(directoryp + "SP_Ms_MT_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    # values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax"],MTdf_pre_tot["BlackHole"]), bins=20)
    # values = percentage(systemSize, values)
    # df = pd.DataFrame({
    # "BlackHole": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Semax", columns="BlackHole", values="frequency")
    # f, ax = plt.subplots(figsize=(36, 32))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax, xticklabels=xlabl, yticklabels=ylabl)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    # ax.set_xlabel("Semimajor Axis inital [AU]",  fontsize=45, labelpad=25)
    # ax.set_ylabel("$M_p$  [$M_{\odot}$]",  fontsize=45, labelpad=25)
    # ax.set_title(f"Black Hole ($M_p$)  - Star ($M_s$) Binary Systems(MTO) N:{len(df['BlackHole'])}",  fontsize=40, pad=30)
    # f.savefig(directoryp + "MT_All_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()

    #%%
    
    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax_firstmsk"],MTdf_pre_tot["BlackHole_firstmsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Star($M_s$) Binaries before first MT N:{np.sum(first)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_preFirstMT_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["Semax_lastmsk"],MTdf_pst_tot["BlackHole_lastmsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(36, 32))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Star($M_s$) Binaries after last MT N:{np.sum(last)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_postLastMT_" + current_time + ".png",   bbox_inches='tight')
    plt.close()


    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax_mainpremsk"],MTdf_pre_tot["BlackHole_mainpremsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries before first MT N:{np.sum(mainsequence_pre)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_preFirstMT_" + current_time + ".png",   bbox_inches='tight')    
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["Semax_mainpstmsk"],MTdf_pst_tot["BlackHole_mainpstmsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["OP_mainpstmsk"],MTdf_pst_tot["BlackHole_mainpstmsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_OP_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax_mainpreTlimmsk"],MTdf_pre_tot["BlackHole_mainpreTlimmsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries before first MT N:{np.sum(mainsequence_pre_Tlim)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_preFirstMT_Tlim_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["Semax_mainpstTlimmsk"],MTdf_pst_tot["BlackHole_mainpstTlimmsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_Tlim)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_Tlim_" + current_time + ".png",bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["OP_mainpstTlimmsk"],MTdf_pst_tot["BlackHole_mainpstTlimmsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})

    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_Tlim)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_Tlim_OP_" + current_time + ".png",bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["Semax_mainpstalimmasudamsk"],MTdf_pst_tot["BlackHole_mainpstalimmasudamsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    print('xbins:', xbins)
    print('ybins:', ybins)
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})
    values = values*searchability_index
    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_alim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_alim_masuda_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["Semax_mainpstTlimmasudamsk"],MTdf_pst_tot["BlackHole_mainpstTlimmasudamsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    print('xbins:', xbins)
    print('ybins:', ybins)
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})
    values = values*searchability_index
    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_Tlim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_Tlim_masuda_" + current_time + ".png",   bbox_inches='tight')
    plt.close()



    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["OP_mainpstalimmasudamsk"],MTdf_pst_tot["BlackHole_mainpstalimmasudamsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})
    values = values*searchability_index
    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_alim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_alim_masuda_OP" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(*mylog(MTdf_pst_tot["OP_mainpstTlimmasudamsk"],MTdf_pst_tot["BlackHole_mainpstTlimmasudamsk"]), bins=20)
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(10**x) for x in xbins]    
    ybins = ['{:.1f}'.format(10**y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})
    values = values*searchability_index
    df = pd.DataFrame(values, columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_alim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_Tlim_masuda_OP" + current_time + ".png",   bbox_inches='tight')
    plt.close()


# Define the bin edges for the histogram
    bh_mass_bins = [5, 9.1, 16.6, 30.2, 54.9, 100]
    orbital_period_bins = [0.3, 0.8, 1.9, 4.8, 11.9, 30.0]

    # Define the 5x5 occurrence rate data (from the figure)
    # occurrence_rates = np.array([
    #     [1.3e-6, 1.8e-6, 1.9e-6, 6.7e-7, 8.7e-9],
    #     [2.3e-6, 3.1e-6, 3.1e-6, 7.9e-7, 3.9e-10],
    #     [4.7e-6, 7.9e-6, 1.9e-6, 1.5e-9, 2.2e-10],
    #     [8.4e-6, 1.8e-5, 1.2e-5, 2.7e-6, 2.0e-9],
    #     [8.2e-6, 1.5e-5, 8.9e-6, 1.3e-8, 2.5e-9]
    # ])
    occurrence_rates = np.array([
        [8.2e-6, 1.5e-5, 8.9e-6, 1.3e-8, 2.5e-9],
        [8.4e-6, 1.8e-5, 1.2e-5, 2.7e-6, 2.0e-9],
        [4.7e-6, 7.9e-6, 1.9e-6, 1.5e-9, 2.2e-10],
        [2.3e-6, 3.1e-6, 3.1e-6, 7.9e-7, 3.9e-10],
        [1.3e-6, 1.8e-6, 1.9e-6, 6.7e-7, 8.7e-9]
    ])   

    # Create the DataFrame
    # df0 = pd.DataFrame(occurrence_rates.T, columns=orbital_period_bins[:-1], index=bh_mass_bins[:-1])
    df0 = pd.DataFrame({
    "BlackHole": np.tile(bh_mass_bins[1:], len(orbital_period_bins)-1),
    "Orbper": np.repeat(orbital_period_bins[1:], len(bh_mass_bins)-1),
    'frequency': occurrence_rates.T.flatten()})

    print('df0:', df0)
    vals0 = df0.pivot(index="BlackHole", columns="Orbper", values="frequency")
    print('vals0:', vals0)
    # Rename index and columns for clarity
    vals_masuda = vals0

    values, xbins, ybins = np.histogram2d(MTdf_pst_tot["OP_mainpstTlimmasudamsk"],MTdf_pst_tot["BlackHole_mainpstTlimmasudamsk"], bins=[(orbital_period_bins),(bh_mass_bins)]) 
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(x) for x in xbins]    
    ybins = ['{:.1f}'.format(y) for y in ybins]

    values = values*searchability_index
    print('values:', values)
    print('x bins:', xbins)
    print('y bins:', ybins)
    df = pd.DataFrame({
    "BlackHole": np.tile(ybins[1:], len(xbins)-1),
    "Orbper": np.repeat(xbins[1:], len(ybins)-1),
    'frequency': values.T.flatten()})
    df = df.pivot(index="BlackHole", columns="Orbper", values="frequency")
    values = df.values
    values0 = vals_masuda.values
    print('values:', values)

    df = pd.DataFrame(np.abs(values-values0), columns=xbins[1:], index=ybins[1:])
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)

    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_alim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_Tlim_masuda_diff_" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    values, xbins, ybins = np.histogram2d(MTdf_pst_tot["OP_mainpstTlimmasudamsk"],MTdf_pst_tot["BlackHole_mainpstTlimmasudamsk"],  bins=[(orbital_period_bins[1:]),(bh_mass_bins[1:])])     
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(x) for x in xbins]    
    ybins = ['{:.1f}'.format(y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})
    values = values*searchability_index
    df = pd.DataFrame(values, columns=xbins, index=ybins)
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_alim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_Tlim_masuda_OP_searchability" + current_time + ".png",   bbox_inches='tight')
    plt.close()
    values, xbins, ybins = np.histogram2d(MTdf_pst_tot["OP_mainpstTlimmasudamsk"],MTdf_pst_tot["BlackHole_mainpstTlimmasudamsk"], bins=[(orbital_period_bins[1:]),(bh_mass_bins[1:])])     
    values = percentage(systemSize, values)
    xbins = ['{:.2f}'.format(x) for x in xbins]    
    ybins = ['{:.1f}'.format(y) for y in ybins]
    # df = pd.DataFrame({
    # "BlackHole": np.tile(ybins[:-1], len(xbins)-1),
    # "Orbper": np.repeat(xbins[:-1], len(ybins)-1),
    # 'frequency': values.T.flatten()})
    values = values
    df = pd.DataFrame(values, columns=xbins, index=ybins)
    print('df:', df)
    vals = df
    # print('values:', values)
    print('vals:', vals)
    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)
    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_alim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "MT_MainSequence_postLastMT_Tlim_masuda_OP_noSearchability" + current_time + ".png",   bbox_inches='tight')
    plt.close()
    print('values:', values)
    print('vals:', vals)
    print('vals_masuda:', vals_masuda)


    # values, xbins, ybins = np.histogram2d(*mylog(MTdf_pre_tot["Semax_MT"],MTdf_pre_tot["Companion"]), bins=20)
    # values = percentage(systemSize, values)
    # df = pd.DataFrame({
    # "Companion": np.repeat(xbins[:-1], len(ybins)-1),
    # "Semax": np.tile(ybins[:-1], len(xbins)-1),
    # 'frequency': values.T.flatten()})
    # vals = df.pivot(index="Companion", columns="Semax", values="frequency")
    # f, ax = plt.subplots(figsize=(18, 16))
    # xlabl = ["{:.2f}".format(10**x) for x in vals.columns]
    # ylabl = ["{:.1f}".format(10**y) for y in vals.index]
    # hmapplot = sns.heatmap(vals, annot=True,fmt=".0f", linewidths=.1, ax=ax)
    # cbar = hmapplot.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=40)
    # ax.invert_yaxis()
    # ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    # ax.set_xlabel("Semimajor Axis (au)",  fontsize=45, labelpad=25)
    # ax.set_ylabel("$M_s$ (solar mass)",  fontsize=45, labelpad=25)
    # #ax.set_title("Black Hole (BH mass)  - Star ($M_s$) Binary Systems (MT)",  fontsize=30, pad=30)
    # f.savefig(directoryp + "_MT_Ms_" + current_time + ".png",   bbox_inches='tight')
    # plt.close()


    #####

    f, ax = plt.subplots(figsize=(18, 16))

    hmapplot = sns.heatmap(vals_masuda, annot=True,fmt=".1e", linewidths=.1, ax=ax)
    cbar = hmapplot.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    ax.invert_yaxis()
    ax.set_xticks([x+0.5 for x in ax.get_xticks()],labels=[i for i in ax.get_xticklabels()])
    ax.set_yticks([y+0.5 for y in ax.get_yticks()],labels=[i for i in ax.get_yticklabels()])
    ax.tick_params(axis='both', which='major',labelsize=30, rotation=45)
    # ax.set_xticks(np.log10(xbins), xlabl)
    # ax.set_yticks(np.log10(ybins), ylabl)
    ax.set_xlabel("Orbital period (days)",  fontsize=45, labelpad=25)
    ax.set_ylabel("BH mass (solar mass)",  fontsize=45, labelpad=25)

    #ax.set_title(f"Black Hole(BH mass)-Main Sequence Star($M_s$) Binaries after last MT N:{np.sum(mainsequence_pst_alim_masuda)}",  fontsize=30, pad=30)
    f.savefig(directoryp + "masuda" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    f, ax = plt.subplots(figsize=(18, 16))
    plt.hist(MTdf_pst_tot["OP_mainpstTlimmasudamsk"], bins=100, histtype='stepfilled', color='yellow', label='T<30 days', log=True)
    plt.hist(MTdf_pst_tot["OP_mainpstTlimmsk"], bins=100, histtype='step', color='r', label='T<70 days', log=True)
    plt.legend(fontsize=30)
    f.savefig(directoryp + "hist_1" + current_time + ".png",   bbox_inches='tight')
    plt.close()

    f, ax = plt.subplots(figsize=(18, 16))
    plt.hist(MTdf_pst_tot["OP_mainpstTlimmsk"], bins=100, histtype='stepfilled', color='cyan', label='T<70 days', log=True)
    plt.hist(MTdf_pst_tot["OP_mainpstmsk"], bins=100, histtype='step', color='k', label='no limit', log=True)
    plt.legend(fontsize=30)
    f.savefig(directoryp + "hist_2" + current_time + ".png",   bbox_inches='tight')
    plt.close()
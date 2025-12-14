#### PARENT SCRIPT: Stage_03_detailed_evolution.py
#----> this version plots population properties.
import os, sys
import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
import logging            # for reading the COMPAS data
import time as t                      # for finding computation time
import datetime as dt
import h5py as h5                # for reading the COMPAS data
import matplotlib.pyplot as plt  # for plotting
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from src.Utils import calculations as calc   # functions from calculations.py
from src.Utils import utils as utils      # functions from utils.py

now = dt.datetime.now()
# Choose the modee to process
mode = 'WD_Enabled_Detailed'

# Import COMPAS specific scripts
compasRootDir = os.environ['COMPAS_ROOT_DIR']
sys.path.append(compasRootDir + '/postProcessing/PythonScripts')
print("COMPAS_ROOT_DIR:", compasRootDir)
print(sys.path)

pathToData = os.path.join(compasRootDir, "luxetenebrae/runs", mode)
pathToFiles = os.path.join(compasRootDir, "luxetenebrae/files", mode)
pathToLogs = os.path.join(compasRootDir, "luxetenebrae/logs", mode)
pathToPlots = os.path.join(compasRootDir, "luxetenebrae/plots", mode)
# Create output directory if it doesn't exist
os.makedirs(pathToFiles, exist_ok=True)
os.makedirs(pathToLogs, exist_ok=True)
os.makedirs(pathToPlots, exist_ok=True)

# Get the script name
script_name = os.path.basename(__file__)
# Configure logging
log_filename = f"{pathToLogs}/{script_name}_{now.strftime('%m.%d')}_{now.strftime('%H%M')}.log"
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect stdout and stderr to the log file
sys.stdout = utils.StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
sys.stderr = utils.StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
# Log the start of the script with the script name
logging.info(f'Script {script_name} started')

# Displays Time
s = dt.datetime.now()
starttime = t.ctime(t.time())
start = t.process_time()
start_time = s.strftime("%d%m%y") + "_" + s.strftime('%H%M')
print("Start time :", start_time)

out = [f for f in os.listdir(pathToFiles) if "stage_02DE_outputs" in f]
print('Output files found:', out)
detailed_output = []
SP_output = []
MT_output = []
CE_output = []
for f in out:
    print("Reading file: ", pathToFiles + "/" + f)

    data = fits.open(pathToFiles + "/" + f)
    print("Number of tables in the file:", len(data))
    detailed_output.append(data[1])
    SP_output.append(data[2])

    if len(data) > 3:
        MT_output.append(data[3])

    if len(data) > 4:
        CE_output.append(data[4])



# Keeps corresponding numerical values. 
# SPs
c = dt.datetime.now()
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')

i = 0

merger_flags_manual = []
merger_flags_mthist = []
merger_flags = []
bhms_final_systems = []
bhbh_final_systems = []

semimajoraxis_average = []
lifetime_average = []
eccentricity_average = []

j = 0
k = 0

primary_bh = 0
secondary_bh = 0
# Choose an output hdf5 file to work with

ps = open(f'{pathToFiles}/survivors.txt', 'r')    
survivors = ps.read().splitlines()


# out = [f for f in os.listdir(pathToFiles) if ".fits" in f]
# data_outputs = []
# for f in out:
#     # print("Reading file: ", pathToData + '/Files/' + mode + "/" + f)

#     data = fits.open(pathToFiles + "/" + f)
#     data_outputs.append(data[1])

runs= [x[0] for x in os.walk(pathToData) if "Detailed_Output" in x[1]]
print('runs:', runs)
data_outputs2 = []
for run in runs[:-1]:
    out2 = [f for f in os.listdir(run) if ".h5" in f]
    for f in out2:
        print("Reading file: ", run  + "/" + f)
        data = h5.File( run + "/" + f, 'r')
        data_outputs2.append(data)

# Keeps corresponding numerical values. 
# SPs
c = dt.datetime.now()
current_time = c.strftime("%d%m%y") + "_" + c.strftime('%H%M')





# f = open(pathToData + '/Files/' + mode + str(s.strftime("%m.%d")) +  "_general_outputs.txt", "a")

# f.writelines(["\n","\n Run :", start_time])
# print(f'Data outputs :, {data_outputs}')
primary_mass = []
secondary_mass = []
survivor_primary_mass = []
survivor_secondary_mass = []
for Data in data_outputs2:
    Data = h5.File(Data.filename, 'r')
    SPs = Data['BSE_System_Parameters']   
    massZamsSP1 = SPs['Mass@ZAMS(1)'][()] 
    massZamsSP2 = SPs['Mass@ZAMS(2)'][()]
    primary_mass = np.concatenate((primary_mass, massZamsSP1))
    secondary_mass = np.concatenate((secondary_mass, massZamsSP2))
    seedsSP = SPs['SEED'][()]
    seedsSP_mask = np.isin(seedsSP, survivors)  
    seeds_matched = seedsSP[seedsSP_mask]  
    survivor_primary_mass = np.concatenate((massZamsSP1[seedsSP_mask], survivor_primary_mass))
    survivor_secondary_mass = np.concatenate((massZamsSP2[seedsSP_mask], survivor_secondary_mass))     
    # Build the dictionary

for Dat in detailed_output:
    
    Data = Dat.data
    print("Batch (of " + str(len(detailed_output)) + " sys.) " + str(i) +  " start time :", current_time)
    semimajoraxis = Data['SemiMajorAxis']
    # print(f'Length of semimajor axis array: {len(semimajoraxis)}')
    # print(f'Minimum semimajor axis: {np.min(semimajoraxis)}')

    if np.any(semimajoraxis < 0):
        print(f"Negative semimajor axis found in batch {i}, skipping this batch")
        k += 1
        continue
    seed = Data['Seed'][0]
    print('seed', seed)
    stellar_type_1 = Data['Stellar_Type_1']
    stellar_type_2 = Data['Stellar_Type_2']
    radius_1 = (Data['Radius(1)']*const.R_sun).to(u.au).value
    radius_2 = (Data['Radius(2)']*const.R_sun).to(u.au).value
    time = Data['Time']
    roche_lobe_1 = (Data['RocheLobe(1)']*const.R_sun).to(u.au).value
    roche_lobe_2 = (Data['RocheLobe(2)']*const.R_sun).to(u.au).value
    mass_1 = Data['Mass(1)']
    mass_2 = Data['Mass(2)']
    mass_zams_1 = Data['Mass@ZAMS(1)']
    mass_zams_2 = Data['Mass@ZAMS(2)']
    mt_history = Data['MT_History']
    eccentricity = Data['Eccentricity']
    # rocheR_periapsis = calculate_periapsis()

    bhms_mask = np.array([utils.is_ms_bh_pair(st1, st2) for st1, st2 in zip(stellar_type_1, stellar_type_2)])

    print(f"Number of MS + BH states in this system: {np.sum(bhms_mask)}")
    if np.sum(bhms_mask) == 0:
        print(f"No BH-MS simultaneous emergence in this batch, skipping")
        j += 1
        continue

    if utils.is_ms_bh_pair(stellar_type_1[-1], stellar_type_2[-1]):
        bhms_final_systems.append(seed)

    if utils.is_bh_bh_pair(stellar_type_1[-1], stellar_type_2[-1]):
        bhbh_final_systems.append(seed)

    time_filtered = time[bhms_mask]   # Assuming 'Time' is the time evolution data
    stype1_filtered = stellar_type_1[bhms_mask] 
    stype2_filtered = stellar_type_2[bhms_mask] 
    radius1_filtered = radius_1[bhms_mask] 
    radius2_filtered = radius_2[bhms_mask] 
    semimajoraxis_filtered = semimajoraxis[bhms_mask] 
    eccentricity_filtered = eccentricity[bhms_mask] 

    bhms_lifetime = time_filtered[-1] - time_filtered[0]
    lifetime_average.append(bhms_lifetime) 

    semaj_average = np.mean(semimajoraxis_filtered)
    semimajoraxis_average.append(semaj_average)

    ecc_average = np.mean(eccentricity_filtered)
    eccentricity_average.append(ecc_average)

    primary_mask = np.array([utils.which_bh(st1, st2) for st1, st2 in zip(stype1_filtered, stype2_filtered)])
    # primary_mass_filtered = mass_zams_1[primary_mask == 1]
    # secondary_mass_filtered = mass_zams_2[primary_mask == 0]
    if np.sum(primary_mask) < len(primary_mask):
        secondary_bh += 1
    else:
        primary_bh += 1

    i = i + 1

# After processing all systems, analyze semimajor axis distribution by percentiles
print("\n=== SEMIMAJOR AXIS POPULATION ANALYSIS ===")
semimajoraxis_array = np.array(semimajoraxis_average)

# Calculate key percentiles
percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
percentile_values = np.percentile(semimajoraxis_array, percentiles)

print(f"Total number of systems analyzed: {len(semimajoraxis_array)}")
print(f"Semimajor axis statistics:")
print(f"  Mean: {np.mean(semimajoraxis_array):.3e} AU")
print(f"  Median: {np.median(semimajoraxis_array):.3e} AU")
print(f"  Standard deviation: {np.std(semimajoraxis_array):.3e} AU")
print(f"  Min: {np.min(semimajoraxis_array):.3e} AU")
print(f"  Max: {np.max(semimajoraxis_array):.3e} AU")

print(f"\nPercentile analysis:")
for p, val in zip(percentiles, percentile_values):
    print(f"  {p:2d}th percentile: {val:.3e} AU")

# Analyze population distribution in different ranges
ranges = [
    (0, percentile_values[2]),  # Bottom 10%
    (percentile_values[2], percentile_values[3]),  # 10th-25th percentile
    (percentile_values[3], percentile_values[4]),  # 25th-50th percentile (Q1-Q2)
    (percentile_values[4], percentile_values[5]),  # 50th-75th percentile (Q2-Q3)
    (percentile_values[5], percentile_values[7]),  # 75th-95th percentile
    (percentile_values[7], np.max(semimajoraxis_array))  # Top 5%
]

range_labels = ['Bottom 10%', '10th-25th', '25th-50th (Q1-Q2)', 
                '50th-75th (Q2-Q3)', '75th-95th', 'Top 5%']

print(f"\nPopulation distribution analysis:")
# Convert to arrays for easier analysis
lifetime_array = np.array(lifetime_average)
eccentricity_array = np.array(eccentricity_average)

for i, (range_label, (min_val, max_val)) in enumerate(zip(range_labels, ranges)):
    mask = (semimajoraxis_array >= min_val) & (semimajoraxis_array <= max_val)
    count = np.sum(mask)
    percentage = (count / len(semimajoraxis_array)) * 100
    
    if count > 0:
        # Semimajor axis statistics
        avg_sa = np.mean(semimajoraxis_array[mask])
        std_sa = np.std(semimajoraxis_array[mask])
        
        # Lifetime statistics
        avg_lifetime = np.mean(lifetime_array[mask])
        std_lifetime = np.std(lifetime_array[mask])
        median_lifetime = np.median(lifetime_array[mask])
        
        # Eccentricity statistics
        avg_ecc = np.mean(eccentricity_array[mask])
        std_ecc = np.std(eccentricity_array[mask])
        median_ecc = np.median(eccentricity_array[mask])
        
        print(f"  {range_label}: {count} systems ({percentage:.1f}%)")
        print(f"    Semimajor Axis: {avg_sa:.3e} ± {std_sa:.3e} AU")
        print(f"    Lifetime: {avg_lifetime:.3e} ± {std_lifetime:.3e} Myr (median: {median_lifetime:.3e})")
        print(f"    Eccentricity: {avg_ecc:.3f} ± {std_ecc:.3f} (median: {median_ecc:.3f})")
    else:
        print(f"  {range_label}: {count} systems ({percentage:.1f}%) - No data")
    print()

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. Histogram with percentile lines
axes[0,0].hist(semimajoraxis_average, bins=50, color='blue', alpha=0.7, density=True)
for p, val in zip([10, 25, 50, 75, 90], percentile_values[2:7]):
    axes[0,0].axvline(val, color='red', linestyle='--', alpha=0.7, 
                     label=f'{p}th percentile' if p in [10, 50, 90] else '')
axes[0,0].set_xlabel('Semimajor Axis (AU)')
axes[0,0].set_ylabel('Density')
axes[0,0].set_yscale('log')
axes[0,0].set_xscale('log')
axes[0,0].set_title('Semimajor Axis Distribution with Percentiles')
axes[0,0].legend()

# 2. Cumulative distribution
sorted_sa = np.sort(semimajoraxis_array)
cumulative = np.arange(1, len(sorted_sa) + 1) / len(sorted_sa) * 100
axes[0,1].plot(sorted_sa, cumulative, 'b-', linewidth=2)
axes[0,1].set_xlabel('Semimajor Axis (AU)')
axes[0,1].set_ylabel('Cumulative Percentage')
axes[0,1].set_xscale('log')
axes[0,1].grid(True, alpha=0.3)
axes[0,1].set_title('Cumulative Distribution of Semimajor Axis')

# 3. Box plot showing quartiles
box_data = [semimajoraxis_array[
    (semimajoraxis_array >= ranges[i][0]) & (semimajoraxis_array <= ranges[i][1])
] for i in range(len(ranges))]
bp = axes[1,0].boxplot(box_data, labels=range_labels, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('lightblue')
axes[1,0].set_ylabel('Semimajor Axis (AU)')
axes[1,0].set_yscale('log')
axes[1,0].set_title('Distribution by Population Ranges')
axes[1,0].tick_params(axis='x', rotation=45)

# 4. Population count by range
range_counts = [np.sum((semimajoraxis_array >= r[0]) & (semimajoraxis_array <= r[1])) 
                for r in ranges]
axes[1,1].bar(range_labels, range_counts, color='skyblue', alpha=0.7)
axes[1,1].set_ylabel('Number of Systems')
axes[1,1].set_title('System Count by Population Range')
axes[1,1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig(pathToPlots + 'semimajor_axis_population_analysis.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# Create additional comprehensive analysis plots for lifetime and eccentricity by semimajor axis percentiles
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Convert to arrays for easier analysis
lifetime_array = np.array(lifetime_average)
eccentricity_array = np.array(eccentricity_average)

# 1. Lifetime vs Semimajor Axis scatter plot with percentile regions
axes[0,0].scatter(semimajoraxis_array, lifetime_array, alpha=0.6, s=20, c='blue')
for p, val in zip([10, 25, 50, 75, 90], percentile_values[2:7]):
    axes[0,0].axvline(val, color='red', linestyle='--', alpha=0.7)
axes[0,0].set_xlabel('Semimajor Axis (AU)')
axes[0,0].set_ylabel('Lifetime (Myr)')
axes[0,0].set_xscale('log')
axes[0,0].set_yscale('log')
axes[0,0].set_title('Lifetime vs Semimajor Axis')
axes[0,0].grid(True, alpha=0.3)

# 2. Eccentricity vs Semimajor Axis scatter plot with percentile regions
axes[0,1].scatter(semimajoraxis_array, eccentricity_array, alpha=0.6, s=20, c='green')
for p, val in zip([10, 25, 50, 75, 90], percentile_values[2:7]):
    axes[0,1].axvline(val, color='red', linestyle='--', alpha=0.7)
axes[0,1].set_xlabel('Semimajor Axis (AU)')
axes[0,1].set_ylabel('Eccentricity')
axes[0,1].set_xscale('log')
axes[0,1].set_title('Eccentricity vs Semimajor Axis')
axes[0,1].grid(True, alpha=0.3)

# 3. 3D scatter plot (Semimajor Axis vs Lifetime vs Eccentricity)
ax_3d = fig.add_subplot(2, 3, 3, projection='3d')
ax_3d.scatter(np.log10(semimajoraxis_array), np.log10(lifetime_array), eccentricity_array, 
              alpha=0.6, s=20, c='purple')
ax_3d.set_xlabel('log₁₀(Semimajor Axis) [AU]')
ax_3d.set_ylabel('log₁₀(Lifetime) [Myr]')
ax_3d.set_zlabel('Eccentricity')
ax_3d.set_title('3D Parameter Space')

# 4. Box plots for lifetime by semimajor axis percentile ranges
lifetime_box_data = [lifetime_array[
    (semimajoraxis_array >= ranges[i][0]) & (semimajoraxis_array <= ranges[i][1])
] for i in range(len(ranges))]
bp_lifetime = axes[1,0].boxplot(lifetime_box_data, labels=range_labels, patch_artist=True)
for patch in bp_lifetime['boxes']:
    patch.set_facecolor('lightcoral')
axes[1,0].set_ylabel('Lifetime (Myr)')
axes[1,0].set_yscale('log')
axes[1,0].set_title('Lifetime Distribution by SA Percentile Range')
axes[1,0].tick_params(axis='x', rotation=45)

# 5. Box plots for eccentricity by semimajor axis percentile ranges
ecc_box_data = [eccentricity_array[
    (semimajoraxis_array >= ranges[i][0]) & (semimajoraxis_array <= ranges[i][1])
] for i in range(len(ranges))]
bp_ecc = axes[1,1].boxplot(ecc_box_data, labels=range_labels, patch_artist=True)
for patch in bp_ecc['boxes']:
    patch.set_facecolor('lightgreen')
axes[1,1].set_ylabel('Eccentricity')
axes[1,1].set_title('Eccentricity Distribution by SA Percentile Range')
axes[1,1].tick_params(axis='x', rotation=45)

# 6. Correlation matrix visualization
correlation_data = np.column_stack([
    np.log10(semimajoraxis_array),
    np.log10(lifetime_array),
    eccentricity_array
])
correlation_matrix = np.corrcoef(correlation_data.T)
im = axes[1,2].imshow(correlation_matrix, cmap='coolwarm', aspect='auto', vmin=-1, vmax=1)
axes[1,2].set_xticks([0, 1, 2])
axes[1,2].set_yticks([0, 1, 2])
axes[1,2].set_xticklabels(['log₁₀(SA)', 'log₁₀(Lifetime)', 'Eccentricity'])
axes[1,2].set_yticklabels(['log₁₀(SA)', 'log₁₀(Lifetime)', 'Eccentricity'])
axes[1,2].set_title('Parameter Correlation Matrix')

# Add correlation values as text
for i in range(3):
    for j in range(3):
        axes[1,2].text(j, i, f'{correlation_matrix[i,j]:.2f}', 
                      ha='center', va='center', fontsize=12, 
                      color='white' if abs(correlation_matrix[i,j]) > 0.5 else 'black')

plt.colorbar(im, ax=axes[1,2], shrink=0.8)

plt.tight_layout()
plt.savefig(pathToPlots + 'comprehensive_parameter_analysis.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# Summary table for each percentile range
print(f"\n=== DETAILED PERCENTILE RANGE SUMMARY ===")
print(f"{'Range':<15} {'Count':<8} {'SA (AU)':<15} {'Lifetime (Myr)':<20} {'Eccentricity':<15}")
print(f"{'='*15} {'='*8} {'='*15} {'='*20} {'='*15}")

for i, (range_label, (min_val, max_val)) in enumerate(zip(range_labels, ranges)):
    mask = (semimajoraxis_array >= min_val) & (semimajoraxis_array <= max_val)
    count = np.sum(mask)
    
    if count > 0:
        sa_stats = f"{np.median(semimajoraxis_array[mask]):.2e}"
        lifetime_stats = f"{np.median(lifetime_array[mask]):.2e}"
        ecc_stats = f"{np.median(eccentricity_array[mask]):.3f}"
    else:
        sa_stats = "N/A"
        lifetime_stats = "N/A"
        ecc_stats = "N/A"
    
    print(f"{range_label:<15} {count:<8} {sa_stats:<15} {lifetime_stats:<20} {ecc_stats:<15}")
    
fig, ax = plt.subplots(3,1, figsize=(3.5, 7))
ax[0].hist(semimajoraxis_average, bins=50, color='blue', alpha=0.7, label='Semimajor Axis')
ax[0].set_xlabel('Semimajor Axis (AU)')
ax[0].set_ylabel('Number of Systems')
ax[0].set_yscale('log')
ax[1].hist(lifetime_average, bins=50, color='green', alpha=0.7, label='Lifetime')
ax[1].set_xlabel('Lifetime (Myr)')
ax[1].set_ylabel('Number of Systems')
ax[1].set_yscale('log')
ax[2].hist(eccentricity_average, bins=50, color='red', alpha=0.7, label='Eccentricity')
ax[2].set_xlabel('Eccentricity')
ax[2].set_ylabel('Number of Systems')
ax[2].set_yscale('log')
#plt.suptitle('BH-MS Systems Average Histogram', fontsize=16)
plt.tight_layout()
plt.savefig(pathToPlots + 'bhms_properties_histogram.png')
plt.close(fig)

fig, ax = plt.subplots(figsize=(3.5, 3.5))
# counts, bin_edges1 = np.histogram(primary_mass, bins=50)
# counts, bin_edges2 = np.histogram(secondary_mass, bins=50)
# counts1, bin_edges = np.histogram(survivor_primary_mass, bins=50)
# counts2, bin_edges = np.histogram(survivor_secondary_mass, bins=50)
# survivor_primary_mass_density = density(counts1, bins=bin_edges1)
# survivor_secondary_mass_density = density(counts2, bins=bin_edges2)
ax.hist(primary_mass, bins=50, color='blue', alpha=0.7, label='Primary Mass', histtype='step')
ax.hist(survivor_primary_mass, bins=50, color='blue', alpha=0.7, label='BHMS Primary Mass', histtype='step', linestyle='--')
ax.hist(secondary_mass, bins=50, color='orange', alpha=0.7, label='Secondary Mass', histtype='step')
ax.hist(survivor_secondary_mass, bins=50, color='orange', alpha=0.7, label='BHMS Secondary Mass', histtype='step', linestyle='--')
ax.set_xlabel('Mass (M☉)')
ax.set_ylabel('Number of Systems')
ax.set_yscale('log')
ax.legend()
#ax.set_title('Initial Mass Distribution of BH-MS Systems')
plt.tight_layout()
plt.savefig(pathToPlots + 'bhms_mass_distribution_initial.png')
plt.close(fig)

fig, ax = plt.subplots(figsize=(3.5, 3.5))
ax.hist(primary_mass, bins=50, color='blue', alpha=0.7, label='Primary Mass', density=True)
ax.hist(secondary_mass, bins=50, color='orange', alpha=0.7, label='Secondary Mass', density=True)
kroupa_fit = utils.kroupa_imf_normalized(5, 150)
ax.plot(np.linspace(5, 150, 1000), kroupa_fit, color='black', linestyle='--', label='Kroupa IMF')
ax.set_xlabel('Mass (M☉)')
ax.set_ylabel('Probability Density')
ax.set_yscale('log')
ax.legend()
plt.tight_layout()
plt.savefig(pathToPlots + 'bhms_mass_distribution_density.png')
plt.close()

# print(f'Systems with merger (R1+R2>SA): {(merger_flags_manual)}')
# print(f'Systems with merger (MT hist): {(merger_flags_mthist)}')
# print(f'Systems with merger (flag): {(merger_flags)}')
# print(f'Systems with BH-MS final state:{(bhms_final_systems)}')



# print(f'Number of systems with merger (R1+R2>SA): {len(merger_flags_manual)}')
# print(f'Number of systems with merger (MT hist): {len(merger_flags_mthist)}')
# print(f'Number of systems with merger (flag): {np.sum(merger_flags)}')
print(f'Number of systems with BH-MS final state:{len(bhms_final_systems)}')
print(f'Number of systems with BH-BH final state:{len(bhbh_final_systems)}')

print(f'Number of systems with primary BH: {primary_bh}')
print(f'Number of systems with secondary BH: {secondary_bh}')
print(f'Total number of systems with BHs: {i}')
print(f'Total number of systems with BH-MS pairs: {len(lifetime_average)}')

print(f'Skipped systems due to absence of BH-MS pairs: {j}')
print(f'Skipped systems due to negative semimajor axis: {k}')

print(f'Minimum primary mass: {np.min(survivor_primary_mass)} M☉')
print(f'Maximum primary mass: {np.max(survivor_primary_mass)} M☉')

# Displays Time
e = dt.datetime.now()
current_time = e.strftime("%d%m%y") + "_" + e.strftime('%H%M')
duration = e - s
print(" \n Finish time :", current_time)
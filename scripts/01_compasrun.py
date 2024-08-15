#DESCRIPTION: This script calls COMPAS via terminal and runs it with given program options. 
#MUST be executed in the directory where COMPAS environment is built.

import os
import time as t
import os.path
import numpy as np

#Choose initial mass function (IMF) : ['KROUPA', 'SALPETER', 'UNIFORM, 'POWERLAW']
IMF = 'KROUPA'
#If POWERLAW is chosen, providing a steepness value alpha, add option --initial-mass-power=alpha
alpha = '-2.3'
#Number of systems per run is n, number of runs is N. Total number of systems is n * N.
n = '360000'
N = 50

starttime = t.ctime(t.time())
start = t.process_time()

for i in range(0,N):

    os.system("COMPAS --mode=BSE --number-of-systems=" + n + " --output-path=/data/a.saricaoglu/Runs/" + IMF + " --initial-mass-function=" + IMF)

endtime = t.ctime(t.time())
duration = t.process_time() - start

if not os.path.exists("/data/a.saricaoglu/Runs/logfiles/" + IMF + "/" +  str(c.strftime("%m.%d"))): 
    os.makedirs("/data/a.saricaoglu/Runs/logfiles/" + IMF + "/" +  str(c.strftime("%m.%d"))) 
  
#Save execution logs to file
f = ["Run :" + str(i) + ";  Started : " + str(starttime), ";  Ended: " + str(endtime) + ";  Duration : " + str(duration)]
np.savetxt("/data/a.saricaoglu/Runs/logfiles/" + IMF + "_" + str(c.strftime("%m.%d")) + ".txt", f)

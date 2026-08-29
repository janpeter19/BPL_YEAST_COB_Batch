# setup applicateion data BPL_YEAST_COB_Batch 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-08-24 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import matplotlib.pyplot as plt 
from pyfmi import load_fmu

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   flag_vendor = 'JM'
   flag_type = 'CS'
   fmu_model ='BPL_YEAST_COB_Batch_windows_jm_cs.fmu'        
   model = load_fmu(fmu_model, log_level=0)  
elif platform.system() == 'Linux':  
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-compiled OpenModelica') 
      if flag_type in ['CS','cs']:         
         fmu_model ='BPL_YEAST_COB.Batch_linux_om_cs.fmu'    
         model = load_fmu(fmu_model, log_level=0) 
      if flag_type in ['ME','me']:         
         fmu_model ='BPL_YEAST_COB.Batch_linux_om_me.fmu'    
         model = load_fmu(fmu_model, log_level=0)
   else:    
      print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = model.simulate_options()
   opts_std['silent_mode'] = True
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary' 
   opts_fast = model.simulate_options()   
   opts_fast['silent_mode'] = True
   opts_fast['ncp'] = 10 
   opts_fast['result_handling'] = 'memory'     
elif flag_type in ['ME', 'me']:
   opts_std = model.simulate_options()
   opts_std["CVode_options"]["verbosity"] = 50 
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary' 
   opts_fast = model.simulate_options()   
   opts_fast["CVode_options"]["verbosity"] = 50 
   opts_fast['ncp'] = 10 
   opts_fast['result_handling'] = 'memory'  
else:    
   print('There is no FMU for this platform')
  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
   MSL_usage = model.get('MSL.usage')[0]
   MSL_version = model.get('MSL.version')[0]
   BPL_version = model.get('BPL.version')[0]
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '4.1.0 - used components: none' 
   MSL_version = '4.1.0'
   BPL_version = 'Bioprocess Library version 2.3.2' 
else:    
   print('There is no FMU for this platform')

# Simulation time
simulationTime = 12.0
prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = model.get_states_list()
stateValue.update(timeDiscreteStates)

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture']

# Provide process diagram on disk
fmu_process_diagram ='BPL_YEAST_COB_Batch_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, ax, lines
#------------------------------------------------------------------------------------------------------------------
   
# Create dictionaries parValue and parLocation
parValue = {}
parValue['V_start'] = 4.5
parValue['VX_start'] = 1.0
parValue['VG_start'] = 10.0
parValue['VE_start'] = 0.0

parValue['mum'] = 0
parValue['qGr'] = 0
parValue['qEr'] = 0
parValue['qO2'] = 0

parLocation = {}
parLocation['V_start'] = 'bioreactor.V_start'
parLocation['VX_start'] = 'bioreactor.m_start[1]'
parLocation['VG_start'] = 'bioreactor.m_start[2]'
parLocation['VE_start'] = 'bioreactor.m_start[3]'

parLocation['mum'] = 'bioreactor.culture.mum'
parLocation['qGr'] = 'bioreactor.culture.qGr'
parLocation['qEr'] = 'bioreactor.culture.qEr'
parLocation['qO2'] = 'bioreactor.culture.qO2'

# Extra only for describe()
parLocation['mu'] = 'bioreactor.culture.mu'

# Parameter value check - especially for hysteresis to avoid runtime error
parCheck = []
parCheck.append("parValue['V_start'] > 0")
parCheck.append("parValue['VX_start'] >= 0")
parCheck.append("parValue['VG_start'] >= 0")

# Create list of diagrams to be plotted by simu()
diagrams = []

# Create an empty list axes to be defined in newplot() and plotted by simu() or show()
ax = []

# Create list of pens for the diagrams
lines = ['-','--',':','-.']
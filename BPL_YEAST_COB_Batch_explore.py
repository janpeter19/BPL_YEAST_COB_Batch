# setup application functions BPL_YEAST_COB_batch, dependent on previous import of functions from fmu_explore 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-08-24 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Specific application functions: newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

# Define standard plots
def newplot(title='Batch cultivation', plotType='TimeSeries'):
   """ Standard plot window 
       title = '' """
    
   # Reset pens
   resetPen()
   
   # Plot diagram 
   if plotType == 'TimeSeries':
         
      ax1 = plt.subplot(3,1,1)
      ax2 = plt.subplot(3,1,2)
      ax3 = plt.subplot(3,1,3)
      
      ax.clear()
      ax.append(ax1)
      ax.append(ax2)
      ax.append(ax3)
    
      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('G [g/L]')
    
      ax[1].grid()
      ax[1].set_ylabel('E [1/h]')
      
      ax[2].grid()
      ax[2].set_ylabel('X [g/L]')      
      ax[2].set_xlabel('Time [h]') 
      
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['bioreactor.c[2]'], label='G', color='b', linestyle=linetype)") 
      diagrams.append("ax[1].plot(t,sim_res['bioreactor.c[3]'], label='E', color='b', linestyle=linetype)") 
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.c[1]'], label='X', color='b', linestyle=linetype)") 
      
      
   elif plotType == 'TimeSeries2':
         
      ax1 = plt.subplot(3,1,1)
      ax2 = plt.subplot(3,1,2)
      ax3 = plt.subplot(3,1,3)

      ax.clear()
      ax.append(ax1)
      ax.append(ax2)
      ax.append(ax3)
    
      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('X [g/L]')
    
      ax[1].grid()
      ax[1].set_ylabel('mu [1/h]')
      
      ax[2].grid()
      ax[2].set_ylabel('G, E [g/L]')      
      ax[2].set_xlabel('Time [h]') 
      
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['bioreactor.c[1]'], label='X', color='r', linestyle=linetype)") 
      diagrams.append("ax[1].step(t,sim_res['bioreactor.culture.mu'], label='mu', color='r', linestyle=linetype)")     
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.c[2]'], label='G', color='b', linestyle=linetype)") 
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.c[3]'], label='E', color='g', linestyle=linetype)") 
      diagrams.append("ax[2].legend(['G','E'])")
 
   elif plotType == 'Extended':
         
      # Transfer of argument to global variable       
      plt.figure()
      ax11 = plt.subplot(2,2,1)
      ax12 = plt.subplot(2,2,2)
      ax21 = plt.subplot(2,2,3)
      ax22 = plt.subplot(2,2,4)
    
      ax11.set_title(title)
      ax11.grid()
      ax11.set_ylabel('G and E [g/L]')
    
      ax21.grid()
      ax21.set_ylabel('X [g/L]')
      ax21.set_xlabel('Time [h]')
      
      ax12.grid()
      ax12.set_ylabel('mu [1/h]')  
      
      ax22.grid()
      ax22.set_ylabel('qO2 [mole/(g*h)]')  
      ax22.set_xlabel('Time [h]')     
      
      diagrams.clear()
      diagrams.append("ax11.plot(t,sim_res['bioreactor.c[2]'], label='G', color='b', linestyle='-')") 
      diagrams.append("ax11.plot(t,sim_res['bioreactor.c[3]'], label='E', color='r', linestyle='-')") 
      diagrams.append("ax21.plot(t,sim_res['bioreactor.c[1]'], label='X', color='b', linestyle='-')")       
      diagrams.append("ax12.plot(t,sim_res['bioreactor.culture.mu'], label='mu', color='b', linestyle='-')")
      
      diagrams.append("ax22.plot(t,sim_res['bioreactor.culture.qO2'], 'b-')")
      
   elif diagram == 'PhasePlane':
      plt.figure()
      ax1 = plt.subplot(1,1,1)
                
      ax1.set_title(title)
      ax1.grid()
      ax1.set_ylabel('G')
      ax1.set_xlabel('X')
                   
   else:
      print("Plot window type not correct") 

# Define and extend describe for the current application
def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""
        
   if name == 'culture':
      print('Saccharomyces cerevisae - default parameters for strain H1022')        
     
   elif name in ['broth', 'liquidphase', 'media']:
      X = model.get('liquidphase.X')[0]; 
      X_description = model.get_variable_description('liquidphase.X'); 
      X_mw = model.get('liquidphase.mw[1]')[0]
      
      G = model.get('liquidphase.G')[0]; 
      G_description = model.get_variable_description('liquidphase.G'); 
      G_mw = model.get('liquidphase.mw[2]')[0]
      
      E = model.get('liquidphase.E')[0]; 
      E_description = model.get_variable_description('liquidphase.E'); 
      E_mw = model.get('liquidphase.mw[3]')[0]

      print('Reactor broth substances included in the model')
      print()
      print(X_description, '  index = ', X, 'molecular weight = ', X_mw, 'Da')
      print(G_description, 'index = ', G, 'molecular weight = ', G_mw, 'Da')
      print(E_description, 'index = ', E, 'molecular weight = ', E_mw, 'Da')
     
   elif name in ['parts']:
      describe_parts(component_list_minimum)

   elif name in ['MSL']:
      describe_MSL()

   else:
      describe_general(name, decimals)
      
#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

FMU_explore_info()
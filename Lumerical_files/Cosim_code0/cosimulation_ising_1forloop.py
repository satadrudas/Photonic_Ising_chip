
import sys, os, random, pdb
import numpy as np

## Uncomment the following if you are using Linux
sys.path.append("/opt/lumerical/interconnect/api/python") # linux


import lumapi

# Loading initial parameters
######################################################################
alpha=0.5
beta=0.5
N_iter=50 # number of iterations
N_runs=1 # for averaging the results for a given hyperparameters

alpha_min=0
alpha_max=3
alpha_step=0.1

beta_min=0
beta_max=3
beta_step=0.01

std=0.05

seed = np.random.randint(100000)
np.random.seed(seed)






# The block below loads the J and N
########################################################################
J_file='/home/satadrudas/lumerical_files/CIM/Maxcut_instances/s_100.txt'

if(J_file):
    f = open(J_file,"r")
    N, number_of_edges = [int(i) for i in f.readline().split()]
    J = np.zeros([N,N])
    lines = f.readlines()
    for line in lines:
        l = line.split()
        r,c,w = int(l[0])-1, int(l[1])-1, float(l[2])
        J[r][c] = w
        J[c][r] = w        
    f.close()
J=beta*(-J)
np.fill_diagonal(J,alpha)

J_flat=J.flatten()
J_max=max(alpha,beta)
J_theta=np.arccos(J_flat/(J_max+1e-15)) ##dividing to normalize J for mzm







# Loading Initial voltages of the MZM, sampling size and other necessary parameters
#####################################################################################

mzm_V_pi=1
mzm_freq=10e9

# Time steps needed
#
# VVM + readout per iteration: N*(N+1)
# Final spin state readouts: N
# Total timesteps: (((N*(N+1))*N_iter)*N_runs)+N

ising_timesteps=(((N*(N+1))*N_iter)*N_runs)
spin_readout_timesteps=2*N
total_timesteps= ising_timesteps+spin_readout_timesteps



ising_time=(1/mzm_freq)*ising_timesteps # in seconds
spin_readout_time=(1/mzm_freq)*spin_readout_timesteps # in seconds
time_window= ising_time + spin_readout_time# total time in seconds

#Nsamples=8192*int(time_window/(5.12e-09)) ## 5.12e-09 seconds have 8192 samples, how many samples does time_window seconds have?

samples_num=100 #samples per data
Nsamples=samples_num*total_timesteps #it means that every (1/mzm_freq) seconds has 1000 samples

homodyne_gating_freq=mzm_freq/N # in Hz

final_spin_values=np.zeros(N)
x_in=np.zeros(N)
x_out=np.zeros(N)
x_out_flat=np.tile(x_out,N)

theta1=x_out_flat-(np.pi/2)+np.random.normal(scale=std,size=(len(x_out_flat)))
theta2=J_theta


mzm_1_v1=2*mzm_V_pi+(mzm_V_pi*theta1/np.pi)
mzm_1_v2=2*mzm_V_pi-(mzm_V_pi*theta1/np.pi)

mzm_2_v1=2*mzm_V_pi+(mzm_V_pi*theta2/np.pi) 
mzm_2_v2=2*mzm_V_pi-(mzm_V_pi*theta2/np.pi) 






# Initializing the lumapi simulation parameters
###############################################################################
h = lumapi.open("interconnect")
lumapi.evalScript(h,''' load("ising.icp"); ''');#if it doesn't open, make sure that the the startup directory of a net INTERCONNECT is corrected

# intialize simply by pushing data to Lumerical script workspace and using same script as Lumerical script example
lumapi.putDouble(h,"time_window",time_window)
lumapi.putDouble(h,"Nsamples",Nsamples)
lumapi.putDouble(h,"mzm_V_pi",mzm_V_pi)

lumapi.evalScript(h,'''
    # turn off history and redrawing
    redrawoff;
    historyoff;

    # time vectors
    t = linspace(0,time_window,Nsamples+1);
    t = t(1:end-1);
    dt = t(2)-t(1);

    # restore design mode in case we are in analysis mode
    switchtodesign;

    # Set up simulation to match time window and number of samples
    # The time windows should match because the INTERCONNECT simulation
    # will stop once the simulation time exceeds the time window.
    # The number of samples defines the INTERCONNECT time step, dt, by dt = time_window/(Nsamples+1).
    # The time steps do NOT have to match, although in this example they do. Indeed,
    # the time step of an external simulator can be variable
    setnamed("::Root Element","time window",time_window);
    setnamed("::Root Element","number of samples",Nsamples);

    # set number of threads to 2. More can run faster but testing is required and is machine dependent
    setnamed("::Root Element","number of threads",2);

    # set the desired attenuation
    setnamed("MZM_1","pi dc voltage",mzm_V_pi);
    setnamed("MZM_2","pi dc voltage",mzm_V_pi);

    # intialize the simulation
    runinitialize;

    # initialize struct used for passing co-sim data
    electronics_control = struct;


    # End initialization
    #####################################################
    
    
    ''')

time = lumapi.getVar(h,"t") # extracting time from lumapi
dt=lumapi.getVar(h,"dt") # extracting dt from lumapi

#the voltages have to be sent to INTERCONNECT for the cosim struct to get those values....we can send the entire vector itself
lumapi.putMatrix(h,"mzm_1_v1",mzm_1_v1)
lumapi.putMatrix(h,"mzm_1_v2",mzm_1_v2)
lumapi.putMatrix(h,"mzm_2_v1",mzm_1_v1) #we dont have to update these since these are just the J_flat, except when we want to read
lumapi.putMatrix(h,"mzm_2_v2",mzm_2_v2) #we dont have to update these since these are just the J_flat, except when we want to read




#The loop is as follows:

# while loop, total timesteps=(N(N+1)*N_iter)+N
#
#loop N_iter times:
#	loop N time:
# 		VVM
# 		1 redout 
#		reset after N readouts for the next itertion
#
#
#loop N times:
#	1 spin modulation in the mzm1, mzm2 set to 1 all the time now
#	1 spin state readout

counter_mzm=0 # the index for the mzm voltage vectors, ranges from 0 to N**2
counter_per_iter=0 # counts the elements as the MAC happens, ranges from 0 to N+1
count_iter=0 # keeps track of the iteratin, ranges from 0 to N
counter2=0
final_read_counter=0
read_flag=0 #bool


for t in time:
	lumapi.putDouble(h,"tval",t)
	lumapi.putDouble(h,"counter_mzm",counter_mzm)

	if(count_iter==N_iter):

		if (read_flag):
			
			readout_sample_increament_count=0

			lumapi.evalScript(h,'''
			
			#read the integrator
			electronics_control.index=6;
			integrator= getvalue("COSIM_1","port",electronics_control);
			
			
			
			#reset port_5=1
			electronics_control.index=5;
			electronics_control.time=tval;
			electronics_control.value=1;
			setvalue("COSIM_1","port",electronics_control);
			
			
			runstep;
			
			#get value from the integrator as struct in the variable 'integrator'

			
			
			''')

			if(sample_counter%samples_num==0):
				readout_sample_increament_count=readout_sample_increament_count+1
				
				if (readout_sample_increament_count==1):
					final_spin_values[counter2]=lumapi.getVar(h,"integrator")
					counter2=counter2+1

					sample_counter=sample_counter+1

					if (counter2==N):
						break

				else:
					read_flag=0 #exit

					lumapi.evalScript('''
					
					
					#reset port_5=0
					electronics_control.index=5;
					electronics_control.time=tval;
					electronics_control.value=0;
					setvalue("COSIM_1","port",electronics_control);
					
					
					#set all the mzm values
					
					electronics_control.index=1;
					electronics_control.time=tval;
					electronics_control.value=mzm_1_v1(counter_mzm);
					setvalue("COSIM_1","port",electronics_control);
					
					electronics_control.index=2;
					electronics_control.time=tval;
					electronics_control.value=mzm_1_v2(counter_mzm);
					setvalue("COSIM_1","port",electronics_control);
					
					electronics_control.index=3;
					electronics_control.time=tval;
					electronics_control.value=mzm_2_v1(counter_mzm);
					setvalue("COSIM_1","port",electronics_control);
					
					electronics_control.index=4;
					electronics_control.time=tval;
					electronics_control.value=mzm_2_v2(counter_mzm);
					setvalue("COSIM_1","port",electronics_control);	
				
					
					
					runstep;		
					
					''')

					sample_counter=sample_counter+1



			else:
				sample_counter=sample_counter+1
				continue


		
		else:

			lumapi.evalScript('''
			
			
			#reset port_5=0
			electronics_control.index=5;
			electronics_control.time=tval;
			electronics_control.value=0;
			setvalue("COSIM_1","port",electronics_control);
			
			
			#set all the mzm values
			
			electronics_control.index=1;
			electronics_control.time=tval;
			electronics_control.value=mzm_1_v1(counter_mzm);
			setvalue("COSIM_1","port",electronics_control);
			
			electronics_control.index=2;
			electronics_control.time=tval;
			electronics_control.value=mzm_1_v2(counter_mzm);
			setvalue("COSIM_1","port",electronics_control);
			
			electronics_control.index=3;
			electronics_control.time=tval;
			electronics_control.value=mzm_2_v1(counter_mzm);
			setvalue("COSIM_1","port",electronics_control);
			
			electronics_control.index=4;
			electronics_control.time=tval;
			electronics_control.value=mzm_2_v2(counter_mzm);
			setvalue("COSIM_1","port",electronics_control);	
		
			
			
			runstep;		
			
			''')

			sample_counter=sample_counter+1

			if (sample_counter%samples_num==0):
				counter_mzm=counter_mzm+1 
				read_flag=1





	if (counter_per_iter%N==0): # just for reading and reset after every vector-vector multiplication,  

		readout_sample_increament_count=0
		lumapi.evalScript(h,'''
		
		#read the integrator
		electronics_control.index=6;
		integrator= getvalue("COSIM_1","port",electronics_control);
		
		
		
		#reset port_5=1
		electronics_control.index=5;
		electronics_control.time=tval;
		electronics_control.value=1;
		setvalue("COSIM_1","port",electronics_control);
		
		
		runstep;
		
		#get value from the integrator as struct in the variable 'integrator'

		
		
		''')
		
		if (sample_counter%samples_num==0): # this condition is carried on from last iteration
			readout_sample_increament_count=readout_sample_increament_count+1

			if (readout_sample_increament_count==1):
				x_in[counter2]=lumapi.getVar(h,"integrator")
				counter2=counter2+1

				sample_counter=sample_counter+1

				if (counter2==N):
					x_out=x_in #update the x_out with the new spin states
					x_out_flat=np.tile(x_out,N)
					## show I update the mzm_1 output string here itself??? ...yes.
					theta1=x_out_flat-(np.pi/2)+np.random.normal(scale=std,size=(len(x_out_flat)))
					mzm_1_v1=2*mzm_V_pi+(mzm_V_pi*theta1/np.pi)
					mzm_1_v2=2*mzm_V_pi-(mzm_V_pi*theta1/np.pi)
					lumapi.putMatrix(h,"mzm_1_v1",mzm_1_v1)
					lumapi.putMatrix(h,"mzm_1_v2",mzm_1_v2)

					counter2=0 # reset the recording counter
					counter_mzm=0 #ranges from 0 to N**2
					count_iter=count_iter+1

					if (count_iter==N_iter):

						x_final=x_out[0:N] 
						theta1=x_final-(np.pi/2)+np.random.normal(scale=std,size=(len(x_final)))
						mzm_1_v1=2*mzm_V_pi+(mzm_V_pi*theta1/np.pi)
						mzm_1_v2=2*mzm_V_pi-(mzm_V_pi*theta1/np.pi)
						
						mzm_2_v1=np.zeros(N)
						mzm_2_v2=np.zeros(N)
						
						#below dont worry about the size, they get resized on new assignment
						lumapi.putMatrix(h,"mzm_1_v1",mzm_1_v1)
						lumapi.putMatrix(h,"mzm_1_v2",mzm_1_v2)
						lumapi.putMatrix(h,"mzm_2_v1",mzm_1_v1) 
						lumapi.putMatrix(h,"mzm_2_v2",mzm_2_v2)

						counter_mzm=0
						sample_counter=0
						counter2=0
						final_read_counter=0
						read_flag=0 #bool

						continue # in the case of an big loop, it will be continue
				
				

			else: # this else condition is the final condition of this loop, once this satisfies then the control exits from the loop, the gets into the outermost else condition
				counter_per_iter=0 #counter_per_iter+1 # basically becomes N+1, .....N+1 end is basically 0 begining

				lumapi.evalScript(h,'''

				
				#reset port_5=0
				electronics_control.index=5;
				electronics_control.time=tval;
				electronics_control.value=0;
				setvalue("COSIM_1","port",electronics_control);
				
				
				#set all the mzm values
				
				electronics_control.index=1;
				electronics_control.time=tval;
				electronics_control.value=mzm_1_v1(counter_mzm);
				setvalue("COSIM_1","port",electronics_control);
				
				electronics_control.index=2;
				electronics_control.time=tval;
				electronics_control.value=mzm_1_v2(counter_mzm);
				setvalue("COSIM_1","port",electronics_control);
				
				electronics_control.index=3;
				electronics_control.time=tval;
				electronics_control.value=mzm_2_v1(counter_mzm);
				setvalue("COSIM_1","port",electronics_control);
				
				electronics_control.index=4;
				electronics_control.time=tval;
				electronics_control.value=mzm_2_v2(counter_mzm);
				setvalue("COSIM_1","port",electronics_control);	
				
				runstep;	
				
				
				''')	

				sample_counter=sample_counter+1



		else:
			sample_counter=sample_counter+1
			
			continue





			
			
	else:
		lumapi.evalScript(h,'''

		
		#reset port_5=0
		electronics_control.index=5;
		electronics_control.time=tval;
		electronics_control.value=0;
		setvalue("COSIM_1","port",electronics_control);
		
		
		#set all the mzm values
		
		electronics_control.index=1;
		electronics_control.time=tval;
		electronics_control.value=mzm_1_v1(counter_mzm);
		setvalue("COSIM_1","port",electronics_control);
		
		electronics_control.index=2;
		electronics_control.time=tval;
		electronics_control.value=mzm_1_v2(counter_mzm);
		setvalue("COSIM_1","port",electronics_control);
		
		electronics_control.index=3;
		electronics_control.time=tval;
		electronics_control.value=mzm_2_v1(counter_mzm);
		setvalue("COSIM_1","port",electronics_control);
		
		electronics_control.index=4;
		electronics_control.time=tval;
		electronics_control.value=mzm_2_v2(counter_mzm);
		setvalue("COSIM_1","port",electronics_control);	
		
		runstep;	
		
		
		''')	

		sample_counter=sample_counter+1

		if (sample_counter%samples_num==0):
			counter_mzm=counter_mzm+1 # remember, the count_mzm only gets updated here till it become N**2 after which it gets reset to 0 
			counter_per_iter=counter_per_iter+1


#final solution is in the "final_spin_values" array
		
		



# def mzm_modulation_step(time_stamp,mzm_index):

# 	lumapi.putDouble(h,"tval",time_stamp)
#	lumapi.putDouble(h,"counter_mzm",mzm_index)

# 	lumapi.evalScript(h,'''

	
# 	#reset port_5=0
# 	electronics_control.index=5;
# 	electronics_control.time=tval;
# 	electronics_control.value=0;
# 	setvalue("COSIM_1","port",electronics_control);
	
	
# 	#set all the mzm values
	
# 	electronics_control.index=1;
# 	electronics_control.time=tval;
# 	electronics_control.value=mzm_1_v1(counter_mzm);
# 	setvalue("COSIM_1","port",electronics_control);
	
# 	electronics_control.index=2;
# 	electronics_control.time=tval;
# 	electronics_control.value=mzm_1_v2(counter_mzm);
# 	setvalue("COSIM_1","port",electronics_control);
	
# 	electronics_control.index=3;
# 	electronics_control.time=tval;
# 	electronics_control.value=mzm_2_v1(counter_mzm);
# 	setvalue("COSIM_1","port",electronics_control);
	
# 	electronics_control.index=4;
# 	electronics_control.time=tval;
# 	electronics_control.value=mzm_2_v2(counter_mzm);
# 	setvalue("COSIM_1","port",electronics_control);	
	
# 	runstep;	
	
	
# 	''')		
	

# def vvm_readout_step(time_stamp):

# 	lumapi.putDouble(h,"tval",time_stamp)
		
# 	lumapi.evalScript(h,'''
	
# 		#read the integrator
# 		electronics_control.index=6;
# 		integrator= getvalue("COSIM_1","port",electronics_control);
		
		
		
# 		#reset port_5=1
# 		electronics_control.index=5;
# 		electronics_control.time=tval;
# 		electronics_control.value=1;
# 		setvalue("COSIM_1","port",electronics_control);
		
		
# 		runstep;
		
# 		#get value from the integrator as struct in the variable 'integrator'

		
		
# 		''')
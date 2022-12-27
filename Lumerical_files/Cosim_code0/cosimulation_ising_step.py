
import sys, os, random, pdb
import numpy as np

## Uncomment the following if you are using Linux
sys.path.append("/opt/lumerical/interconnect/api/python") # linux


import lumapi

######################################################################
# LOADING INITIAL PARAMETERS
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

N_iters=20   
# The block above loads the J and N
###########################################################################

mzm_V_pi=1
mzm_freq=10e9
time_window=(1/mzm_freq)*((N+1)*(N+1))*N_iters*N_runs
Nsamples=8192*int(time_window/(5.12e-09)) ## or instead Nsamples=time_window*5, 5 samples per signal ....5.12e-09 s has 8192 samples

homodyne_gating_freq=mzm_freq/N



x_out=np.zeros(N)
x_out_flat=np.tile(x_out,N)

theta1=x_out_flat-(np.pi/2)+np.random.normal(scale=std,size=(len(x_out_flat)))
theta2=J_theta


mzm_1_v1=2*mzm_V_pi+(mzm_V_pi*theta1/np.pi)
mzm_1_v2=2*mzm_V_pi-(mzm_V_pi*theta1/np.pi)

mzm_2_v1=2*mzm_V_pi+(mzm_V_pi*theta2/np.pi) 
mzm_2_v2=2*mzm_V_pi-(mzm_V_pi*theta2/np.pi) 


###############################################################################
h = lumapi.open("interconnect")
lumapi.evalScript(h,''' load("ising.icp"); ''');#if it doesn't open, make sure that the the startup directory of a net INTERCONNECT is corrected

# intialize simply by pushing data to Lumerical script workspace and using same script as Lumerical script example
lumapi.putDouble(h,"time_window",time_window);
lumapi.putDouble(h,"Nsamples",Nsamples);
lumapi.putDouble(h,"mzm_V_pi",mzm_V_pi);

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

time = lumapi.getVar(h,"t")
dt=lumapi.getVar(h,"dt")

#the voltages have to be sent to INTERCONNECT for the cosim struct to get those values....we can send the entire vector itself
lumapi.putMatrix(h,"mzm_1_v1",mzm_1_v1)
lumapi.putMatrix(h,"mzm_1_v2",mzm_1_v2)
lumapi.putMatrix(h,"mzm_2_v1",mzm_1_v1) #we dont have to update these since these are just the J_flat, except when we want to read
lumapi.putMatrix(h,"mzm_2_v2",mzm_2_v2) #we dont have to update these since these are just the J_flat, except when we want to read


counter=0
counter2=0
count_iter=0

for t in time:###add another loop for the samples

	lumapi.putDouble(h,"tval",t)
	lumapi.putDouble(h,"counter",counter+1)
	
	
	if (count_iter==N_iters): ### here no need samples for now....jus do for a singe samples
		
		x_final=x_out[0:N] ### i dont think we can do just reading of x_out fro the final state because it is alpha*x_in+beta*J*rest, the final spin value is cosine of that
		# can we just take cos of the it digitally instead of the mzi
		theta1=x_final-(np.pi/2)+np.random.normal(scale=std,size=(len(x_final)))
		mzm_1_v1=2*mzm_V_pi+(mzm_V_pi*theta1/np.pi)
		mzm_1_v2=2*mzm_V_pi-(mzm_V_pi*theta1/np.pi)
		
		mzi_2_v1=np.zeros(N)
		mzi_2_v2=np.zeros(N)
		
		#below dont worry about the size, they get resized on new assignment
		lumapi.putMatrix(h,"mzm_1_v1",mzm_1_v1)
		lumapi.putMatrix(h,"mzm_1_v2",mzm_1_v2)
		lumapi.putMatrix(h,"mzm_2_v1",mzm_1_v1) 
		lumapi.putMatrix(h,"mzm_2_v2",mzm_2_v2)
				
		index=0
		
		
		time_stamp=t
		
		for i in range(N):
			lumapi.putDouble(h,"tval",time_stamp)
		
		
			lumapi.evalScript('''
			
			
			#reset port_5=0
			electronics_control.index=5;
			electronics_control.time=tval;
			electronics_control.value=0;
			setvalue("COSIM_1","port",electronics_control);
			
			
			#set all the mzm values
			
			electronics_control.index=1;
			electronics_control.time=tval;
			electronics_control.value=mzm_1_v1(index);
			setvalue("COSIM_1","port",electronics_control);
			
			electronics_control.index=2;
			electronics_control.time=tval;
			electronics_control.value=mzm_1_v2(index);
			setvalue("COSIM_1","port",electronics_control);
			
			electronics_control.index=3;
			electronics_control.time=tval;
			electronics_control.value=mzm_2_v1(index);
			setvalue("COSIM_1","port",electronics_control);
			
			electronics_control.index=4;
			electronics_control.time=tval;
			electronics_control.value=mzm_2_v2(index);
			setvalue("COSIM_1","port",electronics_control);	
			
			runstep;	

			#read the integrator
			electronics_control.index=6;
			integrator= getvalue("COSIM_1","port",electronics_control);
			
			
			
			#reset port_5=1
			electronics_control.index=5;
			electronics_control.time=tval;
			electronics_control.value=1;
			setvalue("COSIM_1","port",electronics_control);
			
			
			runstep;		
			
			
			
			''')
			
			time_stamp=time_stamp+dt
			
			
	## no need samples here also or should I?????? just for consistency
	if (counter%N==0): # just for reading and reset,  
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
		x_in[counter2]=lumapi.getVar(h,"integrator")
		counter2=counter2+1

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
			counter=0 #ranges from 0 to N**2
			count_iter=count_iter+1
			
			
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
		electronics_control.value=mzm_1_v1(counter);
		setvalue("COSIM_1","port",electronics_control);
		
		electronics_control.index=2;
		electronics_control.time=tval;
		electronics_control.value=mzm_1_v2(counter);
		setvalue("COSIM_1","port",electronics_control);
		
		electronics_control.index=3;
		electronics_control.time=tval;
		electronics_control.value=mzm_2_v1(counter);
		setvalue("COSIM_1","port",electronics_control);
		
		electronics_control.index=4;
		electronics_control.time=tval;
		electronics_control.value=mzm_2_v2(counter);
		setvalue("COSIM_1","port",electronics_control);	
		
		runstep;	
		
		
		''')	
		
	counter=counter+1
	






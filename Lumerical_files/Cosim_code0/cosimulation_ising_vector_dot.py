
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

    
# The block above loads the J and N
###########################################################################

mzm_V_pi=1
mzm_freq=10e9
time_window=(1/mzm_freq)*((N+1)*(N+1))*N_iters*N_runs #+1 for reading slot
Nsamples=8192*int(time_window/(5.12e-09)) ## or instead Nsamples=time_window*5, 5 samples per signal ....5.12e-09 s has 8192 samples

homodyne_gating_freq=mzm_freq/N



x_out=np.zeros(N)
x_out_flat=np.tile(x_out,N)

theta1=x_out_flat-(np.pi/2)+np.random.normal(scale=std,size=(len(self.x_out_flat)))
theta2=self.J_theta


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
lumapi.putMatrix(h,"mzm_2_v1",mzm_1_v1 #we dont have to update these since these are just the J_flat, except when we want to read
lumapi.putMatrix(h,"mzm_2_v2",mzm_2_v2) #we dont have to update these since these are just the J_flat, except when we want to read


############################################
#do you think it will be a good idea to repeat the elements of the mzm_1_v1 to match the time????????????????????????????????????????????????????
#or just do it in the loop?????????????????
###########################################
counter=0
counter2=0


'''
do a while loop instead
while(t<=time[-1])
	t=t+dt
'''

t=time[0]
while(t<=time[-1]): ###add another loop for the samples

	if (counter%gating_count==0 and counter>0): # just for reading and reset
		lumapi.evalScript(h,'''
		reset port=1
		
		
		runstep; ### is it needed to just read the port???? yes i guess because of the reset integrator
		
		#get value from the integrator as struct in the variable 'integrator'

		
		
		''')
		x_in[counter2]=lumapi.getVar(h,"integrator")
		counter2=counter2+1
		
		if (counter2==N):
			x_out=x_in #update the x_out with the new spin states
			x_out_flat=np.tile(x_out,N)
			## show I update the mzm_1 output string here itself??? ...yes.
			theta1=x_out_flat-(np.pi/2)+np.random.normal(scale=std,size=(len(self.x_out_flat)))
			mzm_1_v1=2*mzm_V_pi+(mzm_V_pi*theta1/np.pi)
			mzm_1_v2=2*mzm_V_pi-(mzm_V_pi*theta1/np.pi)
			lumapi.putMatrix(h,"mzm_1_v1",mzm_1_v1)
			lumapi.putMatrix(h,"mzm_1_v2",mzm_1_v2)
			
			
			
			counter2=0 # reset the recording counter
			
		t=t+dt
			
			
	else:
		lumapi.evalScript(h,'''
		reset port=0
		#we already put th mzm_volts using putMatrix befor the while loop
		####give input to the mzm as co-sim struct (N*N) long vector
		

		#now just loop in time for 1 vector-vector dot product
		# for ()
		
		
		#pull back time variable from the INTERCONNECT??????
		
		
		
		runstep;
		''')		
		
		
	t=t+dt
	counter=counter+1
	if(counter=np.square(N)):
		counter=0

##how to update t... t=t+v.time???


##make anothe loop for readout














counter=0
counter2=0




t=time[0]
while(t<=time[-1]): ###add another loop for the samples

	if (counter%gating_count==0 and counter>0): # just for reading and reset
		lumapi.evalScript(h,'''
		reset port=1
		
		
		runstep; ### is it needed to just read the port???? yes i guess because of the reset integrator
		
		#get value from the integrator as struct in the variable 'integrator'

		
		
		''')
		x_in[counter2]=lumapi.getVar(h,"integrator")
		counter2=counter2+1
		
		if (counter2==N):
			x_out=x_in #update the x_out with the new spin states
			x_out_flat=np.tile(x_out,N)
			## show I update the mzm_1 output string here itself??? ...yes.
			theta1=x_out_flat-(np.pi/2)+np.random.normal(scale=std,size=(len(self.x_out_flat)))
			mzm_1_v1=2*mzm_V_pi+(mzm_V_pi*theta1/np.pi)
			mzm_1_v2=2*mzm_V_pi-(mzm_V_pi*theta1/np.pi)
			lumapi.putMatrix(h,"mzm_1_v1",mzm_1_v1)
			lumapi.putMatrix(h,"mzm_1_v2",mzm_1_v2)
			
			
			
			counter2=0 # reset the recording counter
			
			
	else:
		lumapi.evalScript(h,'''
		reset port=0
		#give input to the mzm as co-sim struct (N*N) long vector
		
		#we already put th mzm_volts using putMatrix...
		#now just loop in time for 1 vector-vector dot product
		
		pull back time variable from the INTERCONNECT??????
		
		
		
		runstep;
		''')		
		
		
	t=t+dt
	counter=counter+1
	if(counter=np.square(N)):
		counter=0












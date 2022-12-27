#### All that is left is how to retrieve the vector data and send the processed vector data...use getv(), putv()

'''
PORT1 : mzi_1_v1 
PORT2 : mzi_1_v2
PORT3 : mzi_2_v1
PORT4 : mzi_2_v2
PORT5 : 1 at a freq of homodyne_gating_freq
PORT6 : homodyne input at a freq homodyne_gating_freq

'''



'''
import sys, os, random, pdb
import numpy

## Uncomment the following if you are using Linux
sys.path.append("/opt/lumerical/interconnect/api/python") # linux


import lumapi

h = lumapi.open("interconnect")



time_window = 5e-9
MZI_freq=10GHz
homodyne_gating_freq=MZI_freq/N

alpha=0.5
beta=0.5


#THE  J AND N VALUE WILL BE RETRIEVED FROM THE J_FILE, SO THE ASSIGNED VALUE DOESNT MATTER
#######################################
J_file=
J=
N=100 #number of spins

J=beta*(-J)
np.fill_diagonal(J,alpha)
J_flat=J.flatten()
J_max=max(alpha,beta)
J_theta=np.arccos(J_flat/(J_max+1e-15)) ##dividing to normalize J for MZI

########################################

N_iter=50 # number of iterations
N_runs=1 # for averaging the results for a given hyperparameters

alpha_min=0
alpha_max=3
alpha_step=0.1

beta_min=0
beta_max=3
beta_step=0.01

std=0.05
MZI_V_pi=1


#for tunable element H, ....NOT NEEDED SINCE WE HAVE ALREADY GIVEN PHASES
#v1=2*V_pi+(V_pi/np.pi)*(-np.pi/4)
#v2=2*V_pi+(V_pi/np.pi)*(-3*np.pi/4)

# we must load the file after starting INTERCONNECT
lumapi.evalScript(h,''' load("ising.icp"); ''');


#remember that the dot product will have a propportional constant, we have to figure that out
x_in=np.zeros(N) # takes input from homodyne
x_in=x_in*J_max/2 # multiply by J_max since we divided it for J_theta, and divide by 2 because the the homodyne product is having factor 2
x_out=x_in
x_out_flat=np.tile(x_out,N)



theta1=x_out_flat-(np.pi/2)+np.random.normal(scale=std,size=(len(self.x_out_flat)))
theta2=self.J_theta


mzi_1_v1=2*MZI_V_pi+(MZI_V_pi*theta1/np.pi)
mzi_1_v2=2*MZI_V_pi-(MZI_V_pi*theta1/np.pi)

mzi_2_v1=2*MZI_V_pi+(MZI_V_pi*theta2/np.pi)
mzi_2_v2=2*MZI_V_pi-(MZI_V_pi*theta2/np.pi)



'''









































import sys, os, random, pdb
import numpy

## Uncomment the following if you are using Linux
sys.path.append("/opt/lumerical/interconnect/api/python") # linux


import lumapi

h = lumapi.open("interconnect")

#####################################################
# Define user parameters
#
# the target speed: high, low or variable
target_speed = "variable" # "high", "low", "variable"
time_window = 5e-9
Nsamples = 2000
# set the attenuation
attenuation = 10 #dB
# monitor signal at which to flip from low speed operation to high speed
conversion_threshold = 0.0005;
seed = 123456 # seed for generating a pseudo-random bit sequence
# end user parameters
#####################################################

# we must load the file after starting INTERCONNECT
lumapi.evalScript(h,''' load("Ising.icp"); ''');

# intialize simply by pushing data to Lumerical script workspace and using same script as Lumerical script example
lumapi.putDouble(h,"time_window",time_window);
lumapi.putDouble(h,"Nsamples",Nsamples);
lumapi.putDouble(h,"attenuation",attenuation);
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
    setnamed("ATT_1","attenuation",attenuation);

    # intialize the simulation
    runinitialize;

    # initialize struct used for passing co-sim data
    v = struct;


    # End initialization
    #####################################################
    
    
    ''')

# pull the time vector and other allocated vectors back to Python
time = lumapi.getVar(h,"t")
signal = 0.*time
monitor = 0.*time
output = 0.*time
speed = 0.*time

# initialize value of monitor minimum value over bit length
monitor_min = 0;

# define the low and high speed bit rates at 10 GB/s and 25 GB/s
bit_rate = [10e9,25e9]

# reset the random number seed
random.seed(seed)

# intial values of t0 and bit_length to ensure we recalculate correct values in first step
t0 = 0;
bit_length = 0;


counter = 0;
for t in time:
    if t-t0 >= bit_length: # bit has just ended, move on to next bit
        bit = random.random() > 0.5 # set bit
        if target_speed == "high":
            bit_rate_current = bit_rate[1]
        elif target_speed == "low":
            bit_rate_current = bit_rate[0]
        elif target_speed == "variable":
            if monitor_min > conversion_threshold:
                bit_rate_current = bit_rate[1]
            else:
                bit_rate_current = bit_rate[0]
        else:
            print("Error, target_speed must be 'high', 'low' or 'variable'")
            pdb.set_trace()
            lumapi.close(h)
            exit

        bit_length = 1/bit_rate_current;
        t0 = t0 + bit_length;
        bit_rate_current*1e-9;
        monitor_min = 0; # reset monitor min for next bit
    
    speed[counter] = bit_rate_current

    # set signal
    signal[counter] = 0.5*(1+bit);

    lumapi.putDouble(h,"tval",t)
    lumapi.putDouble(h,"s",signal[counter])

    lumapi.evalScript(h,'''
    # push signal into INTERCONNECT simulation
    v.index = 1;
    v.time =  tval;
    v.value = s;
    setvalue("COSIM_1","port",v);
    
    # run INTERCONNECT until maximum time pushed
    runstep;
    
    # read signals from INTERCONNECT simulation at current time
    v.index = 2;
    v2 = getvalue("COSIM_1","port",v);
    v.index = 3;
    v3 = getvalue("COSIM_1","port",v);
    ''');
    monitor[counter] = lumapi.getVar(h,"v2")
    output[counter] = lumapi.getVar(h,"v3")
    
    # update monitor min
    if monitor[counter] > monitor_min:
        monitor_min = monitor[counter]
    
    counter = counter + 1

# run finalize of INTERCONNECT
lumapi.evalScript(h,"runfinalize;")

# push data back to INTERCONNECT for plotting
lumapi.putMatrix(h,"signal",signal)
lumapi.putMatrix(h,"monitor",monitor)
lumapi.putMatrix(h,"speed",speed)
lumapi.putMatrix(h,"output",output)
lumapi.putString(h,"target_speed",target_speed)
lumapi.evalScript(h,'''
# plot recorded data
plot(t*1e9,signal,"time (ns)","signal");
legend("speed = " + target_speed);
setplot("y1 min",-0.1);
setplot("y1 max",1.1);
plot(t*1e9,monitor,"time (ns)","monitor signal");
legend("speed = " + target_speed);
plot(t*1e9,output,"time (ns)","PD output");
legend("speed = " + target_speed);
plot(t*1e9,speed*1e-9,"time (ns)","speed (GB/s)");
legend("speed");

BER_high = getresult("EYE_1","measurement/BER");
BER_low = getresult("EYE_2","measurement/BER");
?"Operating at speed: " + target_speed;
?"Estimated high speed BER: " + num2str(BER_high);
?"Estimated low speed BER: " + num2str(BER_low);
''')

pdb.set_trace()

lumapi.close(h)

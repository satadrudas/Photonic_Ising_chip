import sys, os, random, pdb
import numpy as np

## Uncomment the following if you are using Linux
sys.path.append("/opt/lumerical/v232/python/bin/python3") # linux
sys.path.append("/opt/lumerical/v232/api/python/lumapi.py") 

sys.path.append("/opt/lumerical/v232/python/bin") # linux
sys.path.append("/opt/lumerical/v232/api/python") 

sys.path.append(os.path.dirname(__file__)) #Current directory


import lumapi

h = lumapi.open("interconnect")

input_data=np.array([0.0,0.1,0.2, -0.3,0.4, 0.5,-0.6, 0.7,0.8, -0.9, 1.0 ])
#input_data=np.square(input_data)

# here you don't have to bother with squaring as it doesn't work like the amplitude modulator, but you control the optical amplitude as 'usual'
v_pi= 4 
mzm_freq=1e9
time_per_input=1/mzm_freq
samples_per_input = 50
time_per_sample=time_per_input / samples_per_input
iteration=5
Nsamples = iteration*(len(input_data)+1) *samples_per_input #2000
time_window = Nsamples * time_per_sample

integrator_data=np.zeros(iteration)
dot_product=np.zeros(iteration)

theta=np.arccos(input_data)
mzm_input1=v_pi*theta/np.pi
mzm_input2=-v_pi*theta/np.pi

lumapi.evalScript(h,''' load("1_modulation_and_reading_mzm.icp"); ''');

lumapi.putDouble(h,"time_window",time_window);
lumapi.putDouble(h,"Nsamples",Nsamples);
lumapi.putMatrix(h,"mzm_input1", mzm_input1);
lumapi.putMatrix(h,"mzm_input2", mzm_input2);
lumapi.putDouble(h,"v_pi",v_pi);


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
    
    # the v_pi value
    setnamed("MZM_1", "pi dc voltage", v_pi);

    # set number of threads to 2. More can run faster but testing is required and is machine dependent
    setnamed("::Root Element","number of threads",2);


    # intialize the simulation
    runinitialize;


    # initialize struct used for passing co-sim data
    v = struct;

    # End initialization

    #####################################################
  
    ''')

# pull the time vector and other allocated vectors back to Python
time = lumapi.getVar(h,"t")
mzm_v1 = 0.*time
mzm_v2 = 0.*time
output = 0.*time

counter=0;
input_index=0;

for t in time:

    mzm_v1[counter]=mzm_input1[input_index]
    mzm_v2[counter]=mzm_input2[input_index]
        
    lumapi.putDouble(h,"tval",t)
    lumapi.putDouble(h,"mzm_v1",mzm_v1[counter])
    lumapi.putDouble(h,"mzm_v2",mzm_v2[counter])
    
    lumapi.evalScript(h,'''
    # push signal into INTERCONNECT simulation
    v.index = 1;
    v.time =  tval;
    v.value = mzm_v1;
    setvalue("COSIM_1","port",v);

    v.index = 2;
    v.time =  tval;
    v.value = mzm_v2;
    setvalue("COSIM_1","port",v);    

    # run INTERCONNECT until maximum time pushed
    runstep;
    

    # read signals from INTERCONNECT simulation at current time
    v.index = 3;

    v2 = getvalue("COSIM_1","port",v);
    
    ''')


    output[counter] = lumapi.getVar(h,"v2")
        
    counter=counter+1
    
    if counter%samples_per_input==0:
        input_index = input_index +1
        
        if input_index == len(input_data):
            input_index=0
lumapi.evalScript(h,"runfinalize;")

lumapi.putMatrix(h,"mzm_v1",mzm_v1)
lumapi.putMatrix(h,"mzm_v2",mzm_v2)
lumapi.putMatrix(h,"output",output)

lumapi.evalScript(h,'''
# plot recorded data
plot(t*1e9,mzm_v1,"time (ns)","mzm_v1");
legend("mzm_v1");

plot(t*1e9,mzm_v2,"time (ns)","mzm_v2");
legend("mzm_v2");

plot(t*1e9,output,"time (ns)","output");
legend("output");


''')

pdb.set_trace()

lumapi.close(h)

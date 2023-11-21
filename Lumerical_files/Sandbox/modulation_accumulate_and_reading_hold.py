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


vecotr_size = 7
#input_data=np.random.uniform(0,1,vecotr_size)#np.array([0.1,0.2, 0.3, 0.5, 0.7,0.8, 0.9 ])
input_data=np.array([0.1,0.2, 0.3,0.4, 0.5, 0.6,0.7,0.8, 0.9, 1.0 ])

# the AM seems to modulate the power in propotion to the signal strength and not the "optical amplitude", to do that, take square root of the power
input_data = np.square(input_data)

time_per_data=0.2e-9
samples_per_data = 10
iteration=5
Nsamples = iteration*(len(input_data)+1) *samples_per_data #2000
time_window = Nsamples * time_per_data

integrator_data=np.zeros(iteration)


lumapi.evalScript(h,''' load("2_modulation_accumulate_and_reading.icp"); ''');

lumapi.putDouble(h,"time_window",time_window);
lumapi.putDouble(h,"Nsamples",Nsamples);
lumapi.putMatrix(h,"input_data", input_data);

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
pd_output = 0.*time
integrator_output = 0.*time

counter=0;
counter2=0;
input_index=0;# 0 to len(input_data)+1, the +1 is for the reset step
reset_flag=0;
integrator_index=0;


for t in time:

    lumapi.putDouble(h,"tval",t)

    # for resetting the integrator
    if reset_flag:
        signal[counter]=0
        lumapi.putDouble(h,"s",signal[counter]);
    
        lumapi.evalScript(h,'''
        # push signal into INTERCONNECT simulation
        v.index = 3;
        v.time =  tval;
        v.value = 1;   
        setvalue("COSIM_1","port",v);
        
        
        # push signal into INTERCONNECT simulation
        v.index = 1;
        v.time =  tval;
        v.value = s;
        setvalue("COSIM_1","port",v);
        
    
        # run INTERCONNECT until maximum time pushed
        runstep;
        
    
        # read signals from INTERCONNECT simulation at current time
        #v.index = 2;    
        #v2 = getvalue("COSIM_1","port",v);
        
        ''')
        

        
    else:
    
            
        signal[counter]=input_data[input_index];
            
        lumapi.putDouble(h,"s",signal[counter]);
        
        lumapi.evalScript(h,'''
        # push signal into INTERCONNECT simulation
        v.index = 1;
        v.time =  tval;
        v.value = s;
        setvalue("COSIM_1","port",v);
        
        v.index = 3;
        v.time =  tval;
        v.value = 0;   
        setvalue("COSIM_1","port",v);
        
    
        # run INTERCONNECT until maximum time pushed
        runstep;
        
    
        # read signals from INTERCONNECT simulation at current time
        v.index = 2;    
        v_pd = getvalue("COSIM_1","port",v);

        v.index = 4;    
        v_integrator = getvalue("COSIM_1","port",v);        
        ''')
    
    
        pd_output[counter] = lumapi.getVar(h,"v_pd")
        integrator_output[counter] = lumapi.getVar(h,"v_integrator")
            
    counter=counter+1
    
    if counter%samples_per_data==0:
     
        input_index = input_index +1
    
        if input_index == len(input_data):
            integrator_data[integrator_index]=integrator_output[counter-1]
            integrator_index=integrator_index+1
            if integrator_index==len(integrator_data):
                integrator_index=0     
                       
            reset_flag=1
   
            #input_index=0
        if input_index == len(input_data)+1:
            input_index=0
            reset_flag=0
            # update for the next vector
            #input_data=np.random.uniform(0,1,vecotr_size) # if you want every iteration to have a different vector to accumulate then uncomment this line
            
        
                
print(integrator_data)
lumapi.evalScript(h,"runfinalize;")

lumapi.putMatrix(h,'accumulated_data',integrator_data )

lumapi.putMatrix(h,"signal",signal)
lumapi.putMatrix(h,"pd_output",pd_output)
lumapi.putMatrix(h,"integrator_output",integrator_output)

lumapi.evalScript(h,'''
# plot recorded data
plot(t*1e9,signal,"time (ns)","signal");
legend("signal");

plot(t*1e9,pd_output,"time (ns)","pd_output");
legend("pd_output");

plot(t*1e9,integrator_output,"time (ns)","integrator_output");
legend("integrator_output");


''')

pdb.set_trace()

lumapi.close(h)

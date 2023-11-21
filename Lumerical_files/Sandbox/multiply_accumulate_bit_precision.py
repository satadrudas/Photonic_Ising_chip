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


vector_size = 50
#input_data1=np.random.uniform(-1,1,vector_size)#np.array([0.1,0.2, 0.3, 0.5, 0.7,0.8, 0.9 ])
#input_data2=np.ones(vector_size)

#input_data1=np.array([0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, 0.0])
#input_data2=np.array([ 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])

#input_data1=np.array([ 1.0,1.0, -1.0, -1.0])
#input_data2=np.array([ 1.0,-1.0, 1.0, -1.0])

input_data1=np.random.uniform(-10000,10000,vector_size)#np.array([0.1,0.2, 0.3, 0.5, 0.7,0.8, 0.9 ])
input_data2=np.random.uniform(-10000,10000,vector_size)#np.array([0.1,0.2, 0.3, 0.5, 0.7,0.8, 0.9 ])
#input_data2=np.ones(vector_size)


#input_data1=np.insert(input_data1, 0, 0.0)
input_data1=np.insert(input_data1, len(input_data1), 0.0)

#input_data2=np.insert(input_data2, 0, -1.0 if input_data2[0]<0 else 1.0)
input_data2=np.insert(input_data2, len(input_data2),  1.0)

result=np.dot(input_data1, input_data2)

###!!!!!!!!!careful here, input_data is CHANGING VALUES AFTER FACTORING

max_factor = np.max([np.max(np.abs(input_data1)), np.max(np.abs(input_data2))])
input_data1=input_data1/max_factor
input_data2=input_data2/max_factor


# the AM seems to modulate the power in propotion to the signal strength and not the "optical amplitude", to do that, take square root of the power
#input_data1 = np.square(input_data1)

v_pi= 4 
mzm_freq=1e9
time_per_input=1/mzm_freq
samples_per_input = 3
time_per_sample=time_per_input / samples_per_input
iteration=1
Nsamples = iteration*(len(input_data1)+1) *samples_per_input #2000
time_window = Nsamples * time_per_sample

integrator_data=np.zeros(iteration)
dot_product=np.zeros(iteration)

theta1=np.arccos(input_data1)
mzm1_input1=v_pi*theta1/np.pi
mzm1_input2=-v_pi*theta1/np.pi

theta2=np.arccos(input_data2)
mzm2_input1=v_pi*theta2/np.pi
mzm2_input2=-v_pi*theta2/np.pi

lumapi.evalScript(h,''' load("3_multiply_accumulate.icp"); ''');

lumapi.putDouble(h,"time_window",time_window);
lumapi.putDouble(h,"Nsamples",Nsamples);
lumapi.putMatrix(h,"mzm1_input1", mzm1_input1);
lumapi.putMatrix(h,"mzm1_input2", mzm1_input2);
lumapi.putMatrix(h,"mzm2_input1", mzm2_input1);
lumapi.putMatrix(h,"mzm2_input2", mzm2_input2);
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
    setnamed("MZM_1", "bias voltage 1", 0);
    setnamed("MZM_1", "bias voltage 2", 0);
    setnamed("MZM_1", "extinction ratio", 1e+13);
    setnamed("MZM_1", "insertion loss", 0);
    
    setnamed("MZM_2", "pi dc voltage", v_pi);
    setnamed("MZM_2", "bias voltage 1", 0);
    setnamed("MZM_2", "bias voltage 2", 0);
    setnamed("MZM_2", "extinction ratio", 1e+13);
    setnamed("MZM_2", "insertion loss", 0);

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
mzm1_v1 = 0.*time
mzm1_v2 = 0.*time
mzm2_v1 = 0.*time
mzm2_v2 = 0.*time
pd_output = 0.*time
integrator_output = 0.*time
dot_product_data = 0.*time

counter=0;
input_index=0;# 0 to len(input_data1)+1, the +1 is for the reset step
reset_flag=0;
integrator_index=0;


for t in time:

    lumapi.putDouble(h,"tval",t)

    # for resetting the integrator
    if reset_flag:
        mzm1_v1[counter]=np.arccos(0.0)*v_pi/np.pi
        mzm1_v2[counter]=-np.arccos(0.0)*v_pi/np.pi        
        mzm2_v1[counter]=np.arccos(0.0)*v_pi/np.pi
        mzm2_v2[counter]=-np.arccos(0.0)*v_pi/np.pi

        
    else:
    
            
        mzm1_v1[counter]=mzm1_input1[input_index]
        mzm1_v2[counter]=mzm1_input2[input_index]
        mzm2_v1[counter]=mzm2_input1[input_index]
        mzm2_v2[counter]=mzm2_input2[input_index]
                
    lumapi.putDouble(h,"mzm1_v1",mzm1_v1[counter])
    lumapi.putDouble(h,"mzm1_v2",mzm1_v2[counter])
    lumapi.putDouble(h,"mzm2_v1",mzm2_v1[counter])
    lumapi.putDouble(h,"mzm2_v2",mzm2_v2[counter])
    lumapi.putDouble(h,"reset_flag", reset_flag)

    lumapi.evalScript(h,'''
    # push signal into INTERCONNECT simulation
    v.index = 1;
    v.time =  tval;
    v.value = mzm1_v1;
    setvalue("COSIM_1","port",v);

    v.index = 2;
    v.time =  tval;
    v.value = mzm1_v2;
    setvalue("COSIM_1","port",v);  
    
    v.index = 3;
    v.time =  tval;
    v.value = mzm2_v1;
    setvalue("COSIM_1","port",v);

    v.index = 4;
    v.time =  tval;
    v.value = mzm2_v2;
    setvalue("COSIM_1","port",v);  
    
    v.index = 6;
    v.time =  tval;
    v.value = reset_flag;   
    setvalue("COSIM_1","port",v);
    

    # run INTERCONNECT until maximum time pushed
    runstep;
    

    # read signals from INTERCONNECT simulation at current time
    v.index = 5;    
    v_pd = getvalue("COSIM_1","port",v);

    v.index = 7;    
    v_integrator = getvalue("COSIM_1","port",v);        
    ''')


    pd_output[counter] = lumapi.getVar(h,"v_pd")
    integrator_output[counter] = lumapi.getVar(h,"v_integrator")*(max_factor**2)
    dot_product_data[counter]=integrator_output[counter]/(time_per_input * 0.001)# the 0.001 is to factor out the milliwatt to 1

        
    counter=counter+1
    
    if counter%samples_per_input==0:
     
        input_index = input_index +1
    
        if input_index == len(input_data1):
            reset_flag=1

            integrator_data[integrator_index]=integrator_output[counter-1]
            dot_product[integrator_index]=dot_product_data[counter-1]
            integrator_index=integrator_index+1
            
            if integrator_index==len(integrator_data): # for the next vvm
                integrator_index=0     
                       
   
            #input_index=0
        if input_index == len(input_data1)+1:
            input_index=0
            reset_flag=0
            # update for the next vector
            #input_data1=np.random.uniform(0,1,vecotr_size) # if you want every iteration to have a different vector to accumulate then uncomment this line

###############################################################
      
print("Dot product (actual): "+str(result)+"\n") # comment it when not using      
print("Integrated values (Simulation): "+str(integrator_data))
print("Dot product values (Simulation): "+str(dot_product)+"\n")

tol=1e-3
if np.absolute(result-dot_product)<tol : print("Correct value") 
else: print("Wrong value. NOTE: Sometimes the result says wrong value due to some floating point errors.")

print("\n\n")
###############################################################3

lumapi.evalScript(h,"runfinalize;")

lumapi.putMatrix(h,'accumulated_data',integrator_data )

lumapi.putMatrix(h,"mzm1_v1",mzm1_v1)
lumapi.putMatrix(h,"mzm1_v2",mzm1_v2)
lumapi.putMatrix(h,"mzm2_v1",mzm2_v1)
lumapi.putMatrix(h,"mzm2_v2",mzm2_v2)
lumapi.putMatrix(h,"pd_output",pd_output)
lumapi.putMatrix(h,"integrator_output",integrator_output)
lumapi.putMatrix(h,"dot_product_data",dot_product_data)

lumapi.evalScript(h,'''
# plot recorded data
plot(t*1e9,mzm1_v1,"time (ns)","mzm1_v1");
legend("mzm1_v1");

plot(t*1e9,mzm1_v2,"time (ns)","mzm1_v2");
legend("mzm1_v2");

plot(t*1e9,mzm2_v1,"time (ns)","mzm2_v1");
legend("mzm2_v1");

plot(t*1e9,mzm2_v2,"time (ns)","mzm2_v2");
legend("mzm2_v2");

plot(t*1e9,pd_output,"time (ns)","pd_output");
legend("pd_output");

plot(t*1e9,integrator_output,"time (ns)","integrator_output");
legend("integrator_output");

plot(t*1e9,dot_product_data,"time (ns)","dot_product_data");
legend("dot_product_data");

''')

pdb.set_trace()

lumapi.close(h)

import sys, os, random, pdb
import numpy as np
import matplotlib.pyplot as plt

## Uncomment the following if you are using Linux
sys.path.append("/opt/lumerical/v232/python/bin/python3") # linux
sys.path.append("/opt/lumerical/v232/api/python/lumapi.py") 

sys.path.append("/opt/lumerical/v232/python/bin") # linux
sys.path.append("/opt/lumerical/v232/api/python") 

sys.path.append(os.path.dirname(__file__)) #Current directory


import lumapi

h = lumapi.open("interconnect")


vector_size = 500


v_pi= 4 
mzm_freq=10e9
time_per_input=1/mzm_freq
samples_per_input = 10
time_per_sample=time_per_input / samples_per_input
iteration=1
Nsamples = iteration*vector_size *samples_per_input #2000
time_window = Nsamples * time_per_sample
norm_limiter=0.7
dac_bit_precision=8
dac_v_max=v_pi/2#/np.pi * np.arcsin(norm_limiter)

# normalization is not really needed here, but norm limiter is needed. or instead of norm limiter jor generate from random wit limit
def normalize(input, limiter=norm_limiter):
    normalizing_factor = np.max(np.abs(input))/limiter
    normalized_input = input/normalizing_factor
    return normalized_input, normalizing_factor


def mzm_voltages(input, theta_flag=0, quantization=1, dac_precision=dac_bit_precision):
    '''
    Returns the neccessary voltages of a
    dual drive MZM for the given input
    
    If theta_flag is set, then the theta will be set 
    to the input and no arcsin will be taken
    
    quantization flag is set for dac precision
    '''
    
    if theta_flag:
        theta=input#*1.1
    else:
        theta=np.arcsin(input)
        
    mzm_v1=v_pi*theta/np.pi
    mzm_v2=-v_pi*theta/np.pi   
    
    if quantization:
            
        allowed_dac_volts = np.linspace(-dac_v_max,dac_v_max, num=2**dac_precision-1, endpoint=True)
        bins = allowed_dac_volts + np.abs(allowed_dac_volts[0]-allowed_dac_volts[1])/2
        
        mzm_v1_bin_index = np.digitize(mzm_v1,bins)
        mzm_v2_bin_index = np.digitize(mzm_v2,bins)
        
        for i,j in zip(range(0,len(mzm_v1)),mzm_v1_bin_index):
            mzm_v1[i]=allowed_dac_volts[j]
            
        for i,j in zip(range(0,len(mzm_v2)),mzm_v2_bin_index):
            mzm_v2[i]=allowed_dac_volts[j]
         
    return mzm_v1, mzm_v2  
    
    
###################################################################
input_data1=np.random.uniform(-1,1,vector_size)
input_data2=np.random.uniform(-1,1,vector_size)


result=np.multiply(input_data1, input_data2)

#normalizing the inputs
input_data1, input_data1_normalizing_factor = normalize(input_data1)
input_data2, input_data2_normalizing_factor = normalize(input_data2)


optical_result = np.zeros(len(input_data1))

mzm1_input1, mzm1_input2 = mzm_voltages(input_data1)
mzm2_input1, mzm2_input2 = mzm_voltages(input_data2)


lumapi.evalScript(h,''' load("multiplication_error_test.icp"); ''');

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
    setnamed("MZM_1", "bias voltage 1", -v_pi/2);
    setnamed("MZM_1", "bias voltage 2", v_pi/2);
    setnamed("MZM_1", "extinction ratio", 1e+13);
    setnamed("MZM_1", "insertion loss", 0);
    
    setnamed("MZM_2", "pi dc voltage", v_pi);
    setnamed("MZM_2", "bias voltage 1", -v_pi/2);
    setnamed("MZM_2", "bias voltage 2", v_pi/2);
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


counter=0;
input_index=0;# 0 to len(input_data1)+1, the +1 is for the reset step


for t in time:

    lumapi.putDouble(h,"tval",t)
            
    mzm1_v1[counter]=mzm1_input1[input_index]
    mzm1_v2[counter]=mzm1_input2[input_index]
    mzm2_v1[counter]=mzm2_input1[input_index]
    mzm2_v2[counter]=mzm2_input2[input_index]
                
    lumapi.putDouble(h,"mzm1_v1",mzm1_v1[counter])
    lumapi.putDouble(h,"mzm1_v2",mzm1_v2[counter])
    lumapi.putDouble(h,"mzm2_v1",mzm2_v1[counter])
    lumapi.putDouble(h,"mzm2_v2",mzm2_v2[counter])

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
    

    # run INTERCONNECT until maximum time pushed
    runstep;
    

    # read signals from INTERCONNECT simulation at current time
    v.index = 5;    
    v_pd = getvalue("COSIM_1","port",v);

     
    ''')


    pd_output[counter] = lumapi.getVar(h,"v_pd")
    
    counter=counter+1
    
    if counter%samples_per_input==0:
    
        optical_result[input_index]=pd_output[counter-1]*input_data1_normalizing_factor*input_data2_normalizing_factor/0.001
     
        input_index = input_index +1
        
        if input_index == vector_size:
            break
        
    
error = result - optical_result
###############################################################
      
print("Correct result: "+str(result)+"\n") # comment it when not using      

print("Optical result: "+str(optical_result)+"\n") # comment it when not using      
print("Error: "+str(error)+"\n") # comment it when not using      
print("Error mean:"+str(np.mean(error))+'\n')
print("Error standard deviation:"+str(np.std(error))+'\n')

print("\n\n")
###############################################################3

lumapi.evalScript(h,"runfinalize;")


lumapi.putMatrix(h,"mzm1_v1",mzm1_v1)
lumapi.putMatrix(h,"mzm1_v2",mzm1_v2)
lumapi.putMatrix(h,"mzm2_v1",mzm2_v1)
lumapi.putMatrix(h,"mzm2_v2",mzm2_v2)
lumapi.putMatrix(h,"pd_output",pd_output)
lumapi.putMatrix(h,"result",result)
lumapi.putMatrix(h,"optical_result",optical_result)
lumapi.putMatrix(h,"error",error)

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

plot(result,optical_result,"y", "y_hat", "floating point result vs optical result", "plot type=point, marker style=x");
legend("optical result");

#histc(error, 12);

''')


plt.hist(error, bins='auto')
plt.show()

plt.plot(result,optical_result)
plt.show() 


pdb.set_trace()

lumapi.close(h)

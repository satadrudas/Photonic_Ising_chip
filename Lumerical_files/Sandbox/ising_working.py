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
###################################################################################
def normalize(input):
    normalizing_factor = np.max(np.abs(input))
    normalized_input = input/normalizing_factor
    return normalized_input, normalizing_factor
    
def vvm_preprocess(input, insert_value=0.0):
    '''
    inserts a 0.0 at the end of the vector 
    to get proper integration
    '''
    input=np.insert(input, len(input), insert_value)    
    return input
    
def mzm_voltages(input, theta_flag=0):
    '''
    Returns the neccessary voltages of a
    dual drive MZM for the given input
    
    If theta_flag is set, then the theta will be set 
    to the input and no arccos will be taken
    '''
    if theta_flag:
        theta=input
    else:
        theta=np.arccos(input)
        
    mzm_v1=v_pi*theta/np.pi
    mzm_v2=-v_pi*theta/np.pi    
    return mzm_v1, mzm_v2  


def count_positive_values(arr):
    count = 0
    for value in arr:
        if value > 0:
            count += 1
    return count

def J_hyperparameters(J_mat,alpha, beta):
    J_mat = beta*(-J_mat)#changing to antiferromagnetic
    np.fill_diagonal(J_mat,alpha)
    return J_mat


    
###################################################################################

N_spins = 30
alpha = 1.5
beta=0.0

v_pi= 4
mzm_freq=1e9
time_per_input=1/mzm_freq
samples_per_input = 3
time_per_sample=time_per_input / samples_per_input
N_iterations=30
N_samples = N_iterations*(N_spins+2)*N_spins *(samples_per_input) #+2 because 1 for the integrating part and 1 for the extra 0 we add during reprocessing
time_window = N_samples * time_per_sample
sig = 0.1

integrator_data=np.zeros(N_spins)
dot_product=np.zeros(N_spins)
spin_evolution = np.zeros(shape=(N_iterations+1,N_spins))
    


lumapi.evalScript(h,''' load("3_multiply_accumulate.icp"); ''');

lumapi.putDouble(h,"time_window",time_window);
lumapi.putDouble(h,"N_samples",N_samples);
lumapi.putDouble(h,"v_pi",v_pi);

lumapi.evalScript(h,'''
    # turn off history and redrawing
    redrawoff;
    historyoff;

    # time vectors
    t = linspace(0,time_window,N_samples+1);

    t = t(1:end-1);
    dt = t(2)-t(1);

    # restore design mode in case we are in analysis mode
    switchtodesign;

    # Set up simulation to match time window and number of samples
    # The time windows should match because the INTERCONNECT simulation
    # will stop once the simulation time exceeds the time window.
    # The number of samples defines the INTERCONNECT time step, dt, by dt = time_window/(N_samples+1).
    # The time steps do NOT have to match, although in this example they do. Indeed,
    # the time step of an external simulator can be variable

    setnamed("::Root Element","time window",time_window);
    setnamed("::Root Element","number of samples",N_samples);
    
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
iteration_counter = 0


J = np.diag(np.ones(N_spins),0)#np.random.uniform(-1,1,(N_spins, N_spins))
J_matrix = J_hyperparameters(J,alpha, beta)

spins = np.zeros(N_spins)

sig = 0.05

spins = spins-np.pi/2 + np.random.normal(0,sig,N_spins) 

result=np.dot(J_matrix, np.cos(spins))

J_matrix_normalized, J_matrix_normalization_factor = normalize(J_matrix)

input_data1_normalized=J_matrix_normalized[integrator_index]
#spins_normalized, spins_normalization_factor = normalize(spins)

input_data1 = vvm_preprocess(input_data1_normalized, 0.0)
spins = vvm_preprocess(spins, 0.0) # actually it doenst matter what vlue you put for the preprocessor cuz, the J already has a 0.0, so the product is anyway going to be 0
   
mzm1_input1,mzm1_input2 = mzm_voltages(input_data1, 0)
mzm2_input1,mzm2_input2 = mzm_voltages(spins, 1)



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
    integrator_output[counter] = lumapi.getVar(h,"v_integrator")*J_matrix_normalization_factor#*spins_normalization_factor
    dot_product_data[counter]=integrator_output[counter]/(time_per_input * 0.001)# the 0.001 is to factor out the milliwatt to 1

        
    counter=counter+1
    
    if counter%samples_per_input==0:
     
        input_index = input_index +1
    
        if input_index == len(input_data1):
            reset_flag=1

            integrator_data[integrator_index] = integrator_output[counter-1]
            dot_product[integrator_index] = dot_product_data[counter-1]
            
            integrator_index = integrator_index+1
          
            if integrator_index==len(integrator_data): # for the next mvm
                
                integrator_index=0   
                iteration_counter = iteration_counter + 1  
                
                print(np.array(dot_product))
                spin_evolution[iteration_counter]=np.array(dot_product)
                
                spins = dot_product -np.pi/2 + np.random.normal(0,sig,N_spins) 

                spins = vvm_preprocess(spins)# actually it doenst matter what vlue you put for the preprocessor cuz, the J already has a 0.0, so the product is anyway going to be 0
                mzm2_input1,mzm2_input2 = mzm_voltages(spins, 1)

    
                       
   
            #input_index=0
        if input_index == len(input_data1)+1:
            input_index=0
            reset_flag=0

            #update for the next VVM
            #spin vector stays the same, only J[i] -> J[i+1]
            
            input_data1_normalized=J_matrix_normalized[integrator_index]            
            input_data1 = vvm_preprocess(input_data1_normalized)                
            mzm1_input1,mzm1_input2 = mzm_voltages(input_data1)
            
            
###############################################################
print("Dot product (actual): "+str(result)+"\n") # comment it when not using      
print("Integrated values (Simulation): "+str(integrator_data)+"\n")
print("Dot product values (Simulation): "+str(dot_product)+"\n")
#print("Error in calculation :"+str(np.absolute(result-dot_product))+"\n")


p_num = count_positive_values(dot_product)
n_num = len(dot_product)-p_num

print("Spin positive :"+str(p_num/len(dot_product))+" Spin negative :"+str(n_num/len(dot_product) )+"\n")


#tol=[1e-3]
#if np.absolute(result-dot_product).all()<tol : 
#    print("Correct value") 
#else: 
#    print("Wrong value. NOTE: Sometimes the result says wrong value due to some floating point errors.")

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
lumapi.putMatrix(h,"dot_product",dot_product)##

test=np.random.uniform(-1,1,(N_spins, 17))
lumapi.putMatrix(h,"test",test)##

spin_evolution=np.array(spin_evolution)/alpha
print(spin_evolution)#,spin_evolution.shape )
lumapi.putMatrix(h,"spin_evolution",spin_evolution)##
lumapi.putDouble(h,"N_spins",N_spins)
lumapi.putDouble(h,"N_iterations",N_iterations)


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

x=linspace(1,length(dot_product),length(dot_product));
plot(x,dot_product,"Spins", "Spin value", "Final Spin value", "plot type=bar");
legend("Final spin values");

#x=linspace(1,5,5);
#plot(x,test,"Spins", "Spin value", "Final Spin value", "plot type=line");
#legend("Test");

x=linspace(0,size(spin_evolution,1)-1,size(spin_evolution,1));
plot(x,spin_evolution,"Iteration", "Spin value", " Spin evolution (uncloupled)", "plot type=line");
#legend("Test");

''')

pdb.set_trace()

lumapi.close(h)

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

N_spins = 16
alpha = 0.4
beta=0.5

c_integrator=100e-12 # Farad
v_pi= 4
mzm_freq=1e9
time_per_input=1/mzm_freq
samples_per_input = 3
time_per_sample=time_per_input / samples_per_input
N_iterations=20
N_samples = N_iterations*(N_spins+2)*N_spins *(samples_per_input) #+2 because 1 for the integrating part and 1 for the extra 0 we add during reprocessing
time_window = N_samples * time_per_sample
sig = 0.1

integrator_data = np.zeros(N_spins)
dot_product = np.zeros(N_spins)
hamiltonian_evolution = np.zeros(N_iterations+1)

spin_evolution = np.zeros(shape=(N_iterations+1,N_spins))
noise = np.random.normal(0,sig,(N_iterations+1,N_spins)) 

###################################################################################
def square_lattice_coupling_matrix_generator(n):
    j=np.zeros((n,n))
    sq_len = int(np.sqrt(n))

    for i in range(0,sq_len):
        for k in range(i*sq_len, (i+1)*sq_len-1):
            j[k, k+1]=1
            j[k+1,k]=1 

    for i in range(0,sq_len-1):
        for k in range(i*sq_len, (i+1)*sq_len):
            j[k, k+sq_len]=1
            j[k+sq_len,k]=1

    return j

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

def cut_value(final_spins, J_mat):
    
    value = 0
    final_spins = np.sign(final_spins)
    for i in range(0,N_spins):
        for j in range(0, N_spins):
            if final_spins[i]*final_spins[j]<0:
                value = value + J_mat[i,j]
                    
    return value/2
    

def hamiltonian(spin_value,J_mat):
    H = 0
    spin_value = np.sign(spin_value)
    J_mat=-J_mat ####!!!!!!!!!!! REMEMBER ANTIFERRO, REMOVE THIS WHEN YOU HAVE CORRECTED THE J BEFORE ITSELF
    for row,i in enumerate(spin_value):
        for col,j in enumerate(spin_value):
            H = H -J_mat[row][col]*i*j
    return H/2 # refer the poor man paper
    
###################################################################################


    


lumapi.evalScript(h,''' load("ising.icp"); ''');

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


#J = np.random.uniform(0,1,(N_spins, N_spins))
#J=(J+np.transpose(J))/2
#J = np.diag(np.ones(N_spins),0)

#########################################
#J_file="/home/satadrudas/Photonics/Photonic_Ising_chip/Lumerical_files/Sandbox/Maxcut_instances/s_100.txt"
#f = open(J_file,"r")
#N_spins, number_of_edges = [int(i) for i in f.readline().split()]
#J = np.zeros([N_spins,N_spins])
#lines = f.readlines()
#for line in lines:
#    l = line.split()
#    r,c,w = int(l[0])-1, int(l[1])-1, float(l[2])
#    J[r][c] = w
#    J[c][r] = w        
#f.close()
#J=np.array(J)

#########################################
J = square_lattice_coupling_matrix_generator(N_spins)
J_matrix = J_hyperparameters(J,alpha, beta)

spins = np.zeros(N_spins)


spins = spins-np.pi/2 + noise[0]

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
    integrator_output[counter] = lumapi.getVar(h,"v_integrator")/c_integrator # voltage across the capacitor, , no need quantization here
    dot_product_data[counter]=integrator_output[counter]*J_matrix_normalization_factor*c_integrator/(time_per_input * 0.001)# the 0.001 is to factor out the milliwatt to 1

        
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
                
                #print(np.array(dot_product))
                #spin_evolution[iteration_counter]=np.array(dot_product)
                
                spins = dot_product -np.pi/2 + noise[iteration_counter]
                spin_evolution[iteration_counter]=np.cos(np.array(spins))
                print(spin_evolution[iteration_counter])
                hamiltonian_evolution[iteration_counter] = hamiltonian(spin_evolution[iteration_counter], J)
                print("current hamiltonian: "+str(hamiltonian_evolution[iteration_counter])+"\n\n")

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

final_spins_value = spin_evolution[-1]

print("Dot product (actual): "+str(result)+"\n") # comment it when not using      
print("Integrated values (Simulation): "+str(integrator_data)+"\n")
print("Dot product values (Simulation): "+str(dot_product)+"\n")
print("Hamiltonian Evolution: "+str(hamiltonian_evolution)+"\n")
print("Cut value (alpha = "+str(alpha)+", beta = "+str(beta)+") : "+str(cut_value(spin_evolution[-1], J)))

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
lumapi.putMatrix(h,"final_spins_value",final_spins_value)##
lumapi.putMatrix(h,"hamiltonian_evolution",hamiltonian_evolution)##


spin_evolution=np.array(spin_evolution)
print(spin_evolution)#,spin_evolution.shape )
lumapi.putMatrix(h,"spin_evolution",spin_evolution)##



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

#x=linspace(1,length(dot_product),length(dot_product));
#plot(x,dot_product,"Spins", "Spin value", "Final Spin value", "plot type=bar");
#legend("Final spin values");

x=linspace(1,length(final_spins_value),length(final_spins_value));
plot(x,final_spins_value,"Spins", "Spin value", "Final Spin value", "plot type=bar");
legend("Final spin values");

x=linspace(0,size(spin_evolution,1)-1,size(spin_evolution,1));
plot(x,spin_evolution,"Iterations", "Spin value", " Spin evolution", "plot type=line");
legend("Spins");
#change the title to uncoupled spins when beta=0.0

x=linspace(0,length(hamiltonian_evolution)-1,length(hamiltonian_evolution));
plot(x,hamiltonian_evolution,"Iterations", "Ising Energy", "Ising Energy", "plot type=line");
legend("Ising Energy");

''')

# plotting the spin evoluting...WORKS
#x=np.linspace(0,N_iterations, N_iterations+1)
#plt.plot(x,spin_evolution)
#plt.show()

# for square lattice you can also visualize the checkered board pattern
plt.matshow(final_spins_value.reshape(int(np.sqrt(N_spins)), int(np.sqrt(N_spins))))
plt.colorbar()
plt.show()

pdb.set_trace()

lumapi.close(h)

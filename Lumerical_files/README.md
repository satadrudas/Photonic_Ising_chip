I used Lumerical INTERCONNECT to do some preliminary simulations.

The circuit is driven by a Lumerical INTERCONNECT python co-simulation module.

![alt text](Ising.png)

NOTE:
In the code, the following path seems to work:

```
## Uncomment the following if you are using Linux
sys.path.append("/opt/lumerical/v232/python/bin/python3") # linux
sys.path.append("/opt/lumerical/v232/api/python/lumapi.py") 

sys.path.append("/opt/lumerical/v232/python/bin") # linux
sys.path.append("/opt/lumerical/v232/api/python") 

sys.path.append(os.path.dirname(__file__)) #Current directory
```

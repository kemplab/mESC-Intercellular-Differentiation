import Model_initialParameters as mod
from CountDirectory import newDirect
import sys
import getpass


if(__name__ == '__main__'):
    length=1
    user=getpass.getuser()
    
    ##Set save path for simulations. If you create a new, empty folder 
    ##then everytime you run this code it will make a new Simulation folder
    path="D:\\Final FP Sims\\Tests"

    inh_val=[1.0]
    per_inh=[1.0]
    inh_time=1000
    inh_time1=1000
    inhib_type=None
    ##Current settings don't cause any inhibition.
    ##Possible settings include 'AC' which inhibits production rates and
    ## 'GJ which inhibits transport. Unless you want on/off though you can obtain
    ## this effect by changing other parameters

    
    for j in inh_val:
        
            for i in range(length):
                x=newDirect(path)
                mod.main(x,path,j,inhib_type,inh_time,inh_time1)

import numpy as np

# stress data
stress_low_t = np.array([7.79,3.03,1.63,1.04,.688])
stress_high_t = np.array([7.37,2.77,1.37,.808,.515])
stress_low_r = np.array([7.79,7.1,7.16,7.41,7.37])
stress_high_r = np.array([.688,.599,.547,.527,.515])

def avg_percent_change(arr1,arr2):
    change = arr2 - arr1
    change = change/arr1
    return np.average(change)

stress_vary_r = avg_percent_change(stress_low_r,stress_high_r)
print(stress_vary_r)
#-0.9219841613262995
stress_vary_t = avg_percent_change(stress_low_t,stress_high_t)
print(stress_vary_t)
#-0.15475269415119092

# deflection data
def_low_t = np.array([.819,.218,.084,.0396,.0212])
def_high_t = np.array([.778,.207,.0794,.0374,.0201])
def_low_r = np.array([.819,.811,.801,.79,.778])
def_high_r = np.array([.0212,.02098,.0207,.0204,.0201])

def_vary_r = avg_percent_change(def_low_r,def_high_r)
print(def_vary_r)
#-0.9741489039866066
def_vary_t = avg_percent_change(def_low_t,def_high_t)
print(def_vary_t)
#-0.05254480368553418

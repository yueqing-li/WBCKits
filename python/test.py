import numpy as np

import pyWBCKits as WBCKits

A1 = np.array([[1., 0., 0., 0., 0., 0., 0.],
               [0., 1., 0., 0., 0., 0., 0.],
               [0., 0., 1., 0., 0., 0., 0.]])
b1 = np.array([1., 1., 1.])

A2 = np.array([[0., 0., 0., 1., 0., 0., 0.],
               [0., 0., 0., 0., 1., 0., 0.],
               [0., 0., 0., 0., 0., 1., 0.]])
b2 = np.array([0., 0., 0.])
w2 = 3.*np.ones(3)

A3 = np.array([[0., 0., 0., 0., 0., 1., 0.],
               [0., 0., 0., 0., 0., 0., 1.]])
b3 = 5.*np.array([1., 1.])
w3 = 0.5*np.ones(2)

C1 = np.eye(7)
lb1 = -3.0*np.ones(7)
ub1 =  3.0*np.ones(7)

C2 = np.array([[1., 0., 1., 0., 0., 0., 0.],
               [0., 1., 0., 1., 0., 0., 0.]])
lb2 = np.array([-3., -2.])
ub2 = np.array([3., 2.])

####################
# Options TEST
####################

Options = WBCKits.PyOptions()
print("Default settings:")
Options.Print()
Options.nwsr_max_ = 200
Options.cpu_time_max_ = 1e-2 # second
Options.hqp_algorithm_ = WBCKits.HQP_Origin
Options.decompose_method_ = WBCKits.Decompose_QR
Options.decompose_threshold_ = 1e-3
Options.get_slack_variable_ = True
print("New settings:")
Options.Print()

#####################
# WQP TEST
#####################
wqp = WBCKits.WQP()
print(dir(wqp))
# set options
wqp.SetOptions(Options)
# add tasks
wqp.AddTask(A1,b1)
wqp.AddTask(A2,b2,w2)
wqp.AddTask(A3,b3,w3)

# add contraints
wqp.AddConstraint(C1,lb1,ub1);
wqp.AddConstraint(C2,lb2, ub2);

# solve
wqp.SolveWBC()

result = wqp.GetResult()
print("time of solving QP is ", wqp.GetCostOfQP(), ' s.')
print("wqp result: ", result)

#####################
# RHP TEST
#####################
rhp = WBCKits.RHP()
print(dir(rhp))
# set options
Options.enable_regulation_ = True
Options.weight_regulation_ = 1e-2
rhp.SetOptions(Options)
# add tasks
rhp.AddTask(A1,b1)
rhp.AddTask(A2,b2,w2)
rhp.AddTask(A3,b3,w3)

# add contraints
rhp.AddConstraint(C1,lb1,ub1, hard_constraint=False);
rhp.AddConstraint(C2,lb2, ub2, hard_constraint=False);

# set priority
phi = np.array([[1., 0. , 0.],
                [1., 0.5, 0.],
                [1., 0.5, 0.5]])
chi = [1, 2]
rhp.SetPriority(phi, chi)

# solve
rhp.SolveWBC()

result = rhp.GetResult()
print("rhp result A: ", result)

rhp.Reset()
# add tasks
rhp.AddTask(A1,b1)
rhp.AddTask(A2,b2,w2)
rhp.AddTask(A3,b3,w3)

# add contraints
rhp.AddConstraint(C1,lb1,ub1, hard_constraint=True); # default
rhp.AddConstraint(C2,lb2, ub2, hard_constraint=False);

# set priority
phi = np.array([[1., 1. , 1.],
                [1., 1. , 1.],
                [1., 1. , 1.]])
chi = [1, 1]
rhp.SetPriority(phi, chi)

# solve
rhp.SolveWBC()

result = rhp.GetResult()
print("rhp result B: ", result)

#####################
# HQP TEST
#####################
hqp = WBCKits.HQP()
print(dir(hqp))

# add tasks
hqp.AddTask(A1,b1)
hqp.AddTask(A2,b2,w2)
hqp.AddTask(A3,b3,w3)

# add contraints
hqp.AddConstraint(C1,lb1,ub1, hard_constraint=False);
hqp.AddConstraint(C2,lb2, ub2, hard_constraint=False);

# set priority
phi = [2, 3, 1]
chi = [1, 2]
hqp.SetPriority(phi, chi)

# solve
hqp.SolveWBC()

result = hqp.GetResult()
print("hqp result A: ", result)

hqp.Reset()
# add tasks
hqp.AddTask(A1,b1)
hqp.AddTask(A2,b2,w2)
hqp.AddTask(A3,b3,w3)

# add contraints
hqp.AddConstraint(C1,lb1,ub1, hard_constraint=True); # default
hqp.AddConstraint(C2,lb2, ub2, hard_constraint=False);

# set priority
phi = [1, 1, 1]
chi = [1, 1]
hqp.SetPriority(phi, chi)

# solve
if hqp.SolveWBC():
    print("HQP solved!")

result = hqp.GetResult()
print("hqp result B: ", result)
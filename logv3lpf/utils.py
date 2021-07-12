import opendssdirect as dss
import numpy as np
import pandas as pd
from cmath import sqrt,exp,pi
import scipy
from scipy.sparse import lil_matrix
from scipy.linalg import block_diag
from numpy.linalg import inv
import warnings
from copy import deepcopy

def get_transformer_incidence_matrices(case):
    # Implemented matrices : Ungrounded Wye, Grounded Wye, Delta

    Awye_bus1 = lil_matrix((4,12))
    Awye_bus2 = lil_matrix((4,12))
    AwyeG_bus1 = lil_matrix((3,12))
    AwyeG_bus2 = lil_matrix((3,12))
    Adelta_bus1 = lil_matrix((3,12))
    Adelta_bus2 = lil_matrix((3,12))

    Awye_bus1[[0,1,2,3,3,3],[0,4,8,1,5,9]] = 1
    AwyeG_bus1[[0,1,2],[0,4,8]] = 1

    Awye_bus2[[0,1,2,3,3,3],[2,6,10,3,7,11]] = 1
    AwyeG_bus2[[0,1,2],[2,6,10]] = 1

    Adelta_bus1[[0,0,1,1,2,2],[0,9,1,4,5,8]] = 1
    Adelta_bus2[[0,0,1,1,2,2],[2,11,3,6,7,10]] = 1

    Aconn = {"Wye_bus1":Awye_bus1,"Wye_bus2":Awye_bus2,"WyeG_bus1":AwyeG_bus1,"WyeG_bus2":AwyeG_bus2,"Delta_bus1":Adelta_bus1,"Delta_bus2":Adelta_bus2}

    return Aconn


def get_transformer_yprim(case,i,pu):

    kVs = deepcopy(case.transformers.kVs[i])
    kVspu = [1]*len(kVs)
    conns = case.transformers.Conn[i]
    phases =  deepcopy(case.transformers.phases[i])
    # idx = np.hstack(phases)!=0
    idx = np.array([True,False,True,False,False,True])
    buses = case.transformers.buses[i]
    windings = case.transformers.windings[i]
    Rs = case.transformers.Rs[i]
    noloadloss = float(case.transformers.noloadloss[i])
    imag = float(case.transformers.imag[i])
    Xs = case.transformers.Xs[i]


    for phase in phases:
        if 0 in phase:
            phase.remove(0)
    nl = min([len(x) for x in phases])

    for j in range(0,len(kVs)):
        if (not conns[j] == "Delta") & (nl>=2):
            kVs[j] /= sqrt(3)
            kVspu[j] = sqrt(3)

    kvar = case.transformers.kVAs[i][0]
    kVbases = [case.base[bus]["kVBase"] for bus in buses]

    if windings==2:
        Bt = np.array([[1],[-1]])
        if pu:
            Nt = np.array([[kVspu[0],0],[-kVspu[0],0],[0,kVspu[1]],[0,-kVspu[1]]])
            zsc = nl*(Rs[0]+Rs[1]+1j*Xs[0])
        else:
            Nt = np.array([[1/(kVs[0]*1000),0],[-1/(kVs[0]*1000),0],[0,1/(kVs[1]*1000)],[0,-1/(kVs[1]*1000)]])
            zsc = (Rs[0]+Rs[1]+1j*Xs[0])/(kvar*1000/nl)

        Aconn = get_transformer_incidence_matrices(case)
        Abus1 = Aconn["{:}_bus1".format(conns[0])]
        Abus2 = Aconn["{:}_bus2".format(conns[1])]
        if nl<3:
            Abus1 = Abus1[:nl,:4*nl]
            Abus2 = Abus2[:nl,:4*nl]

        A = scipy.sparse.vstack([Abus1,Abus2])
            
        Taubus1 = np.kron(np.eye(Abus1.shape[0],dtype=int), 1/case.transformers.taps[i][0]) ## Inverse of the tap
        Taubus2 = np.kron(np.eye(Abus2.shape[0],dtype=int), 1/case.transformers.taps[i][1]) ## Inverse of the tap
        Tau = block_diag(Taubus1,Taubus2)

        B = np.kron(np.eye(nl,dtype=int),Bt)
        N = np.kron(np.eye(nl,dtype=int),Nt)
        z = np.kron(np.eye(nl,dtype=int),zsc)

        Yprim = A@N@B@inv(z)@B.T@N.T@A.T

    elif windings==3:
        Bt = np.array([[1,1],[-1,0],[0,-1]])
        if pu:
            Nt = np.array([[kVspu[0],0,0],[-kVspu[0],0,0],[0,kVspu[1],0],[0,-kVspu[1],0],[0,0,kVspu[2]],[0,0,-kVspu[2]]])
            zsc = nl*np.array([[Rs[0] + Rs[1] + 1j*Xs[0],Rs[0]+0.5j*(Xs[0]+Xs[1]-Xs[2])],[Rs[0]+0.5j*(Xs[0]+Xs[1]-Xs[2]),Rs[0] + Rs[2]+1j*Xs[1]]])
        else:
            Nt = np.array([[1/(kVs[0]*1000),0,0],[-1/(kVs[0]*1000),0,0],[0,1/(kVs[1]*1000),0],[0,-1/(kVs[1]*1000),0],[0,0,1/(kVs[2]*1000)],[0,0,-1/(kVs[2]*1000)]])
            zsc = np.array([[Rs[0] + Rs[1] + 1j*Xs[0],Rs[0]+0.5j*(Xs[0]+Xs[1]-Xs[2])],[Rs[0]+0.5j*(Xs[0]+Xs[1]-Xs[2]),Rs[0] + Rs[2]+1j*Xs[1]]])/(kvar*1000/nl)

        Yprim = (Nt@Bt@inv(zsc)@Bt.T@Nt.T)[idx,:][:,idx]
        # yoc = noloadloss/100 - 1j*imag/100
        # Yprim[-2:,-2:] += np.array([[yoc,-yoc],[-yoc,yoc]])


    # if np.linalg.matrix_rank(Yprim)<Yprim.shape[0]:
    #     Yprim += np.random.rand(Yprim.shape[0],Yprim.shape[0])/1000000

    # return Tau@Yprim@Tau
    return Yprim

def create_df(case):
    results = pd.DataFrame()
    for algo in ["openDSS","logv3lpf"]:
        vm = []
        va = []
        phase = []
        buses = []
        for bus in case.results[algo]["vm"]:
            phases = case.bus_phases[bus]
            missing_phases = np.setdiff1d([1,2,3],phases)
            vm.append(case.results[algo]["vm"][bus])
            va.append(case.results[algo]["va"][bus])
            phase.append(phases)
            for missing_phase in missing_phases:
                vm.append(np.nan)
                va.append(np.nan)
                phase.append(missing_phase)
            buses.append([bus]*3)
        vm = np.hstack(vm)
        va = np.hstack(va)
        phase = np.hstack(phase)
        buses = np.hstack(buses)
        df = pd.DataFrame(np.vstack([buses,vm,va,phase]).T,columns = ["bus","vm","va","phase"])
        df["algorithm"] = algo
        results = results.append(df)
    results.vm = results.vm.astype(float)
    results.va = results.va.astype(float)

    return results



    
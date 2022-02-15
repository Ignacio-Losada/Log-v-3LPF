import numpy as np
import pandas as pd
import opendssdirect as dss
from cmath import exp,pi,sqrt,log
from itertools import product
import scipy
from copy import deepcopy
from scipy.sparse import csr_matrix,lil_matrix
from scipy.sparse.linalg import spsolve as solve
from scipy.linalg import block_diag
from scipy.sparse import block_diag as sparse_block_diag
from scipy.sparse import coo_matrix as sparse_matrix
from logv3lpf.utils import get_transformer_yprim
from scipy.linalg import pinv,inv
from math import log

def get_line_matrices(case):

    Nbrl,Nn = case.Nbrl,case.Nn

    # Define incidence matrix
    E = csr_matrix((Nn,2*Nbrl))

    Yu,Yθ,y = [],[],[]
    Yprims = []
    counter = 0

    for i in range(0,len(case.lines)):

        fbus,tbus = case.lines.fbus[i],case.lines.tbus[i]
        zbase = case.base[fbus]["ZBase"]
        fphases = case.lines.fphase[i]
        tphases = case.lines.tphase[i]
        idx_fphases = [case.bus_phases[fbus].index(phase) for phase in fphases]
        idx_tphases = [case.bus_phases[tbus].index(phase) for phase in tphases]
        if len(fphases) == len(tphases):
            nl = len(fphases)
        Δ3_n = case.base[fbus]["Δ3"]
        Δ3_m = case.base[tbus]["Δ3"]

        fidx = [case.node_from_name["{:}.{:}".format(fbus,j)] for j in fphases]
        tidx = [case.node_from_name["{:}.{:}".format(tbus,j)] for j in tphases]

        E[fidx,[j+counter for j in range(0,nl)]] = 1
        counter+=nl
        E[tidx,[j+counter for j in range(0,nl)]] = 1
        counter+=nl

        Δ3nphase = np.diag(Δ3_n[idx_fphases])
        Δ3mphase = np.diag(Δ3_m[idx_tphases])
        Δ3phase_block = block_diag(Δ3nphase, Δ3mphase)


        Ynm_n = case.lines.ynm[i]*zbase
        Ynm_m = -case.lines.ynm[i]*zbase
        Ymn_n = -case.lines.ynm[i]*zbase
        Ymn_m = case.lines.ynm[i]*zbase
        Ynm_s = case.lines.ysh[i]*zbase
        Ymn_s = case.lines.ysh[i]*zbase

        Yprim = np.block([[(Ynm_n+Ynm_s/2),Ynm_m],[Ymn_n,(Ymn_m+Ymn_s/2)]])
        Yprims.append(Yprim)
        Ytilde = Δ3phase_block@Yprim.conj()@Δ3phase_block.conj()
        ytilde = Ytilde.sum(axis=0)

        # Y.append(np.block([Ytilde + np.diag(ytilde), -1j*Ytilde+1j*np.diag(ytilde)]))
        Yu.append(Ytilde + np.diag(ytilde))
        Yθ.append(-1j*Ytilde+1j*np.diag(ytilde))
        y.append(ytilde)

    # Y = csr_matrix(block_diag(*Y))
    if len(y) != 0:
        y = np.hstack(y).reshape(-1,1)
    else:
        y = np.array([])

    case.lines["Yprim"] = Yprims
    case.Yul = Yu
    case.Yθl = Yθ
    case.yl = y
    case.El = E


def get_transformer_matrices(case):

    Nbrt,Nbrr,Nn = case.Nbrt,case.Nbrr,case.Nn

    # Define incidence matrix
    Et = csr_matrix((Nn,Nbrt))

    Yprims = []
    Yprim_idx = []
    Yut,Yθt,yt = [],[],[]
    counter = 0

    # Transformers are modeled as power delivery elements (lines). However, Yprim is not symmetric
    for transformer in case.transformer_names:
        i = case.transformers[case.transformers.name == transformer].index.values[0]
        Yprim_idx.append(i)
        buses = case.transformers.buses[i]
        phases = deepcopy(case.transformers.phases[i])
        [phase.remove(0) for phase in phases if 0 in phase]

        idx_phases = [[case.bus_phases[buses[j]].index(phase) for phase in phases[j]] for j in range(0,len(phases))]
        Δ3s = [np.diag(case.base[buses[j]]["Δ3"][idx_phases[j]]) for j in range(0,len(idx_phases))]

        idx = np.hstack([[case.node_from_name["{:}.{:}".format(buses[j],phase)] for phase in phases[j]] for j in range(0,len(buses))])

        Et[idx,[j+counter for j in range(0,len(idx))]] = 1
        counter += len(idx)

        Δ3phase_block = block_diag(*Δ3s)

        # Yprim,Tau = get_alt_transformer_yprim(case,i)
        Tau,Yprim = get_transformer_yprim(case,i,pu=True) ##BEWARE THIS YPRIM DOES NOT INCLUDE LOSSES
        Yprims.append(Yprim)

        Ytilde = Δ3phase_block@Yprim.conj()@Δ3phase_block.conj()
        ytilde = Ytilde.sum(axis=0)

        Yut.append(Ytilde + np.diag(ytilde))
        Yθt.append(-1j*Ytilde+1j*np.diag(ytilde))
        yt.append(ytilde)

    # Yt = csr_matrix(block_diag(*Yt))
    if len(yt) != 0:
        yt = np.hstack(yt).reshape(-1,1)
    else:
        yt = np.array([])

    case.Yut = Yut
    case.Yθt = Yθt
    case.yt = yt
    case.Et = Et


def get_regulator_matrices(case):

    Nbrt,Nbrr,Nn = case.Nbrt,case.Nbrr,case.Nn

    Er = csr_matrix((Nn,Nbrr))  

    Yprims = []
    Yprim_idx = []
    Yur,Yθr,yr = [],[],[]
    counter = 0
        
    for regulator in case.regulator_names:
        i = case.transformers[case.transformers.name == regulator].index.values[0]
        Yprim_idx.append(i)
        buses = case.transformers.buses[i]
        phases = deepcopy(case.transformers.phases[i])
        [phase.remove(0) for phase in phases if 0 in phase]

        idx_phases = [[case.bus_phases[buses[j]].index(phase) for phase in phases[j]] for j in range(0,len(phases))]
        Δ3s = [np.diag(case.base[buses[j]]["Δ3"][idx_phases[j]]) for j in range(0,len(idx_phases))]


        idx = np.hstack([[case.node_from_name["{:}.{:}".format(buses[j],phase)] for phase in phases[j]] for j in range(0,len(buses))])

        Er[idx,[j+counter for j in range(0,len(idx))]] = 1
        counter += len(idx)

        Δ3phase_block = block_diag(*Δ3s)

        # Yprim,Tau = get_alt_transformer_yprim(case,i)
        # Tau = np.diag(1/np.array(case.transformers.taps[i]))
        # Yprim = Tau@get_transformer_yprim(case,i,pu=True)@Tau ##BEWARE THIS YPRIM DOES NOT INCLUDE LOSSES
        Tau,Yprim = get_transformer_yprim(case,i,pu=True) ##BEWARE THIS YPRIM DOES NOT INCLUDE LOSSES
        Yprims.append(Yprim)
        Ytilde = Δ3phase_block@Tau@Yprim.conj()@Tau@Δ3phase_block
        ytilde = Ytilde.sum(axis=0)

        Yur.append(Ytilde + np.diag(ytilde))
        Yθr.append(-1j*Ytilde+1j*np.diag(ytilde))
        yr.append(ytilde)

    # Yr = csr_matrix(block_diag(*Yr))
    if len(yr) != 0:
        yr = np.hstack(yr).reshape(-1,1)
        # case.transformers["Yprim"] = [Yprims[Yprim_idx.index(j)] for j in range(0,len(Yprim_idx))] 
    else:
        yr = np.array([])

    case.Yur = Yur
    case.Yθr = Yθr
    case.yr = yr
    case.Er = Er


def calculate_base_matrices(case):
    refvm = case.refvm
    refva = case.refva
    refbus = case.refbus
    Nn = case.Nn

    if len(case.yl) == 0:
        Ylbus = csr_matrix((2*Nn,2*Nn))
    else:
        Ylbus = sparse_block_diag((case.El,case.El))@scipy.sparse.bmat([[sparse_block_diag(case.Yul).real,sparse_block_diag(case.Yθl).real],[sparse_block_diag(case.Yul).imag,sparse_block_diag(case.Yθl).imag]])@sparse_block_diag((case.El,case.El)).T
    if len(case.yt) == 0:  
        Ytbus = csr_matrix((2*Nn,2*Nn))
    else:
        Ytbus = sparse_block_diag((case.Et,case.Et))@scipy.sparse.bmat([[sparse_block_diag(case.Yut).real,sparse_block_diag(case.Yθt).real],[sparse_block_diag(case.Yut).imag,sparse_block_diag(case.Yθt).imag]])@sparse_block_diag((case.Et,case.Et)).T
    if len(case.yr) == 0:
        Yrbus = csr_matrix((2*Nn,2*Nn))
    else:
        Yrbus = sparse_block_diag((case.Er,case.Er))@scipy.sparse.bmat([[sparse_block_diag(case.Yur).real,sparse_block_diag(case.Yθr).real],[sparse_block_diag(case.Yur).imag,sparse_block_diag(case.Yθr).imag]])@sparse_block_diag((case.Er,case.Er)).T
    A = Ylbus + Ytbus + Yrbus

    y = np.zeros((Nn,1), dtype=np.complex)
    if len(case.yl) != 0:
        y += case.El@case.yl
    if len(case.yr) != 0:
        y += case.Er@case.yr
    if len(case.yt) != 0:
        y += case.Et@case.yt
    sshunt = np.vstack([y.real,y.imag])

    #Gets masks for the admittance and incidence matrices (gets rid of the reference bus)
    ref_bus_phases_real = [case.node_from_name["{:}.{:}".format(refbus,i)] for i in case.bus_phases[refbus]]
    ref_bus_phases_imag = [el+case.Nn for el in ref_bus_phases_real]
    ref_bus_phases = ref_bus_phases_real+ref_bus_phases_imag
    mask = np.setdiff1d(list(range(0, 2*case.Nn)),ref_bus_phases)

    # Calculate A^-1
    Ainv = inv(A[mask,:][:,mask].toarray()) ## This is what takes most of the time when parsing the files (~2 min in the IEEE8500). Everything else takes no more than 10 seconds even for a large case
    Adrop = A[:,ref_bus_phases]
    sdrop = Adrop@np.array([log(refvm)]*3 + [np.deg2rad(refva)]*3).reshape(-1,1)

    E = scipy.sparse.bmat([[case.EloadYl,None],[None,case.EloadYl]])
    U = csr_matrix(E)[mask,:]
    V = csr_matrix(E.T)[:,mask]

    k = U.shape[1]
    I = csr_matrix(np.eye(k,dtype=float))


    AinvU = Ainv@U
    VAinv = V@Ainv
    VAinvU = VAinv@U
    
    case.mask = mask
    case.Ylbus = Ylbus
    case.Ytbus = Ytbus
    case.Yrbus = Yrbus
    case.A = A
    case.I = I
    case.U = U
    case.V = V
    case.E = E
    case.sshunt = sshunt[mask]
    case.Adrop = Adrop
    case.sdrop = sdrop[mask]
    case.Ainv = Ainv
    case.AinvU = AinvU
    case.VAinv = VAinv
    case.VAinvU = VAinvU


def get_loads(case):
    Ω = case.Ω
    Π = case.Π

    s = []
    Lu = []
    Lθ = []

    for i in range(len(case.loads)):
        phases = np.array(case.loads.phases[i])-1
        bus = case.loads.bus[i]
        nphases = case.loads.nphases[i]
        model = case.loads.model[i]
        kV = case.loads.kV[i]
        kW = case.loads.kW[i]
        kvar = case.loads.kvar[i]
        isdelta = case.loads.isdelta[i]


        idx_phases = [case.bus_phases[bus].index(phase) for phase in case.loads.phases[i]]
        Δ3phase = np.diag(case.base[bus]["Δ3"][idx_phases])

        if model == 1:
            base = case.base[bus]["VABase"]
            p = np.repeat((kW*1000+1j*kvar*1000)/(nphases*base),nphases).reshape(-1,1)
            if isdelta:
                if nphases == 1:
                    Ωphase = Ω[phases,phases[0]].reshape(-1,1)
                    s.append(Ωphase@p)
                elif nphases == 3:
                    s.append(Ω@p)
            else:
                s.append(p)

        elif model == 2:
            base = case.base[bus]["ZBase"]
            if (nphases == 1) | (nphases == 3):
                z = np.repeat((kV*1000)**2 / np.conj((kW*1000+1j*kvar*1000)),nphases)
            elif nphases == 2:
                z = np.repeat((2/3)*(kV*1000)**2 / np.conj((kW*1000+1j*kvar*1000)),nphases)
            z /= base

            if isdelta:
                # Δ3phase = np.diag(case.base[bus]["Δ3"][phases,phases])
                # Delta loads cannot be 2 phase because otherwise the buses they are connected to are ambiguous
                if nphases == 1:
                    Πphase = Π[phases,phases[0]].reshape(-1,1)
                elif nphases == 3:
                    Πphase = Π
                
                X = Δ3phase.conj()@Πphase@np.diag((1/z).conj())@Πphase.T
                s.append(np.sum(Δ3phase@X.T,axis=1).reshape(-1,1))

                Lu.append((Δ3phase@X.T+np.diag(np.sum(Δ3phase@X.T,axis=1))))
                Lθ.append(1j*(np.diag(np.sum(Δ3phase@X.T,axis=1)) - Δ3phase@X.T))           
            else:
                # Δ3phase = np.diag(case.base[bus]["Δ3"][phases,phases])
                s.append((1/z).conj().reshape(-1,1))
                
                X = Δ3phase.conj()@np.diag((1/z).conj())
                load = Δ3phase@X.T+np.diag(np.sum(Δ3phase@X.T,axis=1))
                Lu.append(load)
                Lθ.append(np.zeros(load.shape))




        elif model == 5:
            base = case.base[bus]["IBase"]
            shift = case.base[bus]["Shift"]
            power_angle = np.angle(kW+1j*kvar,deg=True)
            phaseidx = [case.bus_phases[bus].index(phase) for phase in np.array(case.loads.phases[i])]
            # Δ3phase = np.diag(case.base[bus]["Δ3"][phases,phases])
            if (nphases == 1) | isdelta:
                I = np.repeat(abs(kW+1j*kvar)/(nphases*kV),nphases)
            else:
                I = np.repeat(abs(kW+1j*kvar)/(sqrt(3)*nphases*kV),nphases)
            I /= base

            if isdelta:
                va = shift + pi/6
                # Delta loads cannot be 2 phase because otherwise the buses they are connected to are ambiguous
                if nphases == 1:
                    I = (I*np.exp( 1j * (va[phaseidx[0]]-power_angle) )).reshape(-1,1)
                    Πphase = Π[phases,phases[0]].reshape(-1,1)
                elif nphases == 3:
                    I = (I*np.exp( 1j * (va[phaseidx]-power_angle) )).reshape(-1,1)
                    Πphase = Π

                load = Πphase@I.conj()
                sload = Δ3phase@load
                
                s.append(sload.reshape(-1,1))
                Lu.append(np.diag(sload.T[0]))
                Lθ.append(np.diag(1j*sload.T[0]))


            else:
                va = shift
                I = (I*np.exp( 1j * (2*pi/360)*(va[phaseidx]-power_angle))).reshape(-1,1)
                sload = Δ3phase@I.conj()
                s.append(sload)
                Lu.append(np.diag(sload.T[0]))
                Lθ.append(np.diag(1j*sload.T[0]))
        else:
            print("Load {:} is model {:} which is not supported. Try model=1,model=2 or model=5".format(case.loads.name[i],case.loads.model[i]))



    for i in range(len(case.capacitors)):
        phases = np.array(case.capacitors.phases[i])-1
        bus = case.capacitors.bus[i]
        nphases = case.capacitors.nphases[i]
        kV = case.capacitors.kV[i]
        kvar = case.capacitors.kvar[i]
        isdelta = case.capacitors.isdelta[i]
        base = case.base[bus]["ZBase"]

        idx_phases = [case.bus_phases[bus].index(phase) for phase in case.capacitors.phases[i]]
        Δ3phase = np.diag(case.base[bus]["Δ3"][idx_phases])

        if (nphases == 1) | (nphases == 3):
            z = np.repeat((kV*1000)**2 / np.conj((-1j*kvar*1000)),nphases)
        elif nphases == 2:
            z = np.repeat((2/3)*(kV*1000)**2 / np.conj((-1j*kvar*1000)),nphases)
        z /= base

        if isdelta:
            # Δ3phase = np.diag(case.base[bus]["Δ3"][phases,phases])
            # Delta loads cannot be 2 phase because otherwise the buses they are connected to are ambiguous
            if nphases == 1:
                Πphase = Π[phases,phases[0]].reshape(-1,1)
            elif nphases == 3:
                Πphase = Π
            
            X = Δ3phase.conj()@Πphase@np.diag((1/z).conj())@Πphase.T
            s.append(np.sum(Δ3phase@X.T,axis=1).reshape(-1,1))

            Lu.append((Δ3phase@X.T+np.diag(np.sum(Δ3phase@X.T,axis=1))))
            Lθ.append(1j*(np.diag(np.sum(Δ3phase@X.T,axis=1)) - Δ3phase@X.T))           

        else:
            # Δ3phase = np.diag(case.base[bus]["Δ3"][phases,phases])
            s.append((1/z).conj().reshape(-1,1))
            
            X = Δ3phase.conj()@np.diag((1/z).conj())
            load = Δ3phase@X.T+np.diag(np.sum(Δ3phase@X.T,axis=1))
            Lu.append(load)
            Lθ.append(np.zeros(load.shape))

    case.s = -np.vstack(s)
    if len(Lu)!=0:
        Lu = sparse_block_diag(Lu)
    if len(Lθ)!=0:
        Lθ = sparse_block_diag(Lθ)

    case.Lu = Lu
    case.Lθ = Lθ

def rank_k_correction_solve(case,process_solution):

    get_regulator_matrices(case)
    get_loads(case)

    mask = case.mask

    Lu = case.Lu
    Lθ = case.Lθ

    ## Get s vector
    sload = case.EloadS@case.s
    sbus = np.vstack([sload.real,sload.imag])[mask] - case.sdrop - case.sshunt

    if type(Lu) == list:
        Ybusinv = case.Ainv
    else:
        L = scipy.sparse.bmat([[Lu.real,Lθ.real],[Lu.imag,Lθ.imag]])

        ## Applying Woodbury to (A + UBV)^-1 when B is singular 
        ## See https://ecommons.cornell.edu/bitstream/handle/1813/32749/BU-647-M.pdf;jsessionid=8284B7F90591F0956C8438C1EFF7C732?sequence=1
        ## A^-1 - A^-1U(I + BVA^-1U)^1BVA^-1
        Ybusinv = case.Ainv-case.AinvU@inv(case.I+L@case.VAinvU)@L@case.VAinv

    sol = Ybusinv@sbus
    vm = np.exp(sol[:int(len(sol)/2)])
    va = sol[int(len(sol)/2):]

    case.vm = vm
    case.va = va

    if process_solution:
        process_logv3lpf_solution(case)

def update_logv3lpf(case,process_solution):

    calculate_base_matrices(case)

    get_loads(case)

    mask = case.mask

    Lu = case.Lu
    Lθ = case.Lθ

    ## Get s vector
    sload = case.EloadS@case.s
    sbus = np.vstack([sload.real,sload.imag])[mask] - case.sdrop - case.sshunt

    if type(Lu) == list:
        Ybusinv = case.Ainv
    else:
        L = scipy.sparse.bmat([[Lu.real,Lθ.real],[Lu.imag,Lθ.imag]])

        ## Applying Woodbury to (A + UBV)^-1 when B is singular 
        ## See https://ecommons.cornell.edu/bitstream/handle/1813/32749/BU-647-M.pdf;jsessionid=8284B7F90591F0956C8438C1EFF7C732?sequence=1
        ## A^-1 - A^-1U(I + BVA^-1U)^1BVA^-1
        Ybusinv = case.Ainv-case.AinvU@inv(case.I+L@case.VAinvU)@L@case.VAinv

    sol = Ybusinv@sbus
    vm = np.exp(sol[:int(len(sol)/2)])
    va = sol[int(len(sol)/2):]

    case.vm = vm
    case.va = va

    if process_solution:
        process_logv3lpf_solution(case)

def check_logv3lpf_performance(case):
    k = 0
    if len(case.loads[case.loads.model!=1].phases)!=0:
        k += len(np.hstack(case.loads[case.loads.model!=1].phases))
    if len(case.capacitors.phases)!=0:
        k += len(np.hstack(case.capacitors.phases))
    n = case.Nn*2
    for name in case.regcontrols.name:
        for phases in case.transformers[case.transformers.name==name].phases:
            el = np.array(phases)
            k += len(el[el!=0])
    k*=2
    O_KLU = int((2/3)*(n**3)+(2*n)**2)
    O_logv = 3*(k)**3 + (2*n*(k**2)) + n**2 + k
    print("Nodes = ",case.Nn)
    print("k = ",k)
    print("Complexity of KLU decomposition, {:.2e} FLOPS ".format(O_KLU))
    print("Complexity of LogV algorithm, {:.2e} FLOPS ".format(O_logv))
    print("Ratio {:.2f}".format(O_KLU/O_logv))


def process_logv3lpf_solution(case):
    dssvm = case.bus_phases.copy()
    dssva = case.bus_phases.copy()

    ref_bus_phases_real = [case.node_from_name["{:}.{:}".format(case.sourcebus,i)] for i in case.bus_phases[case.sourcebus]]
    vm = np.zeros(case.Nn)
    vm[ref_bus_phases_real] = [case.refvm]*len(ref_bus_phases_real)
    vm[case.mask[:int(len(case.mask)/2)]] = case.vm.T[0]

    for name in dssvm.keys():
        idx = [case.node_from_name["{:}.{:}".format(name,phase)] for phase in case.bus_phases[name]]
        dssvm[name] = vm[idx]

    va = np.zeros(case.Nn)
    va[ref_bus_phases_real] = [0,-2*pi/3,2*pi/3]
    angles = case.va.T[0] +np.hstack([case.base[bus]["Shift"] for bus in np.unique([i.split(".")[0] for i in case.node_from_idx.values()]) if bus!=case.sourcebus])

    va[case.mask[:int(len(case.mask)/2)]] = angles

    for name in dssva.keys():
        idx = [case.node_from_name["{:}.{:}".format(name,phase)] for phase in case.bus_phases[name]]
        degrees = np.rad2deg(va[idx])
        for i in range(0,len(degrees)):
            if degrees[i]>180:
                degrees[i] = -(360-degrees[i])

        dssva[name] = degrees


    if not hasattr(case,"results"):
        case.results = {}

    case.results["logv3lpf"] = {}
    case.results["logv3lpf"]["vm"] = dssvm
    case.results["logv3lpf"]["va"] = dssva
    # self.results["openDSS"]["regtap"] = tap


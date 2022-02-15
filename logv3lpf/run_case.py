import logv3lpf.linpf as linpf
import opendssdirect as dss
from logv3lpf.DSSParser import DSScase
import numpy as np
from logv3lpf.math_functions import polar_to_cartesian
import datetime

class case(object):
    def __init__(self, file,sourcebus,refvm,refva):
        DSScase.__init__(self,file,sourcebus,refvm,refva)
        DSScase.parse_DSS_file(self)
        print("Calculating Line Matrices")
        linpf.get_line_matrices(self)
        print("Calculating Transformer Matrices")
        linpf.get_transformer_matrices(self)
        print("Calculating Regulator Matrices")
        linpf.get_regulator_matrices(self)
        print("Calculating Base Matrices")
        linpf.calculate_base_matrices(self)
        case.sourcebus = sourcebus
        linpf.check_logv3lpf_performance(self)
    
    def run_case(self,algorithm):
        if algorithm=="logv3lpf":
            print("Solving base case logv3lpf...")
            case.run_logv3lpf(self)
        elif algorithm=="NR":
            case.run_newton_raphson(self)
        elif algorithm=="GS":
            case.run_gauss_seidel(self)
        elif algorithm=="FBS":
            case.run_forward_backward_sweep(self)
        elif algorithm=="openDSS":
            print("Solving DSS case in OpenDSS...")
            dss.Solution.Solve()
            print("Storing solution from OpenDSS...")
            DSScase.process_openDSS_solution(self)
        else:
            raise Error('Algorithm not implemented. Options are logv3lpf, NR, GS, FBS. See documentation for details')
    
    def run_logv3lpf(self):
        linpf.rank_k_correction_solve(self,True)
    def update_logv3lpf(self):
        linpf.update_logv3lpf(self,True)

    def run_newton_raphson(self):
        return False
    def run_gauss_seidel(self):
        return False
    def run_forward_backward_sweep(self):
        return False

    def set_load_kvar(self,name,value):
        try:
            self.loads.loc[self.loads.name==name,"kvar"] = value
        except ValueError:
            print("Load {:} does not exist!".format(name))

    def set_load_kW(self,name,value):
        try:
            self.loads.loc[self.loads.name==name,"kW"] = value
        except ValueError:
            print("Load {:} does not exist!".format(name))

    def set_capacitor_kvar(self,name,value):
        try:
            self.capacitors.loc[self.capacitors.name==name,"kvar"] = value
        except ValueError:
            print("Capacitor {:} does not exist!".format(name))

    def set_transformer_tap(self,name,value,bus):
        try:
            self.transformers.loc[self.transformers.name==name, bus] = value
            if name in list(self.regcontrols.name):
                self.regcontrols.loc[self.regcontrols.name==name,"tap"] = value
        except ValueError:
            print("Transformer {:} does not exist or {:} is not a valid bus!".format(name.bus))

    def get_line_current(self,name,pu):
        idx = self.lines[self.lines.name == name].index[0]
        fbus,tbus = self.lines.fbus[idx],self.lines.tbus[idx]
        fphase,tphase = self.lines.fphase[idx],self.lines.tphase[idx]
        fidx = [self.bus_phases[fbus].index(i) for i in fphase]
        tidx = [self.bus_phases[tbus].index(i) for i in tphase]
        un = np.log(np.array([self.results["logv3lpf"]["vm"][fbus][i] for i in fidx]))
        um = np.log(np.array([self.results["logv3lpf"]["vm"][tbus][i] for i in tidx]))
        θn = np.log(np.array([self.results["logv3lpf"]["va"][fbus][i]-case.base[fbus][i] for i in fidx]))
        θm = np.log(np.array([self.results["logv3lpf"]["va"][tbus][i]-case.base[tbus][i] for i in tidx]))
        u = np.hstack([un,um]).reshape(-1,1)
        θ = np.hstack([θn,θm]).reshape(-1,1)
        y = case.yl[case.lines["nphases"][:idx].sum()*2:case.lines["nphases"][:idx].sum()*2+case.lines["nphases"][idx]*2]

        s = case.Yul[idx]@u + case.Yθl[idx]@θ + y

        if pu:
            return s
        else:
            scaling = np.diag([self.base[fbus]["IBase"]]*len(fidx)+[self.base[tbus]["IBase"]]*len(tidx))
            return scaling@I


        



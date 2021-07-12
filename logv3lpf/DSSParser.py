import opendssdirect as dss
from numpy.linalg import inv
from math import pi,sqrt
import pandas as pd
import numpy as np
from opendssdirect.utils import Iterator
from cmath import exp
import types
import inspect
import warnings
import networkx as nx
from copy import deepcopy
from scipy.sparse import csr_matrix


class DSScase(object):
    def __init__(self, file,sourcebus,refvm,refva):
        self.file = file
        self.refbus = sourcebus
        self.refvm = refvm
        self.refva = refva

    def parse_DSS_file(self):
        print("Compiling DSS file in OpenDSS...")
        dss.run_command(r'{:}'.format(self.file))
        print("Parsing DSS file in LogV3LPF format...")
        DSScase.get_capacitors_from_dss(self)
        DSScase.get_transformers_from_dss(self)
        DSScase.get_loads_from_dss(self)
        DSScase.get_lines_from_dss(self)
        DSScase.get_regcontrols_from_dss(self)
        DSScase.get_buses_from_dss(self)
        DSScase.get_constants(self)
        DSScase.create_load_model(self)
        print("Done parsing")


    def get_capacitors_from_dss(self):
        capacitor_names = dss.Capacitors.AllNames()
        buses = []
        nphases = []
        phases = []
        kV = []
        kvar = []
        conn = []
        for name in capacitor_names:
            dss.Capacitors.Name(name)
            kV.append(dss.Capacitors.kV())
            kvar.append(dss.Capacitors.kvar())
            conn.append(dss.Capacitors.IsDelta())
            busesstr = dss.CktElement.BusNames()[0].split(".")
            buses.append(busesstr[0])

            if (len(busesstr)-1)==0:
                nphases.append(3)
                phases.append([1,2,3])
            else:
                nphases.append(len(busesstr[1:]))
                phase = [int(i) for i in busesstr[1:]]
                if 4 in phase:
                    print("Capacitor {:} is ungrounded. Modeling is not accurate".format(name))
                phases.append(phase)

        capacitors_dict = {"name":capacitor_names,"nphases":nphases,"phases":phases,"bus":buses,"kV":kV,"kvar":kvar,"isdelta":conn}
        capacitors = pd.DataFrame(capacitors_dict)
        self.capacitors = capacitors

    def get_transformers_from_dss(self):
        transformer_names = dss.Transformers.AllNames()
        kV = []
        kVA = []
        windings = []
        Conn = []
        R = []
        X = []
        tap = []
        max_tap = []
        min_tap = []
        buses = []
        phases = []
        nphases = []
        isgrounded = []
        noloadloss = []
        imag = []

        for name in transformer_names:
            dss.Transformers.Name(name)
            nwindings = dss.Transformers.NumWindings()
            windings.append(nwindings)
            nbuses = len(dss.CktElement.BusNames())
            busnames = []
            busphases = []
            for i in range(nbuses):
                busstr = dss.CktElement.BusNames()[i].split(".")
                busisgrounded = []
                if (len(busstr)-1) == 0:
                    phase = [1,2,3]
                else:
                    phase = [int(i) for i in busstr[1:]]
                busnames.append(busstr[0])
                busphases.append(phase)
            phases.append(busphases)
            # u, ind = np.unique(busnames, return_index=True)
            # buses.append(u[np.argsort(ind)].tolist())
            buses.append(busnames)
            if nwindings == 2:
                X.append([dss.Transformers.Xhl()])
            else:
                X.append([dss.Transformers.Xhl(),dss.Transformers.Xht(),dss.Transformers.Xlt()])

            dss.Text.Command ('? transformer.{0}.%noloadloss'.format(name))
            noloadloss.append(dss.Text.Result())
            dss.Text.Command ('? transformer.{0}.%imag'.format(name))
            imag.append(dss.Text.Result())

            kVs = []
            kVAs = []
            DeltaConns = []
            Rs = []
            taps = []

            for i in range(nwindings):
                dss.Transformers.Wdg(i+1)
                
                taps.append(dss.Transformers.Tap())
                kVs.append(dss.Transformers.kV())
                kVAs.append(dss.Transformers.kVA())
                DeltaConns.append(dss.Transformers.IsDelta())
                Rs.append(dss.Transformers.R())
            tap.append(taps)
            max_tap.append(dss.Transformers.MaxTap())
            min_tap.append(dss.Transformers.MinTap())
            kV.append(kVs)
            kVA.append(kVAs)
            R.append(Rs)

            Conns = []
            for i in range(0,len(DeltaConns)):
                if DeltaConns[i] == True:
                    Conns.append("Delta")
                else:
                    if 4 in busphases[i]:
                        Conns.append("Wye")
                    else:
                        Conns.append("WyeG")

            Conn.append(Conns)
        if len(buses) == 0:
            transformers_dict = {"name":[],"windings":[], "buses":[], "phases":[],"Conn":[],"max_tap":[],"min_tap":[],"Rs":[],"Xs":[],"kVs":[],"kVAs":[],"taps":[],"%noloadloss":[],"%imag":[]}

        else:
            transformers_dict = {"name":transformer_names,"windings":windings, "buses":buses, "phases":phases,"Conn":Conn,"max_tap":max_tap,"min_tap":min_tap,"Rs":R,"Xs":X,"kVs":kV,"kVAs":kVA,"taps":tap,"noloadloss":noloadloss,"imag":imag}
        self.transformers = pd.DataFrame(transformers_dict)

    def get_regcontrols_from_dss(self):
        regcontrol_names = dss.RegControls.AllNames()
        ct = []
        nt = []
        r = []
        x = []
        d = []
        vreg = []
        tap = []
        band = []
        maxtap = []
        transformers = []

        for name in regcontrol_names:
            dss.RegControls.Name(name)

            ct.append(dss.RegControls.CTPrimary())
            r.append(dss.dss.RegControls.ForwardR())
            x.append(dss.RegControls.ForwardX())
            nt.append(dss.RegControls.PTRatio())
            vreg.append(dss.RegControls.ForwardVreg())
            d.append(dss.RegControls.Delay())
            tap.append(1+dss.RegControls.TapNumber()*0.00625)
            band.append(dss.RegControls.ForwardBand())
            maxtap.append(dss.RegControls.MaxTapChange())
            transformers.append(dss.RegControls.Transformer())
        regcontrol_dict = {"name":regcontrol_names,"transformer":transformers,"vreg":vreg,"band":band,"tap":tap,"ct":ct,"nt":nt,"R":r,"X":x,"delay":d,"max_tap":maxtap}
        self.regcontrols =  pd.DataFrame(regcontrol_dict)


    def get_loads_from_dss(self):
        load_names = dss.Loads.AllNames()
        model = []
        kV = []
        kW = []
        kvar = []
        conn = []
        buses = []
        nphases = []
        phases = []

        for name in load_names:
            dss.Loads.Name(name)

            model.append(dss.Loads.Model())
            kV.append(dss.Loads.kV())
            kW.append(dss.Loads.kW())
            kvar.append(dss.Loads.kvar())
            conn.append(dss.Loads.IsDelta())
            busstr = dss.CktElement.BusNames()[0].split(".")
            phase = [int(i) for i in busstr[1:]]
            if phase == []:
                phase = [1,2,3]
            if 4 in phase:
                print("Load {:} is ungrounded. Modeling is not accurate".format(name))
            buses.append(busstr[0])
            nphases.append(dss.Loads.Phases())
            phases.append(phase)

        loads_dict = {"name":load_names,"nphases":nphases,"phases":phases,"model":model,"bus":np.vstack(buses)[:,0],"kV":kV,"kW":kW,"kvar":kvar,"isdelta":conn}
        loads = pd.DataFrame(loads_dict)
        
        self.loads = loads

    def get_lines_from_dss(self):
        CapConversionFactor=1j*(1e-9)*2*pi*dss.Solution.Frequency()

        line_names = dss.Lines.AllNames()
        buses = []
        switches = []
        phases = []
        fphases = []
        tphases = []
        ynms = []
        yshs = []

        for name in line_names:
            dss.Lines.Name(name)

            Bus1str = dss.Lines.Bus1().split(".")
            Bus2str = dss.Lines.Bus2().split(".")
            Bus1 = Bus1str[0]
            Bus2 = Bus2str[0]

            if (len(Bus1str)-1)==0:
                fphase = [1,2,3]
            else:
                fphase = [int(i) for i in Bus1str[1:]]

            if (len(Bus2str)-1)==0:
                tphase = [1,2,3]
            else:
                tphase = [int(i) for i in Bus2str[1:]]

            fphases.append(fphase)
            tphases.append(tphase)

            buses.append([Bus1,Bus2])
            switches.append(dss.Lines.IsSwitch())
            nphases = dss.Lines.Phases()
            phases.append(nphases)

            ynm = inv((np.array(dss.Lines.RMatrix()).reshape(nphases,nphases)+1j*np.array(dss.Lines.XMatrix()).reshape(nphases,nphases))*dss.Lines.Length())
            ynms.append(ynm)

            ysh = np.array(dss.Lines.CMatrix()).reshape(nphases,nphases)*CapConversionFactor*dss.Lines.Length()
            yshs.append(ysh)

        if len(buses) == 0:
            lines_dict = {"name":[],"nphases":[],"fbus":[],"fphase":[],"tbus":[],"tphase":[],"ynm":[],"ysh":[],"isswitch":[]}
        else:
            lines_dict = {"name":line_names,"nphases":phases,"fbus":np.vstack(buses)[:,0],"fphase":fphases,"tbus":np.vstack(buses)[:,1],"tphase":tphases,"ynm":ynms,"ysh":yshs,"isswitch":switches}
        lines = pd.DataFrame(lines_dict)
        self.lines = lines

    def get_buses_from_dss(self):
        buses = []
        buses.append([bus for buses in self.transformers.buses for bus in buses])
        buses.append([bus_name for bus_name in self.lines.fbus])
        buses.append([bus_name for bus_name in self.lines.tbus])
        buses = list(np.unique(np.hstack(buses)))

        bus_phases = {buses[i]:[] for i in range(0,len(buses))}
        self.bus_phases = bus_phases

        for i in range(0,len(self.lines)):
            fbus = self.lines.fbus[i]
            fphase = self.lines.fphase[i]
            tbus = self.lines.tbus[i]
            tphase = self.lines.tphase[i]

            fbus_add = list(set(fphase) - set(self.bus_phases[fbus]))
            tbus_add = list(set(tphase) - set(self.bus_phases[tbus]))

            [self.bus_phases[fbus].append(i) for i in fbus_add]
            [self.bus_phases[tbus].append(i) for i in tbus_add]

        for i in range(0,len(self.transformers)):
            tbuses = self.transformers.buses[i]
            phases = self.transformers.phases[i]

            for j in range(0,len(tbuses)):
                bus = tbuses[j]
                for phase in phases[j]:
                    if phase == 0:
                        continue
                    else:
                        if phase not in bus_phases[bus]:
                            bus_phases[bus].append(phase)
        
        for i in range(0,len(self.capacitors)):
            bus = self.capacitors.bus[i]
            phase = self.capacitors.phases[i]

            bus_add = list(set(phase) - set(self.bus_phases[bus]))

            [self.bus_phases[bus].append(i) for i in bus_add]

        for i in range(0,len(self.loads)):
            bus = self.loads.bus[i]
            phase = self.loads.phases[i]

            bus_add = list(set(phase) - set(self.bus_phases[bus]))

            [self.bus_phases[bus].append(i) for i in bus_add]
        for bus in self.bus_phases:
            self.bus_phases[bus].sort()

        VAbase = 100e6
        kvbase = []
        for name in buses:
            dss.Circuit.SetActiveBus(name)
            kvbase.append(dss.Bus.kVBase()*1000)

        ### Voltage bases are line to neutral and VABase is single phase. This helps since I am solving for LN voltages. Also
        ### there is no need to divide by 3 to have the constant power loads in pu 
        busbase = {}
        angles = np.array([0,-2*pi/3,2*pi/3])
        for i in range(0,len(buses)):
            idx = np.array(self.bus_phases[buses[i]])-1
            busbase[buses[i]] = {"kVBase":kvbase[i],"VABase":VAbase,"ZBase":(kvbase[i]**2)/VAbase,"IBase":VAbase/(kvbase[i]),"Shift":angles[idx]}

        node_from_name = {}
        node_from_idx = {}
        counter = 0
        for i in range(0,len(buses)):
            for j in bus_phases[buses[i]]:
                node_from_name["{:}.{:}".format(buses[i],j)] = counter
                node_from_idx[counter] = "{:}.{:}".format(buses[i],j)
                counter+=1

        self.node_from_idx = node_from_idx
        self.node_from_name = node_from_name
        self.base = busbase
        DSScase.get_phase_shift(self)
        DSScase.get_delta_matrix(self)

    def get_phase_shift(self):
        G=nx.DiGraph()
        for i in range(0,len(self.lines)):
            G.add_edge(self.lines.fbus[i],self.lines.tbus[i])
        for i in range(0,len(self.transformers)):
            u, ind = np.unique(self.transformers.buses[i], return_index=True)
            buses = u[np.argsort(ind)]
            G.add_edge(buses[0],buses[1])
        for i in range(len(self.transformers)):
            buses = self.transformers.buses[i]
            conn = self.transformers.Conn[i]
            phases = self.transformers.phases[i]
            nwindings = self.transformers.windings[i]
            name = self.transformers.name[i]
            if nwindings == 2:
                if (conn[0]=="Delta") & (conn[1]=="Delta"):
                    # When both sides are delta or wye/wyeG, there is no phase shift
                    continue
                elif ((conn[0]=="Delta") | (conn[1]=="Delta")):
                    # When one side is delta and the other is wye, there is a 30 degree phase shift
                    # Often times, convention is Wye-delta adds -30 degree shift. Wye-delta +30 degree shift
                    # However, openDSS sets phase shift to -30 degree unless otherwise specified (By default, LeadLag = Lag)
                    print("Transformer {:} is delta-wye or wye-delta. LV side bus and buses downstream lag 30 degrees w.r.t. HV side".format(name))
                    downstream = [buses[1]] + [n for n in nx.traversal.bfs_tree(G, buses[1]) if n != buses[1]]
                    for bus in downstream:
                        self.base[bus]["Shift"] = self.base[bus]["Shift"] - pi/6 # -30 degree shift in radians
                        # self.base[bus]["Shift"] -= 0
            elif nwindings == 3: # This is hardcoded, and will work for the IEEE 8500. Revise if you want to generalize the calculation for nwindings>2
                pbus_phases = self.bus_phases[buses[0]]
                pbus_angle = self.base[buses[0]]["Shift"]
                pphase = phases[0][0]
                angle = pbus_angle[pbus_phases.index(pphase)]
                if (conn[0] == "WyeG"):
                    self.base[buses[1]]["Shift"] = np.array([angle,angle+pi])
                elif (conn[0] == "Delta"):
                    self.base[buses[1]]["Shift"] = np.array([angle+pi/6,angle+pi+pi/6])
                downstream = [buses[1]] + [n for n in nx.traversal.bfs_tree(G, buses[1]) if n != buses[1]]
                for bus in downstream:
                        self.base[bus]["Shift"] = self.base[buses[1]]["Shift"] # -30 degree shift in radians


    def get_delta_matrix(self):
        for bus in self.base:
            shift = self.base[bus]["Shift"]
            Δ3 = np.exp(1j*shift)
            self.base[bus]["Δ3"] = Δ3

    def get_constants(self):
        Nn = len(self.node_from_idx)
        Nb = len(self.bus_phases)
        Nbrl = sum(self.lines.nphases)

        Ω = (sqrt(3)/3)*np.array([[exp(-1j*pi/6),0,exp(1j*pi/6)],[exp(1j*pi/6),exp(-1j*pi/6),0],[0,exp(1j*pi/6),exp(-1j*pi/6)]])
        Π = np.array([[1,0,-1],[-1,1,0],[0,-1,1]])

        transformer_names = list(np.setdiff1d(self.transformers.name.values, self.regcontrols.transformer.values))
        regulator_names = list(self.regcontrols.transformer.values)

        tnphases = 0
        rnphases = 0
        for transformer in transformer_names:
            i = self.transformers[self.transformers.name == transformer].index.values[0]
            phases = deepcopy(self.transformers.phases[i])
            [phase.remove(0) for phase in phases if 0 in phase]
            tnphases += len(np.hstack(phases))
        for regulator in regulator_names:
            i = self.transformers[self.transformers.name == regulator].index.values[0]
            phases = deepcopy(self.transformers.phases[i])
            [phase.remove(0) for phase in phases if 0 in phase]
            rnphases += len(np.hstack(phases))
        Nbrt = tnphases
        Nbrr = rnphases
        nloadphases = 0
        for i in self.loads.phases:
            nloadphases += len(i)
        for i in self.capacitors.phases:
            nloadphases += len(i)

        self.Nbrl = Nbrl
        self.Nbrt = Nbrt
        self.Nbrr = Nbrr
        self.Nloadphases = nloadphases
        self.transformer_names = transformer_names
        self.regulator_names = regulator_names
        self.Nn = Nn
        self.Nb = Nb
        self.Ω = Ω
        self.Π = Π


    def create_load_model(self):
        Nn = self.Nn
        EloadS = csr_matrix((Nn,self.Nloadphases))
        counter = 0
        for i in range(0,len(self.loads)):
            bus = self.loads.bus[i]
            idxs = [self.node_from_name["{:}.{:}".format(bus,phase)] for phase in self.loads.phases[i]]
            for idx in idxs:
                EloadS[idx,counter] = 1
                counter += 1
        for i in range(0,len(self.capacitors)):
            bus = self.capacitors.bus[i]
            idxs = [self.node_from_name["{:}.{:}".format(bus,phase)] for phase in self.capacitors.phases[i]]
            for idx in idxs:
                EloadS[idx,counter] = 1
                counter += 1

        NloadYl = self.Nloadphases-sum([len(el) for el in self.loads[self.loads.model==1].phases])
        EloadYl = csr_matrix((Nn,NloadYl))
        counter = 0
        for i in range(0,len(self.loads)):
            if (self.loads.model[i] == 1):
                continue
            else:
                bus = self.loads.bus[i]
                idxs = [self.node_from_name["{:}.{:}".format(bus,phase)] for phase in self.loads.phases[i]]
                for idx in idxs:
                    EloadYl[idx,counter] = 1
                    counter += 1
        for i in range(0,len(self.capacitors)):
            bus = self.capacitors.bus[i]
            idxs = [self.node_from_name["{:}.{:}".format(bus,phase)] for phase in self.capacitors.phases[i]]
            for idx in idxs:
                EloadYl[idx,counter] = 1
                counter += 1  


        self.EloadYl = EloadYl
        self.EloadS = EloadS


    def process_openDSS_solution(self):
        dssvm = self.bus_phases.copy()
        dssva = self.bus_phases.copy()

        for name in dssvm.keys():
            dss.Circuit.SetActiveBus(name)

            x=dss.Bus.PuVoltage()
            a=dss.Bus.VMagAngle()

            y=[]
            ya=[]

            for i in range(0,len(x),2):
                y.append(dss.CmathLib.cabs(x[i],x[i+1]))
                ya.append(a[i+1])
            if len(dssvm[name])!=len(y):
                print("Bus {:}, phase mismatch".format(name))
            dssvm[name] = y
            dssva[name] = ya

        tap = {}
        for name in dss.RegControls.AllNames():
            dss.RegControls.Name(name)
            tap[name] = 1+dss.RegControls.TapNumber()*0.00625


        if not hasattr(self,"results"):
            self.results = {}
        self.results["openDSS"] = {}
        self.results["openDSS"]["vm"] = dssvm
        self.results["openDSS"]["va"] = dssva
        self.results["openDSS"]["regtap"] = tap



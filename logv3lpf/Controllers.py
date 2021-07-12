from logv3lpf.math_functions import polar_to_cartesian

def RegControl(Vreg,band,Vn,Nt,In,Ct,R,X):
    # vr = (Vn/Nt) - ((R+1j*X)/Ct)*In
    vr = abs((Vn/Nt)) - abs(((R+1j*X)/Ct)*In)
    
    # if abs(vr) >= Vreg + band/2:
    if vr >= Vreg + band/2:
        return -1
    # elif abs(vr) <= Vreg - band/2:
    elif vr <= Vreg - band/2:
        return 1
    else:
        return 0
    # return min(abs(Vreg - band/2 - abs(vr)),abs(abs(vr) - (Vreg + band/2)))/0.75

def apply_regcontrols(case):
    tap_changes = []
    for i in range(len(case.regcontrols)):
        # Get regulator constants
        name = case.regcontrols.name[i]
        idx = case.transformers[case.transformers.name == name].index[0]
        Vreg,band,Nt = case.regcontrols.vreg[i], case.regcontrols.band[i], case.regcontrols.nt[i]
        Ct,R,X = case.regcontrols.ct[i], case.regcontrols.R[i], case.regcontrols.X[i]

        # Get transformer info
        buses = case.transformers.buses[idx]
        phases = case.transformers.phases[idx]
        monitoring_bus,tap_bus = buses[0],buses[1]
        taps = case.transformers.taps[idx]
        tap_phase = phases[1][0]
        idx_monitor = case.bus_phases[monitoring_bus].index(tap_phase)
        idx_tap = case.bus_phases[tap_bus].index(tap_phase)

        # Get bases
        kVBase = case.base[monitoring_bus]["kVBase"]
        ZBase = case.base[tap_bus]["ZBase"]


        vn,van = case.results["logv3lpf"]["vm"][monitoring_bus][idx_monitor], case.results["logv3lpf"]["va"][monitoring_bus][idx_monitor]
        vm,vam = case.results["logv3lpf"]["vm"][tap_bus][idx_tap], case.results["logv3lpf"]["va"][tap_bus][idx_tap]
        Vn = polar_to_cartesian(vn,van)*kVBase
        Vm = polar_to_cartesian(vm,vam)*kVBase

        In = case.get_transformer_current(name,False)[0]
        tap_change = RegControl(Vreg,band,Vn,Nt,In,Ct,R,X)

        case.transformers.at[idx,"taps"] = [1,taps[1] + tap_change*0.00625]

        tap_changes.append(tap_change)

    return tap_changes

def apply_controls(case):
    taps = apply_regcontrols(case)

    return taps

    





    




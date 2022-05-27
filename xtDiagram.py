import pandas as pa
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

import ShockCondensedMatter

class TP

def main():
    metal_data = pa.read_excel('/Users/dihiaidrici/Desktop/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='MetalData', header=1, index_col=0)
    speed_data = pa.read_excel('/Users/dihiaidrici/Desktop/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='Speed', header=None, index_col=0)
    position_initial = pa.read_excel('/Users/dihiaidrici/Desktop/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='Position')

    # remember that indexing begins at 0, Rows: Al[0], Mg[1], Ti[2], Zr[3], SS[4], Al2O3[5]
    Al = metal_data.iloc[0]
    Mg = metal_data.iloc[1]
    Ti = metal_data.iloc[2]
    Zr = metal_data.iloc[3]
    SS = metal_data.iloc[4]
    Al2O3 = metal_data.iloc[5]

    Zr_Al2O3 = ShockCondensedMatter(Zr, Al2O3)

    # impactS = speed_data.iloc[0]
    # downA = velocity_dependence(Al, Al2O3, impactS)
    # downM = velocity_dependence(Mg, Al2O3, impactS)
    # downT = velocity_dependence(Ti, Al2O3, impactS)
    # downZ = velocity_dependence(Zr, Al2O3, impactS)
    # downAS = velocity_dependence(Al, SS, impactS)
    # downMS = velocity_dependence(Mg, SS, impactS)
    # downTS = velocity_dependence(Ti, SS, impactS)
    # downZS = velocity_dependence(Zr, SS, impactS)
    # figures(impactS, downA, downM, downT, downZ, downAS, downMS, downTS, downZS)

    Vip = 1 #km/s
    ViT = 0

    down2 = impedance_matching(Zr, Al2O3, Vip, ViT)
    downZ2 = down2[0:5, 0]  # V2, PfL, TfL, rhofL, VsL
    downT2 = down2[5:9, 0]  # PfR, rhofR, VsR
    [x1, t1, V1, x2, t2, V2] = [float(position_initial.xp), 0, Vip, 0, 0, downZ2[4]]
    [xf, tf] = xt_diagram(x1, t1, V1, x2, t2, V2)
    [x1T, t1T, V1T, x1S, t1S, V1S] = [0, 0, downT2[2], float(position_initial.xT), 0, 0]
    [xf_TS, tf_TS] = xt_diagram(x1T, t1T, V1T, x1S, t1S, V1S)

    plt.xlim(-4*10**-6, 15*10**-6), plt.ylim(0, 3*10**-6)
    plt.plot([x1, xf], [t1, tf])
    plt.plot([x2, xf], [t2, tf])
    plt.plot([x1T, xf_TS], [t1T, tf_TS])
    plt.plot([x1S, xf_TS], [t1S, tf_TS])
    plt.show()


    # downA = targetProjectile(Al, Al2O3, impactS)
    # downM = targetProjectile(Mg, Al2O3, impactS)
    # downT = targetProjectile(Ti, Al2O3, impactS)


def xt_diagram(x1, t1, V1, x2, t2, V2):
    tf = ((x1 - x2) + V2*t2 - V1*t1)/(V2-V1)
    xf = x1 + V1*(tf-t1)
    return xf,tf

def impedance_matching(L, R, ViL, ViR):
    down = np.ones((9, 1))
    V0 = 0.1
    V2 = fsolve(pressure_matching, V0, args=(L, R, ViL, ViR), xtol=1e-6)
    down[0, 0] = float(V2)
    down[1:9, 0] = downstream_prop(L, R, V2, ViL, ViR)
    return down


def pressure_matching(V, L, R, ViL, ViR):
    # V: interface speed, L: left material, R: Right material,
    PfL = L.density*(ViL-V)*(L.C01 + (ViL-V)*L.s1) + L.Pi
    PfR = R.density*(V-ViR)*(R.C01 + (V-ViR)*R.s1) + R.Pi
    return PfL - PfR


def downstream_prop(L, R, V, ViL, ViR):
    PfL = float(L.density*(ViL-V)*(L.C01 + (ViL-V)*L.s1) + L.Pi)
    VsL = float(ViL - (L.C01 + (ViL-V)*L.s1))
    rhofL = float(L.density*(VsL-ViL)/(VsL-V))
    delta_energyL = float((1/L.density - 1/rhofL)*(PfL + L.Pi)/2)
    # delta_energy = (PfL*V - L.Pi*ViL)/(L.density*(VsL-ViL)) - (V**2 - ViL**2)/2
    TfL = float(delta_energyL*10**6/L.specHC + L.Ti)

    PfR = float(R.density*(V-ViR)*(R.C01 + (V-ViR)*R.s1) + R.Pi)
    VsR = float(ViR + R.C01 + R.s1*(V-ViR))
    rhofR = float(R.density*(VsR-ViR)/(VsR-V))
    delta_energyR = float((1/R.density - 1/rhofR)*(PfR + R.Pi)/2)
    TfR = float((delta_energyR*10**6/R.specHC)+ R.Ti)

    if PfL == PfR:
        return PfL, TfL, rhofL, VsL, PfR, TfR, rhofR, VsR
    else:
        print("There is a problem")


def velocity_dependence(metal, target, impactS):
    down = np.ones((5, len(impactS)-1))
    for x in range(1, len(impactS)):
        ViL = impactS.iloc[x]  # km/s
        ViR = impactS.iloc[0]
        V0 = 0.1
        V2 = fsolve(impedance_matching, V0, args=(metal, target, ViL, ViR), xtol=1e-6)
        down[0, x-1] = float(V2)
        # Pressure, Temperature, density, shock speed
        down[1:5, x-1] = downstream_prop(metal, target, V2, ViL, ViR)
    return down


def figures(impactS, downA, downM, downT, downZ, downAS, downMS, downTS, downZS):
    plt.figure()
    plt.subplot(211)
    plt.plot(impactS[1:len(impactS)], downA[2, :], 'r--', impactS[1:len(impactS)], downM[2, :], 'b--', impactS[1:len(impactS)], downT[2, :], 'g--', impactS[1:len(impactS)], downZ[2, :], 'm--')
    plt.ylabel('Temperature rise, K')
    plt.xlabel('Impact Speed, km/s')
    plt.subplot(212)
    plt.plot(impactS[1:len(impactS)], downA[1, :], 'r--', impactS[1:len(impactS)], downM[1, :], 'b--', impactS[1:len(impactS)], downT[1, :], 'g--', impactS[1:len(impactS)], downZ[1, :], 'm--')
    plt.ylabel('Pressure rise, GPa')
    plt.xlabel('Impact Speed, km/s')
    plt.savefig('/Users/dihiaidrici/Desktop/EverythingPresentation/shockCondTP', format='eps')
    plt.show()

    plt.figure()
    plt.subplot(211)
    plt.plot(impactS[1:len(impactS)], downA[2, :], 'r--', impactS[1:len(impactS)], downAS[2, :], 'b--')
    plt.ylabel('Temperature rise, K')
    plt.xlabel('Impact Speed, km/s')
    plt.subplot(212)
    plt.plot(impactS[1:len(impactS)], downA[1, :], 'r--', impactS[1:len(impactS)], downAS[1, :], 'b--')
    plt.ylabel('Pressure rise, GPa')
    plt.xlabel('Impact Speed, km/s')
    #plt.savefig('/Users/dihiaidrici/Desktop/EverythingPresentation/shockCondTP', format='eps')
    plt.show()


if __name__ == '__main__':
    main()

import numpy as np
from scipy.optimize import fsolve
from sympy import *
import matplotlib.pyplot as plt


class ShockCMatter:
    def __init__(self, metalData, targetData, steelData, Vip, ViT, xProj, xTar, xTarb, xSteel):
        self.M_i = metalData  # M for metal
        self.T_i = targetData  # T for target
        self.S_i = steelData  # S for steel
        self.xp = xProj
        self.xT = xTar
        self.xTS = xTarb
        self.xS = xSteel
        self.Vip = Vip
        self.ViT = ViT

    def shock_interactions(self):
        # impedance_matching(L, rhoL, PL, TL,ViL, R, rhoR, PR, TR, ViR) START WITH LEFT and then RIGHT
        # 1.Impedance matching projectile-target impact
        state1 = ShockCMatter.impedance_matching(self.M_i, self.M_i.density, self.M_i.Pi, self.M_i.Ti, self.Vip,
                                                 self.T_i, self.T_i.density, self.T_i.Pi, self.T_i.Ti, self.ViT)
        V1 = float(state1[0])  # interface speed
        state1p = state1[1:6]  # rhofL (g/cm^3),PfL (GPa), TfL (K), delta_energyL (MJ/kg) , VsL (km/s)
        state1T = state1[6:11]  # rhofR, PfR, TfR, delta_energyR, VsR

        # Obtain Projectile temperature using Cp = Cp(T). For internal energy we need Cv, but NIST only has Cp, but for solids Cv ~ Cp.
        ti = self.M_i.Ti/1000
        t0 = np.array([state1p[2]/1000])  # guess
        t = fsolve(ShockCMatter.changingCp_temp, t0, args=(ti, self.M_i.A, self.M_i.B, self.M_i.C, self.M_i.D, self.M_i.E,
                                                           self.M_i.F, self.M_i.H, state1p[3]*10**6, self.M_i.M), xtol=1e-6)[0]
        Tfp = t * 1000  # K
        tt = np.linspace(0.01, 5, 100)
        y = state1p[3]*10**6 - (self.M_i.A*(tt-ti) + self.M_i.B*(tt**2 - ti**2)/2 + self.M_i.C*(tt**3 - ti**3)/3 +
                                self.M_i.D*(tt**4 - ti**4)/4 - self.M_i.E*(1/tt - 1/ti))*1000/self.M_i.M
        plt.plot(tt,y, 'r')
        plt.show()
        print(state1p)
        print(Tfp)

        # 2.Shock in projectile - back of the projectile (Vsp_1 and Vi)
        state2 = ShockCMatter.shock_exit(self.M_i, state1p[0], state1p[1], state1p[2], V1)
        V2 = float(state2[0])
        state2p = state2[1:4]  # rhofR, PfR, VsR

        # 3.Impedance matching Alumina reflected shock - steel target (Idealization) (VsT and 0)
        state3 = ShockCMatter.impedance_matching(self.T_i, state1T[0], state1T[1], state1T[2], V1,
                                                 self.S_i, self.S_i.density, self.S_i.Pi, self.S_i.Ti, self.ViT)
        V3 = float(state3[0])
        state3T = state3[1:6]  # Alumina
        state3S = state3[6:11]  # Steel

        ShockCMatter.position_graph(self, V1, state1p, state1T, V2, state2p, V3, state3T, state3S)


    @staticmethod
    def impedance_matching(L, rhoL, PL, TL, ViL, R, rhoR, PR, TR, ViR):
        newState = np.ones(11)  # downstreamProperties
        V0 = np.array([0])  # guess
        V2 = fsolve(ShockCMatter.pressure_matching, V0, args=(L, rhoL, PL, ViL, R, rhoR, PR, ViR), xtol=1e-6)
        newState[0] = V2[0]
        newState[1:11] = ShockCMatter.downstream_prop(L, rhoL, PL, TL, ViL, R, rhoR, PR, TR, ViR, V2[0])
        return newState

    @staticmethod
    def pressure_matching(V, L, rhoL, PL, ViL, R, rhoR, PR, ViR):
        # V: interface speed, L: left material, R: Right material,
        PfL = rhoL*(ViL - V)*(L.C01 + (ViL - V)*L.s1) + PL  # GPa
        PfR = rhoR*(V - ViR)*(R.C01 + (V - ViR)*R.s1) + PR  # GPa
        return PfL - PfR

    @staticmethod
    def changingCp_temp(t, ti, A, B, C, D, E, F, H, delta_E, M):
        LHS = delta_E
        RHS = (A*(t-ti) + B*(t**2 - ti**2)/2 + C*(t**3 - ti**3)/3 + D*(t**4 - ti**4)/4 - E*(1/t - 1/ti))*1000/M
        # RHS = (A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H)*1000/M  # Nist equation, F and H replace ti
        return LHS - RHS

    @staticmethod
    def downstream_prop(L, rhoL, PL, TL, ViL, R, rhoR, PR, TR, ViR, V):
        PfL = rhoL*(ViL - V)*(L.C01 + (ViL - V)*L.s1) + PL   # GPa
        VsL = ViL - (L.C01 + L.s1*(ViL - V))  # km/s
        rhofL = rhoL*(VsL - ViL)/(VsL - V)  # g/cm^3
        delta_energyL = (1/rhoL - 1/rhofL)*(PfL + PL)/2  # MJ/kg
        TfL = (delta_energyL*10**6/L.specHC) + TL  # K

        PfR = rhoR*(V - ViR)*(R.C01 + (V - ViR)*R.s1) + PR
        VsR = ViR + R.C01 + R.s1*(V - ViR)
        rhofR = rhoR*(VsR - ViR)/(VsR - V)
        delta_energyR = (1/rhoR - 1/rhofR)*(PfR + PR)/2
        TfR = (delta_energyR*10**6/R.specHC) + TR
        if int(PfL) == int(PfR):
            return np.array([rhofL, PfL, TfL, delta_energyL, VsL, rhofR, PfR, TfR, delta_energyR, VsR])
        else:
            print("There is a problem")

    @staticmethod
    def shock_exit(R, rho1, P1, T1, V1):
        P2 = R.Pi
        V = symbols('V')
        b = solve((P2 - P1) - rho1*(V - V1)*(R.C01 + R.s1*(V - V1)), V)
        V2 = b[1]
        PfR = float(rho1*(V2 - V1)*(R.C01 + (V2 - V1)*R.s1) + P1)  # GPa
        VsR = float(V1 + R.C01 + R.s1*(V2 - V1))
        rhofR = float(rho1*(VsR - V1)/(VsR - V2))   # g/cm^3
        newState = [V2, rhofR, PfR, VsR]
        return newState

    @staticmethod
    def xt_diagram(x1, t1, V1, x2, t2, V2):
        tf = ((x1 - x2) + V2*t2 - V1*t1)/(V2-V1)
        xf = x1 + V1*(tf-t1)
        return xf, tf

    def position_graph(self, V1, state1p, state1T, V2, state2p, V3, state3T, state3S):
        # Graphs and positions
        # 1. back of the projectile (b), projectile shock (s)
        x1b, t1b, V1b, x1_sp, t1_sp, V1_sp = self.xp, 0, self.Vip, self.xT, 0, state1p[4]
        [x1, t1] = ShockCMatter.xt_diagram(x1b, t1b, V1b, x1_sp, t1_sp, V1_sp)
        # 2. Alumina shock and back of target (where the steel begins)
        x1_sT, t1_sT, V1_sT, x1S, t1S, V1S = self.xT, 0, state1T[4], self.xTS, 0, 0
        [x2, t2] = ShockCMatter.xt_diagram(x1_sT, t1_sT, V1_sT, x1S, t1S, V1S)
        # 3. Alumina reflected shock and interface (I1)
        x3_I1, t3_I1, V3_I1, x3_sT, t3_sT, V3_sT = self.xT, 0, V1, x2, t2, state3T[4]
        [x3, t3] = ShockCMatter.xt_diagram(x3_I1, t3_I1, V3_I1, x3_sT, t3_sT, V3_sT)
        # 4. Reflected expansion wave projectile (E) and interface (I1)
        x4_E, t4_E, V4_E, x4_I1, t4_I1, V4_I1 = x1, t1, state2p[2], 0, 0, V1
        [x4, t4] = ShockCMatter.xt_diagram(x4_E, t4_E, V4_E, x4_I1, t4_I1, V4_I1)
        # 5. back of the projectile
        x5 = x1 + V2*(t4-t1)

        # plt.xlim(x1b, x1S),
        # plt.ylim(0, 3*10**-6)
        plt.plot([x1_sp, x1], [t1_sp, t1])
        plt.plot([x1b, x1], [t1b, t1])
        plt.plot([x1_sT, x2], [t1_sT, t2])
        plt.plot([x1S, x2], [t1S, t2])
        plt.plot([x3_sT, x3], [t3_sT, t3])
        plt.plot([x4_E, x4], [t4_E, t4])
        plt.plot([x4_I1, x4], [t4_I1, t4])
        plt.plot([x1, x5], [t1, t3])
        plt.xlabel('Position, km')
        plt.ylabel('time, s')
        plt.show()

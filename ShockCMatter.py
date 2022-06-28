import numpy as np
from scipy.optimize import fsolve
from sympy import *
import math
import matplotlib.pyplot as plt
from tabulate import tabulate


class ShockCMatter:
    def __init__(self, metalData, targetData, steelData, Vip, ViT, xProj, xTar, xTarb, xSteel):
        # initialize initial state
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
        # impedance_matching(L, rhoL, PL, TL,ViL, R, rhoR, PR, TR, ViR) start with LEFT and then RIGHT
        # 1.Impedance matching projectile-target impact
        state1 = ShockCMatter.impedance_matching(self.M_i, self.M_i.density, self.M_i.Pi, self.M_i.Ti, self.Vip,
                                                 self.T_i, self.T_i.density, self.T_i.Pi, self.T_i.Ti, self.ViT)
        V1 = float(state1[0])  # interface speed
        state1p = state1[1:6]  # rhofL (g/cm^3),PfL (GPa), TfL (K), delta_energyL (MJ/kg) , VsL (km/s)
        state1T = state1[6:11]  # rhofR, PfR, TfR, delta_energyR, VsR

        # Projectile Temperature with Cp = Cp(T). Internal energy requires Cv, NIST only has Cp: For solids Cv ~ Cp.
        ti = self.M_i.Ti/1000  # NIST equation requires temperature/1000
        t0 = np.array([state1p[2]/1000])  # guess
        t = fsolve(ShockCMatter.changingCp_temp, t0, args=(ti, self.M_i.A, self.M_i.B, self.M_i.C, self.M_i.D, self.M_i.E,
                                                           self.M_i.F, self.M_i.H, state1p[3]*10**6, self.M_i.M), xtol=1e-6)[0]
        Tfp = t*1000  # K Final temperature
        # Plot the deltaE - Nist equation to validate the temperature.
        tt = np.linspace(0.29815, 2, 50)
        y = state1p[3]*10**6 - (self.M_i.A*(tt-ti) + self.M_i.B*(tt**2 - ti**2)/2 + self.M_i.C*(tt**3 - ti**3)/3 +
                                self.M_i.D*(tt**4 - ti**4)/4 - self.M_i.E*(1/tt - 1/ti))*1000/self.M_i.M
        plt.plot(tt, y, 'r')
        plt.show()

        # 2.Expansion wave back of the projectile
        V2, Vs_lead = ShockCMatter.expansion_projectile(self.M_i, state1p, V1, 'RRW')
        Vs_trail = self.M_i.CL
        print(V2, Vs_lead)
        delta_E_Raleigh, delta_E_Hugoniot, delta_E = ShockCMatter.delta_energy(self.M_i, self.Vip, state1p)
        t_shock_expansion = fsolve(ShockCMatter.changingCp_temp, t0, args=(ti, self.M_i.A, self.M_i.B, self.M_i.C, self.M_i.D, self.M_i.E,
                                                           self.M_i.F, self.M_i.H, delta_E*10**6, self.M_i.M), xtol=1e-6)[0]
        T_shock_expansion = t_shock_expansion*1000  # K Final temperature
        print(delta_E_Raleigh, delta_E_Hugoniot, delta_E)
        print(Tf)

        # 3.Impedance matching Alumina reflected shock - steel target (Idealization) (VsT and 0)
        # state3 = ShockCMatter.impedance_matching(self.T_i, state1T[0], state1T[1], state1T[2], V1,
        #                                         self.S_i, self.S_i.density, self.S_i.Pi, self.S_i.Ti, self.ViT)
        # V3 = float(state3[0])
        # state3T = state3[1:6]  # Alumina
        # state3S = state3[6:11]  # Steel

        # ShockCMatter.position_graph(self, V1, state1p, state1T, V2, state2p, V3, state3T, state3S)

        # Plot the P-u diagram
        V = np.linspace(-0.5, 1, 20)
        PLp = self.M_i.density*(self.Vip - V)*(self.M_i.C01 + (self.Vip - V)*self.M_i.s1) + self.M_i.Pi  # GPa LRW
        PRp = self.M_i.density*(V - V2)*(self.M_i.C01 + (V - V2)*self.M_i.s1) + self.M_i.Pi  # GPa RRW
        PRT = self.T_i.density*(V - self.ViT)*(self.T_i.C01 + (V - self.ViT)*self.T_i.s1) + self.T_i.Pi  # GPa
        plt.plot(V, PLp, 'r')
        plt.plot(V, PRp, 'b')
        plt.plot(V, PRT, 'g')
        plt.show()

    @staticmethod
    def expansion_projectile(state_i, state_f, Vf, side):
        V0 = 0
        V2 = fsolve(ShockCMatter.expansion_wave, V0, args=(state_i, state_f, Vf, side), xtol=1e-6)
        # first Wave speed traveling through final shocked state of the projectile.
        Vs_lagrange_lead = state_i.C01 + 2*state_i.s1*(Vf - V2)
        Vs_lead = Vs_lagrange_lead + Vf
        return V2, Vs_lead

    @staticmethod
    def expansion_wave(V0, state_i, state_f, Vf, side):
        # Find V0 of the reflection
        P0 = state_i.Pi
        if side == 'LRW':  # LRW: P  = rho0*(V0 - V)*(C0 + s*(V0 - V)) + P0
            return state_i.density*(V0 - Vf)*(state_i.C01 + state_i.s1*(V0 - Vf)) + P0 - state_f[1]
        elif side == 'RRW':  # RRW: P  = rho0*(V - V0)*(C0 + s*(V - V0)) + P0
            return state_i.density*(Vf - V0)*(state_i.C01 + state_i.s1*(Vf - V0)) + P0 - state_f[1]

    @staticmethod
    def impedance_matching(L, rhoL, PL, TL, ViL, R, rhoR, PR, TR, ViR):
        newState = np.ones(11)  # downstreamProperties
        V0 = np.array([0])  # guess
        V_interface = fsolve(ShockCMatter.pressure_matching, V0, args=(L, rhoL, PL, ViL, R, rhoR, PR, ViR), xtol=1e-6)
        newState[0] = V_interface[0]
        newState[1:11] = ShockCMatter.downstream_prop(L, rhoL, PL, TL, ViL, R, rhoR, PR, TR, ViR, V_interface[0])
        return newState

    @staticmethod
    def pressure_matching(V, L, rhoL, PL, ViL, R, rhoR, PR, ViR):
        # V: interface speed, L: left material, R: Right material,
        PfL = rhoL*(ViL - V)*(L.C01 + (ViL - V)*L.s1) + PL  # GPa
        PfR = rhoR*(V - ViR)*(R.C01 + (V - ViR)*R.s1) + PR  # GPa
        return PfL - PfR

    @staticmethod
    def downstream_prop(L, rhoL, PL, TL, ViL, R, rhoR, PR, TR, ViR, V):
        PfL = rhoL*(ViL - V)*(L.C01 + (ViL - V)*L.s1) + PL   # GPa
        VsL = ViL - (L.C01 + L.s1*(ViL - V))  # km/s
        rhofL = rhoL*(VsL - ViL)/(VsL - V)  # g/cm^3
        delta_energyL = (1/rhoL - 1/rhofL)*(PfL + PL)/2  # MJ/kg
        TfL = (delta_energyL*10**6/L.specHC) + TL  #K

        PfR = rhoR*(V - ViR)*(R.C01 + (V - ViR)*R.s1) + PR
        VsR = ViR + (R.C01 + R.s1*(V - ViR))
        rhofR = rhoR*(VsR - ViR)/(VsR - V)
        delta_energyR = (1/rhoR - 1/rhofR)*(PfR + PR)/2
        TfR = (delta_energyR*10**6/R.specHC) + TR
        if int(PfL) == int(PfR):
            return np.array([rhofL, PfL, TfL, delta_energyL, VsL, rhofR, PfR, TfR, delta_energyR, VsR])
        else:
            print("There is a problem")

    @staticmethod
    def changingCp_temp(t, ti, A, B, C, D, E, F, H, delta_E, M):
        LHS = delta_E  #J/kg
        # Multiply by 1000 because of NIST equations, H = KJ/mol, M = kg/mol
        RHS = (A*(t-ti) + B*(t**2 - ti**2)/2 + C*(t**3 - ti**3)/3 + D*(t**4 - ti**4)/4 - E*(1/t - 1/ti))*1000/M
        # RHS = (A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H)*1000/M  # Nist equation, F and H replace ti
        return LHS - RHS

    @staticmethod
    def delta_energy(state_i, Vi, state_f):
        sv0 = 1/state_i.density
        sv = 1/state_f[0]
        Vs = state_f[4]   # It is the shock compression speed because Rayleigh line.
        # Delta_E_Raleigh is the same as the delta energy computed in downstream prop
        delta_E_Raleigh = -(Vs - Vi)**2*(sv - sv0)/sv0 + (Vs - Vi)**2*(sv**2 - sv0**2)/(2*sv0**2) - state_i.Pi*(sv - sv0)
        delta_E_Hugoniot = (state_i.C01/state_i.s1)**2*(math.log((sv0 - state_i.s1*(sv0 - sv))/sv0) +
                                                        state_i.s1*(sv0 - sv)/(sv0 - state_i.s1*(sv0 - sv))) - state_i.Pi*(sv - sv0)
        delta_E = delta_E_Raleigh - delta_E_Hugoniot
        return delta_E_Raleigh, delta_E_Hugoniot, delta_E

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
        plt.ylim(0, 3*10**-6)
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

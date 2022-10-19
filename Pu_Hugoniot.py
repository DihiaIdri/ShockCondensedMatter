PfL = rhoL * (ViL - V) * (L.C01 + (ViL - V) * L.s1) + PL  # GPa


PfR = rhoR * (V - ViR) * (R.C01 + (V - ViR) * R.s1) + PR

# state2 = ShockCMatter.shock_exit(self.M_i, state1p[0], state1p[1], state1p[2], V1)
# V2 = float(state2[0])
# state2p = state2[1:4]  # rhofR, PfR, VsR

# 3.Impedance matching Alumina reflected shock - steel target (Idealization) (VsT and 0)
# state3 = ShockCMatter.impedance_matching(self.T_i, state1T[0], state1T[1], state1T[2], V1,
#                                         self.S_i, self.S_i.density, self.S_i.Pi, self.S_i.Ti, self.ViT)
# V3 = float(state3[0])
# state3T = state3[1:6]  # Alumina
# state3S = state3[6:11]  # Steel

# ShockCMatter.position_graph(self, V1, state1p, state1T, V2, state2p, V3, state3T, state3S)

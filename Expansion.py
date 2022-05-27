from sympy import *

class Expansion:

    @staticmethod
    def expansion(char, metal, state1,Ui):
        Pf = metal.Pi
        rhoi = state1[0]
        Pi = state1[1]
        Ti = state1[2]
        U = symbols('U')
        if char == 'L':
            U = solve((Pf - Pi - rhoi*(Ui - U)*(metal.C01 + (Ui - U)*metal.s1)), U)
            U = U[1][0]
            Us = Ui - (metal.C01 + metal.s1*(Ui - U))
        elif char == 'R':
            U = solve((Pf - Pi - rhoi*(U - Ui)*(metal.C01 + (U - Ui)*metal.s1)), U)
            U = U[1][0]
            Us = float(Ui + metal.C01 + metal.s1*(U - Ui))

        rhof = float(rhoi*(Us - Ui)/(Us - U))
        # Incredibly wrong
        # delta_energy = float((1/rhoi - 1/rhof)*(Pf + Pi)/2)
        # Tf = float((delta_energy*10**6/metal.specHC) + Ti)
        print(np.array([rhof, Pf, Tf, Us, U]))
        return np.array([rhof, Pf, Tf, Us, U])

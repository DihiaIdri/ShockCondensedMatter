import pandas as pa
from ShockCMatter import ShockCMatter
from Expansion import Expansion


def main():
    metal_data = pa.read_excel('/Users/dihiaidrici/Desktop/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='MetalData', header=1, index_col=0)
    position_initial = pa.read_excel('/Users/dihiaidrici/Desktop/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='Position')

    # remember that indexing begins at 0, Rows: Al[0], Mg[1], Ti[2], Zr[3], SS[4], Al2O3[5]
    Al = metal_data.iloc[0]
    Mg = metal_data.iloc[1]
    Ti = metal_data.iloc[2]
    Zr = metal_data.iloc[3]
    SS = metal_data.iloc[4]
    Al2O3 = metal_data.iloc[5]

    Vip = 1  # km/s
    ViT = 0
    U1p = 0
    xProj = position_initial.xp
    xTar = position_initial.xT
    xTarb = position_initial.xTS
    xSteel = position_initial.xS

    # Al_Al2O3 = ShockCMatter(Al, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel)
    # Al_state1 = Al_Al2O3.shock_interactions()
    # Mg_Al2O3 = ShockCMatter(Mg, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel)
    # Mg_state1 = Mg_Al2O3.shock_interactions()
    # Ti_Al2O3 = ShockCMatter(Ti, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel)
    # Ti_state1 =  Ti_Al2O3.shock_interactions()
    Zr_Al2O3 = ShockCMatter(Zr, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel)
    Zr_Al2O3.shock_interactions()

    # Here I am looking at the temperature, speed of the fragments once expansion occurs.
    # my calculations are either bad or I am missing something from the physics.
    # Zr_Al2O3_jetting = Expansion.expansion('R', Zr, Zr_state1, U1p)
    # print(Zr_Al2O3_jetting)


if __name__ == '__main__':
    main()
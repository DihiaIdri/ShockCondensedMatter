import numpy as np
import pandas as pa
from ShockCMatter import ShockCMatter
from graphForComparison import graphForComparison
from tabulate import tabulate


def main():
    # Extracting data from excel Sheet
    metal_data = pa.read_excel('/Users/dihiaidrici/Desktop/SecondPaper/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='MetalData', header=1, index_col=0)
    initial_xv = pa.read_excel('/Users/dihiaidrici/Desktop/SecondPaper/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='initial_xv', header=1)
    speedRange = pa.read_excel('/Users/dihiaidrici/Desktop/SecondPaper/EverythingSecondPaper/ShockInMetalProperties.xlsx', sheet_name='SpeedRange', header=0, index_col=0)

    # remember that indexing begins at 0, Rows: Al[0], Mg[1], Ti[2], Zr[3], SS[4], Al2O3[5]
    Al = metal_data.iloc[0]
    Mg = metal_data.iloc[1]
    Ti = metal_data.iloc[2]
    Zr = metal_data.iloc[3]
    SS = metal_data.iloc[4]
    Al2O3 = metal_data.iloc[5]
    ZrNbCuNiAl = metal_data.iloc[6]

    data_xv = initial_xv.iloc[0]
    Vip = data_xv.Vip  # km/s
    ViT = data_xv.ViT
    xProj = data_xv.xp
    xTar = data_xv.xT
    xTarb = data_xv.xTS
    xSteel = data_xv.xS

    # Calculations
    Al_Al2O3 = ShockCMatter(Al, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel, 'Al')
    Al_state1 = Al_Al2O3.shock_interactions()
    # Mg_Al2O3 = ShockCMatter(Mg, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel, 'Mg')
    # Mg_state1 = Mg_Al2O3.shock_interactions()
    # Ti_Al2O3 = ShockCMatter(Ti, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel, 'Ti')
    # Ti_state1 =  Ti_Al2O3.shock_interactions()
    Zr_Al2O3 = ShockCMatter(Zr, Al2O3, SS, Vip, ViT, xProj, xTar, xTarb, xSteel, 'Zr')
    Zr_state1 = Zr_Al2O3.shock_interactions()
    # ZrNbCuNiAl_SS = ShockCMatter(ZrNbCuNiAl, SS, SS, Vip, ViT, xProj, xTar, xTarb, xSteel)
    # ZrNbCuNiAl_state1 = ZrNbCuNiAl_SS.shock_interactions()

    # Tabulate Data
    header1 = ["Metal", "u, km/s", "Us, km/s", "\u03C1, g/cm\u00b3", "P, GPa", "\u0394e Ra, MJ/kg", "T, K",
               "\u0394e Hu, MJ/kg", "T, K"]
    data = [Al_state1, Zr_state1]
    print(tabulate(data, headers=header1, tablefmt="heavy_outline"))

    # Plot comparison Graphs
    T_Al = np.zeros(len(speedRange.columns))
    P_Al = np.zeros(len(speedRange.columns))
    T_Zr = np.zeros(len(speedRange.columns))
    P_Zr = np.zeros(len(speedRange.columns))
    for i in range(len(speedRange.columns)):
        u1 = speedRange.columns[i]
        Al_Al2O3 = ShockCMatter(Al, Al2O3, SS, u1, ViT, xProj, xTar, xTarb, xSteel, 'Al')
        Al_state1 = Al_Al2O3.shock_interactions()
        T_Al[i] = Al_state1[6]  # T = 6, T_SE = 8
        P_Al[i] = Al_state1[4]  # T = 6, T_SE = 8

        Zr_Al2O3 = ShockCMatter(Zr, Al2O3, SS, u1, ViT, xProj, xTar, xTarb, xSteel, 'Zr')
        Zr_state1 = Zr_Al2O3.shock_interactions()
        T_Zr[i] = Zr_state1[6]
        P_Zr[i] = Zr_state1[4]

    Al_graphs = graphForComparison(T_Al, P_Al, T_Zr, P_Zr, speedRange, 'Al', 'Zr')
    Al_graphs.temperatureVsSpeed()
    Al_graphs.pressureVsSpeed()
    Al_graphs.temperatureVsSpeedComparison()
    Al_graphs.pressureVsSpeedComparison()

    Zr_graphs = graphForComparison(T_Zr, P_Zr, T_Al, P_Al, speedRange, 'Zr', 'Al')
    Zr_graphs.temperatureVsSpeed()
    Zr_graphs.pressureVsSpeed()


if __name__ == '__main__':
    main()
import matplotlib.pyplot as plt


class graphForComparison:

    def __init__(self, T, P, Tcomp, Pcomp, speedRange, name, namecomp):
        # initialize initial state
        self.speed = speedRange.columns  # km/s
        self.T = T  # K
        self.P = P  # GPa
        self.Tcomp = Tcomp  # temperature and pressure for comparison
        self.Pcomp = Pcomp
        self.name = name
        self.namecomp = namecomp

    def temperatureVsSpeed(self):
        plt.plot(self.speed, self.T, label=self.name)
        plt.xlabel('Impact speed, km/s')
        plt.ylabel('Temperature behind shock, K')
        plt.show()

    def pressureVsSpeed(self):
        plt.plot(self.speed, self.P, label=self.name)
        plt.xlabel('Impact speed, km/s')
        plt.ylabel('Pressure behind shock, GPa')
        plt.show()

    def temperatureVsSpeedComparison(self):
        plt.plot(self.speed, self.T, "-b", label=self.name)
        plt.plot(self.speed, self.Tcomp, "-r", label=self.namecomp)
        plt.xlabel('Impact speed, km/s')
        plt.ylabel('Temperature behind shock, K')
        plt.legend(loc="upper left")
        plt.show()

    def pressureVsSpeedComparison(self):
        plt.plot(self.speed, self.P, "-b", label=self.name)
        plt.plot(self.speed, self.Pcomp, "-r", label=self.namecomp)
        plt.xlabel('Impact speed, km/s')
        plt.ylabel('Pressure behind shock, GPa')
        plt.legend(loc="upper left")
        plt.show()

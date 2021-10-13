import math
from _datetime import *
import numpy as np
from itertools import product
from matplotlib import pyplot as plt
from matplotlib.pyplot import rc
import matplotlib.patches as mpatches


def emptySatelliteDictionary():
    numberOfSatelliteWithoutOne = 31
    for key in range(1, numberOfSatelliteWithoutOne + 1):
        Satellite.satelliteDictionary[key] = []


class Satellite:
    satelliteDictionary = {}
    mask = 0

    def __init__(self, satelliteNr, selDate):
        self.satelliteNr = satelliteNr
        self.selDate = selDate
        self.satelliteCoords = np.zeros([0, 3])
        self.satellitesParameters = []

    @staticmethod
    def openFile():
        emptySatelliteDictionary()

        with open("almanac.yuma.week0098.503808.txt", 'r', encoding='utf-8') as openFile:
            numberSatellite = 1
            for line in openFile:
                parametrSatellite = ''.join(line.split()).split(':')
                if not isEmptyLine(parametrSatellite):
                    Satellite.satelliteDictionary[numberSatellite].append(parametrSatellite)
                else:
                    numberSatellite += 1

    def parametersOfAlmanach(self):
        satelliteName = self.satelliteDictionary[self.satelliteNr][0][0]
        idSatellite = int(self.satelliteDictionary[self.satelliteNr][1][1])
        health = float(self.satelliteDictionary[self.satelliteNr][2][1])
        eccentricity = float(self.satelliteDictionary[self.satelliteNr][3][1])
        timeOfApplicability = float(self.satelliteDictionary[self.satelliteNr][4][1])
        orbitalInclination = float(self.satelliteDictionary[self.satelliteNr][5][1])
        RoRA = float(self.satelliteDictionary[self.satelliteNr][6][1])
        aSqrt = float(self.satelliteDictionary[self.satelliteNr][7][1])
        RAaToA = float(self.satelliteDictionary[self.satelliteNr][8][1])
        argumentOfPerigee = float(self.satelliteDictionary[self.satelliteNr][9][1])
        meanAnomaly = float(self.satelliteDictionary[self.satelliteNr][10][1])
        Af0 = float(self.satelliteDictionary[self.satelliteNr][11][1])
        Af1 = float(self.satelliteDictionary[self.satelliteNr][12][1])
        gpsWeek = int(self.satelliteDictionary[self.satelliteNr][13][1])
        self.satellitesParameters = [satelliteName, idSatellite, health, eccentricity, timeOfApplicability,
                                     orbitalInclination, RoRA, aSqrt, RAaToA, argumentOfPerigee, meanAnomaly,
                                     Af0, Af1, gpsWeek]

    def calculateTimeAlmanach(self):
        year, month, day, hour, minute, second = setVariableFromDate(self.selDate)
        numberOfDays = date.toordinal(date(year, month, day)) - date.toordinal(date(1980, 1, 6))
        numberOfWeeks = numberOfDays // 7
        daysInWeeks = numberOfDays % 7
        numberOfWeeks = numberOfWeeks - 2048
        tow = daysInWeeks * 86400 + hour * 3600 + minute * 60 + second
        timeA = numberOfWeeks * 604800 + tow
        return timeA

    # -Czas jaki upłynął-
    def getDeltaTime(self):
        self.parametersOfAlmanach()
        toA = self.satellitesParameters[4]
        gpsWeek = self.satellitesParameters[13]
        timeA = self.calculateTimeAlmanach()
        toa = gpsWeek * 604800 + toA
        tK = timeA - toa
        return tK

    # -Średnia prędkość kątowa-
    def getMeanAngleVelocity__n(self):
        self.parametersOfAlmanach()
        ni = 3.986004415 * 10 ** 14
        A = self.satellitesParameters[7] ** 2
        n = math.sqrt(ni / A ** 3)
        return n

    # -Poprawiona anomalia średnia na epokę tK-    _Mk_
    def getCorrectedMeanAnomaly(self):
        self.parametersOfAlmanach()
        meanAnomaly = self.satellitesParameters[10]
        tK = self.getDeltaTime()
        n = self.getMeanAngleVelocity__n()
        mK = meanAnomaly + n * tK
        mK %= 2 * math.pi
        return mK

    # -Anomalia mimośrodowa (Równanie Kepplera)-   _Ek_
    def getEccentricAnomaly(self):
        self.parametersOfAlmanach()
        eccentricity = self.satellitesParameters[3]
        mK = self.getCorrectedMeanAnomaly()
        Ek_Old = mK
        Ek_new = mK + eccentricity * math.sin(Ek_Old)
        while abs(Ek_new - Ek_Old) > 10 ** (-15):
            Ek_Old = Ek_new
            Ek_new = mK + eccentricity * math.sin(Ek_Old)
        return Ek_new

    # -Anomalia prawdziwa-   _Vk_
    def getRealAnomaly(self):
        self.parametersOfAlmanach()
        eccentricity = self.satellitesParameters[3]
        Ek_new = self.getEccentricAnomaly()
        vK = math.atan2(math.sqrt(1 - eccentricity ** 2) * math.sin(Ek_new), math.cos(Ek_new) - eccentricity)
        if vK < 0:
            vK += 2 * math.pi
        return vK

    # -Argument szerokości-
    def getWidthArgument(self):
        self.parametersOfAlmanach()
        argumentOfPerigee = self.satellitesParameters[9]
        vK = self.getRealAnomaly()
        fiK = vK + argumentOfPerigee
        return fiK

    # -Promień orbity-
    def getOrbitRadius(self):
        self.parametersOfAlmanach()
        A = self.satellitesParameters[7] ** 2
        eccentricity = self.satellitesParameters[3]
        Ek_new = self.getEccentricAnomaly()
        rk = A * (1 - eccentricity * math.cos(Ek_new))
        return rk

    # -Pozycja satelity w układzie orbity-
    def getSatPosInTheOrbit(self):
        self.parametersOfAlmanach()
        fiK = self.getWidthArgument()
        rk = self.getOrbitRadius()
        xk = rk * math.cos(fiK)
        yk = rk * math.sin(fiK)
        return xk, yk

    # -Poprawiona długość węzła-
    def getCorrLenOfTheAscendingNode(self):
        self.parametersOfAlmanach()
        wE = 7.2921151467 * 10 ** (-5)
        toA = self.satellitesParameters[4]
        RoRA = self.satellitesParameters[6]
        RAaToA = self.satellitesParameters[8]
        tK = self.getDeltaTime()
        OmegaK = RAaToA + (RoRA - wE) * tK - wE * toA
        return OmegaK

    # -Pozycja satelity w układzie geocentrycznym ECEF-
    def posOfTheSatInGeocentricSys(self):
        self.openFile()
        self.parametersOfAlmanach()
        orbitalInclination = self.satellitesParameters[5]
        OmegaK = self.getCorrLenOfTheAscendingNode()
        xk, yk = self.getSatPosInTheOrbit()
        Xk = xk * math.cos(OmegaK) - yk * math.cos(orbitalInclination) * math.sin(OmegaK)
        Yk = xk * math.sin(OmegaK) + yk * math.cos(orbitalInclination) * math.cos(OmegaK)
        Zk = yk * math.sin(orbitalInclination)
        return np.array([Xk, Yk, Zk])

    def getParametrsForPosition(self):
        self.parametersOfAlmanach()
        dT = self.getDeltaTime()
        n = self.getMeanAngleVelocity__n()
        mK = self.getCorrectedMeanAnomaly()
        Ek_new = self.getEccentricAnomaly()
        vK = self.getRealAnomaly()
        fiK = self.getWidthArgument()
        rk = self.getOrbitRadius()
        xk, yk = self.getSatPosInTheOrbit()
        OmegaK = self.getCorrLenOfTheAscendingNode()
        Xk, Yk, Zk = self.posOfTheSatInGeocentricSys()
        return ["dT: " + str(dT), "n: " + str(n), "mK: " + str(mK),
                "Ek: " + str(Ek_new), "vK: " + str(vK), "fiK: " + str(fiK),
                "rK: " + str(rk), "xk: " + str(xk), "yk: " + str(yk),
                "OmegaK: " + str(OmegaK), "Xk: " + str(Xk),
                "Yk: " + str(Yk), "Zk: " + str(Zk)]

    def tableOfSatellitesCoords(self):
        satelliteCoords = self.posOfTheSatInGeocentricSys()
        self.satelliteCoords = np.append(self.satelliteCoords, satelliteCoords)

    def allSatellitesCoord(self):
        self.openFile()
        for satelliteNr in range(1, 32):
            self.satelliteNr = satelliteNr
            self.parametersOfAlmanach()
            self.tableOfSatellitesCoords()
        return self.satelliteCoords


class elipsPoint:
    def __init__(self, *args):
        if len(args) > 0:
            self.fi = np.deg2rad(args[0])
            self.Lambda = np.deg2rad(args[1])
            self.h = 100
        elif len(args) == 0:
            self.fi = np.deg2rad(52.0)
            self.Lambda = np.deg2rad(21.0)
            self.h = 100

    def getFi(self):
        return self.fi

    def getLambda(self):
        return self.Lambda

    def getH(self):
        return self.h


class gps:
    def __init__(self, point, satellite):
        self.point = point
        self.satellite = satellite

    def calculateXYZs(self):
        e2 = 0.00669438002290
        a = 6378137
        fi = self.point.getFi()
        Lambda = self.point.getLambda()
        h = self.point.getH()
        N = a / (math.sqrt(1 - (e2 * (math.sin(fi)) ** 2)))
        X = (N + h) * math.cos(fi) * math.cos(Lambda)
        Y = (N + h) * math.cos(fi) * math.sin(Lambda)
        Z = (N * (1 - e2) + h) * math.sin(fi)
        return np.array([X, Y, Z])

    def vectorOfSatelliteGPS(self):
        XYZs = self.calculateXYZs()
        XYZr = self.satellite.posOfTheSatInGeocentricSys()
        Xrs = np.transpose(np.array(XYZr - XYZs))
        return Xrs

    def NEU(self):
        fi = self.point.getFi()
        Lambda = self.point.getLambda()
        return np.array([[-math.sin(fi) * math.cos(Lambda), -math.sin(Lambda), math.cos(fi) * math.cos(Lambda)],
                         [-math.sin(fi) * math.sin(Lambda), math.cos(Lambda), math.cos(fi) * math.sin(Lambda)],
                         [math.cos(fi), 0, math.sin(fi)]])

    def vectorNEU(self):
        return np.dot(np.transpose(self.NEU()), self.vectorOfSatelliteGPS())

    def azymut(self):
        XrsNEU = self.vectorNEU()
        azymut = np.rad2deg(math.atan2(XrsNEU[1], XrsNEU[0]))
        if azymut > 0:
            azymut = azymut
        else:
            azymut = azymut + np.rad2deg(2 * math.pi)
        return azymut

    def elevation(self):
        XrsNEU = self.vectorNEU()
        elevation = np.rad2deg(math.asin(XrsNEU[2] / math.sqrt(XrsNEU[0] ** 2 + XrsNEU[1] ** 2 + XrsNEU[2] ** 2)))
        return elevation

    def vectorLengthOfSatellite(self):
        XYZs = self.calculateXYZs()
        XYZr = self.satellite.posOfTheSatInGeocentricSys()
        lengthVec = math.sqrt((float(XYZr[0]) - float(XYZs[0])) ** 2 + (float(XYZr[1]) - float(XYZs[1])) ** 2 +
                              (float(XYZr[2]) - float(XYZs[2])) ** 2)
        return lengthVec


def plotElevation(elevation, day):
    dayGPS = date(day[0], day[1], day[2])
    anotherDayGPS = dayGPS + timedelta(days=1)
    fig, ax = plt.subplots()
    plt.style.use('ggplot')
    fig.set_size_inches(17, 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set(ylim=(0, 100))
    plt.subplots_adjust(bottom=0.15)
    hours = [str(i) + ":" + (str(j) if j > 0 else str(j) + "0") for i in range(0, 25) for j in range(0, 51, 10) if
             i < 24 or j <= 0]
    indSatellitetoChange = 9
    for i, elem in enumerate(elevation):
        if i > indSatellitetoChange:
            ax.plot(hours, elem, lw=1.7, alpha=0.7, label=f"G{i + 2}")
        else:
            ax.plot(hours, elem, lw=1.7, alpha=0.7, label=f"G0{i + 1}" if i != indSatellitetoChange else f"G{i + 1}")
    labels = [str(hours[i]) for i in range(0, len(hours), 12)]
    labels[-1] = anotherDayGPS.strftime("%d-%m-%Y")
    plt.title(f'Day: {dayGPS.strftime("%d-%m-%Y")}')
    plt.xticks(labels)
    ax.set_xlabel("GPS time [h]")
    ax.set_ylabel("elevation [o]")
    ax.legend(bbox_to_anchor=(1.1, 1.18), loc='upper right')
    plt.show()
    # fig.savefig("plotsSatellite1-3_.pdf", bbox_inches="tight", transparent=True)


def getHoursPerDay():
    hours = [i for i in range(0, 25)]
    minutes = [i for i in range(0, 51, 10)]
    return [(hours, minutes) for hours, minutes in product(hours, minutes) if hours < 24 or minutes <= 0]


def plotGPS():
    elev = []
    gpsCoord = elipsPoint()
    gpsDay = [2021, 3, 1, 0, 00, 00]
    for i in range(1, 32):
        elev.append([])
        for hour, minutes in getHoursPerDay():
            selectedDate = [gpsDay[0], gpsDay[1], gpsDay[2], hour, minutes, gpsDay[5]]
            satellite = Satellite(i, selectedDate)
            recv = gps(gpsCoord, satellite)
            elev[i - 1].append(recv.elevation())
    plotElevation(elev, gpsDay)


def Skyplot(satPosition):
    # sat_positions - [PRN, el, az] w stopniach
    rc('grid', color='gray', linewidth=1, linestyle='--')
    fontSize = 20
    rc('xtick', labelsize=fontSize)
    rc('ytick', labelsize=fontSize)
    rc('font', size=fontSize)
    # define colors

    green = '#467821'
    # blue = '#348ABD'
    # red = '#A60628'
    # orange = '#E24A33'
    # purple = '#7A68A6'

    # start ploting
    fig = plt.figure(figsize=(8, 6))
    plt.subplots_adjust(bottom=0.1,
                        top=0.85,
                        left=0.1,
                        right=0.74)
    ax = fig.add_subplot(polar=True)  # define a polar type of coordinates
    ax.set_theta_zero_location('N')  # ustawienie kierunku północy na górze wykresu
    ax.set_theta_direction(-1)  # ustawienie kierunku przyrostu azymutu w prawo

    PG = 0  # zliczanie satelitów GPS

    for (PRN, el, az) in satPosition:
        PG += 1
        # show sat number
        ax.annotate(PRN,
                    xy=(np.deg2rad(az[int(len(az) / 2)]), 90 - el[int(len(az) / 2)]),
                    bbox=dict(boxstyle="round", fc=green, alpha=0.5),
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    color='k')
    Gps = mpatches.Patch(color=green, label='{:02.0f}  GPS'.format(PG))
    plt.legend(handles=[Gps], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # axis ticks descriptions

    plt.title(f'Skyplot')
    ax.set_yticks(range(0, 90 + 10, 10))  # Define the yticks
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)

    # saving and showing plot
    # plt.savefig('satellite_skyplot.pdf')
    plt.show()  # wyświetleni


def satPositions():
    sat_positions = []
    for i in range(1, 32):
        Satellite(i, [2021, 3, 1, 0, 00, 00]).parametersOfAlmanach()
        GPS = gps(elipsPoint(), Satellite(i, [2021, 3, 1, 0, 00, 00]))
        el = GPS.elevation()
        az = GPS.azymut()
        if el > 0:
            sat_positions.append(["G" + str([i]), [el], [az]])
    return sat_positions


def plotSkyplot():
    satPosition = satPositions()
    Skyplot(satPosition)


def isEmptyLine(line):
    return line == ['']


def setVariableFromDate(data):
    return data[0], data[1], data[2], data[3], data[4], data[5]


def gpsVisibleSatellite(gpsCoord, selectedDate):
    visibleSatellite = []
    for i in range(1, 32):
        satellite = Satellite(i, selectedDate)
        recv = gps(gpsCoord, satellite)
        if recv.elevation() > satellite.mask:
            visibleSatellite.append(satellite.satelliteNr)
    return visibleSatellite


def GPSDay():
    gpsday = [2021, 3, 1, 0, 00, 00]
    return gpsday


def getQ(gpsCoord, hour):
    matrixA = np.zeros([0, 4])
    for nrSatellite in gpsVisibleSatellite(gpsCoord, hour):
        satellite = Satellite(nrSatellite, hour)
        recv = gps(gpsCoord, satellite)
        satelliteXYZ = satellite.posOfTheSatInGeocentricSys()
        gpsXYZ = recv.calculateXYZs()
        lengthVector = recv.vectorLengthOfSatellite()
        matrixA = np.append(matrixA, np.array(
            [-(satelliteXYZ[0] - gpsXYZ[0]) / lengthVector, -(satelliteXYZ[1] - gpsXYZ[1]) / lengthVector,
             -(satelliteXYZ[2] - gpsXYZ[2]) / lengthVector, 1]))
    matrixA = np.reshape(matrixA, (matrixA.size // 4, 4))
    matrixQ = np.linalg.inv(np.dot(matrixA.T, matrixA))
    return matrixQ


def GDOP(gpsCoord, hour):
    matrixQ = getQ(gpsCoord, hour)
    return math.sqrt(matrixQ.diagonal()[0] + matrixQ.diagonal()[1] + matrixQ.diagonal()[2] + matrixQ.diagonal()[3])


def PDOP(gpsCoord, hour):
    matrixQ = getQ(gpsCoord, hour)
    return math.sqrt(matrixQ.diagonal()[0] + matrixQ.diagonal()[1] + matrixQ.diagonal()[2])


def TDOP(gpsCoord, hour):
    matrixQ = getQ(gpsCoord, hour)
    return math.sqrt(matrixQ.diagonal()[3])


def QNEU(gpsCoord, hour):
    gpsday = GPSDay()
    satellite = Satellite(None, gpsday)
    recv = gps(gpsCoord, satellite)
    rNEU = recv.NEU()
    qXYZ = getQ(gpsCoord, hour)[0:3, 0:3]
    qNEU = np.dot(np.dot(rNEU.T, qXYZ), rNEU)
    return qNEU


def HDOP(gpsCoord, hour):
    matrixQNEU = QNEU(gpsCoord, hour)
    return math.sqrt(matrixQNEU.diagonal()[0] + matrixQNEU.diagonal()[1])


def VDOP(gpsCoord, hour):
    matrixQNEU = QNEU(gpsCoord, hour)
    return math.sqrt(matrixQNEU.diagonal()[2])


def DOP(gpsCoord, hour):
    gdop = GDOP(gpsCoord, hour)
    pdop = PDOP(gpsCoord, hour)
    tdop = TDOP(gpsCoord, hour)
    return [gdop, pdop, tdop]


def getTimeInSeconds(data):
    dDay = date.toordinal(date(data[0], data[1], data[2])) - date.toordinal(date(1980, 1, 6))
    week = dDay // 7
    day = dDay % 7
    week = week - 2048
    tow = day * 86400 + data[3] * 3600 + data[4] * 60 + data[5]
    timeK = week * 604800 + tow
    return timeK


"""

def plotDOP():
    # legend = ['Geometryczne osłabienie precyzji', 'Pozycyjne osłabienie precyzji', 'Osłabienie precyzji w czasie']

    def get_DOP(Q):
        gdop = math.sqrt(Q[0][0] + Q[1][1] + Q[2][2] + Q[3][3])
        pdop = math.sqrt(Q[0][0] + Q[1][1] + Q[2][2])
        tdop = math.sqrt(Q[3][3])

        return [gdop, pdop, tdop]

    gpsStart = getTimeInSeconds([2021, 3, 1, 0, 0, 00])
    gpsEnd = getTimeInSeconds([2021, 3, 2, 0, 0, 00])
    plt.figure(figsize=(10, 5))
    ax = plt.gca()
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')
    timeInSeconds = []
    DOPtab = []
    for sec in range(gpsStart, gpsEnd, 600):
        timeInSeconds.append(sec)
        visible = []
        for i in range(1, 32):
            GPS = gps(elipsPoint(), Satellite(i, [2021, 3, 1, 0, 00, 00]))
            el = GPS.elevation()
            if el > 0:
                visible.append(i)
        matrix_A = np.zeros([0, 4])
        for i in range(1, 32):
            GPS = gps(elipsPoint(), Satellite(i, [2021, 3, 1, 0, 00, 00]))
            xs, ys, zs = Satellite(i, [2021, 3, 1, 0, 00, 00]).posOfTheSatInGeocentricSys()
            x0, y0, z0 = GPS.calculateXYZs()
            p0s = math.sqrt(((xs - x0) ** 2) + ((ys - y0) ** 2) + ((zs - z0) ** 2))
            matrix_A = np.append(matrix_A, np.array([-(xs - x0) / p0s, -(ys - y0) / p0s, -(zs - z0) / p0s, 1]))
        matrix_A = np.reshape(matrix_A, (matrix_A.size // 4, 4))
        matrix_Q = np.linalg.inv(np.dot(matrix_A.T, matrix_A))
        dop = get_DOP(matrix_Q)
        DOPtab.append(dop)

        plt.plot(DOPtab)

    # plt.legend(legend)
    plt.title("DOP's", color='black')
    plt.xlabel("Time", color='white')
    plt.ylabel("Values", color='white')
    plt.show()



def plotVisible():
    for i in range(1, 32):
        Satellite(i, [2021, 3, 1, 0, 00, 00]).parametersOfAlmanach()
        alm_data = Satellite(i, [2021, 3, 1, 0, 00, 00]).satellitesParameters
        visibleSatellite = []
        gpsStart = getTimeInSeconds([2021, 3, 1, 0, 0, 00])
        gpsEnd = getTimeInSeconds([2021, 3, 2, 0, 0, 00])
        plt.figure(figsize=(10, 5), facecolor='black')
        for sec in range(gpsStart, gpsEnd, 600):
            visibleSatellite.append(0)
        for idx, value in enumerate(alm_data):
            time_in_sec = []
            i = 0
            for sec in range(gpsStart, gpsEnd, 600):
                sec = Satellite(i, [2021, 3, 1, 0, 00, 00]).calculateTimeAlmanach()
                GPS = gps(elipsPoint(), Satellite(i, [2021, 3, 1, 0, 00, 00]))
                el = GPS.elevation()
                time_in_sec.append(sec)

            if el > 0:
                visibleSatellite[i] += 1
            i += 1

        ax = plt.gca()
        ax.set_facecolor('black')
        ax.spines['left'].set_color('white')
        ax.spines['bottom'].set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        plt.title("Number of Visible Satellites", color='white')
        plt.xlabel("Time", color='white')
        plt.ylabel("Number of Satellites", color='white')
        plt.plot(time_in_sec, visibleSatellite)
        plt.show()
"""

if __name__ == "__main__":
    Satellite.openFile()
    # print(Satellite.satelliteDictionary)
    # print(Satellite(1, [2021, 3, 1, 0, 00, 00]).posOfTheSatInGeocentricSys())
    # Satellite(1, [2021, 3, 1, 0, 00, 00]).parametersOfAlmanach()
    # print(Satellite(1, [2021, 3, 1, 0, 00, 00]).satellitesParameters)
    # print(Satellite(1, [2021, 3, 1, 0, 00, 00]).getParametrsForPosition())
    # print(Satellite(None, [2021, 3, 1, 0, 00, 00]).allSatellitesCoord())
    # print(gps(elipsPoint(), Satellite(1, [2021, 3, 1, 0, 00, 00])).vectorOfSatelliteGPS())
    # print(gps(elipsPoint(), Satellite(1, [2021, 3, 1, 0, 00, 00])).NEU())
    # print(gps(elipsPoint(), Satellite(1, [2021, 3, 1, 0, 00, 00])).vectorNEU())
    # print(gps(elipsPoint(), Satellite(1, [2021, 3, 1, 0, 00, 00])).azymutAndElevation())
    # print(gps(elipsPoint(), Satellite(1, [2021, 3, 1, 0, 00, 00])).vectorLengthOfSatellite())
    # print(gpsVisibleSatellite(elipsPoint(), [2021, 3, 1, 0, 00, 00]))
    # print(getQ(elipsPoint(), [2021, 3, 1, 0, 00, 00]))
    # print(QNEU(elipsPoint(), [2021, 3, 1, 0, 00, 00]))
    # print(VDOP(elipsPoint(), [2021, 3, 1, 0, 00, 00]))
    plotGPS()
    plotSkyplot()
    # plotDOP()
    # plotVisible()

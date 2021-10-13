import math
from datetime import *
import numpy as np
from itertools import product

np.set_printoptions(suppress=True)

class Satellite:
    def __init__(self, nrSatellite, selectedDate):
        self.satelliteDict = {}
        self.nrSatellite = nrSatellite
        self.parametrOfSatellite = []
        self.coordSatellite = np.zeros([0, 3])
        self.selectedDate = selectedDate

    def getNrSatellite(self):
        return self.nrSatellite

    def getDate(self):
        return self.selectedDate

    def openFile(self):
        fileName = 'D:\PyCharm 2020.1.2\projects\lmanac.yuma.week0098.503808.txt'
        # fileName = tkinter.filedialog.askopenfilename(filetypes=(("txt file", "*.txt"),
        # ("all files", "*.*")))
        self.emptySatelliteDictionary()

        with open(fileName, 'r', encoding='utf-8') as openFile:
            numberSatellite = 1
            for line in openFile:
                parametrSatellite = ''.join(line.split()).split(':')
                if not isEmptyLine(parametrSatellite):
                    self.satelliteDict[numberSatellite].append(parametrSatellite)
                else:
                    numberSatellite += 1

    def emptySatelliteDictionary(self):
        numberOfSatelliteWithoutOne = 31
        for key in range(1, numberOfSatelliteWithoutOne + 1):
            self.satelliteDict[key] = []

    def setTableOfParametr(self):
        nameSatellite = self.satelliteDict[self.nrSatellite][0][0]
        idSatellite = int(self.satelliteDict[self.nrSatellite][1][1])
        health = float(self.satelliteDict[self.nrSatellite][2][1])
        eParametr = float(self.satelliteDict[self.nrSatellite][3][1])
        toA = float(self.satelliteDict[self.nrSatellite][4][1])
        orbitIncl = float(self.satelliteDict[self.nrSatellite][5][1])
        RoRA = float(self.satelliteDict[self.nrSatellite][6][1])
        aSqrt = float(self.satelliteDict[self.nrSatellite][7][1])
        RAW = float(self.satelliteDict[self.nrSatellite][8][1])
        perigeeParametr = float(self.satelliteDict[self.nrSatellite][9][1])
        mAnomalia = float(self.satelliteDict[self.nrSatellite][10][1])
        Af0 = float(self.satelliteDict[self.nrSatellite][11][1])
        Af1 = float(self.satelliteDict[self.nrSatellite][12][1])
        gps_week = int(self.satelliteDict[self.nrSatellite][13][1])
        self.parametrOfSatellite = [nameSatellite, idSatellite, health, eParametr, toA, orbitIncl, RoRA, aSqrt, RAW,
                                    perigeeParametr, mAnomalia, Af0, Af1, gps_week]

    def getParametr(self):
        deltaTime = self.getDeltaTime()
        nParametr = self.getNParametr()
        Mk = self.getMk()
        Ek = self.getEk()
        vk = self.getVk()
        fiK = self.getFiK()
        rK = self.getRadiusK()
        xk, yk = self.getPositionSatellite()
        omegaK = self.getOmegaK()
        Xk, Yk, Zk = self.getCoordSatellite()
        return [deltaTime, nParametr, Mk, Ek, vk, fiK, rK, xk, yk, omegaK, Xk, Yk, Zk]

    def getDeltaTime(self):
        toA = self.parametrOfSatellite[4]
        gps_week = self.parametrOfSatellite[13]
        tSelectDay = self.calcTimeOfAlmanach()
        toa = gps_week * 604800 + toA
        return tSelectDay - toa

    def calcTimeOfAlmanach(self):
        year, month, day, hour, minute, second = setVariableFromDate(self.selectedDate)
        daysFromStartGPStoDate = date.toordinal(date(year, month, day)) - date.toordinal(date(1980, 1, 6))
        weekFromStartGPStoDate = daysFromStartGPStoDate // 7
        numberOfDayInWeek = daysFromStartGPStoDate % 7
        towParametr = numberOfDayInWeek * 86400 + hour * 3600 + minute * 60 + second
        weekFromStartGPStoDate -= 2048
        return weekFromStartGPStoDate * 604800 + towParametr

    def getNParametr(self):
        U = 3.986004415 * 10 ** 14
        a = self.parametrOfSatellite[7] ** 2
        return math.sqrt(U / a ** 3)

    def getMk(self):
        mAnomalia = self.parametrOfSatellite[10]
        deltaTime = self.getDeltaTime()
        nParametr = self.getNParametr()
        mK = mAnomalia + nParametr * deltaTime
        mK %= 2 * math.pi
        return mK

    def getEk(self):
        e = self.parametrOfSatellite[3]
        Mk = self.getMk()
        previousE = Mk
        nextE = Mk + e * math.sin(previousE)
        while abs((previousE - nextE)) > 10 ** (-15):
            previousE = nextE
            nextE = Mk + e * math.sin(previousE)
        return nextE

    def getVk(self):
        e = self.parametrOfSatellite[3]
        Ek = self.getEk()
        vk = math.atan2(math.sqrt(1 - e ** 2) * math.sin(Ek), math.cos(Ek) - e)
        return vk if vk > 0 else vk + 2 * math.pi

    def getFiK(self):
        vk = self.getVk()
        perigeeParametr = self.parametrOfSatellite[9]
        return vk + perigeeParametr

    def getRadiusK(self):
        e = self.parametrOfSatellite[3]
        a = self.parametrOfSatellite[7] ** 2
        Ek = self.getEk()
        return a * (1 - e * math.cos(Ek))

    def getPositionSatellite(self):
        rK = self.getRadiusK()
        fiK = self.getFiK()
        xk = rK * math.cos(fiK)
        yk = rK * math.sin(fiK)
        return xk, yk

    def getOmegaK(self):
        wE = 7.2921151467 * 10 ** (-5)
        toA = self.parametrOfSatellite[4]
        RoRA = self.parametrOfSatellite[6]
        RAW = self.parametrOfSatellite[8]
        deltaTime = self.getDeltaTime()
        return RAW + (RoRA - wE) * deltaTime - wE * toA

    def getCoordSatellite(self):
        self.openFile()
        self.setTableOfParametr()
        orbitIncl = self.parametrOfSatellite[5]
        xk, yk = self.getPositionSatellite()
        omegaK = self.getOmegaK()
        Xk = xk * math.cos(omegaK) - yk * math.cos(orbitIncl) * math.sin(omegaK)
        Yk = xk * math.sin(omegaK) + yk * math.cos(orbitIncl) * math.cos(omegaK)
        Zk = yk * math.sin(orbitIncl)
        return np.array([Xk, Yk, Zk])

    def getCoordForAllSatellite(self):
        self.openFile()
        for nrSatellite in range(1, 32):
            self.nrSatellite = nrSatellite
            self.setTableOfParametr()
            self.setTableOfCoordSatellite()
        return self.coordSatellite

    def setTableOfCoordSatellite(self):
        coordSatellite = self.getCoordSatellite()
        self.coordSatellite = np.append(self.coordSatellite, coordSatellite, axis=0)

    def saveTableOfCoordSattelite(self):
        np.savetxt('coords.txt', self.coordSatellite, delimiter='    ', fmt='%.3f')

    def checkAnswer(self):
        answer = self.getParametr()
        print(answer[0], True if answer[0] == 187392.0 else False)
        print(answer[1], True if answer[1] == 0.00014585706906799988 else False)
        print(answer[2], True if answer[2] == 1.0839412390722885 else False)
        print(answer[3], True if answer[3] == 1.093207444795061 else False)
        print(answer[4], True if round(answer[4], 14) == 1.10249610963853 else False)
        print(answer[5], True if round(answer[5], 14) == 1.92220922863853 else False)
        print("%.3f" % answer[6], True if round(answer[6], 3) == 26432597.103 else False)
        print("%.3f" % answer[7], True if round(answer[7], 3) == -9098752.972 else False)
        print("%.3f" % answer[8], True if round(answer[8], 3) == 24817229.579 else False)
        print(answer[9], True if answer[9] == -52.216185427270695 else False)
        print("%.3f" % answer[10], True if round(answer[10], 3) == 16152870.357 else False)
        print("%.3f" % answer[11], True if round(answer[11], 3) == 3347394.887 else False)
        print("%.3f" % answer[12], True if round(answer[12], 3) == 20653375.422 else False)


def isEmptyLine(line):
    return line == ['']


def setVariableFromDate(data):
    return data[0], data[1], data[2], data[3], data[4], data[5]


class ElipsoidPoint:
    def __init__(self, *args):
        if len(args) > 0:
            self.fi = np.deg2rad(args[0])
            self.lamb = np.deg2rad(args[1])
            self.h = np.deg2rad(args[2])
        else:
            self.fi = np.deg2rad(52.0)
            self.lamb = np.deg2rad(21.0)
            self.h = 100.0

    def getFi(self):
        return self.fi

    def getLambda(self):
        return self.lamb

    def getH(self):
        return self.h


class GPS:
    def __init__(self, point, satellite):
        self.point = point
        self.satellite = satellite

    def calcXYZ(self):
        e2 = 0.00669438002290
        N = self.getN()
        fi, lamb, h = self.getCoord()
        X = (N + h) * math.cos(fi) * math.cos(h)
        Y = (N + h) * math.cos(fi) * math.sin(lamb)
        Z = (N * (1 - e2) + h) * math.sin(fi)
        return np.array([X, Y, Z])

    def getN(self):
        aElipsoid = 6378137
        e2 = 0.00669438002290
        return aElipsoid / (math.sqrt(1 - e2 * math.sin(self.point.getFi()) ** 2))

    def getCoord(self):
        return self.point.getFi(), self.point.getLambda(), self.point.getH()

    def getVectorSatelliteGPS(self):
        geocentric_XYZ = self.calcXYZ()
        coordSatellite = self.satellite.getCoordSatellite()
        return np.transpose(np.array(coordSatellite - geocentric_XYZ))

    def getNEU(self):
        fi = self.point.getFi()
        lamb = self.point.getLambda()
        return np.array([[-math.sin(fi) * math.cos(lamb), -math.sin(lamb), math.cos(fi) * math.cos(lamb)],
                         [-math.sin(fi) * math.sin(lamb), math.cos(lamb), math.cos(fi) * math.sin(lamb)],
                         [math.cos(fi), 0, math.sin(fi)]])

    def getVectorNEU(self):
        return np.dot(np.transpose(self.getNEU()), self.getVectorSatelliteGPS())

    def getAzymuth(self):
        NEU = self.getVectorNEU()
        return np.rad2deg(math.atan2(NEU[1], NEU[0])) if np.rad2deg(math.atan2(NEU[1], NEU[0])) > 0 else np.rad2deg(math.atan2(NEU[1], NEU[0])) + 360

    def getElevation(self):
        NEU = self.getVectorNEU()
        return np.rad2deg(math.asin(NEU[2] / (math.sqrt(NEU[0] ** 2 + NEU[1] ** 2 + NEU[2] ** 2))))

    def getCosZ(self):
        NEU = self.getVectorNEU()
        return NEU[2] / (math.sqrt(NEU[0] ** 2 + NEU[1] ** 2 + NEU[2] ** 2))


def plotElevation(elevation, day):
    dayGPS = date(day[0], day[1], day[2])
    anotherDayGPS = dayGPS + timedelta(days=1)
    fig, ax = plt.subplots()
    plt.style.use('ggplot')
    fig.set_size_inches(17, 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set(ylim=(0, 100))
    ax.legend(bbox_to_anchor=(1.15, 1.2))
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
    fig.savefig("plotsSatellite1-3.pdf", bbox_inches="tight", transparent=True)


def coordSatelliteCalculator():
    selectedDate = [2021, 3, 1, 0, 00, 00]
    nrSatellite = None
    satellite = Satellite(nrSatellite, selectedDate)
    print(satellite.getCoordForAllSatellite())


def getHoursPerDay():
    hours = [i for i in range(0, 25)]
    minutes = [i for i in range(0, 51, 10)]
    return [(hours, minutes) for hours, minutes in product(hours, minutes) if hours < 24 or minutes <= 0]


def plotGPS():
    elev = []
    gpsCoord = ElipsoidPoint()
    gpsDay = [2021, 3, 1, 0, 00, 00]
    for i in range(1, 32):
        elev.append([])
        for hour, minutes in getHoursPerDay():
            selectedDate = [gpsDay[0], gpsDay[1], gpsDay[2], hour, minutes, gpsDay[5]]
            satellite = Satellite(i, selectedDate)
            recv = GPS(gpsCoord, satellite)
            elev[i - 1].append(recv.getElevation())
    plotElevation(elev, gpsDay)

if __name__ == "__main__":
    plotGPS()
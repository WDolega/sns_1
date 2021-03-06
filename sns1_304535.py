import sys
import math
import numpy as np
import webbrowser
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc, rcParams, grid
import matplotlib.patches as mpatches
matplotlib.use('Qt5Agg')
from datetime import *
from itertools import product
from random import randint
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QComboBox, QLabel, QLineEdit, QGridLayout, QPushButton, QMessageBox, QFileDialog, \
    QHBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


def linkEvent():
    webbrowser.open('http://www.navcen.uscg.gov/?pageName=currentAlmanac&format=yuma')


def isEmptyLine(line):
    return line == ['']


def setVariableFromDate(data):
    return data[0], data[1], data[2], data[3], data[4], data[5]


def getYear():
    years = [str(i) for i in range(2021, 1979, -1)]
    return years


def getMonth():
    months = [str(i) for i in range(1, 13)]
    return months


def getSat():
    sats = [str(i) for i in range(31, 0, -1)]
    return sats


def gpsVisibleSatellite(gpsCoord, selectedDate, file_path, cutoff, numb_sat):
    visibleSatellite = []
    for i in range(1, numb_sat):
        satellite = Satellite(i, selectedDate, file_path)
        gps = GPS(gpsCoord, satellite, cutoff)
        try:
            if gps.getElevation() > 0:
                visibleSatellite.append(satellite.nrSatellite)
        except TypeError:
            pass
    return visibleSatellite


def getQ(coords, selectedDate, file_path, cutoff, numb_sat):
    matrixA = np.zeros([0, 4])
    gpsCoord = ElipsoidPoint(coords)
    for nrSatellite in gpsVisibleSatellite(gpsCoord, selectedDate, file_path, cutoff, numb_sat):
        satellite = Satellite(nrSatellite, selectedDate, file_path)
        satelliteXYZ = satellite.getCoordSatellite()
        gps = GPS(gpsCoord, satellite, cutoff)
        gpsXYZ = gps.calcXYZ()
        lengthVector = gps.vectorLengthOfSatellite()
        matrixA = np.append(matrixA, np.array(
            [-(satelliteXYZ[0] - gpsXYZ[0]) / lengthVector, -(satelliteXYZ[1] - gpsXYZ[1]) / lengthVector,
             -(satelliteXYZ[2] - gpsXYZ[2]) / lengthVector, 1]))
    matrixA = np.reshape(matrixA, (matrixA.size // 4, 4))
    matrixQ = np.linalg.inv(np.dot(matrixA.T, matrixA))
    return matrixQ


def GDOP(matQ):
    matrixQ = matQ
    return math.sqrt(matrixQ.diagonal()[0] + matrixQ.diagonal()[1] + matrixQ.diagonal()[2] + matrixQ.diagonal()[3])


def PDOP(matQ):
    matrixQ = matQ
    return math.sqrt(matrixQ.diagonal()[0] + matrixQ.diagonal()[1] + matrixQ.diagonal()[2])


def TDOP(matQ):
    matrixQ = matQ
    return math.sqrt(matrixQ.diagonal()[3])


def QNEU(gpsCoord, matQ, selectedDate, file_path, cutoff, numb_sat):
    for nrSatellite in gpsVisibleSatellite(gpsCoord, selectedDate, file_path, cutoff, numb_sat):
        satellite = Satellite(nrSatellite, selectedDate, file_path)
        gps = GPS(gpsCoord, satellite, cutoff)
        rNEU = gps.getNEU()
        qXYZ = matQ[0:3, 0:3]
        qNEU = np.dot(np.dot(rNEU.T, qXYZ), rNEU)
    return qNEU


def HDOP(gpsCoord, matQ, selectedDate, file_path, cutoff, numb_sat):
    matrixQNEU = QNEU(gpsCoord, matQ, selectedDate, file_path, cutoff, numb_sat)
    return math.sqrt(matrixQNEU.diagonal()[0] + matrixQNEU.diagonal()[1])


def VDOP(gpsCoord, matQ, selectedDate, file_path, cutoff, numb_sat):
    matrixQNEU = QNEU(gpsCoord, matQ, selectedDate, file_path, cutoff, numb_sat)
    return math.sqrt(matrixQNEU.diagonal()[2])


def plotDops(coords, y, m, d, file_path, cutoff, numb_sat):
    gdop = []
    pdop = []
    tdop = []
    # hdop = []
    # vdop = []
    gpsDay = [y, m, d, 0, 00, 00]
    for hour, minutes in getHoursPerDay():
        selectedDate = [gpsDay[0], gpsDay[1], gpsDay[2], hour, minutes, gpsDay[5]]
        matrix = getQ(coords, selectedDate, file_path, cutoff, numb_sat)
        gdop.append(GDOP(matrix))
        pdop.append(PDOP(matrix))
        tdop.append(TDOP(matrix))
        # hdop.append(HDOP(coords, matrix, selectedDate, file_path, cutoff, numb_sat))
        # vdop.append(VDOP(coords, matrix, selectedDate, file_path, cutoff, numb_sat))

    dayGPS = date(gpsDay[0], gpsDay[1], gpsDay[2])
    fig, ax = plt.subplots()
    plt.style.use('ggplot')
    fig.set_size_inches(19, 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    hours = [str(i) + ":" + (str(j) if j > 0 else str(j) + "0") for i in range(0, 25) for j in range(0, 51, 10) if
             i < 24 or j <= 0]
    ax.plot(hours, gdop, color='#5fdce3', lw=1.7, alpha=0.7, label='GDOP')
    ax.plot(hours, pdop, color='#59f74a', lw=1.7, alpha=0.7, label='PDOP')
    ax.plot(hours, tdop, color='#a60010', lw=1.7, alpha=0.7, label='TDOP')
    # ax.plot(hours, hdop, color='#ff9500', lw=1.7, alpha=0.7, label='HDOP')
    # ax.plot(hours, vdop, color='#ff19fb', lw=1.7, alpha=0.7, label='VDOP')
    ax.legend()
    x_labels = [str(i) + ':00' for i in range(0, 25, 2)]
    plt.xlim(0, 6)
    plt.ylim(bottom=0)
    plt.title(f'DOP\'s {dayGPS.strftime("%d-%m-%Y")}')
    plt.xticks(x_labels)
    plt.yticks(np.arange(0, max(gdop) + 1, 1.0))
    ax.set_xlabel("GPS time [h]")
    ax.set_ylabel("DOP's values")
    plt.show()


def plotElevation(elevation, day):
    dayGPS = date(day[0], day[1], day[2])
    anotherDayGPS = dayGPS + timedelta(days=1)
    fig, ax = plt.subplots()
    plt.style.use('ggplot')
    fig.set_size_inches(19, 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set(ylim=(0, 90), xlim=(0, 24))
    plt.subplots_adjust(bottom=0.10)
    hours = [str(i) + ":" + (str(j) if j > 0 else str(j) + "0") for i in range(0, 25) for j in range(0, 51, 10) if
             i < 24 or j <= 0]
    indSatellitetoChange = 9
    for i, elem in enumerate(elevation):
        if i > indSatellitetoChange:
            ax.plot(hours, elem, lw=1.7, alpha=0.7, label=f"G{i + 2}")
        else:
            ax.plot(hours, elem, lw=2.7, alpha=0.7, label=f"G0{i + 1}" if i != indSatellitetoChange else f"G{i + 1}")
    x_labels = [str(i) + ':00' for i in range(0, 25, 2)]
    y_labels = [i for i in range(0, 91, 10)]
    plt.title(f'Elevation {dayGPS.strftime("%d-%m-%Y")}')
    plt.xticks(x_labels)
    plt.yticks(y_labels)
    ax.set_xlabel("GPS time [h]")
    ax.set_ylabel("Elevation [??]")
    ax.legend(loc='upper right', bbox_to_anchor=(1.12, 1.10), borderaxespad=0., ncol=2)
    plt.show()


def plotVisibleSat(coords, y, m, d, file_path, cutoff, numb_sat):
    vis_sat = []
    visible = []
    gpsCoord = ElipsoidPoint(coords)
    gpsDay = [y, m, d, 0, 00, 00]
    for hour, minutes in getHoursPerDay():
        p = len(vis_sat)
        for i in range(1, numb_sat):
            selectedDate = [gpsDay[0], gpsDay[1], gpsDay[2], hour, minutes, gpsDay[5]]
            satellite = Satellite(i, selectedDate, file_path)
            gps = GPS(gpsCoord, satellite, cutoff)
            try:
                if gps.getElevation() > cutoff:
                    vis_sat.append(i)
            except TypeError:
                pass
        visible.append(len(vis_sat) - p)
    dayGPS = date(gpsDay[0], gpsDay[1], gpsDay[2])
    fig, ax = plt.subplots()
    plt.style.use('ggplot')
    fig.set_size_inches(19, 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    hours = [str(i) + ":" + (str(j) if j > 0 else str(j) + "0") for i in range(0, 25) for j in range(0, 51, 10) if
             i < 24 or j <= 0]
    ax.plot(hours, visible, color='#a60010', lw=1.7, alpha=0.7)
    x_labels = [str(i) + ':00' for i in range(0, 25, 2)]
    plt.xlim(0, 24)
    plt.ylim(bottom=0)
    plt.title(f'Number of visible satellites {dayGPS.strftime("%d-%m-%Y")}')
    plt.xticks(x_labels)
    plt.yticks(np.arange(0, max(visible) + 1, 1.0))
    plt.fill_between(hours, visible, facecolor='#f7818c')
    ax.set_xlabel("GPS time [h]")
    ax.set_ylabel("Elevation [??]")
    plt.show()


def getHoursPerDay():
    hours = [i for i in range(0, 25)]
    minutes = [i for i in range(0, 51, 10)]
    return [(hours, minutes) for hours, minutes in product(hours, minutes) if hours < 24 or minutes <= 0]


def plotGPS(coords, y, m, d, file_path, cutoff, numb_sat):
    elev = []
    gpsCoord = ElipsoidPoint(coords)
    gpsDay = [y, m, d, 0, 00, 00]
    for i in range(1, numb_sat):
        elev.append([])
        for hour, minutes in getHoursPerDay():
            selectedDate = [gpsDay[0], gpsDay[1], gpsDay[2], hour, minutes, gpsDay[5]]
            satellite = Satellite(i, selectedDate, file_path)
            gps = GPS(gpsCoord, satellite, cutoff)
            elev[i - 1].append(gps.getElevation())

    plotElevation(elev, gpsDay)


def plotSkyplot(coords, y, m, d, h, min, file_path, cutoff, numb_sat):
    sat_positions = []
    gpsCoord = ElipsoidPoint(coords)
    gpsDay = [y, m, d, h, min, 00]
    for i in range(1, numb_sat):
        selectedDate = [gpsDay[0], gpsDay[1], gpsDay[2], gpsDay[3], gpsDay[4], gpsDay[5]]
        satellite = Satellite(i, selectedDate, file_path)
        gps = GPS(gpsCoord, satellite, cutoff)
        el = gps.getElevation()
        az = gps.getAzymuth()
        try:
            if el > 0:
                sat_positions.append(["G" + str(i), [el], [az]])
        except TypeError:
            pass
    # elev = []
    # azy = []
    # for i in range(1, numb_sat):
    #     elev.append([])
    #     azy.append([])
    #     for hour, minutes in getHoursPerDay():
    #         selectedDate = [gpsDay[0], gpsDay[1], gpsDay[2], hour, minutes, gpsDay[5]]
    #         satellite = Satellite(i, selectedDate, file_path)
    #         gps = GPS(gpsCoord, satellite, cutoff)
    #         try:
    #             e = gps.getElevation()
    #             if e > 0:
    #                 elev[i - 1].append(90-e)
    #                 azy[i - 1].append(gps.getAzymuth())
    #         except TypeError:
    #             pass
    # polar_cord = zip(azy, elev)
    # print(len(elev))
    # print(len(azy))
    rc('grid', color='gray', linewidth=1, linestyle='--')
    fontsize = 15
    rc('xtick', labelsize=fontsize)
    rc('ytick', labelsize=fontsize)
    rc('font', size=fontsize)

    green = '#467821'
    blue = '#348ABD'
    red = '#A60628'
    orange = '#E24A33'
    purple = '#7A68A6'

    fig = plt.figure(figsize=(12, 6))
    plt.subplots_adjust(bottom=0.1,
                        top=0.85,
                        left=0.1,
                        right=0.74)
    ax = fig.add_subplot(polar=True)  # define a polar type of coordinates
    ax.set_theta_zero_location('N')  # ustawienie kierunku p????nocy na g??rze wykresu
    ax.set_theta_direction(-1)  # ustawienie kierunku przyrostu azymutu w prawo
    PG = 0  # zliczanie satelit??w GPS

    for (PRN, el, az) in sat_positions:
        PG += 1
        ax.annotate(PRN,
                    xy=(np.deg2rad(az[int(len(az) / 2)]), 90 - el[int(len(az) / 2)]),
                    bbox=dict(boxstyle="round", fc=green, alpha=0.5),
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    color='k')

    # for azimuth, elevation in polar_cord:
    #     plt.plot(azimuth, elevation, linewidth=0.5)
    mask = []
    rads = np.arange(0, (2 * np.pi), 0.01)
    for rad in rads:
        mask.append(90-cutoff)
    plt.plot(rads, mask, color='#a60010', label='Mask')
    gps = mpatches.Patch(color=green, label='{:02.0f}  Visible GPS'.format(PG))
    plt.legend(handles=[gps], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    dayGPS = datetime(gpsDay[0], gpsDay[1], gpsDay[2], gpsDay[3], gpsDay[4], gpsDay[5])
    plt.title(f'Polar plot {dayGPS.strftime("%d-%m-%Y %H:%M:%S")}')
    ax.set_yticks(range(0, 90 + 10, 10))  # Define the yticks
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    plt.show()


class Satellite:
    def __init__(self, nrSatellite, selectedDate, file_path):
        self.satelliteDict = {}
        self.nrSatellite = nrSatellite
        self.parametrOfSatellite = []
        self.coordSatellite = np.zeros([0, 3])
        self.selectedDate = selectedDate
        self.file_path = file_path

    def getNrSatellite(self):
        return self.nrSatellite

    def getDate(self):
        return self.selectedDate

    def openFile(self):
        self.emptySatelliteDictionary()
        with open(self.file_path[0], 'r', encoding='utf-8') as openFile:
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



class ElipsoidPoint:
    def __init__(self, args):
        self.fi = np.deg2rad(args[0])
        self.lamb = np.deg2rad(args[1])
        self.h = args[2]

    def getFi(self):
        return self.fi

    def getLambda(self):
        return self.lamb

    def getH(self):
        return self.h


class GPS:
    def __init__(self, point, satellite, cutoff):
        self.point = point
        self.satellite = satellite
        self.cutoff = cutoff

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
        return np.rad2deg(math.atan2(NEU[1], NEU[0])) if np.rad2deg(math.atan2(NEU[1], NEU[0])) > 0 else \
            np.rad2deg(math.atan2(NEU[1], NEU[0])) + 360

    def getElevation(self):
        NEU = self.getVectorNEU()
        if np.rad2deg(math.asin(NEU[2] / (math.sqrt(NEU[0] ** 2 + NEU[1] ** 2 + NEU[2] ** 2)))) >= self.cutoff:
            return np.rad2deg(math.asin(NEU[2] / (math.sqrt(NEU[0] ** 2 + NEU[1] ** 2 + NEU[2] ** 2))))
        else:
            pass

    def getCosZ(self):
        NEU = self.getVectorNEU()
        return NEU[2] / (math.sqrt(NEU[0] ** 2 + NEU[1] ** 2 + NEU[2] ** 2))

    def vectorLengthOfSatellite(self):
        XYZs = self.calcXYZ()
        XYZr = self.satellite.getCoordSatellite()
        lengthVec = math.sqrt((float(XYZr[0]) - float(XYZs[0])) ** 2 + (float(XYZr[1]) - float(XYZs[1])) ** 2 +
                              (float(XYZr[2]) - float(XYZs[2])) ** 2)
        return lengthVec


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.figure = plt.figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        toolbar = NavigationToolbar(self.canvas, self)
        self.date_edit = QtWidgets.QDateTimeEdit(calendarPopup=True)
        self.date_edit.calendarWidget().setLocale(QtCore.QLocale(QtCore.QLocale.English))
        self.date_edit.setDateTime(QtCore.QDateTime.currentDateTime())
        self.date_edit.setFont(QtGui.QFont('Courier', 13))
        self.date_edit.setStyleSheet("background-color: black")
        self.date_edit.setMaximumWidth(200)
        lab_year = QLabel(' ', self)
        lab_month = QLabel('Month: ', self)
        lab_sat = QLabel('Number of satellites:', self)
        lab_date = QLabel('Date of observation:', self)
        lab_year.setMaximumWidth(55)
        lab_month.setMaximumWidth(70)
        lab_year.setFont(QtGui.QFont('Courier', 13))
        lab_month.setFont(QtGui.QFont('Courier', 13))
        lab_sat.setFont(QtGui.QFont('Courier', 13))
        lab_date.setFont(QtGui.QFont('Courier', 13))
        self.c_year = QComboBox()
        self.c_month = QComboBox()
        self.numb_sat = QComboBox()
        self.c_year.setFont(QtGui.QFont('Courier', 13))
        self.c_month.setFont(QtGui.QFont('Courier', 13))
        self.numb_sat.setFont(QtGui.QFont('Courier', 13))
        self.c_year.setMaximumWidth(80)
        self.c_month.setMaximumWidth(50)
        self.numb_sat.setMaximumWidth(50)
        self.c_year.setStyleSheet("QComboBox{"
                                  "background-color: transparent;\n"
                                  "background-image:url(background.png);}\n"
                                  "QComboBox QAbstractItemView{border: 0px;color:orange}")
        self.c_month.setStyleSheet("QComboBox{color: rgb(255, 149, 0 );\n"
                                   "background-color: transparent;\n"
                                   "background-image:url(background.png);}\n"
                                   "QComboBox QAbstractItemView{border: 0px;color:orange}")
        self.numb_sat.setStyleSheet("QComboBox{background-color: rgb(255, 149, 0); border: 0.5px solid;"
                                    "border-color: white}"
                                    "QComboBox QAbstractItemView{background-color: rgb(255, 149, 0 );"
                                    "border: 0.5px solid; border-color: white}")
        years_list = getYear()
        months_list = getMonth()
        sat_list = getSat()
        self.c_year.addItems(years_list)
        self.c_month.addItems(months_list)
        self.numb_sat.addItems(sat_list)
        self.comb_lay = QHBoxLayout()
        self.comb_lay.setAlignment(QtCore.Qt.AlignLeft)
        self.comb_lay.addWidget(lab_date)
        self.comb_lay.addWidget(self.date_edit)
        self.comb_lay.addWidget(lab_year)
        self.comb_lay.addWidget(lab_sat)
        self.comb_lay.addWidget(self.numb_sat)
        # self.comb_lay.addWidget(lab_year)
        # self.comb_lay.addWidget(self.c_year)
        # self.comb_lay.addWidget(lab_month)
        # self.comb_lay.addWidget(self.c_month)
        # self.comb_lay.addWidget(lab_day)
        # self.comb_lay.addWidget(self.c_day)
        self.latitude = QLineEdit()
        self.longitude = QLineEdit()
        self.height = QLineEdit()
        self.cutoff = QLineEdit()
        self.latitude.setMaximumWidth(200)
        self.longitude.setMaximumWidth(200)
        self.height.setMaximumWidth(200)
        self.cutoff.setMaximumWidth(200)
        self.latitude.setFont(QtGui.QFont('Courier', 13))
        self.longitude.setFont(QtGui.QFont('Courier', 13))
        self.height.setFont(QtGui.QFont('Courier', 13))
        self.cutoff.setFont(QtGui.QFont('Courier', 13))
        self.latitude.setPlaceholderText('From -90 to 90')
        self.longitude.setPlaceholderText('From -180 to 180')
        self.cutoff.setPlaceholderText('From 0 to 90')
        self.latitude.setStyleSheet('background-color: rgb(255, 149, 0 )')
        self.longitude.setStyleSheet('background-color: black')
        self.height.setStyleSheet('background-color: black')
        self.cutoff.setStyleSheet('background-color: rgb(255, 149, 0 )')
        lab1 = QLabel(':Latitude [??]', self)
        lab2 = QLabel(':Longitude [??]', self)
        lab3 = QLabel(':Height [m]', self)
        lab4 = QLabel(':Mask [??]', self)
        lab1.setFont(QtGui.QFont('Courier', 13))
        lab2.setFont(QtGui.QFont('Courier', 13))
        lab3.setFont(QtGui.QFont('Courier', 13))
        lab4.setFont(QtGui.QFont('Courier', 13))
        load_btn = QPushButton('Load almanac', self)
        load2_btn = QPushButton('Download current almanac', self)
        end_btn = QPushButton('&Exit', self)
        rim_btn = QPushButton('Draw elevation', self)
        clear_btn = QPushButton('Draw visible satellites', self)
        skyplot_btn = QPushButton('Draw skyplot', self)
        dop_btn = QPushButton('Draw DOP', self)
        end_btn.resize(end_btn.sizeHint())
        load_btn.setMaximumWidth(700)
        load_btn.setFont(QtGui.QFont('Courier', 13))
        load2_btn.setFont(QtGui.QFont('Courier', 13))
        rim_btn.setFont(QtGui.QFont('Courier', 13))
        clear_btn.setFont(QtGui.QFont('Courier', 13))
        end_btn.setFont(QtGui.QFont('Courier', 13))
        skyplot_btn.setFont(QtGui.QFont('Courier', 13))
        dop_btn.setFont(QtGui.QFont('Courier', 13))
        load_btn.setStyleSheet('QPushButton { background-color: black; color: white; font: bold; '
                               'border: 0.5px solid ; border-color: white}'
                               'QPushButton:pressed { background-color: rgb(255, 149, 0 )}')
        load2_btn.setStyleSheet('QPushButton { background-color: rgb(255, 149, 0 ); color: white; font: bold; '
                                'border: 0.5px solid ; border-color: white}'
                                'QPushButton:pressed { background-color: black}')
        rim_btn.setStyleSheet('QPushButton { background-color: black; color: white; font: bold; '
                              'border: 0.5px solid ; border-color: white}'
                              'QPushButton:pressed { background-color: rgb(255, 149, 0 )}')
        clear_btn.setStyleSheet('QPushButton { background-color: rgb(255, 149, 0 ); color: white; font: bold; '
                                'border: 0.5px solid ; border-color: white}'
                                'QPushButton:pressed { background-color: black}')
        end_btn.setStyleSheet('QPushButton { background-color: rgb(255, 149, 0 ); color: white; font: bold; '
                              'border: 0.5px solid ; border-color: white}'
                              'QPushButton:pressed { background-color: black}')
        skyplot_btn.setStyleSheet('QPushButton { background-color: black; color: white; font: bold; '
                                  'border: 0.5px solid ; border-color: white}'
                                  'QPushButton:pressed { background-color: rgb(255, 149, 0 )}')
        dop_btn.setStyleSheet('QPushButton { background-color: rgb(255, 149, 0 ); color: white; font: bold; '
                              'border: 0.5px solid ; border-color: white}'
                              'QPushButton:pressed { background-color: black}')
        layout = QGridLayout()
        layout.setSpacing(0)
        layout.addLayout(self.comb_lay, 0, 0, 1, 5)
        layout.addWidget(dop_btn, 5, 0, 1, 2)
        layout.addWidget(load_btn, 3, 3, 1, 2)
        layout.addWidget(load2_btn, 3, 0, 1, 2)
        layout.addWidget(rim_btn, 4, 0, 1, 2)
        layout.addWidget(clear_btn, 4, 3, 1, 2)
        layout.addWidget(skyplot_btn, 5, 3, 1, 2)
        layout.addWidget(end_btn, 6, 3, 1, 2 )
        layout.addWidget(self.latitude, 1, 0, 1, 1)
        layout.addWidget(self.longitude, 1, 3, 1, 1)
        layout.addWidget(self.height, 2, 0, 1, 1)
        layout.addWidget(self.cutoff, 2, 3, 1, 1)
        layout.addWidget(lab1, 1, 1, 1, 1)
        layout.addWidget(lab2, 1, 4, 1, 1)
        layout.addWidget(lab3, 2, 1, 1, 1)
        layout.addWidget(lab4, 2, 4, 1, 1)
        # layout.addWidget(self.canvas, 1, 2, 7, 1)
        layout.setColumnStretch(2, 1)
        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        end_btn.clicked.connect(self.end)
        load2_btn.clicked.connect(linkEvent)
        load_btn.clicked.connect(self.work)
        rim_btn.clicked.connect(self.work)
        clear_btn.clicked.connect(self.work)
        skyplot_btn.clicked.connect(self.work)
        dop_btn.clicked.connect(self.work)
        MainWindow.setFont(self, (QtGui.QFont('Courier', 3)))
        self.setGeometry(0, 0, 800, 400)
        self.setStyleSheet('background-color: grey; color: white; font: bold')
        self.setCentralWidget(widget)
        self.setWindowIcon(QIcon('png-transparent-computer-icons-satellite-receiver-television-angle-logo.png'))
        self.setWindowTitle('SNS - Planning')
        self.show()

    def end(self):
        self.close()

    def closeEvent(self, event):

        odp = QMessageBox.question(
            self, 'Communication',
            'Are you sure to exit?    ',
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if odp == QMessageBox.Yes:
            event.accept()

        else:
            event.ignore()

    def work(self):
        send = self.sender()
        global filename
        try:
            if send.text() == 'Load almanac':
                try:
                    filename = QFileDialog.getOpenFileName()
                except OSError:
                    pass
            coords = [float(self.latitude.text()), float(self.longitude.text()), float(self.height.text())]
            date = self.date_edit.dateTime().toString("yyyyMMddHHmm")
            y = int(date[0] + date[1] + date[2] + date[3])
            m = int(date[4] + date[5])
            d = int(date[6] + date[7])
            h = int(date[8] + date[9])
            min = int(date[10] + date[11])
            sat_numb = int(self.numb_sat.currentText()) + 1
            cutoff = float(self.cutoff.text())
            if send.text() == 'Draw elevation':
                try:
                    plotGPS(coords, y, m, d, filename, cutoff, sat_numb)
                except (NameError, FileNotFoundError):
                    QMessageBox.warning(self, 'Error!', '  Load almanac!         ', QMessageBox.Ok)
            if send.text() == 'Draw visible satellites':
                try:
                    plotVisibleSat(coords, y, m, d, filename, cutoff, sat_numb)
                except (NameError, FileNotFoundError):
                    QMessageBox.warning(self, 'Error!', '  Load almanac!         ', QMessageBox.Ok)
            if send.text() == 'Draw skyplot':
                try:
                    plotSkyplot(coords, y, m, d, h, min, filename, cutoff, sat_numb)
                except (NameError, FileNotFoundError):
                    QMessageBox.warning(self, 'Error!', '  Load almanac!         ', QMessageBox.Ok)
            if send.text() == 'Draw DOP':
                try:
                    plotDops(coords, y, m, d, filename, cutoff, sat_numb)
                except (NameError, FileNotFoundError):
                    QMessageBox.warning(self, 'Error!', '  Load almanac!         ', QMessageBox.Ok)
        except ValueError:
            QMessageBox.warning(self, 'Error!', '  Enter correct data!         ', QMessageBox.Ok)


app = QtWidgets.QApplication(sys.argv)
w = MainWindow()
app.exec_()

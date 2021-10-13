import math
import numpy as np
from datetime import date

alm = open("current1.txt", 'r')
alm_lines = alm.readlines()
alm_data = []
for idx, value in enumerate(alm_lines):
    if value[:3] == 'ID:':
        sat_block = alm_lines[idx:idx+13]
        sat_data = []
        for row in sat_block:
            data = row.split(':')
            sat_data.append(float(data[1]))
        alm_data.append(sat_data)

svprn = alm_data[0][0]
health = alm_data[0][1]
e = alm_data[0][2]
toA = alm_data[0][3]
i = alm_data[0][4]
Omega = alm_data[0][5]
sqrtA = alm_data[0][6]
Omega0 = alm_data[0][7]
w = alm_data[0][8]
M0 = alm_data[0][9]
af0 = alm_data[0][10]
af1 = alm_data[0][11]
gps_week = alm_data[0][12]

data = [2021, 3, 1, 0, 0, 0]

dday = date.toordinal(date(data[0], data[1], data[2])) - date.toordinal(date(1980, 1, 6))
week = dday // 7
day = dday % 7
week = week - 2048
tow = day * 86400 + data[3] * 3600 + data[4] * 60 + data[5]
time = week * 604800 + tow

toa = gps_week * 604800 + toA
tK = time - toa

print(tK)

u = 3.986004415 * 10 ** 14
wE = 7.2921151467 * 10**(-5)

a = sqrtA**2
n = math.sqrt(u / a ** 3)

print(n)

Mk = M0 + n * tK
Mk = Mk % math.pi

print(Mk)

E_old = Mk
E_new = Mk + e * math.sin(E_old)
while abs(E_new - E_old) > 10**(-15):
    E_old = E_new
    E_new = Mk + e * math.sin(E_old)

print(E_new)

Vk = math.atan2(math.sqrt(1-e**2)*math.sin(E_new), math.cos(E_new)-e)
if Vk < 0:
    Vk += 360

print(Vk)

w = alm_data[0][8]
FiK = Vk + w

print(FiK)

rk = a * (1 - e * math.cos(E_new))

print(rk)

xk = rk * math.cos(FiK)
yk = rk * math.sin(FiK)

print(xk, yk)

OmegaK = Omega0 + (Omega - wE) * tK - wE * toA

print(OmegaK)

Xk = xk * math.cos(OmegaK) - yk * math.cos(i) * math.sin(OmegaK)
Yk = xk * math.sin(OmegaK) + yk * math.cos(i) * math.cos(OmegaK)
Zk = yk * math.sin(i)

print("Współrzędne satelity: ", Xk, Yk, Zk)

phi = np.deg2rad(52)
lam = np.deg2rad(21)
h = 100
a_ = 6378137
e2 = 0.00669438002290
N = a_ / math.sqrt(1 - e2 * math.sin(phi) ** 2)

Xr = (N + h) * math.cos(phi) * math.cos(lam)
Yr = (N + h) * math.cos(phi) * math.sin(lam)
Zr = (N * (1 - e2) + h) * math.sin(phi)

XYZr = np.array([Xr, Yr, Zr])
XYZo = np.array([Xk - Xr, Yk - Yr, Zk - Zr])
XYZo.transpose()
print(XYZo)

u_neu = np.array([math.cos(phi)*math.cos(lam), math.cos(phi)*math.sin(lam), math.sin(phi)])
n_neu = np.array([-(math.sin(phi)*math.cos(lam)), -(math.sin(phi)*math.sin(lam)), math.cos(phi)])
e_neu = np.array([-(math.sin(lam)), math.cos(lam), 0.0])
neu_v = np.array([[n_neu], [e_neu], [u_neu]])
print("Wektor NEU: ", neu_v)

neu = np.dot(neu_v, XYZo)
print("Współrzędne NEU: ", neu)

Az = math.atan2(neu[1][0], neu[0][0])
Az = np.rad2deg(Az)
if Az < 0 or Az > 360:
    Az += 360
print("Azymut wynosi: ", Az)

z = math.acos(neu[2][0]/math.sqrt(neu[0][0]**2 + neu[1][0]**2 + neu[2][0]**2))
z = np.rad2deg(z)
print("Kąt zenitalny wynosi: ", z)

el = math.asin(neu[2][0]/math.sqrt(neu[0][0]**2 + neu[1][0]**2 + neu[2][0]**2))
el = np.rad2deg(el)
print("Elewacja wynosi: ", el)

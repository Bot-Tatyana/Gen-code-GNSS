from scipy import integrate
from pylab import*
import numpy as np
import matplotlib.pyplot as plt

                                        #ГЕНЕРАЦИЯ СА КОДА
code_length = 1023
def get_code(code_length): #генерация СА кода
    a0 = []
    b0 = []
    a = [1] * 10
    b = [1] * 10
    CA = np.ones(code_length)
    for j in range(code_length):
        a0 = a[9] #выводим
        a1 = (a[2]+a[9])%2 #считаем
        a = a[-1:] + a[0:9] #сдвигаем
        a[0] = a1 #вставляем
        b0 = (b[1]+b[7])%2 #считаем и выводим для PRN7 спутника
        b1 = (b[5]+b[7]+b[8]+b[9])%2 #считаем
        b2 = (b[1]+b[2] + b1)%2 #считаем
        b = b[-1:] + b[0:9] #сдвигаем
        b[0] = b2 #вставляем
        CA[j] = (a0+b0)%2
        CA[np.where(CA < 1)] = -1
    return CA
CA = get_code(code_length) #СА код

                                        #ПОСТОЯННЫЕ
f_0 = 10.23*1e6 #носимальная частота
fCA = f_0/10 #частота СА
f1 = f_0*154 #частота L1
f2 = f_0*120 #частота L2
fs = 1.023*1e11 #
v = 3888.88888 #скорость спутника м/с
c = 3*1e8 #скорость света

digital_signal_S = []
time_S = []

time =[] #ось времени
digital_signal =[] #СА сод
PSK =[] #модулированнный сигнал
carrier_signal =[] #опорный сигнал
multiplier_signal =[] #перемноженный сигнал
phase_time = []

"""def get_signal_S(CA, ti, fc): #сигнал со спутника без физ. возд.
    fc = fc #частота несущей
    t_bit = 1/fCA #длительность бита
    ibit = int(ti/t_bit) #номер бита
    CA_bit_S = CA[ibit%len(CA)] #СА сигнал
    return CA_bit_S"""

'''for ti_S in arange(0, 5/fCA, 1/fs): #сигнал со спутника без физ. возд.
    CA_bit_S = get_signal_S(CA, ti_S, f1)
    time_S.append(ti_S*fCA) #время c масштабом 10^-6
    digital_signal_S.append(CA_bit_S) #СА сигнал'''


                                        #СИГНАЛ ПРИШЕДШИЙ СО СПУТНИКА
def get_signal(ti, CA_code, fc, elevation = 60*pi/180): #сигнал со спутника c физ. возд. модулированный
    dop = 2*v*np.cos(elevation)/c
    fc = fc + fc*dop #частота несущей
    t_bit = 1/(fCA + fCA*dop) #длительность бита
    ibit = int(ti/t_bit) #номер бита
    CA_bit = CA[ibit%len(CA)] #СА сигнал
    signal_mod = sin(2*pi*fc*ti-pi*(CA_bit+1)/2) #модулированный сигнал
    #carrier = sin(2*pi*fc*ti) #опорный сигнал
    return signal_mod

def get_signal_CA(ti, CA_code, fc, elevation = 60*pi/180): #сигнал со спутника c физ. возд. CA
    dop = 2*v*np.cos(elevation)/c
    fc = fc + fc*dop #частота несущей
    t_bit = 1/(fCA + fCA*dop) #длительность бита
    ibit = int(ti/t_bit) #номер бита
    CA_bit = CA[ibit%len(CA)] #СА сигнал
    return CA_bit


"""for ti in arange(0, 5/fCA, 1/fs): #сигнал со спутника без физ. возд.
    signal_mod, CA_bit = get_signal(CA, ti, f1)
    time.append(ti*fCA) #время c масштабом 10^-6
    digital_signal.append(CA_bit) #СА сигнал
    PSK.append(signal_mod) #модулированный сигнал"""

                                        #ДЕМОДУЛЯЦИЯ СИГНАЛА ПРИШЕДШЕГО СО СПУТНИКА

signal_mod_all = []
x_i = []
y_i = []
ca_bit_phase_all = []
ca_bit = []


def detector_phase(get_signal, CA,fc, n):
    position_0 = 0
    for i in range(1, n): #для области интегрирования около примерной смены фазы
        wch = 100 # кол-во периодов сигнала wch = 1/fc * 100 
        chunck_i_position = position_0 + 1/fCA*i - 1/fc*wch/2 #начало chunck
        chunck_length = 1/fc*wch #длина chunck
        length_window_all = 1/fc*40 #длина окон интегрирования
        step = chunck_length/(wch*10) #шаг окон интегрирования
        for j in range(0, int((chunck_length - length_window_all)/step)): #область chanck
            s = chunck_i_position + step*j - (1/fc)/2 #начало 1 окна
            m = chunck_i_position + step*j + length_window_all/2 #конец 1, начало 2 окон
            f = chunck_i_position + step*j + length_window_all + (1/fc)/2 #конец 2 окна
            integ_window_1 = integrate.quad(get_signal, s, m, args=(CA,fc)) #интегрирование 1 окна от s до m
            integ_window_2 = integrate.quad(get_signal, m, f, args=(CA,fc)) #интегрирование 2 окна от m до f
            #print(integ_window_1[0] - integ_window_2[0])
            if integ_window_1[0] * integ_window_2[0] >=0: #условие для фазы
                if integ_window_1[0] > 0 and integ_window_2[0] > 0: #условие 1 или -1
                    ca_bit_phase = 1
                elif integ_window_1[0] < 0 and integ_window_2[0] < 0:
                    ca_bit_phase = -1
                ca_bit_phase_all.append(ca_bit_phase) #масив из 1 и -1
                x_i.append(f - (1/fc)/4)
                y_i.append(get_signal(f - (1/fc)/4, CA, fc))
                break
    position_0 = f
    return  x_i, y_i, ca_bit_phase_all

n = 10
for ti in arange(0, n/fCA, 1/fs):
    signal_mod_all.append(get_signal(ti, CA, f1))
    time.append(ti)
    digital_signal_S.append(get_signal_CA(ti, CA, f1)) #СА сигнал

x_i, y_i, ca_bit_phase_all = detector_phase(get_signal, CA, f1, n+1)

'''def get_CA(x_i, ca_bit_phase_all, fCA):
    for i in range(0, len(ca_bit_phase_all)):
        delta_t = x_i[i + 1] - x_i[i]
        quantity_bit = int(delta_t/(1/fCA))
        ca_bit += sorted(ca_bit_phase_all[i]*quantity_bit)
    return ca_bit
print(get_CA(x_i, ca_bit_phase_all, fCA))'''

print(x_i)

fig = plt.figure()
ax_1 = fig.add_subplot(2, 1, 1) #для сигнала на спутнике
ax_1.plot(time, signal_mod_all)
ax_1.scatter(x_i, y_i,  c='r')
ax_2 = fig.add_subplot(2, 1,2) #для сигнала CA
ax_2.plot(time, digital_signal_S)
#ax_2.scatter(time_X, signal_Y, c='r')
show()
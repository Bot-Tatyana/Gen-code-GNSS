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
fCA = 1.023*1e6 #частота СА
f1 = 1575.42*1e6 #частота L1
f2 = 1227.6*1e6 #частота L2
fs = 1.023*1e11 #
v = 3888.88888 #скорость спутника м/с
c = 3*1e8 #скорость света

time =[] #ось времени
digital_signal =[] #СА сод
PSK =[] #модулированнный сигнал
carrier_signal =[] #опорный сигнал
multiplier_signal =[] #перемноженный сигнал
phase_time = []

                                        #ГЕНЕРАЦИЯ СИГНАЛА ПРИШЕДШЕГО СО СПУТНИКА
def get_signal(CA, ti, fc, elevation = pi/2):
    dop = 2*v*np.cos(elevation)/c
    fc = fc + fc*dop #частота несущей
    t_bit = 1/(fCA + fCA*dop) #длительность бита
    ibit = int(ti/t_bit) #номер бита
    CA_bit = CA[ibit%len(CA)] #СА сигнал
    signal_mod = sin(2*pi*fc*ti-pi*(CA_bit+1)/2) #модулированный сигнал
    carrier = sin(2*pi*fc*ti) #опорный сигнал
    return carrier, signal_mod, CA_bit

for ti in arange(0, 10/fCA, 1/fs):
    carrier, signal_mod, CA_bit = get_signal(CA, ti, f1)
    time.append(ti*fCA) #время c масштабом 10^-6
    carrier_signal.append(carrier) #опорный сигнал
    digital_signal.append(CA_bit) #СА сигнал
    PSK.append(signal_mod) #модулированный сигнал
    multiplier_signal.append(-carrier*signal_mod) #для демодуляции

                                        #ДЕМОДУЛЯЦИЯ СИГНАЛА ПРИШЕДШЕГО СО СПУТНИКА
integ_res = [] #
mul_sig_range = [] #
index_phase_res = [] #
signal_phase = [] #
signal_Y = [] #
time_X = [] #
signal_phase_res = [] #
d = int #
duration_CA = int((1/fCA)/(1/fs)) #длина массива в 1 бите
step = int(duration_CA/20) #шаг интегрирования  2, 4, 8, 12, 24 ...""" #шаги интегрирования

def detector_phase(multiplier_signal):

    a = 0 #
    b = int(len(multiplier_signal)/step) #для всех скачков фазы должен быть
    for i in range(a, b):
        c = int(i*step)
        mul_sig_range = multiplier_signal[c:c+duration_CA] #область интегрирования функции
        d = len(mul_sig_range)
        integ_res_i = integrate.simps(mul_sig_range) #интегрирование
        integ_res.append(integ_res_i) #результат игтегрирования
        if integ_res_i >= -1.0 and integ_res_i <= 1.0:
            index_num_phase = c + int(duration_CA/2)
            index_phase_res.append(index_num_phase) #массив с индексами смены фазы по сигналу (общее)
            signal_phase.append(multiplier_signal[index_num_phase]) #значеие сигнала в смене фазы (не обязательно)
            signal_Y.append(PSK[index_num_phase])
            time_X.append(time[index_num_phase])
            signal_phase_res.append(signal_phase)
    return d, integ_res, index_phase_res, signal_Y, time_X
    #return d, integ_res

d, integ_res, index_phase_res, signal_Y, time_X = detector_phase(multiplier_signal)
#d, integ_res = detector_phase(multiplier_signal)

"""def integ_step(): # результат интегрирования одного шага
    time_step = time[0:step] #
    mul_step_sig = multiplier_signal[0:step]
    res_integ_step = integrate.simps(mul_step_sig, time_step)
    return res_integ_step"""

print(time_X)

fig = plt.figure()
ax_1 = fig.add_subplot(2, 1, 1)
ax_2 = fig.add_subplot(2, 1, 2)
ax_1.plot(time, digital_signal)
ax_2.plot(time, PSK)
ax_2.scatter(time_X, signal_Y, c='r')
plt.show()

"""f, ax = plt.subplots(2, 1)
ax[0].plot(time, digital_signal)
ax[1].plot(time, PSK)
ax.scatter(signal_Y, time_X, c = 'r')
ax[0].axis([0, 4, -1.5, 1.5])
ax[1].axis([0, 4, -1.5, 1.5])
plt.show()"""
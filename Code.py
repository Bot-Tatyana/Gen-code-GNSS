from scipy import integrate
from pylab import*
import numpy as np
import matplotlib.pyplot as plt

                                        #ПОСТОЯННЫЕ
f_0 = 10.23*1e6 #номинальная частота
fCA = f_0/10 #частота СА
f1 = f_0*154 #частота L1
f2 = f_0*120 #частота L2
fs = 1.0*1e11 #
v = 3888.88888 #скорость спутника м/с
c = 3*1e8 #скорость света

#ГЕНЕРАЦИЯ СИГНАЛА
#генерация СА кода
def get_code(code_length, bits): 
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
        b0 = (b[bits[0]]+b[bits[1]])%2 #считаем и выводим для PRN7 спутника
        b1 = (b[5]+b[7]+b[8]+b[9])%2 #считаем
        b2 = (b[1]+b[2] + b1)%2 #считаем
        b = b[-1:] + b[0:9] #сдвигаем
        b[0] = b2 #вставляем
        CA[j] = (a0+b0)%2
        CA[np.where(CA < 1)] = -1
    return CA

#TODO add PRN                                            
def get_code_prn(code_length, prn):
    bits = (2-1,6-1), (3-1,7-1), (4-1,8-1), (5-1,9-1), (1-1,9-1), (2-1,10-1), (1-1,8-1), (2-1,9-1), (3-1,10-1), (2-1,3-1), (3-1,4-1), (5-1,6-1), (-1,7-1), (7-1,8-1), (8-1,9-1), (9-1,10-1), (1-1,4-1), (2-1,5-1), (3-1,6-1), (4-1,7-1), (5-1,8-1), (6-1,9-1), (1-1,3-1), (4-1,6-1), (5-1,7-1), (6-1,8-1), (7-1,9-1), (8-1,10-1), (1-1,6-1), (2-1,7-1), (3-1,8-1), (4-1,9-1) # make for prh
    return get_code(code_length, bits[prn-1]) #возвращает СА код

#угол возвышения спутника
def get_elevation(): #prn, ti
    #MAKE ccalculations
    return 20*np.pi/180 #возвращает угол в радианах

#модулированный сигнал для интегрирования
def get_signal(ti, CA, fc): #сигнал со спутника модулированный для интегрирования
    elevation = get_elevation()
    dop = 2*v*np.cos(elevation)/c
    fc = fc + fc*dop #частота несущей
    t_bit = 1/(fCA + fCA*dop) #длительность бита
    ibit = int(ti/t_bit) #номер бита
    CA_bit = CA[ibit%len(CA)] #СА сигнал
    signal_mod = np.sin(2*np.pi*fc*ti-np.pi*(CA_bit+1)/2) #модулированный сигнал
    return signal_mod

#ДЕМОДУЛЯЦИЯ
def detector_phase(code_length, prn, fc, wch, wlen, step1): # wch=6, wlen=3, step=1/10
    CA = get_code_prn(code_length + 10, prn)
    """
    wch - кол-во периодов несущего сигнала                            
    """
    t = []
    ca_bit_phase_t = []
    position_0 = 0
    for i in range(1, code_length + 10): #для области интегрирования около примерной смены фазы
        chunck_i_position = position_0 + 1/fCA*i - 1/fc*wch/2 #начало области
        chunck_length = 1/fc*wch #длина бласти                                              
        length_window_all = 1/fc*wlen #длина окон интегрирования                             
        step = (1/fc)* step1 #шаг окон интегрирования                                                            
        for j in range(0, int((chunck_length - length_window_all)/step)): #область окон интегрирования
            s = chunck_i_position + step*j #начало 1 окна
            m = chunck_i_position + step*j + length_window_all/2 #конец 1, начало 2 окон
            f = chunck_i_position + step*j + length_window_all #конец 2 окна
            integ_window_1 = integrate.quad(get_signal, s, m, args=(CA,fc)) #интегрирование 1 окна от s до m
            integ_window_2 = integrate.quad(get_signal, m, f, args=(CA,fc)) #интегрирование 2 окна от m до f
            if integ_window_1[0]*integ_window_2[0] > 0: #условие для фазы

                if integ_window_1[0] > 0 and integ_window_2[0] > 0: #условие для знака бита (1 или -1)
                    ca_bit_phase = 1
                elif integ_window_1[0] < 0 and integ_window_2[0] < 0:
                    ca_bit_phase = -1
                
                ca_bit_phase_t.append(ca_bit_phase) #создание массива битов (1 и -1)
                t.append(f - (1/fc)/2) #время, где меняется фаза
                position_0 = f - (1/fCA)*i - (1/fc)/2 #сдвиг начальной позиции
                break
    return  t, ca_bit_phase_t

#восстановление кода
def restore_code(code_length, prn, fc, wch, wlen, step1):
    t, ca_bit_phase_t = detector_phase(code_length, prn, fc, wch, wlen, step1)
    quantity_bit_all = []
    t_i = []
    ca_bit = [] #востановленный код
    def get_CA_bit(t, fCA): #функция определения кол-во битов между фазами
        t = [0] + t #добавление начала времени
        for j in range(0, len(t)-1):
            delta_t = t[j + 1] - t[j] #время между фазами
            quantity_bit = round(delta_t/(1/fCA)) #кол-во битов между фазами
            quantity_bit_all.append(quantity_bit) #массив кол-ва битов между фазами
        return quantity_bit_all
    quantity_bit_all = get_CA_bit(t, fCA) #вызываем функцию определения кол-во битов
    for k in range(0, len(ca_bit_phase_t)-1):
        ca_bit += sorted([ca_bit_phase_t[k]]*quantity_bit_all[k])
    for i in range(0, len(t)-1):
        if t[i] < code_length*1/fCA + 6/fc:
            t_i.append(t[i]) 
    ca_bit = ca_bit[0:code_length]
    return t_i, ca_bit

#битовая ошибка
def error_CA(code_length, prn, fc, wch, wlen, step1):
    t_i, ca_bit = restore_code(code_length, prn, fc, wch, wlen, step1)
    CA = get_code_prn(code_length, prn)
    error_bit = 0
    for i in range(0, len(ca_bit)):
        if CA[i] != ca_bit[i]:
            error_bit +=1
    return error_bit, CA, ca_bit

def error_phase(code_length, prn, fc, wch, wlen, step1):
    #эталонные значения
    def t_bit_signal(fc): 
        elevation = get_elevation()
        dop = 2*v*np.cos(elevation)/c
        fd = fc + fc*dop #частота несущей с доблером
        t_bit = 1/(fCA + fCA*dop) #длительность бита
        return t_bit, fd

    #эталонные значения времени, где произошел сдвиг фаз
    def get_phase(code_length, prn, fc):
        t_bit, fd = t_bit_signal(fc)
        CA = get_code_prn(code_length, prn)
        TCAd = []
        for i in range(0, code_length):
            b = t_bit*i
            a = b - 1/fd/2
            c = b + 1/fd/2
            integ_window_1 = integrate.quad(get_signal, a, b, args=(CA,fc)) #интегрирование 1 окна от a до b
            integ_window_2 = integrate.quad(get_signal, b, c, args=(CA,fc)) #интегрирование 2 окна от b до c
            if integ_window_1[0]*integ_window_2[0] > 0:
                if integ_window_1[0] < 0 and integ_window_2[0] <0:
                    Td = b
                if integ_window_1[0] > 0 and integ_window_2[0] > 0:
                    Td = b
                TCAd.append(Td)
        return TCAd

    #доверительгый итервал для фаз
    def span_phase(code_length, prn, fc):
        span_left = []
        span_right = []
        TCAd = get_phase(code_length, prn, fc)
        t_bit, fd = t_bit_signal(fc)
        for e in range(0, len(TCAd)-1):
            span_left.append(TCAd[e] - 1/fd)
            span_right.append(TCAd[e] + 1/fd)
        return span_left, span_right, TCAd

    span_left, span_right, TCAd = span_phase(code_length, prn, fc)
    t_i, ca_bit = restore_code(code_length, prn, fc, wch, wlen, step1)
    phase = 0
    for j in range(0, len(span_right)-1):
        for i in range(0, len(t_i)-1):
            if t_i[i] > span_left[j]  and t_i[i] < span_right[j] :
                phase +=1
    phase = abs(len(t_i) - phase)
    return phase, t_i, TCAd

#phase - кол-во несовпадений сдвиг фаз
#t_i - точки времени сдвигов фаз restored_code
#TCAd - точки времени сдвигов фаз correct_code
#error_bit - кол-во несовпадений бит
#CA - СА код
#ca_bit - восстановленный СА код
#code_length - длина кода
#prn - номер спутника
#fc - несущая частота
#wch - область интегрирования
#wlen - окна интегрировнаия
#step1 - шаг окон интегрирования

code_length = 1023 
prn=19
fc=f1
wch=6
wlen=3
step1=1/10

#ошибки
def get_tested_phase_detector(code_length, prn, fc, wch, wlen, step1):
    error_bit, CA, ca_bit = error_CA(code_length, prn, fc, wch, wlen, step1)
    phase, t_i_len, TCAd = error_phase(code_length, prn, fc, wch, wlen, step1)
    TCAd_len = len(TCAd)
    return {#'correct_code': CA,
            #'restored_code': ca_bit,
            'error_bit': round(error_bit/code_length, 2),
            'error_phase': round(phase/TCAd_len, 2),
            'up_errors': 0,
            'down_errors': 0,
            'nochange_errors': 0,
            }

print(get_tested_phase_detector(code_length, prn, fc, wch, wlen, step1))
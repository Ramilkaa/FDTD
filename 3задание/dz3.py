# -*- coding: utf-8 -*-

import numpy
import math
from numpy.fft import fft, fftshift
import matplotlib.pyplot as plt
import tools

class Sampler:
    def __init__(self, discrete: float):
        self.discrete = discrete

    def sample(self, x: float) -> int:
        return math.floor(x / self.discrete + 0.5)


class GaussianMod:
        
    '''
    Источник, создающий модулированный гауссов импульс
    '''

    def __init__(self, magnitude, dg, wg, Nl, Sc):
        
        '''
        magnitude - максимальное значение в источнике;
        dg - коэффициент, задающий начальную задержку гауссова импульса;
        wg - коэффициент, задающий ширину гауссова импульса;
        Nl - количество отсчетов на длину волны;
        Sc - число Куранта.
        '''
        self.magnitude = magnitude
        self.dg = dg
        self.wg = wg
        self.Nl = Nl
        self.Sc = Sc

    def getE(self, time):
        return self.magnitude * numpy.sin(2 * numpy.pi * self.Sc * time / self.Nl) * numpy.exp(-((time - self.dg) / self.wg) ** 2)

if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

     # Время расчета в секундах
    maxTime_s = 15e-9

    # Скорость света в вакууме
    c = 299792458.0

    # Дискрет по пространству в м
    dx = 5e-4

    # Размер области моделирования в метрах
    maxSize_m = 1.5

    eps1 = 4

      # Скорость распространения волны
    v=c/numpy.sqrt(eps1)

    # Переход к дискретным отсчетам
    # Дискрет по времени
    dt = dx * Sc / v

    sampler_x = Sampler(dx)
    sampler_t = Sampler(dt)

    # Время расчета в отсчетах
    maxTime = sampler_t.sample(maxTime_s)

    # Размер области моделирования в отсчетах
    maxSize = sampler_x.sample(maxSize_m)

     # Положение источника в метрах
    sourcePos_m = maxSize_m/2
    sourcePos = math.floor(sourcePos_m / dx + 0.5) # Положение источника в отсчетах

    probesPos_m = maxSize_m * 0.75
    # Датчики для регистрации поля
    probesPos = [math.floor( probesPos_m / dx + 0.5)]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]
    
    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)
    eps[:] = eps1

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize - 1)
    
    source = GaussianMod(1, 100.0, 50.0, 60.0, Sc)
    # Ez[1] в предыдущий момент времени
    oldEzLeft = Ez[1]

    # Ez[-2] в предыдущий момент времени
    oldEzRight = Ez[-2]

    # Расчет коэффициентов для граничных условий
    tempLeft = Sc / numpy.sqrt(mu[0] * eps[0])
    koeffABCLeft = (tempLeft - 1) / (tempLeft + 1)

    tempRight = Sc / numpy.sqrt(mu[-1] * eps[-1])
    koeffABCRight = (tempRight - 1) / (tempRight + 1)

    # Параметры отображения поля E
    display_field = Ez
    display_ylabel = 'Ez, В/м'
    display_ymin = -1.1
    display_ymax = 1.1

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве

    display = tools.AnimateFieldDisplay(dx, dt,
                                        maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel)

    display.activate()
    display.drawProbes(probesPos)
    display.drawSources([sourcePos])

    for q in range(maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= Sc / (W0 * mu[sourcePos - 1]) * source.getE(q-0.5)
        #Hy[sourcePos - 1] += (-1/W0)*source.getE(q-0.5)

        # Расчет компоненты поля E
        Hy_shift = Hy[: -1]
        Ez[1:-1] = Ez[1: -1] + (Hy[1:] - Hy_shift) * Sc * W0 / eps[1: -1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += (Sc / (numpy.sqrt(eps[sourcePos] * mu[sourcePos])) *
                          source.getE(q + 1))

        # Граничные условия ABC первой степени
        Ez[0] = oldEzLeft + koeffABCLeft * (Ez[1] - Ez[0])
        oldEzLeft = Ez[1]

        Ez[-1] = oldEzRight + koeffABCRight * (Ez[-2] - Ez[-1])
        oldEzRight = Ez[-2]

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 50 == 0:
            display.updateData(display_field, q)

    display.stop()

    EzSpec = fftshift(numpy.abs(fft(probe.E)))
df = 1.0 / (maxTime * dt)
freq = numpy.arange(-maxTime / 2 , maxTime / 2 , 1)*df
tlist = numpy.arange(0, maxTime * dt, dt)

# Вывод сигнала и спектра
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.set_xlim(3e-9, 8e-9)
ax1.set_xlabel('t, с')
ax1.set_ylabel('Ez, В/м')
ax1.plot(tlist, probe.E/numpy.max(probe.E))
ax1.minorticks_on()
ax1.grid()
ax2.set_xlim(-20e9, 20e9)
ax2.set_ylim(0, 1.1)
ax2.set_xlabel('f, Гц')
ax2.set_ylabel('|S| / |Smax|, б/р')
ax2.plot(freq, EzSpec / numpy.max(EzSpec))
ax2.minorticks_on()
ax2.grid()
plt.show()

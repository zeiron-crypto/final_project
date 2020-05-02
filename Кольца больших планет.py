# Кольца больших планет

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from scipy.integrate import odeint

def rings_generator(xc=0, yc=0, r=0, N=0, V=0, vxc=0, vyc=0):
    """ Функция распределения частиц в кольце объекта, с параметрами:
        xc - координата "х" центра кольца
        yc - координата "у" центра кольца
        vxc - "х" компонента скорости центрального объекта
        vyc - "у" компонента скорости центрального объекта
        r - радиус кольца, относительно центра
        N - количество точек в кольце
        V - скорость движения частиц в кольце, направление скорости касательно
            к радиусу, проведенному к точке кольца из центра (хс, ус)
    """
    # Создание таблицы с параметрами координат и скоростей объектов
    coordinate = np.ndarray(shape=(N+1, 2))
    velocity = np.ndarray(shape=(N+1, 2))

    # Создание цикла для определения координат и скоростей к заданному объекту
    for i in range(0, N+1, 1):

        # Определение ула для объектов
        alpha = 2 * np.pi / (N + 1) * i

        # Определение координат для первого кольца
        x = xc + r * np.sin(alpha)
        y = yc + r * np.cos(alpha)

        # Подставление координат в таблицу
        coordinate[i, 0] = x
        coordinate[i, 1] = y

        # Определение проекций скоростей в разных четвертях тригонометрической окружности
        # Первая четверть
        if x>=0 and y>=0:
            Vx = vxc - V * np.cos(alpha)
            Vy = vyc + V * np.sin(alpha)

        # Вторая четверть
        elif x<0 and y>=0:
            Vx = vxc - V * np.sin(alpha - 3 * np.pi / 2)
            Vy = vyc - V * np.cos(alpha - 3 * np.pi / 2)

        # Третья четверть
        elif x<0 and y<0:
            Vx = vxc - V * np.cos(alpha)
            Vy = vyc + V * np.sin(alpha)

        # Четвёртая четверть
        elif x>=0 and y<0:
            Vx = vxc + V * np.sin(alpha - np.pi / 2)
            Vy = vyc + V * np.cos(alpha - np.pi / 2)
        else:
            print('error')

        # Подставление проекций скоростей  в  таблицу
        velocity[i, 0] = Vx
        velocity[i, 1] = Vy

    return coordinate, velocity

# Определяем переменную величину
seconds_in_year = 365 * 24 * 60 * 60
seconds_in_day = 24 * 60 * 60
years = 2.5
t = np.arange(0, years*seconds_in_year, 5*seconds_in_day)

# Определяем функцию для системы диф. уравнений
def move_func(s, t):
    (x1, v_x1, y1, v_y1,
     x2, v_x2, y2, v_y2,
     x3, v_x3, y3, v_y3,
     x4, v_x4, y4, v_y4,
     x5, v_x5, y5, v_y5,
     x6, v_x6, y6, v_y6) = s

   
    dxdt1 = v_x1
    dv_xdt1 = - G * mc * (x1 - xc) / ((x1 - xc)**2 + (y1 - yc)**2)**1.5
    dydt1 = v_y1
    dv_ydt1 = - G * mc * (y1 - yc) / ((x1 - xc)**2 + (y1 - yc)**2)**1.5

   
    dxdt2 = v_x2
    dv_xdt2 = - G * mc * (x2 - xc) / ((x2 - x1)**2 + (y2 - y1)**2)**1.5

    dydt2 = v_y2
    dv_ydt2 = - G * mc * (y2 - yc) / ((x2 - xc)**2 + (y2 - yc)**2)**1.5


    
    dxdt3 = v_x3
    dv_xdt3 = - G * mc * (x3 - xc) / ((x3 - xc)**2 + (y3 - yc)**2)**1.5

    dydt3 = v_y3
    dv_ydt3 = - G * mc * (y3 - yc) / ((x3 - xc)**2 + (y3 - yc)**2)**1.5
    
   
    
    dxdt4 = v_x4
    dv_xdt4 = - G * mc * (x4 - xc) / ((x4 - xc)**2 + (y4 - yc)**2)**1.5

    dydt4 = v_y4
    dv_ydt4 = - G * mc * (y4 - yc) / ((x4 - xc)**2 + (y4 - yc)**2)**1.5
    
    
    
    dxdt5 = v_x5
    dv_xdt5 = - G * mc * (x5 - xc) / ((x5 - xc)**2 + (y5 - yc)**2)**1.5

    dydt5 = v_y5
    dv_ydt5 = - G * mc * (y5 - yc) / ((x5 - xc)**2 + (y5 - yc)**2)**1.5
    
   
    
    dxdt6 = v_x6
    dv_xdt6 = - G * mc * (x6 - xc) / ((x6 - xc)**2 + (y6 - yc)**2)**1.5

    dydt6 = v_y6
    dv_ydt6 = - G * mc * (y6 - yc) / ((x6 - xc)**2 + (y6 - yc)**2)**1.5

    return (dxdt1, dv_xdt1, dydt1, dv_ydt1,
            dxdt2, dv_xdt2, dydt2, dv_ydt2,
            dxdt3, dv_xdt3, dydt3, dv_ydt3,
            dxdt4, dv_xdt4, dydt4, dv_ydt4,
            dxdt5, dv_xdt5, dydt5, dv_ydt5,
            dxdt6, dv_xdt6, dydt6, dv_ydt6)

# Определяем начальные значения и параметры, входящие в систему диф. уравнений
N = 5
ring = rings_generator(0, 0, 149*10**9, N, 29000)

coor = ring[0]
vel = ring[1]

x10 = coor[0,0]
v_x10 = vel[0,0]
y10 = coor[0,1]
v_y10 = vel[0,1]

x20 = coor[1,0]
v_x20 = vel[1,0]
y20 = coor[1,1]
v_y20 = vel[1,1]

x30 = coor[2,0]
v_x30 = vel[2,0]
y30 = coor[2,1]
v_y30 = vel[2,1]

x40 = coor[3,0]
v_x40 = vel[3,0]
y40 = coor[3,1]
v_y40 = vel[3,1]

x50 = coor[4,0]
v_x50 = vel[4,0]
y50 = coor[4,1]
v_y50 = vel[4,1]

x60 = coor[5,0]
v_x60 = vel[5,0]
y60 = coor[5,1]
v_y60 = vel[5,1]

s0 = (x10, v_x10, y10, v_y10,
      x20, v_x20, y20, v_y20,
      x30, v_x30, y30, v_y30,
      x40, v_x40, y40, v_y40,
      x50, v_x50, y50, v_y50,
      x60, v_x60, y60, v_y60)

mc = 1.9 * 10**(30)
xc = 0
yc = 0

G = 6.67 * 10**(-11)

# Решаем систему диф. уравнений
sol = odeint(move_func, s0, t)

# Строим решение в виде графика и анимируем
fig = plt.figure()

bodys = []

for i in range(0, len(t), 1):
    
    body1, = plt.plot(sol[i, 0], sol[i, 2], 'o', color='y')

    
    body2, = plt.plot(sol[i, 4], sol[i, 6], 'o', color='y')

    
    body3, = plt.plot(sol[i, 8], sol[i, 10], 'o', color='y')
    
    
    body4, = plt.plot(sol[i, 12], sol[i, 14], 'o', color='y')
    
    
    body5, = plt.plot(sol[i, 16], sol[i,18], 'o', color='y')
   
    
    body6, = plt.plot(sol[i, 20], sol[i,22], 'o', color='y')

    bodys.append([ body1, body2, body3, body4, body5, body6])

plt.plot([0], [0], 'o', color='b')

ani = ArtistAnimation(fig, bodys, interval=50)

plt.axis('equal')
ani.save('N_body_G+K2.gif')
plt.show()
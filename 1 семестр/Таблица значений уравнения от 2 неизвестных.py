# 1-(1+x)*y+x/(y^2)-1/(1-y)^2=0
# rosh:

import matplotlib.pyplot as plt

from math import sqrt, cos, pi, tan, atan, sin, asin

eps = 0.0001
aps = {}


def BinSolver_min(x, y):
    f = lambda x, y: (x / (1 / (1 - y) ** 2 + (1 + x) * y - 1)) ** 0.5
    eq1 = f(x, y)
    eq2 = f(x, eq1)
    while abs(eq2 - eq1) / min(abs(eq2), abs(eq1)) > eps:
        eq1 = eq2
        eq2 = f(x, eq2)
    return eq2


def BinSolver_max(x, y):
    f = lambda x, y: 1 - (1 - (1 + x) * y + x / y ** 2) ** (-0.5)
    eq1 = f(x, y)
    eq2 = f(x, eq1)
    while abs(eq2 - eq1) / min(abs(eq2), abs(eq1)) > eps:
        eq1 = eq2
        eq2 = f(x, eq2)
    return eq2


def ROSH_BinSolver_max(r, m, a):
    F = lambda r, m: m / r + 1 / (sqrt(r * r + 1 - 2 * r)) + 1 / 2 * (1 / (1 + m) + r * r * (1 + m) - 2 * r)
    C0 = F(r, m)
    #print(C0)
    f = lambda r, m, a: m / (C0 - (1 / (sqrt(r * r + 1 - 2 * r * cos(a / (180 / pi)))) + 1 / 2 * (1 / (1 + m) + r * r * (1 + m) - 2 * r * cos(a / (180 / pi)))))
    eq1 = f(r, m, a)
    eq2 = f(eq1, m, a)
    while abs(eq2 - eq1) / min(abs(eq2), abs(eq1)) > eps:
        eq1 = eq2
        eq2 = f(eq2, m, a)
    return eq2


def main_lagr(to_pr):
    global aps
    if True:
        m = 1e-9
        ans = BinSolver_min(m, 0.5)
        aps[m] = ans
    for i in range(1, 580):
        j = i / 1000
        ans = BinSolver_min(j, 0.5)
        if to_pr:
            print(str(j).replace('.', ','), str(BinSolver_min(j, 0.5)).replace('.', ','), sep='\t')
        aps[i / 1000] = ans

    for i in range(580, 1001):
        j = i / 1000
        ans = BinSolver_max(j, 0.5)
        aps[i / 1000] = ans
        if to_pr:
            print(str(j).replace('.', ','), str(BinSolver_max(j, 0.5)).replace('.', ','), sep='\t')
    return 0


def ROSH_BinSolver_min(r, m, a):
    F = lambda r, m: m / r + 1 / (sqrt(r * r + 1 - 2 * r)) + 1 / 2 * (1 / (1 + m) + r * r * (1 + m) - 2 * r * (1 + m))
    C0 = F(r, m)
    #print(C0)
    f = lambda r, m, a: m / (C0 - (1 / (sqrt(r * r + 1 - 2 * r * cos(a / (180 / pi)))) + 1 / 2 * (1 / (1 + m) + r * r * (1 + m) - 2 * r * (1 + m) * cos(a / (180 / pi)))))
    eq1 = f(r, m, a)
    eq2 = f(eq1, m, a)
    while abs(eq2 - eq1) / min(abs(eq2), abs(eq1)) > eps:
        eq1 = eq2
        eq2 = f(eq2, m, a)
    return eq2


def ROSH_vert_BinSolver(r, m, a):
    F = lambda r, m: m / r + 1 / (sqrt(r * r + 1 - 2 * r)) + (1 + m) / 2 * (1 / (1 + m) - r) ** 2
    C0 = F(r, m)
    #print(C0)
    f = lambda r, m, a: m / (C0 - (1 / (sqrt(r * r + 1 - 2 * r * cos(a / (180 / pi)))) + (1 + m) / 2 * (1 / (1 + m) - r * cos(a / (180 / pi))) ** 2))
    eq1 = f(r, m, a)
    eq2 = f(eq1, m, a)
    while abs(eq2 - eq1) / min(abs(eq2), abs(eq1)) > eps:
        eq1 = eq2
        eq2 = f(eq2, m, a)
    return eq2


def main_rosh(a, k, m, reverse):
    #print(str(ROSH_BinSolver_max(aps[m] / k, 0.1, a)).replace('.', ','))
    if reverse:
        return ROSH_BinSolver_max((1 - aps[m]) / k, 1 / m, a)
    else:
        return ROSH_BinSolver_max(aps[m] / k, m, a)


def main_rosh_vert(a, k, m, reverse):
    #print(str(ROSH_BinSolver_max(aps[m] / k, 0.1, a)).replace('.', ','))
    if reverse:
        return ROSH_vert_BinSolver((1 - aps[m]) / k, 1 / m, a)
    else:
        return ROSH_vert_BinSolver(aps[m] / k, m, a)



def releaser():
    main_lagr(False)
    for zzz in range(1, 2):
        m = 1e-9

        ax = plt.subplot(111, projection='polar')
        ax.grid(True)
        plt.xticks([])
        plt.yticks([0, 0, 0.5, 1])
        labels = ax.set_yticklabels([], fontsize=0)
        plt.title("m = " + (str(m) + "000")[:5])

        r = [0]     # Не спрашивайте зачем эти 3 строки
        theta = [1]
        ax.plot(theta, r, color='r', linewidth=1)

        #r = [2]  # И про эти тоже
        #theta = [0]
        #ax.plot(theta, r, color='r', linewidth=1)

        k = 1.1
        while k < 4:
            r = []
            for i in range(0, 360):
                rp = main_rosh(i, k, m, False)
                r.append(rp)

            theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

            ax.plot(theta, r, color='b', linewidth=1)
            k *= k ** 0.5 * 1.1

        k = 1
        r = []
        for i in range(0, 360):
            r.append(main_rosh(i, k, m, False))

        theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

        ax.plot(theta, r, color='b', linewidth=3)

        # ogo 1/m

        k = 1
        while k < 4:
            r = []
            theta = []
            for i in range(0, 360):
                rp = main_rosh(i, k, m, True)
                theta_, r_ = converter(i, rp)
                r.append(r_)
                theta.append(theta_)

            ax.plot(theta, r, color='b', linewidth=1)
            k *= k ** 0.5 * 1.1

        k = 1
        r = []
        theta = []
        for i in range(0, 360):
            rp = main_rosh(i, k, m, True)
            theta_, r_ = converter(i, rp)
            r.append(r_)
            theta.append(theta_)

        ax.plot(theta, r, color='b', linewidth=3)

        # --------------------------

        k = 1.1
        while k < 4:
            r = []
            for i in range(0, 360):
                rp = main_rosh_vert(i, k, m, False)
                r.append(rp)

            theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

            ax.plot(theta, r, color='g', linewidth=1)
            k *= k ** 0.5 * 1.1

        k = 1
        r = []
        for i in range(0, 360):
            r.append(main_rosh_vert(i, k, m, False))

        theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

        ax.plot(theta, r, color='g', linewidth=3)

        # ogo 1/m

        k = 1
        while k < 4:
            r = []
            theta = []
            for i in range(0, 360):
                rp = main_rosh_vert(i, k, m, True)
                theta_, r_ = converter(i, rp)
                r.append(r_)
                theta.append(theta_)

            ax.plot(theta, r, color='g', linewidth=1)
            k *= k ** 0.5 * 1.1

        k = 1
        r = []
        theta = []
        for i in range(0, 360):
            rp = main_rosh_vert(i, k, m, True)
            theta_, r_ = converter(i, rp)
            r.append(r_)
            theta.append(theta_)

        ax.plot(theta, r, color='g', linewidth=3)

        # --------------------------

        plt.show()

        # plt.savefig("two_star_compare\\chart{}.png".format(str(m * 1000)))

        plt.close()


def one_star():
    main_lagr(False)
    for zzz in range(1, 32):
        m = zzz ** 2 / 1000
        print(m)

        ax = plt.subplot(111, projection='polar')
        ax.grid(True)
        plt.xticks([])
        plt.yticks([0, 0, 0.5, 1])
        labels = ax.set_yticklabels([], fontsize=0)
        plt.title("m = " + (str(m) + "000")[:5])

        r = [0]  # Не спрашивайте зачем эти 3 строки
        theta = [1]
        ax.plot(theta, r, color='r', linewidth=1)

        k = 1.1
        while k < 4:
            r = []
            for i in range(0, 360):
                rp = main_rosh(i, k, m, False)
                r.append(rp)

            theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

            ax.plot(theta, r, color='b', linewidth=1)
            k *= k ** 0.5 * 1.1

        k = 1
        r = []
        for i in range(0, 360):
            r.append(main_rosh(i, k, m, False))

        theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

        ax.plot(theta, r, color='b', linewidth=3)

        # ------------------------------------

        k = 1.1
        while k < 4:
            r = []
            for i in range(0, 360):
                rp = main_rosh_vert(i, k, m, False)
                r.append(rp)

            theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

            ax.plot(theta, r, color='g', linewidth=1)
            k *= k ** 0.5 * 1.1

        k = 1
        r = []
        for i in range(0, 360):
            r.append(main_rosh_vert(i, k, m, False))

        theta = list(map(lambda x: x / (180 / pi), range(0, 360)))

        ax.plot(theta, r, color='g', linewidth=3)

        # -------------------------

        plt.savefig("one_star_compare\\chart{}.png".format(str(m * 1000)))

        plt.close()


def converter(theta, r):
    theta = theta / (180 / pi)
    # theta_ = atan(tan(theta) * (r * cos(theta)) / (1 - r * cos(theta)))
    r_ = (r * r + 1 - 2 * r * cos(theta)) ** 0.5
    theta_ = r / r_ * sin(theta)
    return theta_, r_


if __name__ == '__main__':
    releaser()
    # one_star()
    # main_lagr()

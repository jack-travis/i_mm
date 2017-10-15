import math

#pyplot needed only for __main__
if __name__ == "__main__":
    from matplotlib import pyplot, rc
    rc('font', size=12)

#unfortunately not easily generalisable
def wind_speeds(p, y_min, y_max, rho, f, N=10):
    """Calculates the geostrophic wind speeds between altitudes $y_min and $y_max.
    Calculation of the pressure gradient uses finite differences, 2nd order.
    :param p: Pressure as a function of altitude
    :type p: function
    :param y_min: Minimum altitude
    :param y_max: Maximum altitude
    :param rho: Mean air density
    :param f: Coriolis parameter
    :param N: Number of points of altitude to calculate at
    :return: Dictionary of altitude and wind speed values
    :rtype: dict
    """
    js = range(N+1)
    D_y = float(y_max - y_min) / N
    ys = [j * D_y + y_min for j in js]
    winds = {}
    for j in js:
        #calculate pressure gradient using finite differences
        y = ys[j]
        d_p = None
        if j == 0:
            d_p = (p(y + D_y) - p(y)) / D_y
        elif j == N:
            d_p = (p(y) - p(y - D_y)) / D_y
        else:
            d_p = (p(y + D_y) - p(y - D_y)) / (2*D_y)
        winds[y] = -(1/(rho*f)) * d_p
    return winds

#calculate u for y in [y_min .. y_max]
L = 2.4e6
pa = 1e5
pb = 200
y_min = 0
y_max = 1e6
rho = 1
f = 1e-4

if __name__ == "__main__":
    N = 10

    p = lambda y : pa + pb*math.cos(y*math.pi/L)

    winds = wind_speeds(p, y_min, y_max, rho, f, N)

    ys = winds.keys()
    wind_calc = winds.values()
    #calculate the wind speed analytically, and thus find the error
    wind_actual = [math.sin(y * math.pi/L) * pb*math.pi / (L*rho*f) for y in ys]
    errors = [abs(c - a) for c,a in zip(wind_calc,wind_actual)]

    #plot values and errors
    pyplot.suptitle("Numerical evaluation of geostrophic wind speed ({0} intervals)".format(N))
    pyplot.subplot(121)
    pyplot.title("Geostrophic wind speed (m s$^{-1}$)")
    pyplot.scatter(ys, wind_calc)
    pyplot.xlabel("Altitude (m)")
    pyplot.subplot(122)
    pyplot.title("Absolute error magnitude")
    pyplot.scatter(ys, errors)
    pyplot.yscale("log")
    pyplot.ylim([10**(int(math.log(min(errors),10))-1),1]) #kludgey
    pyplot.xlabel("Altitude (m)")
    pyplot.subplots_adjust(left=0.05,bottom=0.1,right=0.97,top=0.9,wspace=0.15)
    pyplot.show(block=False)

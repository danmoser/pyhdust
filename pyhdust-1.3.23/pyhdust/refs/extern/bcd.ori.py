import numpy as np
from matplotlib import pyplot as plt
import os


def balmerJump(filename):
    data = np.loadtxt(filename)

    # -- first range:
    w1 = (data[:,1]>=3.5797)*(data[:,1]<=3.6989)*(data[:,2]>=.98)*(data[:,2]<=1.0)

    # -- second range
    w2 = (data[:,1]>=3.53)*(data[:,1]<=3.5658)*(data[:,2]>=.30)*(data[:,2]<=.64)

    # -- linear fit for each range
    c1 = np.polyfit(data[:,1][w1], data[:,3][w1], 1)
    c2 = np.polyfit(data[:,1][w2], data[:,3][w2], 1)

    # -- computing offset at 3700A:
    x0 = np.log10(3700.)
    offset = np.polyval(c1, x0)-np.polyval(c2, x0)
    print 'offset at %4.0fA = %5.3e'%(x0, offset)

    # -- balmer series, using the Rydberg formula
    B = 3645.0682
    m = 2.0 # for Balmer
    n = np.arange(10,18)
    wl = B*(n**2/(n**2-m**2))
    # -- maxima as middle points between minima
    wlMax = 0.5*(wl[1:]+wl[:-1])
    # -- distance between minima:
    wlStep = np.abs(np.diff(wl))
    # -- fit all maxima
    maxima = []
    for i in range(len(wlMax)):
        wl0 = wlMax[i]
        dwl = wlStep[i]/6.

        x = data[:,1][np.abs(data[:,0]-wl0)<dwl]
        y = data[:,3][np.abs(data[:,0]-wl0)<dwl]
        c = np.polyfit(x, y, 2)
        # -- actual position of maximum
        wl1 = -0.5*c[1]/c[0]
        # -- (wl at max, max, poly coef, width of the fit)
        maxima.append((wl1, np.polyval(c, wl1), c, dwl))


    # -- plot data
    plt.figure(0)
    plt.clf()
    plt.plot(data[:,1][data[:,1]>3.5051], data[:,3][data[:,1]>3.5051],
             alpha=0.3, color='k', label=filename)
    plt.plot(data[:,1][w1], data[:,3][w1], '.r', label='range 1', alpha=0.1)
    plt.plot(data[:,1][w2], data[:,3][w2], '.b', label='range 2', alpha=0.1)
    plt.xlabel("log (wavelength A)")
    plt.ylabel("log (flux)")
    plt.title('HD91373')

    # -- plot linear fit:
    x = np.linspace(x0-0.06, x0+0.2, 1000)
    #x = np.linspace(x0-500, x0+1500, 1000)

    plt.plot(x, np.polyval(c1, x), '-r')
    plt.plot(x, np.polyval(c2, x), '-b')
    #plt.legend()

    # -- Balmer Jump
    plt.plot([x0, x0], [np.polyval(c1, x0), np.polyval(c2, x0)], '-g')
    intersect = 0.5*(np.polyval(c1, x)+ np.polyval(c2, x)) # (red line+blue line)/2

    maximaCur = np.interp(x, [m[0] for m in maxima][::-1], [m[1] for m in maxima][::-1])
    plt.plot([m[0] for m in maxima][::-1], [m[1] for m in maxima][::-1],'y')
    #plt.plot(x, maximaCur, '-g')
    wlIntersect = x[np.argmin(np.abs(intersect-maximaCur))]
    print 'intersection = %6.1fA'%(10**wlIntersect)

    # -- position of maxima
    plt.plot(x, intersect, color='orange', linestyle='dashed')

    # -- top of the line polynomial fit:
    for m in maxima:
        x = np.linspace(m[0]-m[3], m[0]+m[3], 1.)
        plt.plot(x, np.polyval(m[2], x), '-', color='orange')
    plt.xlim(3.5, 3.8)
    plt.ylim(-12.9, -11.4)

    return offset, 10**wlIntersect


def analyzeAll():
    filenames = os.listdir('./')
    filenames = filter(lambda x: x.startswith('HD') and not x.endswith('.pdf'), filenames)
    print filenames
    res = {}
    for f in filenames:
        print '*'*5, f, '*'*5
        res[f] = balmerJump(f)
        plt.figure(0)
        #plt.savefig(f+'.png')
        plt.savefig(f+".pdf", format="pdf", dpi=600)
    return res

result = analyzeAll()
print '-'*30
for k in result.keys():
    print k, result[k]
    #np.savetxt('result',(k))

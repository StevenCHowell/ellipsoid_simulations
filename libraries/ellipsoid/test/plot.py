import matplotlib.pylab as plt
import locale
import sys

X,Y=[],[]
for line  in open(sys.argv[1]).readlines():
    words=line.split()
    X.append(locale.atof(words[0]))
    Y.append(locale.atof(words[1]))

plt.plot(X,Y)
plt.xlim([3,16])
plt.ylim([-1.2,0.01])
plt.show()

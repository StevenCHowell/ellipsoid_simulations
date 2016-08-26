import locale,sys
import matplotlib.pylab as plt

f=sys.argv[1]
Time,BoxLength,Density,Pressure=[],[],[],[]
for line in open(f).readlines()[1:]:
    words=line.split()
    try:
        Time.append(locale.atof(words[0]))
        BoxLength.append(locale.atof(words[1]))
        Density.append(locale.atof(words[2]))
        Pressure.append(locale.atof(words[3]))
    except:
        continue

plt.plot(Time,BoxLength,label='Box_Length')
plt.plot(Time,Density,label='Density')
plt.plot(Time,Pressure,label='Pressure')
plt.legend(loc='best')
plt.show()

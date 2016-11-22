import locale
import sys,os,glob
import numpy
import matplotlib.pylab as plt

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
import matplotlib
matplotlib.rc('font', **font)

def readIgorTxt(f):
    flag = False
    para = ''
    for line in open(f).readlines():
        words = line.split()
        if not len(words): continue
        word = words[0]
        if word=='scale' or word=='Scale':
            flag = True
            continue
        elif word=='SLD':
            flag = False
        if flag:
            para += line
        if 'chisq' in line:
            para += '\n'+line
        if 'Sqrt(X' in line:
            para += ' '.join(line.split()[3:])
    para += '\n\n'
    #return unicode(para)
    return para.decode('latin1')

def readexpfile(f):
    Q,I,E=[],[],[]
    for line in open(f).readlines():
        words=line.split()
        Q.append(locale.atof(words[0]))
        I.append(locale.atof(words[1]))
        E.append(locale.atof(words[2]))
    return numpy.array(Q),numpy.array(I)/I[0],numpy.array(E)

def readmodelfile(f):
    Q,I=[],[]
    for line in open(f).readlines():
        words=line.split()
        Q.append(locale.atof(words[0]))
        I.append(locale.atof(words[1]))
    return numpy.array(Q),numpy.array(I)/I[0]


pdb = os.getcwd().split()[-1].split('_')[-1]
Q,Iexp,Eexp = readexpfile(glob.glob('output/*++.txt')[0])
plt.errorbar(Q,Iexp,yerr=Eexp,color='k', fmt='o',label='Experiment (all atom)\n')

for folder in sys.argv[1:]:
    Q,Imodel = readmodelfile(glob.glob(folder+'/*++.txt')[0])
    para = readIgorTxt(folder+'Report.txt')
    plt.plot(Q,Imodel,label='=========='+os.path.split(folder)[0].split('_')[-1]+'==========\n'+para)


plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.xlabel('Q/cm-1')
plt.xlim([min(Q),max(Q)])
plt.ylabel('I(Q)')
plt.title('IGOR curve fitting to all-atom scattering (exp) of '+pdb)
plt.show()

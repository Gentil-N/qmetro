import numpy as np
from matplotlib import pyplot
from matplotlib import pyplot as plt

def read_file(file_name: str):
    file = open(file_name, 'r')
    lines = file.readlines()
    splitted_time = lines[0].split()
    time = [float(t) for t in splitted_time]
    result_list = []
    line_index = 0
    print("Parsing lines...", end='')
    for line in lines[1:]:
        #print(line_index, "...", end='')
        splitted = line.split()
        #for val in splitted:
            #print(val)
            #print(complex(val))
        result_list.append([complex(val) for val in splitted])
    print("Done")
    return [time, result_list]


#data = read_file("diff_equ_test.dat")
#time = np.linspace(data[0][0], data[0][-1], 1000)
#exact = np.exp(time / 2 - np.sin(2 * time) / 4)

#fig, ax = pyplot.subplots(1, 1)
#ax.plot(data[0], np.real(data[1][0]), label='rk4')
#ax.plot(time, exact, label='exact')
#ax.legend()
#ax.set_xlabel('time')
#ax.set_ylabel('???')

data = read_file("corr_ada_2.dat")
#time = np.linspace(data[0][0], data[0][-1], 1000)
#exact = np.exp(np.array(data[0]) / 2 - np.sin(2 * np.array(data[0])) / 4)

fig, ax = pyplot.subplots(1, 1)
ax.plot(data[0], np.real(data[1][0]), label='ada')
ax.plot(data[0], np.real(data[1][3]), label='sz')
#ax.plot(data[0], exact, label='exact')
ax.legend()
ax.set_xlabel('time')
ax.set_ylabel('???')

plt.show()

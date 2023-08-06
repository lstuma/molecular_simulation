import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.use('TkAgg')
plt.style.use("Solarize_Light2")

def view_test(filepaths):
    for filepath in filepaths:
        with open(filepath, 'r') as f:
            label = f.readline()
            data = f.read().split('\n')[:-1]
            data_x = [float(x.split(',')[0]) for x in data]
            data_y = [float(x.split(',')[1]) for x in data]
            plt.plot(data_x, data_y, label=label)
    if filepaths:
        plt.legend(loc="upper left")
        plt.show()

def scatter_test(filepaths):
    for filepath in filepaths:
        with open(filepath, 'r') as f:
            label = f.readline()
            data = f.read().split('\n')[:-1]
            data_x = [float(x.split(',')[0]) for x in data]
            data_y = [float(x.split(',')[1]) for x in data]
            plt.scatter(data_x, data_y, label=label)
    if filepaths:
        plt.legend(loc="upper left")
        plt.show()

def heatmap_test(filepaths):
    for filepath in filepaths:
        with open(filepath, 'r') as f:
            data = [[float(y) for y in x.split(',')[:-1]] for x in f.read().split('\n')][:-1]
            ax = sns.heatmap(data, linewidth=0, cmap="YlGnBu", vmin=-2e+12, vmax=2e+12)
    if filepaths:
        plt.show()


if __name__ == '__main__':
    # tests visualized
    view_test([
        '/home/lstuma/programming/projects/molecular_simulation/simulations/test0',
    ])
    view_test([
        '/home/lstuma/programming/projects/molecular_simulation/simulations/test1_0',
        '/home/lstuma/programming/projects/molecular_simulation/simulations/test1_1',
        '/home/lstuma/programming/projects/molecular_simulation/simulations/test1_2',
        '/home/lstuma/programming/projects/molecular_simulation/simulations/test1_3',
    ]);
    heatmap_test([
    #    '/home/lstuma/programming/projects/molecular_simulation/simulations/test2',
    ])

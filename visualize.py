import numpy as np
import matplotlib.pyplot as plt


def visualizeVelocity(velocity_file, iteratedTime_file, output_file):
    velocities = np.loadtxt(velocity_file, dtype=float)
    iteratedTime = np.loadtxt(iteratedTime_file, dtype=int)

    length = [i * 0.01 for i in range(101)]
    index = [f"t={10 * i}[s]" for i in range(iteratedTime)]

    fig = plt.figure()
    # ... (ここでの設定は元のコードに合わせています) ...

    for j in range(iteratedTime):
        plt.plot(length, velocities[j])

    plt.xlabel("Length[m]")
    plt.ylabel("Velocity[m/s]")

    plt.legend(
        index,
        bbox_to_anchor=(0.5, -0.1),
        loc="upper center",
        borderaxespad=0.5,
        ncol=4,
        fontsize=10,
    )

    plt.rcParams["font.family"] = "sans-serif"  # 使用するフォント
    # x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams["xtick.direction"] = "in"
    # y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams["ytick.direction"] = "in"
    # plt.rcParams["xtick.major.width"] = 1.0  # x軸主目盛り線の線幅
    # plt.rcParams["ytick.major.width"] = 1.0  # y軸主目盛り線の線幅
    plt.rcParams["font.size"] = 12  # フォントの大きさ
    plt.rcParams["axes.linewidth"] = 1.0  # 軸の線幅edge linewidth。囲みの太さ
    # ↓グラフ位置の調整
    plt.subplots_adjust(
        left=0.2, bottom=0.2, right=None, top=0.97, wspace=None, hspace=None
    )

    fig.savefig(output_file)


def visualizeArea(area_file, iteratedTime_file, output_file):
    areas = np.loadtxt(area_file, dtype=float)
    iteratedTime = np.loadtxt(iteratedTime_file, dtype=int)

    length = [i * 0.01 for i in range(101)]
    index = [f"t={10 * i}[s]" for i in range(iteratedTime)]

    fig = plt.figure()
    # ... (ここでの設定は元のコードに合わせています) ...

    for i in range(iteratedTime):
        plt.plot(length, areas[i])

    plt.xlabel("Length[m]")
    plt.ylabel("Area[m^2]")

    plt.legend(
        index,
        bbox_to_anchor=(0.5, -0.1),
        loc="upper center",
        borderaxespad=0.5,
        ncol=4,
        fontsize=10,
    )

    plt.rcParams["font.family"] = "sans-serif"  # 使用するフォント
    # x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams["xtick.direction"] = "in"
    # y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams["ytick.direction"] = "in"
    # plt.rcParams["xtick.major.width"] = 1.0  # x軸主目盛り線の線幅
    # plt.rcParams["ytick.major.width"] = 1.0  # y軸主目盛り線の線幅
    plt.rcParams["font.size"] = 12  # フォントの大きさ
    plt.rcParams["axes.linewidth"] = 1.0  # 軸の線幅edge linewidth。囲みの太さ
    # ↓グラフ位置の調整
    plt.subplots_adjust(
        left=0.2, bottom=0.2, right=None, top=0.97, wspace=None, hspace=None
    )

    fig.savefig(output_file)


def visualizeFlowQuantity(quantity_file, iteratedTime_file, output_file):
    flowQuantities = np.loadtxt(quantity_file, dtype=float)
    iteratedTime = np.loadtxt(iteratedTime_file, dtype=int)

    length = [i * 0.01 for i in range(101)]
    index = [f"t={10 * i}[s]" for i in range(iteratedTime)]

    fig = plt.figure()
    # ... (ここでの設定は元のコードに合わせています) ...

    for j in range(iteratedTime):
        plt.plot(length, flowQuantities[j])

    plt.xlabel("Length[m]")
    plt.ylabel("Quantity[m^3/s]")

    plt.legend(
        index,
        bbox_to_anchor=(0.5, -0.1),
        loc="upper center",
        borderaxespad=0.5,
        ncol=4,
        fontsize=10,
    )

    plt.rcParams["font.family"] = "sans-serif"  # 使用するフォント
    # x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams["xtick.direction"] = "in"
    # y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams["ytick.direction"] = "in"
    # plt.rcParams["xtick.major.width"] = 1.0  # x軸主目盛り線の線幅
    # plt.rcParams["ytick.major.width"] = 1.0  # y軸主目盛り線の線幅
    plt.rcParams["font.size"] = 12  # フォントの大きさ
    plt.rcParams["axes.linewidth"] = 1.0  # 軸の線幅edge linewidth。囲みの太さ
    # ↓グラフ位置の調整
    plt.subplots_adjust(
        left=0.2, bottom=0.2, right=None, top=0.97, wspace=None, hspace=None
    )

    fig.savefig(output_file)


def visualizePressure(file1, file2, file3, output_file):
    # Load the data
    data1 = np.loadtxt(file1, dtype=float)
    data2 = np.loadtxt(file2, dtype=float)
    data3 = np.loadtxt(file3, dtype=float)

    # Generate x values
    x_values = np.arange(len(data1)) / 1200

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(x_values, data1, linestyle="-", label="Data 1")
    ax.plot(x_values, data2, linestyle="--", label="Data 2")
    ax.plot(x_values, data3, linestyle=":", label="Data 3")

    # Set labels and legend
    ax.set_xlabel("Normalized index")
    ax.set_ylabel("Values")
    ax.legend()

    # Save the figure
    fig.savefig(output_file)


# Execute the functions
visualizeVelocity(
    "output/dat/velocity.dat", "output/dat/iteratedTime.dat", "output/png/velocity.png"
)
visualizeArea(
    "output/dat/area.dat", "output/dat/iteratedTime.dat", "output/png/area.png"
)
visualizeFlowQuantity(
    "output/dat/flowQuantity.dat",
    "output/dat/iteratedTime.dat",
    "output/png/flowQuantity.png",
)
visualizePressure(
    "output/dat/pressure1.dat",
    "output/dat/pressure2.dat",
    "output/dat/pressure3.dat",
    "output/png/pressure.png",
)

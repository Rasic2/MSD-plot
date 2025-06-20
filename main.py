import os
from pathlib import Path

import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read, write
from scipy.stats import linregress


def calculate_diffusion(msd_ana: msd.EinsteinMSD, x_value, y_value):
    # 假设后20%的区间为线性
    N = len(msd_ana.times)
    index_fit = int(N * 0.8)
    slope, intercept, r_value, p_value, std_err = linregress(x_value, y_value)
    D = slope / 6  # 3维的情况
    D_value = str(f'Diffusion Coefficient: {D} Å^2/ps')
    return D_value


def data_XDATCAR(work_dir: Path):
    # 检查所需要的文件
    required_files = ['INCAR', 'XDATCAR']
    missing_files = []
    for file in required_files:
        file_path = work_dir / file
        if file_path.exists() == False:
            missing_files.append(file)
    if missing_files:
        raise FileNotFoundError(f"错误: 缺少必需文件: {', '.join(missing_files)}")
    INCAR = work_dir / 'INCAR'
    XDATCAR = work_dir / 'XDATCAR'
    atoms_list = read(XDATCAR, index=':')
    write('traj.xyz', atoms_list)
    # 获取y_value
    u = mda.Universe('traj.xyz')
    H = u.select_atoms('all')
    msd_ana = msd.EinsteinMSD(H, msd_type='xyz', fft=True)
    msd_ana.run()
    y_value = msd_ana.results.timeseries
    # 获取x_value：
    dt_fs = None
    with open(INCAR) as f:
        lines = f.readlines()
        for line in lines:
            if 'POTIM' in line and '=' in line:
                dt_fs = float(line.split('=')[1].strip().split(' ')[0])
                break
        if dt_fs is None:
            raise ValueError("错误：缺少必要参数：POTIM")
    nframes = msd_ana.n_frames
    x_value = np.arange(nframes) * dt_fs / 1000  # 换为 ps
    D_value = calculate_diffusion(msd_ana, x_value, y_value)
    print(D_value)
    return x_value, y_value


def data_dump(work_dir: Path):
    in_files = list(work_dir.glob('*.in'))
    dump_files = list(work_dir.glob('*dump*'))
    if len(in_files) != 1:
        raise FileNotFoundError(f"缺少in文件或in文件不唯一")
    if len(dump_files) != 1:
        raise FileNotFoundError(f"缺少dump文件或dump文件不唯一")
    in_file = in_files[0]
    dump = dump_files[0]
    # 获取y_value
    u = mda.Universe(str(dump), format='LAMMPSDUMP')
    H = u.select_atoms('all')
    msd_ana = msd.EinsteinMSD(H, msd_type='xyz', fft=True)
    msd_ana.run()
    y_value = msd_ana.results.timeseries
    # 获取x_value
    with open(in_file) as f:
        lines = f.readlines()
        for line in lines:
            if 'timestep' in line:
                print(line)
                dt_fs = float(line.split(' ')[-1])
    nframes = msd_ana.n_frames
    x_value = np.arange(nframes) * dt_fs
    D = calculate_diffusion(msd_ana, x_value, y_value)
    print(D)
    return x_value, y_value


def data_xtc(work_dir: Path):
    tpr_files = list(work_dir.glob('*.tpr'))
    xtc_files = list(work_dir.glob('*.xtc'))
    if len(tpr_files) != 1:
        raise FileNotFoundError(f"缺少.tpr文件或不唯一")
    if len(xtc_files) != 1:
        raise FileNotFoundError(f"缺少.xtc文件或不唯一")
    tpr = tpr_files[0]
    xtc = xtc_files[0]
    # 获取y_value
    u = mda.Universe(str(tpr), str(xtc))
    H = u.select_atoms('all')
    msd_ana = msd.EinsteinMSD(H, msd_type='xyz', fft=True)
    msd_ana.run()
    y_value = msd_ana.results.timeseries
    # 获取x_value
    x_value = msd_ana.times
    D_value = calculate_diffusion(msd_ana, x_value, y_value)
    print(D_value)
    return x_value, y_value


# def data_trr(work_dir: Path):
#     trr_files = list(work_dir.glob('*.trr'))
#     gro_files = list(work_dir.glob('*.gro'))
#     if len(trr_files) != 1:
#         raise FileNotFoundError(f"缺少.gro文件或不唯一")
#     if len(gro_files) != 1:
#         raise FileNotFoundError(f"缺少.trr文件或不唯一")
#     trr = trr_files[0]
#     gro = gro_files[0]
#     # 获取y_value
#     u = mda.Universe(gro, trr, )
#     H = u.select_atoms('all')
#     msd_ana = msd.EinsteinMSD(H, msd_type='xyz', fft=True)
#     msd_ana.run()
#     y_value = msd_ana.results.timeseries
#     # 获取x_value
#     x_value = msd_ana.times
#     calculate_diffusion(msd_ana, x_value, y_value)
#     D = calculate_diffusion(msd_ana, x_value, y_value)
#     print(D)
#     return x_value, y_value


def calculate_MSD(work_dir: Path):
    if not work_dir.exists():
        print("work_dir does not exist")
    dir_str = str(work_dir)
    found_file = None
    file_type = None
    for filename in os.listdir(work_dir):
        file_path = work_dir / filename
        if "dump" in filename:
            found_file = file_path
            file_type = "dump"
            break  # 找到第一个匹配文件就停止搜索

        elif filename.endswith(".xtc"):
            found_file = file_path
            file_type = "xtc"
            break

        # elif filename.endswith(".trr"):
        #     found_file = file_path
        #     file_type = "trr"
        #     break

        elif "XDATCAR" in filename:
            found_file = file_path
            file_type = "XDATCAR"
            break
    if found_file is None:
        print(f"Error: No supported trajectory file found in {work_dir}")
        print("Supported types: dump files, .xtc, .trr, XDATCAR")
        return
    if file_type == "dump":
        x_value, y_value = data_dump(work_dir)
    elif file_type == "xtc":
        x_value, y_value = data_xtc(work_dir)
    # elif file_type == "trr":
    #     x_value, y_value = data_trr(work_dir)
    elif file_type == "XDATCAR":
        x_value, y_value = data_XDATCAR(work_dir)
    plt.plot(x_value, y_value)
    plt.xlabel('Time (ps)')
    plt.ylabel('MSD ($Å^2$)')
    plt.show()
    # 线性拟合后段曲线（提取扩散系数）


if __name__ == '__main__':
    # work_dir = Path('/mnt/d/project/MSD-plot/example/vasp/')
    work_dir = Path('/mnt/d/project/MSD-plot/example/lammps/')
    # work_dir = Path('/mnt/d/project/MSD-plot/example/Gromax/xtc')
    # work_dir = Path('/mnt/d/project/MSD-plot/example/Gromax/trr')
    calculate_MSD(work_dir)

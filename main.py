import os
from pathlib import Path

import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing.tests.test_typing import FAIL_DIR


def data_XDATCAR(work_dir: Path):
    # 检查所需要的文件
    required_files = ['POSCAR', 'INCAR', 'XDATCAR']
    missing_files = []
    for file in required_files:
        file_path = work_dir / file
        if not file_path.exists():
            missing_files.append(file)
    if missing_files:
        raise FileNotFoundError(f"错误: 缺少必需文件: {', '.join(missing_files)}")
    POSCAR = work_dir / 'POSCAR'
    INCAR = work_dir / 'INCAR'
    XDATCAR = work_dir / 'XDATCAR'
    # 获取y_value
    u = mda.Universe(str(POSCAR), str(XDATCAR))
    H = u.select_atoms('all')
    msd_ana = msd.EinsteinMSD(H, msd_type='xyz', fft=True)
    msd_ana.run()
    y_value = msd_ana.results.msds
    # 获取x_value：
    dt_fs = None
    with open(INCAR) as f:
        lines = f.readlines()
        for line in lines:
            if 'POTIM' in line and '=' in line:
                dt_fs = float(line.split('=')[1].split('')[0])
                break
        if dt_fs is None:
            raise ValueError("错误：缺少必要参数：POTIM")
    n_frames = len(msd_ana.results.times)
    x_value = np.arange(n_frames) * dt_fs / 1000  # 换为 ps

    return x_value, y_value


def data_dump(work_dir: Path):
    in_files = list(work_dir.glob('*.in'))
    dump_files = list(work_dir.glob('dump'))
    data_files= list(work_dir.glob('data'))
    if len(in_files) != 1:
        raise FileNotFoundError(f"缺少in文件或in文件不唯一")
    if len(dump_files) != 1:
        raise FileNotFoundError(f"缺少dump文件或dump文件不唯一")
    if len(data_files) != 1:
        raise FileNotFoundError(f"缺少data文件或data文件不唯一")
    in_file =in_files[0]
    dump=dump_files[0]
    data_file=data_files[0]
    # 获取y_value
    u = mda.Universe(str(data_file),str(dump))
    H = u.select_atoms('all')
    msd_ana = msd.EinsteinMSD(H, msd_type='xyz', fft=True)
    msd_ana.run()
    y_value = msd_ana.results.msds
    # 获取x_value
    with open(in_file) as f:
        lines = f.readlines()
        for line in lines:
            if 'timestep' in line:
                dt_fs=float(line.split(' ')[1])
    x_value= np.arange(n_frames) * dt_fs
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
    y_value = msd_ana.results.msds
    # 获取x_value
    x_value = msd_ana.times
    return x_value, y_value


def data_trr(work_dir: Path):
    tpr_files = list(work_dir.glob('*.tpr'))
    trr_files = list(work_dir.glob('*.trr'))
    if len(tpr_files) != 1:
        raise FileNotFoundError(f"缺少.tpr文件或不唯一")
    if len(trr_files) != 1:
        raise FileNotFoundError(f"缺少.trr文件或不唯一")
    tpr = tpr_files[0]
    trr = trr_files[0]
    # 获取y_value
    u = mda.Universe(str(tpr), str(trr))
    H = u.select_atoms('all')
    msd_ana = msd.EinsteinMSD(H, msd_type='xyz', fft=True)
    msd_ana.run()
    y_value = msd_ana.results.msds
    # 获取x_value
    x_value = msd_ana.times
    return x_value, y_value


def calculate_MSD(work_dir: Path):
    if not work_dir.exists():
        print("work_dir does not exist")
    dir_str = str(work_dir)
    for filename in os.listdir(dir_str):
        if filename.endswith(".dump"):
            x_lable, y_lable = data_dump()

            elif suffix == ".xtc":
            msd_calc = data_xtc(work_dir)

            elif suffix == ".trr":
            msd_calc = data_trr(work_dir)
            elif file_name == "XDATCAR":
            msd_calc = data_XDATCAR(work_dir)
            else:
            print("unknown filetype")
    plt.plot(x_value, y_value)
    plt.xlabel('Time (ps)')
    plt.ylabel('MSD ($Å^2$)')
    plt.show()
    # 线性拟合后段曲线（提取扩散系数）
    from scipy.stats import linregress
    # 假设后20%的区间为线性
    N = len(msd_ana.times)
    index_fit = int(N * 0.8)
    slope, intercept, r_value, p_value, std_err = linregress(
        msd_ana.times[index_fit:], msd_ana.results.msds[index_fit:]
    )
    D = slope / 6  # 3维的情况
    print(f'Diffusion Coefficient: {D} Å^2/ps')

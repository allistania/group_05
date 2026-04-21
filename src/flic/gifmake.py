#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import argparse
from PIL import Image
import re

def extract_frame_number(filename):
    match = re.search(r'frame_(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def create_all_vars_gif(frame_dir, output_dir, duration=100, slice_y=None):
    # Создаём выходную директорию, если её нет
    os.makedirs(output_dir, exist_ok=True)

    variables = ['rho', 'u', 'v', 'P', 'e_internal']
    var_labels = {'rho': 'Density', 'u': 'Velocity X', 'v': 'Velocity Y', 'P': 'Pressure', 'e_internal': 'Internal Energy'}

    pattern = os.path.join(frame_dir, 'frame_*.csv')
    all_files = glob.glob(pattern)

    # Сортируем файлы по номеру кадра
    file_number_pairs = []
    for f in all_files:
        num = extract_frame_number(os.path.basename(f))
        if num is not None:
            file_number_pairs.append((num, f))
    file_number_pairs.sort()  # сортировка по номеру
    files = [f for num, f in file_number_pairs]

    if not files:
        print(f"No files matching {pattern}")
        return

    print(f"Found {len(files)} frames")
    images = []

    for i, f in enumerate(files):
        print(f"Processing {i+1}/{len(files)}: {os.path.basename(f)}")
        frame_num = extract_frame_number(os.path.basename(f))
        if frame_num is None:
            frame_num = i

        df = pd.read_csv(f, comment='#')

        x = df['x'].values
        y = df['y'].values

        # Определяем размеры сетки
        nx = len(np.unique(x))
        ny = len(np.unique(y))

        # Создаём фигуру с подграфиками: 2 строки, 3 столбца (последний пустой или можно разместить инфо)
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        # Убираем лишний подграфик (axes[1,2])
        fig.delaxes(axes[1,2])

        # Общий заголовок
        fig.suptitle(f'Frame {frame_num}', fontsize=16)

        for idx, var in enumerate(variables):
            row = idx // 3
            col = idx % 3
            ax = axes[row, col]

            z = df[var].values

            if slice_y is None:
                # 2D цветная карта
                Z = z.reshape(ny, nx)
                X = x.reshape(ny, nx)
                Y = y.reshape(ny, nx)

                im = ax.pcolormesh(X, Y, Z, shading='auto')
                plt.colorbar(im, ax=ax)
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_aspect('equal')
            else:
                # 1D срез при постоянном y
                tol = 1e-6
                mask = np.abs(y - slice_y) < tol
                if not np.any(mask):
                    print(f"Warning: y={slice_y} not found, skipping variable {var}")
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                else:
                    x_slice = x[mask]
                    z_slice = z[mask]
                    idx_sort = np.argsort(x_slice)
                    x_slice = x_slice[idx_sort]
                    z_slice = z_slice[idx_sort]
                    ax.plot(x_slice, z_slice, 'b-')
                    ax.set_xlabel('x')
                    ax.grid(True)

            ax.set_title(var_labels.get(var, var))

        # Если slice_y задан, можно добавить информацию о срезе
        if slice_y is not None:
            axes[1,2].text(0.5, 0.5, f'slice at y = {slice_y}', ha='center', va='center', transform=axes[1,2].transAxes)

        plt.tight_layout()
        # Сохраняем временный PNG
        temp_png = os.path.join(output_dir, f"temp_frame_{i:04d}.png")
        plt.savefig(temp_png, dpi=100)
        plt.close()
        images.append(Image.open(temp_png))

    if images:
        output_gif = os.path.join(output_dir, "all_vars.gif")
        images[0].save(output_gif, save_all=True, append_images=images[1:],
                       duration=duration, loop=0)
        print(f"GIF saved as {output_gif}")

        # Удаляем временные PNG
        for png in glob.glob(os.path.join(output_dir, "temp_frame_*.png")):
            os.remove(png)
    else:
        print("No images created")

def main():
    parser = argparse.ArgumentParser(description="Create a GIF with all variables from simulation frames")
    parser.add_argument('frame_dir', nargs='?', default='forGIF/G',
                        help='Folder with frame_*.csv files (default: forGIF/G)')
    parser.add_argument('--output_dir', default='plots', help='Output directory for GIF (default: plots)')
    parser.add_argument('--duration', type=int, default=100, help='Frame duration (ms)')
    parser.add_argument('--slice_y', type=float, default=None, help='y-coordinate for 1D slice (optional)')
    args = parser.parse_args()

    create_all_vars_gif(args.frame_dir, args.output_dir, args.duration, args.slice_y)

if __name__ == "__main__":
    main()
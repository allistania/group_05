# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import glob
import re

# Create directory for plots
os.makedirs('plots', exist_ok=True)

# Find all analytical and numerical solution files
analytical_files = glob.glob("forGIF/tocnoe/tochnoe*.csv")
numerical_files = glob.glob("forGIF/G/hislennoereshG*.csv")

if not analytical_files or not numerical_files:
    print("Error: No files found!")
    if not analytical_files:
        print("No analytical files found in forGIF/tocnoe/ directory!")
    if not numerical_files:
        print("No numerical files found in forGIF/G/ directory!")
    exit()

print(f"Found {len(analytical_files)} analytical files")
print(f"Found {len(numerical_files)} numerical files")

# Функция для извлечения номера из имени файла
def extract_number(filename, pattern):
    match = re.search(pattern, filename)
    if match:
        return int(match.group(1))
    return -1

# Сортируем файлы по номеру
analytical_files_sorted = sorted([f for f in analytical_files if extract_number(f, r'tochnoe(\d+)\.csv') >= 0], 
                                key=lambda x: extract_number(x, r'tochnoe(\d+)\.csv'))

numerical_files_sorted = sorted([f for f in numerical_files if extract_number(f, r'hislennoereshG(\d+)\.csv') >= 0], 
                               key=lambda x: extract_number(x, r'hislennoereshG(\d+)\.csv'))

print(f"After filtering and sorting: {len(analytical_files_sorted)} analytical files, {len(numerical_files_sorted)} numerical files")

if len(analytical_files_sorted) == 0 or len(numerical_files_sorted) == 0:
    print("No valid files found after filtering!")
    exit()

# Используем минимальное количество файлов для синхронизации анимации
num_frames = min(len(analytical_files_sorted), len(numerical_files_sorted))
print(f"Using {num_frames} frames for animation")

# Выводим список файлов для отладки
print("First 5 analytical files:")
for i, f in enumerate(analytical_files_sorted[:5]):
    print(f"  {i}: {f}")

print("First 5 numerical files:")
for i, f in enumerate(numerical_files_sorted[:5]):
    print(f"  {i}: {f}")

# Создаем фигуру и подграфики
fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# Инициализируем пустые линии для анимации
lines = []

# График давления
axs[0,0].grid()
axs[0,0].set_xlabel("x")
axs[0,0].set_ylabel("P")
analytical_line_p, = axs[0,0].plot([], [], color='#086522', label="Analytical solution", linewidth=2)
numerical_line_p, = axs[0,0].plot([], [], '-o', color='#FF4500', markersize=3, label="Numerical solution")
axs[0,0].legend()
lines.append((analytical_line_p, numerical_line_p))

# График скорости
axs[0,1].grid()
axs[0,1].set_xlabel("x")
axs[0,1].set_ylabel("U")
analytical_line_u, = axs[0,1].plot([], [], color='#086522', label="Analytical solution", linewidth=2)
numerical_line_u, = axs[0,1].plot([], [], '-o', color='#FF4500', markersize=3, label="Numerical solution")
axs[0,1].legend()
lines.append((analytical_line_u, numerical_line_u))

# График плотности
axs[1,0].grid()
axs[1,0].set_xlabel("x")
axs[1,0].set_ylabel("rho")
analytical_line_rho, = axs[1,0].plot([], [], color='#086522', label="Analytical solution", linewidth=2)
numerical_line_rho, = axs[1,0].plot([], [], '-o', color='#FF4500', markersize=3, label="Numerical solution")
axs[1,0].legend()
lines.append((analytical_line_rho, numerical_line_rho))

# График внутренней энергии
axs[1,1].grid()
axs[1,1].set_xlabel("x")
axs[1,1].set_ylabel("I")
analytical_line_i, = axs[1,1].plot([], [], color='#086522', label="Analytical solution", linewidth=2)
numerical_line_i, = axs[1,1].plot([], [], '-o', color='#FF4500', markersize=3, label="Numerical solution")
axs[1,1].legend()
lines.append((analytical_line_i, numerical_line_i))

# Инициализируем пределы осей
x_min, x_max = None, None
p_min, p_max = None, None
u_min, u_max = None, None
rho_min, rho_max = None, None
i_min, i_max = None, None

# Функция для чтения данных из файла
def read_data(file_path):
    with open(file_path) as f:
        labels = f.readline()
        data = f.readlines()
    
    data = list(map(lambda x: x.split(','), data))
    
    t = list(map(lambda x: float(x[0]), data))
    xs = list(map(lambda x: float(x[1]), data))
    rhos = list(map(lambda x: float(x[2]), data))
    Ps = list(map(lambda x: float(x[3]), data))
    us = list(map(lambda x: float(x[4]), data))
    Is = list(map(lambda x: float(x[5]), data))
    
    return t, xs, rhos, Ps, us, Is

# Сначала прочитаем все данные для установки правильных пределов осей
print("Calculating axis limits...")
for i in range(num_frames):
    analytical_file = analytical_files_sorted[i]
    numerical_file = numerical_files_sorted[i]
    
    try:
        # Читаем аналитические данные
        t_analytical, xs_analytical, rhos_analytical, Ps_analytical, us_analytical, Is_analytical = read_data(analytical_file)
        
        # Обновляем пределы осей
        if x_min is None:
            x_min = min(xs_analytical)
            x_max = max(xs_analytical)
            p_min = min(Ps_analytical)
            p_max = max(Ps_analytical)
            u_min = min(us_analytical)
            u_max = max(us_analytical)
            rho_min = min(rhos_analytical)
            rho_max = max(rhos_analytical)
            i_min = min(Is_analytical)
            i_max = max(Is_analytical)
        else:
            x_min = min(x_min, min(xs_analytical))
            x_max = max(x_max, max(xs_analytical))
            p_min = min(p_min, min(Ps_analytical))
            p_max = max(p_max, max(Ps_analytical))
            u_min = min(u_min, min(us_analytical))
            u_max = max(u_max, max(us_analytical))
            rho_min = min(rho_min, min(rhos_analytical))
            rho_max = max(rho_max, max(rhos_analytical))
            i_min = min(i_min, min(Is_analytical))
            i_max = max(i_max, max(Is_analytical))
            
    except Exception as e:
        print(f"Error reading files for axis limits: {e}")

# Устанавливаем пределы осей
def set_axis_limits():
    for ax in axs.flat:
        ax.set_xlim(x_min, x_max)
    
    # Пределы для давления
    axs[0,0].set_ylim(p_min * 0.9, p_max * 1.1)
    
    # Пределы для скорости
    axs[0,1].set_ylim(u_min * 1.1 if u_min < 0 else u_min * 0.9, 
                      u_max * 1.1 if u_max > 0 else u_max * 0.9)
    
    # Пределы для плотности
    axs[1,0].set_ylim(rho_min * 0.9, rho_max * 1.1)
    
    # Пределы для внутренней энергии
    axs[1,1].set_ylim(i_min * 0.9, i_max * 1.1)

set_axis_limits()
plt.tight_layout()

# Функция инициализации анимации
def init():
    for analytical_line, numerical_line in lines:
        analytical_line.set_data([], [])
        numerical_line.set_data([], [])
    return [line for pair in lines for line in pair]

# Функция анимации
def animate(frame):
    analytical_file = analytical_files_sorted[frame]
    numerical_file = numerical_files_sorted[frame]
    
    try:
        # Читаем аналитические данные
        t_analytical, xs_analytical, rhos_analytical, Ps_analytical, us_analytical, Is_analytical = read_data(analytical_file)
        
        # Читаем численные данные
        t_numerical, xs_numerical, rhos_numerical, Ps_numerical, us_numerical, Is_numerical = read_data(numerical_file)
        
        # Обновляем заголовок с текущим временем
        fig.suptitle(f"Solve GD equations at t={t_analytical[0]:.3f} s. (Frame {frame+1}/{num_frames})", fontsize=14)
        
        # Обновляем данные для всех подграфиков
        lines[0][0].set_data(xs_analytical, Ps_analytical)  # Аналитическое давление
        lines[0][1].set_data(xs_numerical, Ps_numerical)    # Численное давление
        
        lines[1][0].set_data(xs_analytical, us_analytical)  # Аналитическая скорость
        lines[1][1].set_data(xs_numerical, us_numerical)    # Численная скорость
        
        lines[2][0].set_data(xs_analytical, rhos_analytical)  # Аналитическая плотность
        lines[2][1].set_data(xs_numerical, rhos_numerical)    # Численная плотность
        
        lines[3][0].set_data(xs_analytical, Is_analytical)  # Аналитическая внутренняя энергия
        lines[3][1].set_data(xs_numerical, Is_numerical)    # Численная внутренняя энергия
        
        print(f"Processed frame {frame+1}/{num_frames}: t={t_analytical[0]:.3f}")
        
    except Exception as e:
        print(f"Error processing files: {e}")
        # Возвращаем пустые данные в случае ошибки
        for analytical_line, numerical_line in lines:
            analytical_line.set_data([], [])
            numerical_line.set_data([], [])
    
    return [line for pair in lines for line in pair]

# Создаем анимацию
print("Creating animation...")
anim = animation.FuncAnimation(
    fig, 
    animate, 
    frames=num_frames,
    init_func=init,
    interval=500,  # 400 мс между кадрами
    blit=True,
    repeat=True
)

# Сохраняем как GIF
print("Saving GIF...")
anim.save('plots/animation.gif', writer='pillow', fps=2, dpi=100)
print("GIF created successfully: plots/animation.gif")

# Показываем информацию о использованных файлах
print(f"\nAnimation completed!")
print(f"Total frames: {num_frames}")
print(f"Analytical files: from {os.path.basename(analytical_files_sorted[0])} to {os.path.basename(analytical_files_sorted[num_frames-1])}")
print(f"Numerical files: from {os.path.basename(numerical_files_sorted[0])} to {os.path.basename(numerical_files_sorted[num_frames-1])}")


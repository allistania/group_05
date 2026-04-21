#!/bin/bash
set -e

echo "=== Запуск метода годунова==="

# Функция для ожидания завершения задачи
wait_for_job() {
    local job_id=$1
    local job_name=$2
    
    echo "Ожидание завершения $job_name "
    while squeue -j $job_id 2>/dev/null | grep -q $job_id; do
        sleep 10
    done
    echo "$job_name завершил выполнение"
}

# Запускаем первую задачу
echo "Начало работы кода на С++"
JOB_ID=$(sbatch --parsable run.slurm)
wait_for_job $JOB_ID "кода на С++"

# Запускаем вторую задачу
echo "Начало работы кода на питоне"
JOB_ID1=$(sbatch --parsable run1.slurm)
wait_for_job $JOB_ID1 "кода на питоне"

# Запускаем третью задачу
echo "Начало работы кода на питоне"
JOB_ID2=$(sbatch --parsable run2.slurm)
wait_for_job $JOB_ID2 "кода на питоне"

echo "Конец выполнения всех программ"


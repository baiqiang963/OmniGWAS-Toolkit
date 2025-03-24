

# 检查 result.mat 是否存在
if [ -f "$save_folder/result.mat" ]; then
    echo "检测到 result.mat 已存在，跳过 MATLAB 生成步骤。"
else
    echo "未找到 result.mat，开始运行 MATLAB 生成..."
    cd "$s_pleiofdr" || { echo "无法进入目录 $s_pleiofdr"; exit 1; }
    module load matlab

    # 后台运行 MATLAB 并获取 PID
    matlab -nodisplay -nosplash < runme.m &
    matlab_pid=$!
    echo "MATLAB 进程 PID: $matlab_pid"

    # 进入保存目录并等待文件生成
    cd "$save_folder" || { echo "无法进入目录 $save_folder"; exit 1; }
    while [ ! -f "result.mat" ]; do
        sleep 60
    done
    # 检测到文件后立即终止 MATLAB
    echo "检测到 result.mat 已生成，终止 MATLAB 进程..."
    kill $matlab_pid 2>/dev/null
    wait $matlab_pid 2>/dev/null  # 等待进程终止避免僵尸进程
fi

cd "$save_folder"
module load anaconda/Anaconda2
echo "正在将 result.mat 转换为 CSV..."
python /home/tengli/BQ_GWAS/python_convert-master/fdrmat2csv.py \
    --mat $save_folder/result.mat \
    --ref /home/tengli/BQ_GWAS/pleiofdr-master/9545380.ref

echo "正在执行 Clump 分析..."
python /home/tengli/BQ_GWAS/python_convert-master/sumstats.py clump \
    --clump-field FDR \
    --force \
    --plink /home/tengli/R/x86_64-pc-linux-gnu-library/4.2/plinkbinr/bin/plink_Linux \
    --sumstats $save_folder/result.mat.csv \
    --bfile-chr /home/tengli/BQ_GWAS/pleiofdr-master/ref9545380_bfile/chr \
    --exclude-ranges 6:25119106-33854733 8:7200000-12500000 \
    --clump-p1 0.05 \
    --out result.clump

echo "处理完成！"

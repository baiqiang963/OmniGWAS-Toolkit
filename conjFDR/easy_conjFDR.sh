#!/bin/bash
#========================================================================
#Usage: sh easy_conjFDR trait1 trait1_n trait2 trait2_n save_folder conjfdr/condfdr s_pleiofdr s_python_convert
#========================================================================

cd $5
DIR="$(basename "$1" | sed "s/\..*$//")_$(basename "$3" | sed "s/\..*$//")_$6"
if [ -d "$DIR" ]; then
    echo "文件夹已存在：$DIR"
else
    echo "文件夹不存在，正在创建：$DIR"
    mkdir -p "$DIR"
    echo "文件夹已创建：$DIR"
fi
cd $5/$DIR
if [ -f "$7/traitfolder/$(basename "$1" | sed "s/\..*$//").mat" ]&& [ -f "$7/traitfolder/$(basename "$3" | sed "s/\..*$//").mat" ]; then
  echo "两个转换过后的文件都在$7/traitfolder/中有备份"
else
  echo "其中一个或两个文件不存在。对GWAS执行格式化操作"
  cp $(dirname "$0")/convert_sumstats.sh $5/$DIR
  sed -i "1a\trait1=$1" $5/$DIR/convert_sumstats.sh
  sed -i "2a\trait1_n=$2" $5/$DIR/convert_sumstats.sh
  sed -i "3a\trait2=$3" $5/$DIR/convert_sumstats.sh
  sed -i "4a\trait2_n=$4" $5/$DIR/convert_sumstats.sh
  sed -i "5a\save_folder=$7/traitfolder" $5/$DIR/convert_sumstats.sh
  sed -i "6a\s_python_convert=$8" $5/$DIR/convert_sumstats.sh
  sh convert_sumstats.sh
fi

cp $7/config_default.txt $5/$DIR/config.txt
file1_base=$(basename "$1" | sed "s/\..*$//")
file2_base=$(basename "$3" | sed "s/\..*$//")

echo "File1 base: $file1_base"
echo "File1 name: $file1_name"
echo "File2 base: $file2_base"
echo "File2 name: $file2_name"

sed -i "12s|.*|traitfolder=$7/traitfolder/|" config.txt
sed -i "15s|.*|traitfile1=${file1_base}.mat|" config.txt
sed -i "16s|.*|traitname1=${file1_base}|" config.txt
sed -i "17s|.*|traitfiles={'${file2_base}.mat'}|" config.txt
sed -i "18s|.*|traitnames={'${file2_base}'}|" config.txt
sed -i "21s|.*|outputdir=$5/$DIR|" config.txt
if [ "$6" == "conjfdr" ]; then
    sed -i "25s/.*/stattype=conjfdr/" config.txt
    sed -i "26s/.*/fdrthresh=0.05/" config.txt
fi
if [ "$6" == "condfdr" ]; then
    sed -i "25s/.*/stattype=condfdr/" config.txt
    sed -i "26s/.*/fdrthresh=0.01/" config.txt
fi
if [ -z "$6" ]; then
    sed -i "25s/.*/stattype=conjfdr/" config.txt
    sed -i "26s/.*/fdrthresh=0.05/" config.txt
fi

cp $5/$DIR/config.txt $7/config.txt

# 复制 run_matlab.sh
cp $(dirname "$0")/run_matlab.sh "$5/$DIR"

# 在文件的开头添加一行（插入到第1行）
sed -i "1s|^|save_folder=$5/$DIR\n|" "$5/$DIR/run_matlab.sh"

# 复制 easy_conjFDR2.sh
cp $(dirname "$0")/easy_conjFDR2.sh "$5/$DIR"

# 在第3行后添加 save_folder 变量
sed -i "3a\save_folder=$5/$DIR" "$5/$DIR/easy_conjFDR2.sh"

# 在第4行后添加 trait1_file 变量
sed -i "4a\trait1_file=$(basename "$1" | sed "s/\.[^.]*$//").csv" "$5/$DIR/easy_conjFDR2.sh"

# 在第5行后添加 trait2_file 变量
sed -i "5a\trait2_file=$(basename "$3" | sed "s/\.[^.]*$//").csv" "$5/$DIR/easy_conjFDR2.sh"

sed -i "1a\s_pleiofdr=$7" $5/$DIR/run_matlab.sh
export PATH="/opt/tsce4/torque6/bin:$PATH"
qsub -q fat $5/$DIR/run_matlab.sh

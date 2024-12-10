import os
import zipfile
import pandas as pd
import shutil
import gzip
import scanpy as sc
import subprocess

def process_matrix_zip(zip_file_path, delete_after=False):
    """
    处理指定的.zip文件，解压.zip后查看其中的.gz文件，并查看前十行和后十行内容，
    最终删除解压出来的matrix.mtx, features.tsv和barcodes.tsv文件。
    
    参数：
        zip_file_path (str): 指向.zip文件的路径。
        delete_after (bool): 是否在处理后删除解压出来的文件。
    
    返回：
        None
    """
    # 判断zip文件是否存在
    if not os.path.exists(zip_file_path):
        print(f"错误: {zip_file_path} 不存在！")
        return

    # 获取zip文件所在目录，并切换到该目录
    current_directory = os.path.dirname(zip_file_path)
    os.chdir(current_directory)
    print(f"当前工作目录已切换到: {current_directory}")
    
    # 解压zip文件到当前目录
    with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
        zip_ref.extractall(current_directory)
        print(f"已解压: {zip_file_path} 到当前目录 {current_directory}")
    
    # 获取解压后的文件列表
    files = os.listdir(current_directory)
    
    # 解压所有.gz文件
    extracted_files = {}
    for file in files:
        if file.endswith(".gz"):
            file_path = os.path.join(current_directory, file)
            with gzip.open(file_path, 'rb') as f_in:
                extracted_file_path = os.path.join(current_directory, file[:-3])  # 删除.gz后缀
                with open(extracted_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                extracted_files[file] = extracted_file_path
                print(f"已解压: {file} 到 {extracted_file_path}")
    
    # 读取解压的features.tsv
    if 'features.tsv.gz' in extracted_files:
        features_file = extracted_files['features.tsv.gz']
        features = pd.read_csv(features_file, sep='\t', header=None)
        print("Features.tsv前十行:")
        print(features.head(10))  # 打印前10行
        print("Features.tsv后十行:")
        print(features.tail(10))  # 打印后10行
    
    # 读取解压的barcodes.tsv
    if 'barcodes.tsv.gz' in extracted_files:
        barcodes_file = extracted_files['barcodes.tsv.gz']
        barcodes = pd.read_csv(barcodes_file, sep='\t', header=None)
        print("Barcodes.tsv前十行:")
        print(barcodes.head(10))  # 打印前10行
        print("Barcodes.tsv后十行:")
        print(barcodes.tail(10))  # 打印后10行
    
    # 读取celltypes.csv
    if 'celltypes.csv' in extracted_files:
        celltypes_file = extracted_files['celltypes.csv']
        celltypes = pd.read_csv(celltypes_file)
        print("Celltypes.csv前十行:")
        print(celltypes.head(10))  # 打印前10行
        print("Celltypes.csv后十行:")
        print(celltypes.tail(10))  # 打印后10行
    
        # 读取解压的matrix.mtx
    # 读取matrix.mtx文件，使用Linux命令head和tail读取前后50行
    if 'matrix.mtx.gz' in extracted_files:
        mtx_file = extracted_files['matrix.mtx.gz']
        
        # 使用head命令读取matrix.mtx文件的前50行
        print("Matrix.mtx前50行:")
        head_cmd = f"head -n 50 {mtx_file}"
        head_output = subprocess.check_output(head_cmd, shell=True, text=True)
        print(head_output)

        # 使用tail命令读取matrix.mtx文件的后50行
        print("Matrix.mtx后50行:")
        tail_cmd = f"tail -n 50 {mtx_file}"
        tail_output = subprocess.check_output(tail_cmd, shell=True, text=True)
        print(tail_output)
    
    # 如果delete_after为True，删除解压出来的matrix.mtx, features.tsv, barcodes.tsv文件
    if delete_after:
        for file_name in ['matrix.mtx', 'features.tsv', 'barcodes.tsv']:
            if file_name in extracted_files:
                os.remove(extracted_files[file_name])
                print(f"已删除解压的文件: {extracted_files[file_name]}")

# 示例调用
zip_file_path = "F:/biowork/raw_data/CELLxGENE/HTAPP-783-SMP-4081_scRNA-seq/h5ad_output/matrix.zip"  # 传入zip文件路径
process_matrix_zip(zip_file_path, delete_after=True)

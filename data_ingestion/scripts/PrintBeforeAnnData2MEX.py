import os
import h5py
import scanpy as sc
import pandas as pd
import sys

def process_h5ad_file(file_path, print_all=False):
    # 验证文件路径是否存在
    print(file_path)
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # 读取文件
    obj = h5py.File(file_path)
    print("Keys in the H5AD file:", obj.keys())

    # 用Scanpy读取文件
    ad = sc.read_h5ad(file_path)

    # 提取raw数据
    raw_ad = ad.raw.to_adata()
    print("\n raw_ad = ad.raw.to_adata()")

    # 设置pandas显示选项
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    pd.set_option('display.max_rows', 30)
    pd.set_option('display.max_colwidth', 100)

    # 打印数据   
    if print_all:
        # 遍历 raw_ad.obs 的所有列，并输出每列的数据
        print("\nAll columns in raw_ad.obs:")
        for column in raw_ad.obs.columns:
            print(f"Column: {column}")
            print(raw_ad.obs[column].head(10))  # 打印每列的前10行数据

        # 遍历 raw_ad.var 的所有列，并输出每列的数据
        print("\nAll columns in raw_ad.var:")
        for column in raw_ad.var.columns:
            print(f"Column: {column}")
            print(raw_ad.var[column].head(10))  # 打印每列的前10行数据

    # 判断 'cell_type' 列是否存在并打印
    if 'cell_type' in raw_ad.obs.columns:
        print("\nFirst 10 rows of 'cell_type':")
        print(raw_ad.obs[["cell_type"]].head(10))
    else:
        print("'cell_type' column is not found in raw_ad.obs.")

    # 判断 'barcode_key' 列是否存在并打印
    if 'barcode_key' in raw_ad.obs.columns:
        print("\nFirst 10 rows of 'barcode_key':")
        print(raw_ad.obs[["barcode_key"]].head(10))
    else:
        print("'barcode_key' column is not found in raw_ad.obs.")

    print("\nFirst row of raw_ad.var:")
    print(raw_ad.var.head(10).to_string(index=False))

# 调用函数
if __name__ == "__main__":
    file_path = "F:/biowork/raw_data/CELLxGENE/HTAPP-783-SMP-4081_scRNA-seq/HTAPP-783-SMP-4081_scRNA-seq.h5ad"  # 替换为实际的.h5ad文件路径
    file_path2 = "F:/biowork/raw_data/CELLxGENE/A single-cell and spatially-resolved atlas of human breast cancers/A single-cell and spatially-resolved atlas of human breast cancers.h5ad"  # 替换为实际的.h5ad文件路径

    process_h5ad_file(file_path,True)

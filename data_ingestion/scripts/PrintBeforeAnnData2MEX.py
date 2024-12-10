import os
import h5py
import scanpy as sc
import pandas as pd

def process_h5ad_file(file_path):
    # 验证文件路径是否存在
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # 读取文件
    obj = h5py.File(file_path)
    print("Keys in the H5AD file:", obj.keys())

    # 用Scanpy读取文件
    ad = sc.read_h5ad(file_path)

    # 提取raw数据
    raw_ad = ad.raw.to_adata()

    # 设置pandas显示选项
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_rows', 10)
    pd.set_option('display.max_colwidth', 100)

    # 打印数据
    print("First row of raw_ad.obs:")
    print(raw_ad.obs[["cell_type"]].head(1))

    print("\nFirst row of raw_ad.var:")
    print(raw_ad.var.head(10).to_string(index=False))

# 调用函数
if __name__ == "__main__":
    file_path = "/home/lcr/Desktop/BioWork/WorkData/A single-cell and spatially-resolved atlas of human breast cancers/A single-cell and spatially-resolved atlas of human breast cancers.h5ad"  # 替换为实际的.h5ad文件路径
    process_h5ad_file(file_path)

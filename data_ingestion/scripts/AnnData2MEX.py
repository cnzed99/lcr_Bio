import os
import scanpy as sc
import pandas as pd
import scipy
import scipy.sparse as sparse
import zipfile
import tempfile
import gzip


def h5ad_to_10x(ad,
                gene_id_key="gene_id",
                gene_name_key="gene_name",
                cell_type_key="cell_type",
                output_path="matrix.zip",
                barcode_key=None,
                subsample_rate=0.9):
    """
    将 AnnData 的 h5ad 文件转换为 10x 格式，并生成包含最终矩阵的 zip 文件。

    参数说明：
        ad: 一个 scanpy 的 AnnData 对象。
        gene_id_key: 基因 ID 所在的键名。
        gene_name_key: 基因名称所在的键名。
        cell_type_key: 细胞类型所在的键名。
        barcode_key: 可选，条形码所在的键名。如果未设置，则假设位于 `obs` 的索引中。
        output_path: 保存 zip 文件的路径。
        subsample_rate: 可选，子采样比率（0~1 之间）。
    """
    if subsample_rate:
        # 按指定比率子采样
        sc.pp.subsample(ad, subsample_rate)

    # 处理基因信息
    genes = ad.var.reset_index()[[gene_id_key, gene_name_key]]
    genes["feature_type"] = ["Gene Expression"] * len(genes)

    # 处理条形码信息
    if barcode_key:
        barcodes = ad.obs[[barcode_key]]
    else:
        barcodes = pd.DataFrame(ad.obs.index)

    # 处理细胞类型信息
    celltypes = ad.obs[[cell_type_key]].reset_index()
    celltypes.columns = ["barcode", "annotation"]

     # 确保将数据矩阵转换为整数类型
    ad.X = ad.X.astype('int32')  # 转换为 int32 类型
    print(ad.X.dtype)  # 输出矩阵的数据类型，应该是 int32 或 int64

    with tempfile.TemporaryDirectory() as tmp_dir:
        # 写入压缩稀疏矩阵
        with gzip.open(os.path.join(tmp_dir, "matrix.mtx.gz"), "w") as handle:
            scipy.io.mmwrite(handle, sparse.csc_matrix(ad.X.T))

        # 写入基因、条形码和细胞类型信息
        genes.to_csv(os.path.join(tmp_dir, "features.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")
        barcodes.to_csv(os.path.join(tmp_dir, "barcodes.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")
        celltypes.to_csv(os.path.join(tmp_dir, "celltypes.csv"), index=False)

        # 打包成 zip 文件
        with zipfile.ZipFile(output_path, "w") as zip_handle:
            for file in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz", "celltypes.csv"]:
                zip_handle.write(os.path.join(tmp_dir, file), arcname=file)


def process_h5ad(file_path, output_path):
    """
    处理 h5ad 文件，将其转换为 10x 格式，并保存结果。

    参数说明：
        file_path: 输入 h5ad 文件的路径。
        output_path: 输出 zip 文件的保存路径。
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件 {file_path} 不存在。")

    # 切换工作目录到文件所在目录
    file_dir = os.path.dirname(file_path)
    os.chdir(file_dir)

    # 读取 AnnData 对象
    ad = sc.read_h5ad(os.path.basename(file_path))
    
    # 调用 h5ad_to_10x 函数进行转换
    h5ad_to_10x(ad,
                gene_id_key="feature_name",       # 基因 ID 键名
                gene_name_key="feature_reference",  # 基因名称键名
                cell_type_key="cell_type",     # 细胞类型键名
                output_path=output_path)       # 输出路径


if __name__ == "__main__":
    # 示例调用参数
    input_file = "/home/lcr/Desktop/BioWork/WorkData/A single-cell and spatially-resolved atlas of human breast cancers/A single-cell and spatially-resolved atlas of human breast cancers.h5ad"  # 替换为实际的 h5ad 文件路径
    output_file = "/home/lcr/Desktop/BioWork/WorkData/OutPutData/A single-cell and spatially-resolved atlas of human breast cancers-matrix.zip"  # 替换为实际的输出文件路径
    process_h5ad(input_file, output_file)
    
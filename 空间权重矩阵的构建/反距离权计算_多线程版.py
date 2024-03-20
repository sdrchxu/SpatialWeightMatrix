#在使用使用cdist函数的基础上通过多线程计算缩短运算时间
import geopandas as gpd
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from multiprocessing import Pool,freeze_support

def compute_distances(df):
    coordinates = df[['X', 'Y', 'Z']].values
    distances = cdist(coordinates, df[['X', 'Y', 'Z']].values)
    return distances

if __name__ == '__main__':

    freeze_support()
    # 读取Shapefile文件
    shp_path = r"C:\Users\Administrator\Desktop\空间权重矩阵计算\shp文件\云南省已处理好的火点.shp"
    gdf = gpd.read_file(shp_path)

    # 提取所需的字段（ID, X, Y, Z）
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()

    # 获取CPU核心数量，即线程数
    num_threads = 28

    # 将数据集分成与线程数相等的批次
    batch_size = len(df) // num_threads + (1 if len(df) % num_threads != 0 else 0)
    batches = [df.iloc[i:i+batch_size] for i in range(0, len(df), batch_size)]

    # 使用线程池并行计算距离矩阵
    with Pool(num_threads) as pool:
        results = pool.map(compute_distances, batches)

    # 合并结果
    distances = np.concatenate(results, axis=0)

    # 将距离矩阵保存到DataFrame中
    distances_df = pd.DataFrame(distances)

    # 保存距离矩阵为CSV文件
    output_file_path = r'C:\Users\Administrator\Desktop\空间权重矩阵计算\输出的矩阵\spatial_weights_matrix反距离权.csv'
    distances_df.to_csv(output_file_path)





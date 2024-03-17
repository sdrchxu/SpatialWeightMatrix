import geopandas as gpd
import pandas as pd
import numpy as np
from pysal.lib import weights
from scipy.spatial.distance import cdist
from pysal.explore import esda
from splot.esda import plot_moran
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.sparse import csr_matrix

#多线程版本代码
def inverse_weights(gdf, L0):
    """
    计算反距离权矩阵

    参数：
    ------------
    gdf:         shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
    """
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()
    coordinates = df[['X', 'Y', 'Z']].values
    distances = cdist(coordinates, coordinates)
    spatial_weights = np.where(distances == 0, 0, np.where(distances <= L0, 1 / distances, 0))
    spatial_weights_df = pd.DataFrame(spatial_weights)
    return spatial_weights_df

def gauss_weights(gdf, L0):
    """
    计算高斯权矩阵

    参数：
    ------------
    gdf:    shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
    """
    coordinates = gdf[['X', 'Y', 'Z']].values
    distances = cdist(coordinates, coordinates)
    weights = np.exp(-distances ** 2 / (2 * L0 ** 2))
    spatial_weights_df = pd.DataFrame(weights, index=gdf['ID'], columns=gdf['ID'])
    return spatial_weights_df

def threshold_weights(gdf, L0):
    """
    计算阈值权矩阵

    参数：
    ------------
    gdf:    shapefile文件对象
    L0:          阈值

    返回值：
    -------------
    包含空间权重矩阵的DataFrame
    """
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()
    coordinates = df[['X', 'Y', 'Z']].values
    distances = cdist(coordinates, coordinates)
    spatial_weights = np.where(distances <= L0, 1, 0)
    spatial_weights_df = pd.DataFrame(spatial_weights)
    return spatial_weights_df


def calculate_moran(gdf, field, w, std):
    moran = esda.Moran(gdf[field], w)
    return moran


def global_moran(shp_path, field, output_file, distance_function, threshold=float('inf'), std=True):

    """
    计算全局Moran'I指数，调用该函数后会计算全局Moran'I指数，
    并在指定路径生成txt和png分析结果

    参数：
    -----------
    shp_path:          shapefile文件位置
    field:             输入字段
    distance_function: 空间关系概念化函数，可选threshold/gauss/inverse
    output_file:       输出路径，输出文件名无需写扩展名
    threshold:         距离阈值，留空则为不设置阈值
    std:               标准化方法，True为行标准化，False为不进行标准化

    返回值：
    -------------
    None
    """
    print("Reading Shapefile...")
    gdf = gpd.read_file(shp_path)
    print("Calculating Weights Matrix...")
    if distance_function == 'threshold':
        spatial_weights_df = threshold_weights(gdf, threshold)
    elif distance_function == 'gauss':
        spatial_weights_df = gauss_weights(gdf, threshold)
    elif distance_function == 'inverse':
        spatial_weights_df = inverse_weights(gdf, threshold)

    # Convert DataFrame to numpy array
    sparse_matrix = spatial_weights_df.values

    # Parallelize Moran's I calculation
    print("Calculating Moran'I...")
    w = weights.WSP(csr_matrix(sparse_matrix))
    w = weights.W.from_WSP(w)
    if std:
        w.transform = 'r'

    #n_jobs表示使用全部核心
    results = Parallel(n_jobs=-1)(
        delayed(calculate_moran)(gdf, field, w, std) for _ in range(10))

    moran = max(results, key=lambda x: x.I)

    print("Saving Result Reports...")
    output_txt = output_file + '.txt'
    output_png = output_file + '.png'

    with open(output_txt, 'w') as f:
        f.write('Moran\'s I: {}\n'.format(moran.I))
        f.write('P-value: {}\n'.format(moran.p_sim))
        f.write('Z-score: {}\n'.format(moran.z_sim))

    plot_moran(moran, zstandard=True, figsize=(20, 8))
    plt.savefig(output_png)
    plt.close()
    print("Done!")

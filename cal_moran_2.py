import geopandas as gpd
import pandas as pd
import numpy as np
from pysal.lib import weights
from scipy.sparse import csr_matrix
import scipy
import shapefile
from scipy.spatial.distance import cdist
from pysal.explore import esda
from splot.esda import plot_moran
import matplotlib.pyplot as plt

def inverse_weights(gdf,L0):
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
    # 提取所需的字段（ID, X, Y, Z）
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()


    # 计算距离矩阵
    coordinates = df[['X', 'Y', 'Z']].values
    # coordinates = df[['X', 'Y']].values
    distances = cdist(coordinates, coordinates)
    # 应用反距离权重
    spatial_weights = np.where(distances==0, 0, np.where(distances <= L0,1/distances,0))
    # 将空间权重矩阵保存到DataFrame中
    spatial_weights_df = pd.DataFrame(spatial_weights)
    return spatial_weights_df




def gauss_weights(gdf,L0):
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
    # Converting to numpy array
    coordinates = gdf[['X', 'Y', 'Z']].values
    # Calculating distance matrix
    distances = cdist(coordinates, coordinates)

    # Applying Gaussian function to calculate weights
    weights = np.exp(-distances**2 / (2 * L0**2))

    # Creating DataFrame with weights and setting index and column names
    spatial_weights_df = pd.DataFrame(weights, index=gdf['ID'], columns=gdf['ID'])

    return spatial_weights_df



def threshold_weights(gdf,L0):
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
    # 提取所需的字段（ID, X, Y, Z）
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()
    coordinates = df[['X', 'Y', 'Z']].values

    #计算要素间的距离
    # coordinates = df[['X', 'Y']].values
    distances = cdist(coordinates, coordinates)

    # 对得到的距离进行判断并赋值
    spatial_weights = np.where(distances<=L0, 1, 0)

    # 将空间权重矩阵保存到DataFrame中
    spatial_weights_df = pd.DataFrame(spatial_weights)
    return spatial_weights_df


def global_moran(shp_path,field,output_file,distance_function,threshold=float('inf'),std=True):
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
    sf = shapefile.Reader(shp_path)
    print("Calculating Weights Matrix...")
    if distance_function=='threshold':
        spatial_weights_df = threshold_weights(gdf,threshold)
    elif distance_function=='gauss':
        spatial_weights_df = gauss_weights(sf,threshold)  #高斯权函数需要传入一个Shapefile对象
    elif distance_function=='inverse':
        spatial_weights_df = inverse_weights(gdf,threshold)
    sparse_matrix=scipy.sparse.csr_matrix(spatial_weights_df.values)

    print("Calculating Moran'I...")
    # 创建空间权重对象,WSP 类是pysal中的一个子类，表示"weights spatial",专用于表示处理空间权重矩阵的子类
    wsp = weights.WSP(sparse_matrix)
    #WSP无法直接参与Moran'I指数计算，所以需要将它转为w类
    w=weights.W.from_WSP(wsp)
    if std==True:
        w.transform = 'r' #对空间权重矩阵进行行标准化

    # 如果权重矩阵是对称的，可以使用Sym的子类
    # w = weights.WSP(symmetrize='True', values=spatial_weights_df.values)

    moran = esda.Moran(gdf[field], w)

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

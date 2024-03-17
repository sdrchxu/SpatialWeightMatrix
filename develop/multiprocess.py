import pandas as pd
import geopandas as gpd
from scipy.spatial.distance import cdist
import numpy as np
from multiprocessing import Pool, cpu_count
import time


def process_weights(args):
    i, spatial_weights_length, spatial_weights, gdf = args
    weight_info = []
    for j in range(spatial_weights_length):
        weight = spatial_weights[i, j]
        if weight != 0:
            weight_info.append(f"{gdf.iloc[i]['ID']} {gdf.iloc[j]['ID']} {weight}")
    return weight_info

if __name__ == '__main__':
    time_start = time.time()
    shp_path="F:\\大创数据\\中间产出的数据\\云南省和黄淮海平原已处理好的火点\\黄淮海平原逐月火点\\黄淮海平原已处理好的火点_10月.shp"
    swm_path="D:\\Lenovo\\Desktop\\云南大学\\大创\\程序代码\\空间权重矩阵测试\\swm测试\\黄淮海平原已处理好的火点_10月.swm"
    txt_path="D:\\Lenovo\\Desktop\\云南大学\\大创\\程序代码\\空间权重矩阵测试\\swm测试\\黄淮海平原已处理好的火点_10月.txt"
    print("读取shp文件")
    gdf = gpd.read_file(shp_path)
    df = gdf[['ID', 'X', 'Y', 'Z']].copy()
    coordinates = gdf[['X', 'Y', 'Z']].values
    print("计算要素距离矩阵")
    L0 = 55000
    distances = cdist(coordinates, coordinates)
    print("计算空间权重矩阵")
    spatial_weights = np.where((distances <= L0) & (distances != 0), 1, 0)

    spatial_weights_length = len(spatial_weights)
    print("生成权重列表")
    weight_info = []
    #按CPU核心数分配进程
    with Pool(cpu_count()) as p:
        args = [(i, spatial_weights_length, spatial_weights, gdf) for i in range(spatial_weights_length)]
        weight_info = p.map(process_weights, args)
    #将嵌套列表展开
    weight_info = [item for sublist in weight_info for item in sublist]

    print("生成txt文件")
    with open(txt_path, 'w', encoding="ascii") as f:
        f.write("ID\n")
        for info in weight_info:
            f.write(f"{info}\n")

    time_end=time.time()
    print("运行时间：",time_end-time_start,"秒")





import geopandas as gpd
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

# 读取Shapefile文件
shp_path = r"C:\Users\Administrator\Desktop\空间权重矩阵计算\shp文件\云南省已处理好的火点.shp"
gdf = gpd.read_file(shp_path)

# 提取所需的字段（ID, X, Y, Z）
df = gdf[['ID', 'X', 'Y', 'Z']].copy()

# 设置距离阈值L0
L0 = 100000

# 计算距离矩阵
coordinates = df[['X', 'Y', 'Z']].values
distances = cdist(coordinates, coordinates)

# 应用反距离权重（距离的倒数，加上一个小数值防止除以0错误）
spatial_weights = np.where(distances <= L0, 1 / (distances + 1e-9), 0)

# 将空间权重矩阵保存到DataFrame中
spatial_weights_df = pd.DataFrame(spatial_weights)

# 保存空间权重矩阵为CSV文件
output_file_path = r'C:\Users\Administrator\Desktop\空间权重矩阵计算\输出的矩阵\spatial_weights_matrix反距离权.csv'
spatial_weights_df.to_csv(output_file_path)

import sys
import cal_moran_3


shp_path="F:\大创数据\中间产出的数据\云南省和黄淮海平原已处理好的火点\云南省逐月火点\云南省已处理好的火点_3月.shp"
field="BRIGHT_T31"
output_file="D:\Lenovo\Desktop\云南大学\大创\程序代码\空间权重矩阵测试\云南省3月火点（阈值，多线程）"
cal_moran_3.global_moran(shp_path,field,output_file,distance_function='inverse',threshold=55000,std=True)
import cal_moran_2

shp_path="F:\大创数据\中间产出的数据\云南省和黄淮海平原已处理好的火点\云南省逐月火点\云南省已处理好的火点_1月.shp"
output_file_path="D:\Lenovo\Desktop\云南大学\大创\程序代码\空间权重矩阵测试\高斯权函数结果\云南省1月火点高斯权重.txt"
cal_moran_2.output_gauss_weights(shp_path,L0=55000,elevation=False,output_file_path=output_file_path)
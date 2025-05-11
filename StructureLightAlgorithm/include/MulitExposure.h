#pragma once
//这个文件中将复现不同多重曝光方法

#include "ALLHeads.h"
#include <opencv2/opencv.hpp>


//定义多重曝光方法基类
class HDR_Reconstruction {


public:

	HDR_Reconstruction();

	~HDR_Reconstruction();

	///*
	// * @ brief   生成曝光时间
	// * @ return  返回生成的曝光时间序列
	//*/
	//virtual std::vector<double> GenerateExposure() = 0;

	///*
	// * @ brief   生成自适应条纹图案序列
	// * @ return  返回生成的条纹图案序列
	//*/
	//virtual std::vector<cv::Mat> GenerateAdaptiveStripePattern() = 0;

	///*
	// * @ brief   生成灰度直方图
	// * @ return  返回灰度直方图图案
	//*/
	//virtual cv::Mat GenerationGrayHistogram() = 0;


private:


};


//参考文献：汪锦航,卢荣胜,刘端茂.高动态范围表面自适应条纹投影测量方法[J].光学学报,2021,41(19):145-154.
/*
	本方法是一种基于递归的自适应条纹投影方法。该论文使用的是单目方法
*/
class HDR_AdaptiveStripePattern : public HDR_Reconstruction {

public:

	HDR_AdaptiveStripePattern();


	~HDR_AdaptiveStripePattern();


	//cv::Mat GenerationGrayHistogram() override;







private:




};



//参考文献：冯维,汤少靖,赵晓冬,等.基于自适应条纹的高反光表面三维面形测量方法[J].光学学报,2020,40(05):119-127.
/*
	一种基于图像融合和插值预测的自适应条纹投影方法
*/

class HDR_AdaptiveStripePattern_FW : public HDR_Reconstruction {

public:

	HDR_AdaptiveStripePattern_FW();


	~HDR_AdaptiveStripePattern_FW();


	//cv::Mat GenerationGrayHistogram() override;




private:
 



public:



	/*
	 * @ brief   生成掩膜图像序列：编列图像序列，当前图像的灰度值是整个图像序列中的最大值，并且这个最大值超过阈值
	 *			 将当前掩膜图像像素值置为0，如果当前图像灰度值不是最大值或者最大值未超过阈值，则将掩膜图像像素值置为1
	 * @ param   DiffGrayValuePatterns     相机拍摄的不同灰度值的均匀灰度图像
	 * @ param   MaskImages                生成的掩膜图像序列         
	*/
	void getMaskImages(std::vector<cv::Mat>& DiffGrayValuePatterns, std::vector<cv::Mat>& MaskImages);

	/*
	 * @ brief   求取有效均匀灰度图像序列：将灰度值高于阈值的像素值置为0，其他保持不变
	 * @ param   DiffGrayValuePatterns      相机拍摄的不同灰度值的均匀灰度图像
	 * @ param   GraySequences              生成的有效灰度图像序列
	*/
	void getEffectiveGraySequences(std::vector<cv::Mat>& DiffGrayValuePatterns, std::vector<cv::Mat>& GraySequences);

	/*
	 * @ brief   根据掩膜图像和有效均匀灰度图像序列生成投影图像：将两组序列图像相乘在做累加
	 * @ param   MaskImages             生成的掩膜图像序列
	 * @ param   GraySequences          生成的有效灰度图像序列
	*/
	void generateProjectedImages(std::vector<cv::Mat>& MaskImages, std::vector<cv::Mat>& GraySequences);


	//上面的图像并不是最佳投影灰度值，采用差值预测快速查找法求得最佳投影灰度值


	//进行相机-投影仪标定


	//总结：其实最终就是生成了一个单一最大灰度值的条纹序列，即最大值不是255而是计算出来的某一个值


};


//参考文献：[1]汪运,郭建英,梁浚哲,等.面向高反光表面的结构光面形测量方法[J].中国光学(中英文),2025,18(01):42-52.
/*
	依赖相机响应曲线，计算反应被测表面反射强度的辐照度图像，通过模糊C均值聚类方法，自适应分割不同辐照度区域
	并获得各区域中心的辐照度，利用相机响应曲线预测不同反射区的最优曝光时间
*/
class HDR_CameraResponseCurve : public HDR_Reconstruction {

public:

	/*
	 * @ brief   基于相机响应曲线的高动态范围重建，构造函数
	 * @ param   filepath               存放不同曝光时间图像的文件路径
	 * @ param   Exposures              曝光时间序列
	*/
	HDR_CameraResponseCurve(std::string& filepath,
							std::string& resultpath,
							std::vector<float>& Exposures,
							int width,
							int height);


	~HDR_CameraResponseCurve();
	 
public:




	/*
	 * @ brief   获取相机的相应曲线，传入不同曝光下的图像序列
	 * @ param   DiffExposureImages     不同曝光下的图像序列
	 * @ param   Exposures              曝光时间序列
	 * @ param   ResponseCurveImg       生成相机响应曲线数组，表明照度的对数和相机响应的灰度值之间的关系
	 * PS 针对如何计算相机的响应曲线参考：https://blog.csdn.net/wc781708249/article/details/78490976
	 * https://blog.csdn.net/u013049912/article/details/85157402
	*/
	bool getCameraResponseCurve(std::vector<cv::Mat>& DiffExposureImages, 
								std::vector<float>& Exposures,
								cv::Mat& ResponseCurveImg);


	/*
	 * @ brief   计算辐照度图像，传入单张图像,这个函数可能还不对
	 * @ param   Image           单次曝光下拍摄到的图像
	 * @ param   Exposure        单次的曝光时间        
	 * @ param   gain            图像增益
	*/
	bool CalculateIrradianceImage(cv::Mat& Image, float Exposurem, float gain);




private:
	//存放曝光图像序列的文件
	std::string M_ExposureImagesFilepath;

	//存放相机相应图像结果
	std::string M_CamResponseResultpath;

	//曝光时间序列
	std::vector<float> M_Exposures;

	int M_Height;
	int M_Width;

};
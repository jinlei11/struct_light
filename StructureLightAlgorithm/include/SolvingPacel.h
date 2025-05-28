#pragma once
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <math.h>
#include <filesystem>
#include <string>
#include <algorithm>
#include <functional>
#include "FileManipulation.h"
#include "PhaseErrorAnalyze.h"


#define pi CV_PI

enum CameraSolutioin {
	A5131 = 0,
	A5201 = 1,
	A7300 = 2,
	A7500 = 3,
};


enum CameraPosition {
	Left = 0,
	Right = 1,
};



class SolvingPacel {
public:
	/*
	 * @ brief  相位信息计算类构造函数：针对于单次曝光
	 * @ param  ImagepathLeft             存放左相机拍摄的投影图案的路径
	 * @ param  ImagepathRight            存放右相机拍摄的投影图案的路径
	 * @ param  CalibrationResultPath     存放相机标定结果
	 * @ param  GrayImgsNum               格雷码图像的数量
	 * @ param  PhaseImgsNum              相移图案数量
	 * @ param  width                     图像的宽度
	 * @ param  height                    图像的高度
	 * @ param  AncillaryPatternNum       投影的辅助图片的数量
	 * @ param  isRename                  是否需要重命名图片
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
				 const std::string& CalibrationResultPath,
				 int GrayImgsNum,
				 int PhaseImgsNum,
				 int width,
				 int height,
				 int AncillaryPatternNum,
				 bool isRename);

	/*
	 * @ brief  相位信息计算类构造函数：针对于单次曝光
	 * @ param  imagepath            存放拍摄的投影图案的路径
	 * @ param  GrayImgsNum          格雷码图像的数量
	 * @ param  PhaseImgsNum         相移图案数量
	 * @ param  Form                 相机的型号，由此来计算图像的宽度和高度
	 * @ param  AncillaryPatternNum  投影的辅助图片的数量
 	 * @ param  isRename             是否需要重命名图片
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
		         const std::string& CalibrationResultPath,
				 int NumOfGray,
				 int NumOfPhase,
				 enum CameraSolutioin Form,
				 int AncillaryPatternNum,
				 bool isRename);

	/*
	 * @ brief  相位信息计算类构造函数：针对与多次曝光
	 * @ param  ImagepathLeft             存放左相机拍摄的投影图案的路径
	 * @ param  ImagepathRight            存放右相机拍摄的投影图案的路径
	 * @ param  CalibrationResultPath     存放相机标定结果
	 * @ param  GrayImgsNum               格雷码图像的数量
	 * @ param  PhaseImgsNum              相移图案数量
	 * @ param  width                     图像的宽度
	 * @ param  height                    图像的高度
	 * @ param  ExposureNums              多重曝光的曝光次数
	 * @ param  AncillaryPatternNum       投影的辅助图片的数量
	 * @ param  isRename                  是否需要重命名图片
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
		         const std::string& CalibrationResultPath,
				 int GrayImgsNum,
				 int PhaseImgsNum,
				 int width,
				 int height,
				 int ExposureNums,
				 int AncillaryPatternNum,
		         bool isRename);

	/*
	 * @ brief  相位信息计算类构造函数：针对与多次曝光
	 * @ param  ImagepathLeft             存放左相机拍摄的投影图案的路径
	 * @ param  ImagepathRight            存放右相机拍摄的投影图案的路径
	 * @ param  CalibrationResultPath     存放相机标定结果
	 * @ param  GrayImgsNum               格雷码图像的数量
	 * @ param  PhaseImgsNum              相移图案数量
	 * @ param  Form                      相机的型号，由此来计算图像的宽度和高度
	 * @ param  ExposureNums              多重曝光的曝光次数
	 * @ param  AncillaryPatternNum       投影的辅助图片的数量
	 * @ param  isRename                  是否需要重命名图片
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
		         const std::string& CalibrationResultPath,
				 int NumOfGray,
				 int NumOfPhase,
				 enum CameraSolutioin Form,
				 int ExposureNums,
				 int AncillaryPatternNum,
				 bool isRename);


	/*
	 * @ brief  对单组投影图案进行相位解算
	 * @ param  ProjectionPatterns   单组投影图案
	 * @ param  Phase                解算出的解包裹相位
	 * @ param  FringeOrder          条纹级次图
	 * @ param  regulatory           调制度B的值
	*/
	void SingleDePhase(const std::vector<cv::Mat>& ProjectionPatterns,
									 cv::Mat& Phase,
									 cv::Mat& FringeOrder,
									 int regulatory);
	 

	/*
	 * @ brief  对多组投影图案进行相位解算,该函数直接调用单次解算函数，每组曝光的投影图案直接解相互不干扰
	 * @ param  ProjectionPatterns   单组投影图案
	 * @ param  SovelPacelPhases     解算出的解包裹相位序列
	 * @ param  StripesLevel         条纹级次
	 * @ param  ExposureNums         多重曝光的曝光次数
	 * @ param  regulatory           调制度B的值
	*/
	void MultiDePhase(const std::vector<cv::Mat>& ProjectionPatterns, 
					  std::vector<cv::Mat>& SovelPacelPhases, 
					  std::vector<cv::Mat>& StripesLevel,
					  int ExposureNums, 
					  int regulatory);


	/*
	 * @ brief  对多重曝光的解包裹相位进行融合
	 * @ param  SovelPacelPhases     解算出的解包裹相位序列
	 * @ param  ExposureNums         多重曝光的曝光次数
	 * @ param  IntegrationPhase     融合之后的相位
	 * @ param  GrayThresholds       灰度阈值
	 * @ param  position             表明是左相机还是右相机
	*/
	void PhasesIntegration(std::vector<cv::Mat>& SovelPacelPhases, 
						   int ExposureNums,
						   cv::Mat& IntegrationPhase,
						   int GrayThresholds,
		                   CameraPosition position);


	/*
	 * @ brief  极线矫正函数
	 * @ param  filename                  存放相机标定参数的文件路径
	 * @ param  unwrappedLeft             解算出的左相机的解包裹相位
	 * @ param  unwrappedRight            解算出的左相机的解包裹相位
	 * @ param  jiaozhengL                极线矫正后的左相机解包裹相位
	 * @ param  jiaozhengR                极线矫正后的右相机解包裹相位
	 * @ param  isShowPolarrectification  是否显示极线矫正的结果
	 * @ return  双目标定的Q矩阵
	*/
	cv::Mat readCameraInfoAndPolarrectification(const std::string& filename,
												cv::Mat& unwrappedLeft,
												cv::Mat& unwrappedRight,
												cv::Mat& jiaozhengL,
												cv::Mat& jiaozhengR,
												bool isShowPolarrectification = false);

	/*
	 * @ brief  双目结构光相位恢复
	 * @ param  correction_left      极线矫正后的左相机展开相位
	 * @ param  correction_right     极线矫正后的右相机展开相位
	 * @ param  StripesLevel_Left    左相机条纹级次
	 * @ param  StripesLevel_Right   右相机条纹级次
	 * @ param  GrayThreshold        灰度阈值
	 * @ param  isSaved              是否保存极线校正后的图像
	 * @ return  双目标定的Q矩阵
	*/
	cv::Mat BinocularPhaseRecovery(cv::Mat& correction_left, 
								   cv::Mat& correction_right, 
								   std::vector<cv::Mat>& StripesLevel_Left,
								   std::vector<cv::Mat>& StripesLevel_Right,
								   int GrayThreshold, bool isSaved = false);


	/*
	 * @ brief  双目结构光相位恢复,这个函数无需重复传入格雷码图像，只传入一次曝光下的格雷码图像即可
	 * @ param  correction_left      极线矫正后的左相机展开相位
	 * @ param  correction_right     极线矫正后的右相机展开相位
	 * @ param  Q                     双目标定的Q矩阵
	*/
	bool BinocularPhaseRecovery(cv::Mat& correction_left, cv::Mat& correction_right ,cv::Mat& Q);



	/*
	 * @ brief  读取已经计算完的左右相机极线矫正后的图像以及Q矩阵
	 * @ param  path                 存放图片的路径
	 * @ param  correction_left      极线矫正后的左相机展开相位
	 * @ param  correction_right     极线矫正后的右相机展开相位
	 * @ param  Q                     双目标定的Q矩阵
	*/
	void readRecoveryPhaseImageAndQ(std::string& path, cv::Mat& correction_left, cv::Mat& correction_right, cv::Mat& Q);





	//Binarise the image Change the pixel value to 0 or 255
	std::vector<cv::Mat> changeValue(std::vector<cv::Mat>& grays,
									 cv::Mat& threshold,
									 int width,
									 int height);


	//将float类型转变为double类型了说
	cv::Mat Get_wrappedphase_double();
	cv::Mat getunwtapped_double();
	cv::Mat RemoveBackground_double();


private:
	/*
	 * @ brief  以下函数用于相位展开计算
	 * @ struct    TrigValues                         预计算三角函数值的结构体
	 * @ function  generateQualityMask                生成质量掩码
	 * @ function  decodeGrayCodeInline               解码格雷码获取条纹级次
	 * @ function  unwrapPhaseInline                  相位展开
	*/

	// 预计算三角函数值的结构体
	struct TrigValues {
		std::vector<float> sinVals, cosVals;
		float shiftVal, invPhaseImgsNum;

		TrigValues(int phaseImgsNum) {
			shiftVal = static_cast<float>(CV_2PI) / phaseImgsNum;
			invPhaseImgsNum = 1.0f / phaseImgsNum;
			sinVals.resize(phaseImgsNum);
			cosVals.resize(phaseImgsNum);

			for (int k = 0; k < phaseImgsNum; ++k) {
				const float angle = k * shiftVal;
				sinVals[k] = sin(angle);
				cosVals[k] = cos(angle);
			}
		}
	};

	// 生成质量掩码
	cv::Mat generateQualityMask(const cv::Mat& modulation, const cv::Mat& averageIntensity) {
		cv::Mat phaseError = PhaseErrorAnalyzer::calculateTheoreticalPhaseError(modulation, averageIntensity, 1.0f);
		cv::Mat qualityMask = cv::Mat::zeros(M_Height, M_Width, CV_8UC1);

		cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range) {
			for (int i = range.start; i < range.end; ++i) {
				const float* errorPtr = phaseError.ptr<float>(i);
				uchar* qualityPtr = qualityMask.ptr<uchar>(i);

				for (int j = 0; j < M_Width; ++j) {
					qualityPtr[j] = (errorPtr[j] < 0.1) ? 255 :
						(errorPtr[j] < 0.5) ? 128 : 0;
				}
			}
			});

		return qualityMask;
	}

	// 内联辅助函数，减少函数调用开销
	inline std::pair<int, int> decodeGrayCodeInline(const std::vector<const float*>& grayImgPtrs, 
													int j,
													float threshold) {
		int K1 = 0, tempVal = 0;
		// 解码格雷码获取K1，初始值为零，因为与零异或保持不变
		for (int k = 0; k < M_GrayImgsNum - 1; ++k) {
			tempVal ^= (grayImgPtrs[k][j] > threshold) ? 1 : 0;
			K1 = (K1 << 1) + tempVal;
		}
		// 处理互补格雷码获取K2
		tempVal ^= (grayImgPtrs[M_GrayImgsNum - 1][j] > threshold) ? 1 : 0;
		const int K2 = ((K1 << 1) + tempVal + 1) >> 1; // 使用位移代替除法
		return { K1, K2 };
	}

	inline std::pair<float, int> unwrapPhaseInline(float wrappedPhase, int K1, int K2) {
		if (wrappedPhase <= CV_PI * 0.5f) {
			return { wrappedPhase + CV_2PI * K2, K2 };
		}
		else if (wrappedPhase >= CV_PI * 1.5f) {
			return { wrappedPhase + CV_2PI * (K2 - 1), K2 - 1 };
		}
		else {
			return { wrappedPhase + CV_2PI * K1, K1 };
		}
	}


private:
	//是否使用多重曝光标志位
	bool M_IsMulExposure = false;
	//多重曝光次数
	int M_ExposureNums;

	//投影的辅助图片的数量
	int M_AncillaryPatternNum;

	//左相机图片存放路径
	std::string M_ImagesPathLeft;
	//右相机图片存放路径
	std::string M_ImagesPathRight;
	//相机标定参数存放路径
	std::string M_CalibrationResultPath;

	//格雷码图像的数量
	int M_GrayImgsNum;
	//条纹图像的数量
	int M_PhaseImgsNum;
	//图像的尺寸
	int M_Width;
	int M_Height;

	//存放左相机图像序列
	std::vector<cv::Mat> M_Patterns_Left;
	//存放右相机图像序列
	std::vector<cv::Mat> M_Patterns_Right;


	std::vector<cv::Mat> M_AncillaryPatterns;


	//存放所有图像
	std::vector<cv::Mat> M_Allp;
	//存放二值化后的图像
	std::vector<cv::Mat> M_BinaryP;
	//存放包裹相位
	cv::Mat M_pha;
	cv::Mat M_wrapped;
	//存放K1矩阵，K2矩阵
	cv::Mat M_K1;
	cv::Mat M_K2;
	//存放解包裹图像
	cv::Mat M_unwrapped;
	//置信度
	cv::Mat M_confidence;

	//文件操作
	FileManipulation M_GetImage;
	//是否需要重命名图片
	bool M_isRename;



};

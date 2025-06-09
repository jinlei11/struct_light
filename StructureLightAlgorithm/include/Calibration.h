#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <time.h>
#include <filesystem>
#include <algorithm>
#include <string>
#include <functional>
#include <opencv2/calib3d/calib3d.hpp>
#include "FileManipulation.h"


//Calibration plate models
enum CalibrationPlate {
	Plate_GP050 = 0,
	Plate_doubleCircle = 1,
	Plate_bigCircle = 2,
	Plate_Big_chess = 3,
	Plate_smallCircle = 4,
	Plate_chess811 = 5
};

//极线矫正精度评估
struct EpipolarErrorStats {
	double mean_error;
	double std_error;
	double max_error;
	int num_matches;
};




class CameraCalibration {

public:
	/*
	 * @ brief  相机标定类构造函数
	 * @ prarm  Plate                标定板的具体型号
	 * @ prarm  filepath_left        左相机图像路径
	 * @ param  filepath_Right       右相机图像路径
	 * @ param  ResultPath           标定结果保存路径
	 * @ param  isDrawResult         是否显示标定结果标志位
	*/
	CameraCalibration(enum CalibrationPlate Plate, 
					  const std::string& filepath_left,
					  const std::string& filepath_Right,
					  const std::string& ResultPath,
		              bool isDrawResult);

	/*
	 * @ brief  相机标定类构造函数
	 * @ prarm  board_size           标定板的角点个数，行数和列数
	 * @ param  squareSize           标定板的角点间距
	 * @ prarm  filepath_left        左相机图像路径
	 * @ param  filepath_Right       右相机图像路径
	 * @ param  ResultPath           标定结果保存路径
	 * @ param  CalibrationPlateMode 标定板类型
	 * @ param  isDrawResult         是否显示标定结果标志位
	*/
	CameraCalibration(const cv::Size& board_size,
					  const cv::Size& squareSize, 
					  const std::string& filepath_left,
					  const std::string& filepath_Right,
					  const std::string& ResultPath,
					  int CalibrationPlateMode,
					  bool isDrawResult);
	
	/*
	 * @ brief  获取圆形标定板像素坐标
	 * @ param  IsDrawResult         是否显示标定结果标志位
	*/
	void GetImagePointsCircle(bool IsDrawResult);

	/*
	 * @ brief  获取圆形标定板像素坐标,反转标定板颜色
	 * @ param  IsDrawResult         是否显示标定结果标志位
	*/
	void GetImagePointsCircleInvert(bool IsDrawResult);


	/*
	 * @ brief  获取棋盘格标定板像素坐标
	 * @ param  IsDrawResult         是否显示标定结果标志位
	*/
	void GetImagePointsChess(bool IsDrawResult);

	/*
	 * @ brief  获取世界坐标点
	*/
	void GetObjectPoints();

	/*
	 * @ brief  相机标定
	*/
	cv::Mat MyCameraCalibrationAndSave(); 

	/*
	 * @ brief  相机标定
	*/
	cv::Mat M_MyCameraCalibrationAndSave();




private:
	//文件操作
	FileManipulation M_GetImage;

	//标定图片存放路径：左相机
	std::string M_CalibrationPicsPath_Left;
	//标定图片存放路径：右相机
	std::string M_CalibrationPicsPath_Right;

	//Stores the image acquired by the left camera
	std::vector<cv::Mat>M_imgLs;
	//Storing images acquired by the right camera
	std::vector<cv::Mat>M_imgRs;

	//Store the left camera image coordinate points
	std::vector<std::vector<cv::Point2f>> M_imgsPointsL;
	//Store the right camera image coordinate points
	std::vector<std::vector<cv::Point2f>> M_imgsPointsR;

	//store world coordinate system coordinate points
	std::vector<std::vector<cv::Point3f>> M_objectPoints;

	//Number of corner points of a checkerboard grid or circular calibration board (Width, Height)
	cv::Size M_board_size;

	//Checkerboard or circular calibration board with square dimensions or circle spacing
	cv::Size M_squareSize;

	//标定结果路径
	std::string M_CalibrationResultPath;

	//是否显示角点识别结果
	bool M_IsDrawResult;


	//Distinguishing calibration plate types
	int M_CalibPlateType;

private:

	/*
	 * @ brief  寻找圆形标定板角点
	 * @ prarm  imgs                 圆形标定板图像容器
	 * @ param  imgPoints            存放找到的圆形标定板角点坐标
	 * @ prarm  boardSize            标定板角点数
	 * @ param  blob                 滤波器
	 * @ param  showCorners         是否显示标定结果标志位
	*/
	void ProcessImagePointsCircle(const std::vector<cv::Mat>& imgs,
								 std::vector<std::vector<cv::Point2f>>& imgPoints,
								 const cv::Size& boardSize,
								 const cv::Ptr<cv::FeatureDetector>& blob,
								 bool showCorners);

	/*
	 * @ brief  寻找圆形标定板角点,反转颜色
	 * @ prarm  imgs                 圆形标定板图像容器
	 * @ param  imgPoints            存放找到的圆形标定板角点坐标
	 * @ prarm  boardSize            标定板角点数
	 * @ param  blob                 滤波器
	 * @ param  showCorners          是否显示标定结果标志位
	*/
	void ProcessImagePointsCircleInvert(const std::vector<cv::Mat>& imgs,
										std::vector<std::vector<cv::Point2f>>& imgPoints,
										const cv::Size& boardSize,
										const cv::Ptr<cv::FeatureDetector>& blob,
										bool showCorners);


	/*
	 * @ brief  寻找棋盘格标定板角点
	 * @ prarm  imgs                 棋盘格标定板图像容器
	 * @ param  imgPoints            存放找到的棋盘格标定板角点坐标
	 * @ prarm  boardSize            标定板角点数
	 * @ param  showCorners          是否显示标定结果标志位
	*/
	void ProcessImagePointsChess(const std::vector<cv::Mat>& imgs,
								 std::vector<std::vector<cv::Point2f>>& imgPoints,
								 const cv::Size& boardSize,
								 bool showCorners);


	/**********************************   标定相关变量   *******************************************/

	//MatrixL:  左相机内参矩阵
	//distL:    左相机畸变系数矩阵
	//MatrixR:  右相机内参矩阵
	//distR：   右相机畸变系数矩阵
	//R：       相机之间的旋转矩阵R;
	//T：       相机之间的平移矩阵T
	//E：       相机之间的本征矩阵（旋转+平移）
	//F:        相机之间的基本矩阵  本征矩阵+两个相机之间的内参F:
	//Q:        4x4的视差深度映射矩阵
	//R1        矫正旋转矩阵,将第一个相机坐标系下未矫正的点变换到第一个相机矫正坐标系下，即 R_{左矫正坐标系}{左未矫正坐标系}
	//R2:       矫正旋转矩阵,将第二个相机坐标系下未矫正的点变换到第二个相机矫正坐标系下，即 R_{右矫正坐标系}{右未矫正坐标系}
	//P1:       3x4左相机投影矩阵。将左矫正坐标系下的点投影到左矫正坐标系图像平面坐标系
	//P2:        3x4右相机投影矩阵。将左矫正坐标系下的点投影到右矫正坐标系图像平面坐标系
	cv::Mat M_MatrixL, M_MatrixR, M_DistL, M_DistR;
	std::vector<cv::Mat> M_rvecsL, M_tvecsL, M_rvecsR, M_tvecsR;
	cv::Mat M_R, M_T, M_E, M_F;
	cv::Mat M_R1, M_R2, M_P1, M_P2, M_Q;
	//map11 输出的X坐标重映射参数
	//map12 输出的Y坐标重映射参数
	//map21 输出的X坐标重映射参数
	//map22 输出的Y坐标重映射参数
	cv::Mat M_map11, M_map12, M_map21, M_map22;


	/**********************************   标定相关变量   *******************************************/

	/*
	 * @ brief  相机标定相关函数
	 * @ function   ReadCalibrationImages         读取标定图像
	 * @ function   ExtractFeaturePoints          提取角点并计算世界坐标点
	 * @ function   CalibrateSingleCameras        应用opencv库进行相机标定
	 * @ function   PrintReprojectionError        应用opencve计算重投影误差
	 * @ function   StereoCalibration             双目相机标定
	 * @ function   StereoRectificationAndRemap   极线矫正并可视化结果
	 * @ function   PrepareRectifiedImage         极线矫正
	 * @ function   DrawEpipolarLines             极线矫正后的可视化：绘制线
	 * @ function   DrawEpipolarLines             极线矫正后的可视化：展示或保存图像
	*/

	void ReadCalibrationImages();

	bool ExtractFeaturePoints();

	void CalibrateSingleCameras();

	void PrintReprojectionError(const std::vector<std::vector<cv::Point3f>>& objPts,
								const std::vector<std::vector<cv::Point2f>>& imgPts,
								const std::vector<cv::Mat>& rvecs,
								const std::vector<cv::Mat>& tvecs,
								const cv::Mat& camMat, const cv::Mat& distCoeffs,
								const std::string& camName);

	void StereoCalibration();

	void StereoRectificationAndRemap();

	void DrawEpipolarLines(cv::Mat& img);

	void ShowAndSaveResult(const cv::Mat& result, int index);



	/*
	 * @ brief  极线矫正精度评估函数
	 * @ function   evaluateEpipolarError         读取标定图像
	 * @ function   ExtractFeaturePoints          提取角点并计算世界坐标点
	 * @ function   CalibrateSingleCameras        应用opencv库进行相机标定
	 * @ function   PrintReprojectionError        应用opencve计算重投影误差
	*/

	EpipolarErrorStats evaluateEpipolarError(const cv::Mat& img_left, const cv::Mat& img_right);
};
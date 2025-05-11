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
	Plate_smallCircle = 4
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


};
#pragma once

#include <opencv2/opencv.hpp>
#include <iostream>
#include <opencv2/core/hal/interface.h>
#include <opencv2/core/mat.hpp>
#include <istream>
#include <ostream>
#include <fstream>


class CloudPointsGeneration {
public:
	/*
	 * @ brief  点云数据生成类构造函数
	 * @ param   CloudPointsSavePath       点云数据保存路径
	*/
	CloudPointsGeneration(const std::string& CloudPointsSavePath);


	/*
	 * @ brief  计算视差图：相同相位值在左右相机展开相位图像的像素坐标之差
	 * @ param   UnwrappedLeft       左相机展开相位
	 * @ param   UnwrappedRight      右相机展开相位
	 * @ param   disparity           视差图
	 * @ param   maxDisparity        匹配相位时的最大相位差
	*/
	cv::Mat CalculateDisparity(cv::Mat& UnwrappedLeft, 
							   cv::Mat& UnwrappedRight,  
							   cv::Mat& disparity, 
							   int maxDisparity);


	/*
	 * @ brief  计算深度图：根据视差图，结合相机参数计算该视差值下对应的三维真实深度
	 * @ param   disparity        视差图
	 * @ param   Q                双目相机标定的Q矩阵
	 * @ param   mindep           深度图中阈值的最小值
	 * @ param   maxdep           深度图中阈值的最大值
	 * @ param   minDisparity     视差图中阈值的最小值
	 * @ param   maxDisparity     视差图中阈值的最大值
	*/
	cv::Mat CalculateDepthMap(cv::Mat& Disparity, 
							  cv::Mat& Q, 
							  float mindep, 
							  float maxdep,
							  float minDisparity,
							  float maxDisparity);


	cv::Mat CalculateDepthMap_(cv::Mat& Disparity,
							   cv::Mat& Q,
							   cv::Mat& DepthMap,
							   float mindep,
							   float maxdep,
							   float minDisparity,
							   float maxDisparity);


	/*
	 * @ brief  生成三位点云坐标数据，并保存点云数据位txt文件
	 * @ param   DepthMap          根据视差结果计算出来的真实世界坐标系的Z坐标
	 * @ param   Q                 双目相机标定的Q矩阵
	 * @ param   CloudPointsPath   保存点云文件的路径
	*/
	void SaveCloudPointsToTxt(cv::Mat& DepthMap, 
							  cv::Mat& Q, 
							  std::string& CloudPointsPath);


	~CloudPointsGeneration();

private:
	//Storage Parallax
	cv::Mat M_Disparity;

	//点云保存路径
	std::string M_CloudPointsSavePath;
};
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
	 * @ brief  �������������๹�캯��
	 * @ param   CloudPointsSavePath       �������ݱ���·��
	*/
	CloudPointsGeneration(const std::string& CloudPointsSavePath);


	/*
	 * @ brief  �����Ӳ�ͼ����ͬ��λֵ���������չ����λͼ�����������֮��
	 * @ param   UnwrappedLeft       �����չ����λ
	 * @ param   UnwrappedRight      �����չ����λ
	 * @ param   disparity           �Ӳ�ͼ
	 * @ param   maxDisparity        ƥ����λʱ�������λ��
	*/
	cv::Mat CalculateDisparity(cv::Mat& UnwrappedLeft, 
							   cv::Mat& UnwrappedRight,  
							   cv::Mat& disparity, 
							   int maxDisparity);


	/*
	 * @ brief  �������ͼ�������Ӳ�ͼ������������������Ӳ�ֵ�¶�Ӧ����ά��ʵ���
	 * @ param   disparity        �Ӳ�ͼ
	 * @ param   Q                ˫Ŀ����궨��Q����
	 * @ param   mindep           ���ͼ����ֵ����Сֵ
	 * @ param   maxdep           ���ͼ����ֵ�����ֵ
	 * @ param   minDisparity     �Ӳ�ͼ����ֵ����Сֵ
	 * @ param   maxDisparity     �Ӳ�ͼ����ֵ�����ֵ
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
	 * @ brief  ������λ�����������ݣ��������������λtxt�ļ�
	 * @ param   DepthMap          �����Ӳ��������������ʵ��������ϵ��Z����
	 * @ param   Q                 ˫Ŀ����궨��Q����
	 * @ param   CloudPointsPath   ��������ļ���·��
	*/
	void SaveCloudPointsToTxt(cv::Mat& DepthMap, 
							  cv::Mat& Q, 
							  std::string& CloudPointsPath);


	~CloudPointsGeneration();

private:
	//Storage Parallax
	cv::Mat M_Disparity;

	//���Ʊ���·��
	std::string M_CloudPointsSavePath;
};
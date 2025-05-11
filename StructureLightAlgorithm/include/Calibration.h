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
	 * @ brief  ����궨�๹�캯��
	 * @ prarm  Plate                �궨��ľ����ͺ�
	 * @ prarm  filepath_left        �����ͼ��·��
	 * @ param  filepath_Right       �����ͼ��·��
	 * @ param  ResultPath           �궨�������·��
	 * @ param  isDrawResult         �Ƿ���ʾ�궨�����־λ
	*/
	CameraCalibration(enum CalibrationPlate Plate, 
					  const std::string& filepath_left,
					  const std::string& filepath_Right,
					  const std::string& ResultPath,
		              bool isDrawResult);

	/*
	 * @ brief  ����궨�๹�캯��
	 * @ prarm  board_size           �궨��Ľǵ����������������
	 * @ param  squareSize           �궨��Ľǵ���
	 * @ prarm  filepath_left        �����ͼ��·��
	 * @ param  filepath_Right       �����ͼ��·��
	 * @ param  ResultPath           �궨�������·��
	 * @ param  CalibrationPlateMode �궨������
	 * @ param  isDrawResult         �Ƿ���ʾ�궨�����־λ
	*/
	CameraCalibration(const cv::Size& board_size,
					  const cv::Size& squareSize, 
					  const std::string& filepath_left,
					  const std::string& filepath_Right,
					  const std::string& ResultPath,
					  int CalibrationPlateMode,
					  bool isDrawResult);
	
	/*
	 * @ brief  ��ȡԲ�α궨����������
	 * @ param  IsDrawResult         �Ƿ���ʾ�궨�����־λ
	*/
	void GetImagePointsCircle(bool IsDrawResult);


	/*
	 * @ brief  ��ȡ���̸�궨����������
	 * @ param  IsDrawResult         �Ƿ���ʾ�궨�����־λ
	*/
	void GetImagePointsChess(bool IsDrawResult);

	/*
	 * @ brief  ��ȡ���������
	*/
	void GetObjectPoints();

	/*
	 * @ brief  ����궨
	*/
	cv::Mat MyCameraCalibrationAndSave(); 

	/*
	 * @ brief  ����궨
	*/
	cv::Mat M_MyCameraCalibrationAndSave();




private:
	//�ļ�����
	FileManipulation M_GetImage;

	//�궨ͼƬ���·���������
	std::string M_CalibrationPicsPath_Left;
	//�궨ͼƬ���·���������
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

	//�궨���·��
	std::string M_CalibrationResultPath;

	//�Ƿ���ʾ�ǵ�ʶ����
	bool M_IsDrawResult;


	//Distinguishing calibration plate types
	int M_CalibPlateType;

private:

	/*
	 * @ brief  Ѱ��Բ�α궨��ǵ�
	 * @ prarm  imgs                 Բ�α궨��ͼ������
	 * @ param  imgPoints            ����ҵ���Բ�α궨��ǵ�����
	 * @ prarm  boardSize            �궨��ǵ���
	 * @ param  blob                 �˲���
	 * @ param  showCorners         �Ƿ���ʾ�궨�����־λ
	*/
	void ProcessImagePointsCircle(const std::vector<cv::Mat>& imgs,
								 std::vector<std::vector<cv::Point2f>>& imgPoints,
								 const cv::Size& boardSize,
								 const cv::Ptr<cv::FeatureDetector>& blob,
								 bool showCorners);

	/*
	 * @ brief  Ѱ�����̸�궨��ǵ�
	 * @ prarm  imgs                 ���̸�궨��ͼ������
	 * @ param  imgPoints            ����ҵ������̸�궨��ǵ�����
	 * @ prarm  boardSize            �궨��ǵ���
	 * @ param  showCorners          �Ƿ���ʾ�궨�����־λ
	*/
	void ProcessImagePointsChess(const std::vector<cv::Mat>& imgs,
								 std::vector<std::vector<cv::Point2f>>& imgPoints,
								 const cv::Size& boardSize,
								 bool showCorners);


};
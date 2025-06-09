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

//���߽�����������
struct EpipolarErrorStats {
	double mean_error;
	double std_error;
	double max_error;
	int num_matches;
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
	 * @ brief  ��ȡԲ�α궨����������,��ת�궨����ɫ
	 * @ param  IsDrawResult         �Ƿ���ʾ�궨�����־λ
	*/
	void GetImagePointsCircleInvert(bool IsDrawResult);


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
	 * @ brief  Ѱ��Բ�α궨��ǵ�,��ת��ɫ
	 * @ prarm  imgs                 Բ�α궨��ͼ������
	 * @ param  imgPoints            ����ҵ���Բ�α궨��ǵ�����
	 * @ prarm  boardSize            �궨��ǵ���
	 * @ param  blob                 �˲���
	 * @ param  showCorners          �Ƿ���ʾ�궨�����־λ
	*/
	void ProcessImagePointsCircleInvert(const std::vector<cv::Mat>& imgs,
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


	/**********************************   �궨��ر���   *******************************************/

	//MatrixL:  ������ڲξ���
	//distL:    ���������ϵ������
	//MatrixR:  ������ڲξ���
	//distR��   ���������ϵ������
	//R��       ���֮�����ת����R;
	//T��       ���֮���ƽ�ƾ���T
	//E��       ���֮��ı���������ת+ƽ�ƣ�
	//F:        ���֮��Ļ�������  ��������+�������֮����ڲ�F:
	//Q:        4x4���Ӳ����ӳ�����
	//R1        ������ת����,����һ���������ϵ��δ�����ĵ�任����һ�������������ϵ�£��� R_{���������ϵ}{��δ��������ϵ}
	//R2:       ������ת����,���ڶ����������ϵ��δ�����ĵ�任���ڶ��������������ϵ�£��� R_{�ҽ�������ϵ}{��δ��������ϵ}
	//P1:       3x4�����ͶӰ���󡣽����������ϵ�µĵ�ͶӰ�����������ϵͼ��ƽ������ϵ
	//P2:        3x4�����ͶӰ���󡣽����������ϵ�µĵ�ͶӰ���ҽ�������ϵͼ��ƽ������ϵ
	cv::Mat M_MatrixL, M_MatrixR, M_DistL, M_DistR;
	std::vector<cv::Mat> M_rvecsL, M_tvecsL, M_rvecsR, M_tvecsR;
	cv::Mat M_R, M_T, M_E, M_F;
	cv::Mat M_R1, M_R2, M_P1, M_P2, M_Q;
	//map11 �����X������ӳ�����
	//map12 �����Y������ӳ�����
	//map21 �����X������ӳ�����
	//map22 �����Y������ӳ�����
	cv::Mat M_map11, M_map12, M_map21, M_map22;


	/**********************************   �궨��ر���   *******************************************/

	/*
	 * @ brief  ����궨��غ���
	 * @ function   ReadCalibrationImages         ��ȡ�궨ͼ��
	 * @ function   ExtractFeaturePoints          ��ȡ�ǵ㲢�������������
	 * @ function   CalibrateSingleCameras        Ӧ��opencv���������궨
	 * @ function   PrintReprojectionError        Ӧ��opencve������ͶӰ���
	 * @ function   StereoCalibration             ˫Ŀ����궨
	 * @ function   StereoRectificationAndRemap   ���߽��������ӻ����
	 * @ function   PrepareRectifiedImage         ���߽���
	 * @ function   DrawEpipolarLines             ���߽�����Ŀ��ӻ���������
	 * @ function   DrawEpipolarLines             ���߽�����Ŀ��ӻ���չʾ�򱣴�ͼ��
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
	 * @ brief  ���߽���������������
	 * @ function   evaluateEpipolarError         ��ȡ�궨ͼ��
	 * @ function   ExtractFeaturePoints          ��ȡ�ǵ㲢�������������
	 * @ function   CalibrateSingleCameras        Ӧ��opencv���������궨
	 * @ function   PrintReprojectionError        Ӧ��opencve������ͶӰ���
	*/

	EpipolarErrorStats evaluateEpipolarError(const cv::Mat& img_left, const cv::Mat& img_right);
};
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
	 * @ brief  ��λ��Ϣ�����๹�캯��������ڵ����ع�
	 * @ param  ImagepathLeft             �������������ͶӰͼ����·��
	 * @ param  ImagepathRight            �������������ͶӰͼ����·��
	 * @ param  CalibrationResultPath     �������궨���
	 * @ param  GrayImgsNum               ������ͼ�������
	 * @ param  PhaseImgsNum              ����ͼ������
	 * @ param  width                     ͼ��Ŀ��
	 * @ param  height                    ͼ��ĸ߶�
	 * @ param  AncillaryPatternNum       ͶӰ�ĸ���ͼƬ������
	 * @ param  isRename                  �Ƿ���Ҫ������ͼƬ
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
	 * @ brief  ��λ��Ϣ�����๹�캯��������ڵ����ع�
	 * @ param  imagepath            ��������ͶӰͼ����·��
	 * @ param  GrayImgsNum          ������ͼ�������
	 * @ param  PhaseImgsNum         ����ͼ������
	 * @ param  Form                 ������ͺţ��ɴ�������ͼ��Ŀ�Ⱥ͸߶�
	 * @ param  AncillaryPatternNum  ͶӰ�ĸ���ͼƬ������
 	 * @ param  isRename             �Ƿ���Ҫ������ͼƬ
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
	 * @ brief  ��λ��Ϣ�����๹�캯������������ع�
	 * @ param  ImagepathLeft             �������������ͶӰͼ����·��
	 * @ param  ImagepathRight            �������������ͶӰͼ����·��
	 * @ param  CalibrationResultPath     �������궨���
	 * @ param  GrayImgsNum               ������ͼ�������
	 * @ param  PhaseImgsNum              ����ͼ������
	 * @ param  width                     ͼ��Ŀ��
	 * @ param  height                    ͼ��ĸ߶�
	 * @ param  ExposureNums              �����ع���ع����
	 * @ param  AncillaryPatternNum       ͶӰ�ĸ���ͼƬ������
	 * @ param  isRename                  �Ƿ���Ҫ������ͼƬ
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
	 * @ brief  ��λ��Ϣ�����๹�캯������������ع�
	 * @ param  ImagepathLeft             �������������ͶӰͼ����·��
	 * @ param  ImagepathRight            �������������ͶӰͼ����·��
	 * @ param  CalibrationResultPath     �������궨���
	 * @ param  GrayImgsNum               ������ͼ�������
	 * @ param  PhaseImgsNum              ����ͼ������
	 * @ param  Form                      ������ͺţ��ɴ�������ͼ��Ŀ�Ⱥ͸߶�
	 * @ param  ExposureNums              �����ع���ع����
	 * @ param  AncillaryPatternNum       ͶӰ�ĸ���ͼƬ������
	 * @ param  isRename                  �Ƿ���Ҫ������ͼƬ
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
	 * @ brief  �Ե���ͶӰͼ��������λ����
	 * @ param  ProjectionPatterns   ����ͶӰͼ��
	 * @ param  Phase                ������Ľ������λ
	 * @ param  regulatory           ���ƶ�B��ֵ
	*/
	void SingleDePhase(const std::vector<cv::Mat>& ProjectionPatterns, 
					   cv::Mat& Phase,
					   int regulatory);


	/*
	 * @ brief  �Զ���ͶӰͼ��������λ����,�ú���ֱ�ӵ��õ��ν��㺯����ÿ���ع��ͶӰͼ��ֱ�ӽ��໥������
	 * @ param  ProjectionPatterns   ����ͶӰͼ��
	 * @ param  SovelPacelPhases     ������Ľ������λ����
	 * @ param  ExposureNums         �����ع���ع����
	 * @ param  regulatory           ���ƶ�B��ֵ
	*/
	void MultiDePhase(const std::vector<cv::Mat>& ProjectionPatterns, 
					  std::vector<cv::Mat>& SovelPacelPhases, 
					  int ExposureNums, 
					  int regulatory);


	/*
	 * @ brief  �Զ����ع�Ľ������λ�����ں�
	 * @ param  SovelPacelPhases     ������Ľ������λ����
	 * @ param  ExposureNums         �����ع���ع����
	 * @ param  IntegrationPhase     �ں�֮�����λ
	 * @ param  GrayThresholds       �Ҷ���ֵ
	 * @ param  position             ��������������������
	*/
	void PhasesIntegration(std::vector<cv::Mat>& SovelPacelPhases, 
						   int ExposureNums,
						   cv::Mat& IntegrationPhase,
						   int GrayThresholds,
		                   CameraPosition position);


	/*
	 * @ brief  ���߽�������
	 * @ param  filename                  �������궨�������ļ�·��
	 * @ param  unwrappedLeft             �������������Ľ������λ
	 * @ param  unwrappedRight            �������������Ľ������λ
	 * @ param  jiaozhengL                ���߽������������������λ
	 * @ param  jiaozhengR                ���߽������������������λ
	 * @ param  isShowPolarrectification  �Ƿ���ʾ���߽����Ľ��
	 * @ return  ˫Ŀ�궨��Q����
	*/
	cv::Mat readCameraInfoAndPolarrectification(const std::string& filename,
												cv::Mat& unwrappedLeft,
												cv::Mat& unwrappedRight,
												cv::Mat& jiaozhengL,
												cv::Mat& jiaozhengR,
												bool isShowPolarrectification = false);

	/*
	 * @ brief  ˫Ŀ�ṹ����λ�ָ�
	 * @ param  correction_left      ���߽�����������չ����λ
	 * @ param  correction_right     ���߽�����������չ����λ
	 * @ param  GrayThreshold        �Ҷ���ֵ
	 * @ param  isSaved              �Ƿ񱣴漫��У�����ͼ��
	 * @ return  ˫Ŀ�궨��Q����
	*/
	cv::Mat BinocularPhaseRecovery(cv::Mat& correction_left, cv::Mat& correction_right, int GrayThreshold, bool isSaved = false);


	/*
	 * @ brief  ˫Ŀ�ṹ����λ�ָ�,������������ظ����������ͼ��ֻ����һ���ع��µĸ�����ͼ�񼴿�
	 * @ param  correction_left      ���߽�����������չ����λ
	 * @ param  correction_right     ���߽�����������չ����λ
	 * @ param  Q                     ˫Ŀ�궨��Q����
	*/
	bool BinocularPhaseRecovery(cv::Mat& correction_left, cv::Mat& correction_right ,cv::Mat& Q);



	/*
	 * @ brief  ��ȡ�Ѿ������������������߽������ͼ���Լ�Q����
	 * @ param  path                 ���ͼƬ��·��
	 * @ param  correction_left      ���߽�����������չ����λ
	 * @ param  correction_right     ���߽�����������չ����λ
	 * @ param  Q                     ˫Ŀ�궨��Q����
	*/
	void readRecoveryPhaseImageAndQ(std::string& path, cv::Mat& correction_left, cv::Mat& correction_right, cv::Mat& Q);





	//Binarise the image Change the pixel value to 0 or 255
	std::vector<cv::Mat> changeValue(std::vector<cv::Mat>& grays,
									 cv::Mat& threshold,
									 int width,
									 int height);

	//Calculation of binarised Gray code image using the mean of the phase-shifted image
	std::vector<cv::Mat> binaryzation(const std::string& path);

	//Another version of binarization 
	//requires passing the number of gray code images and the number of phase-shifted images.
	std::vector<cv::Mat> binaryzation(int NumOfGray,int NumOfPhase);

	//Generate and save binarized images
	std::vector<cv::Mat> binaryzation(int NumOfGray,
									  int NumOfPhase,
									  const std::string& path);


	//Here are some functions for solving the phase of a parcel
	//Obtain the K1 and K2 matrices
	void get_K1_K2();

	//Implementation of decoding Gray Code
	std::string jiema(const std::string& gray);

	//Implement binary to decimal conversion
	int twoToten(std::string gray);
	int BinaryToDecimal(const std::string& gray);

	//Calculate parcel phase
	cv::Mat Get_wrappedphase();
	cv::Mat Get_wrappedphase1();

	//Solve the parcel phase
	cv::Mat getunwrapped(cv::Mat& wrapped_pha);
	cv::Mat getunwrapped();


	////����A
	//inline cv::Mat calculateConfidenceMap() {
	//	cv::Mat confidenceMap(M_Width, M_Width, CV_32FC1, cv::Scalar(0.f));
	//	cv::Mat temp;
	//	for (auto& img : M_Phase) {
	//		img.convertTo(temp, CV_32FC1);
	//		confidenceMap += temp / 4.f;
	//	}

	//	return confidenceMap;
	//}

	//Remove background based on confidence level B
	cv::Mat RemoveBackground(cv::Mat Unwrapped, cv::Mat confidence);
	cv::Mat RemoveBackground();




	//Convert image file names in a folder to order from 0
	void renameImagesInFolder(const std::string& folderPath);

	//Make sure to read the image files in numerical order
	bool numericStringCompare(const std::string& s1, const std::string& s2);




	//��float����ת��Ϊdouble������˵
	cv::Mat Get_wrappedphase_double();
	cv::Mat getunwtapped_double();
	cv::Mat RemoveBackground_double();


private:
	//�Ƿ�ʹ�ö����ع��־λ
	bool M_IsMulExposure = false;
	//�����ع����
	int M_ExposureNums;

	//ͶӰ�ĸ���ͼƬ������
	int M_AncillaryPatternNum;

	//�����ͼƬ���·��
	std::string M_ImagesPathLeft;
	//�����ͼƬ���·��
	std::string M_ImagesPathRight;
	//����궨�������·��
	std::string M_CalibrationResultPath;

	//������ͼ�������
	int M_GrayImgsNum;
	//����ͼ�������
	int M_PhaseImgsNum;
	//ͼ��ĳߴ�
	int M_Width;
	int M_Height;

	//��������ͼ������
	std::vector<cv::Mat> M_Patterns_Left;
	//��������ͼ������
	std::vector<cv::Mat> M_Patterns_Right;


	std::vector<cv::Mat> M_AncillaryPatterns;


	//�������ͼ��
	std::vector<cv::Mat> M_Allp;
	//��Ŷ�ֵ�����ͼ��
	std::vector<cv::Mat> M_BinaryP;
	//��Ű�����λ
	cv::Mat M_pha;
	cv::Mat M_wrapped;
	//���K1����K2����
	cv::Mat M_K1;
	cv::Mat M_K2;
	//��Ž����ͼ��
	cv::Mat M_unwrapped;
	//���Ŷ�
	cv::Mat M_confidence;

	//�ļ�����
	FileManipulation M_GetImage;
	//�Ƿ���Ҫ������ͼƬ
	bool M_isRename;



};

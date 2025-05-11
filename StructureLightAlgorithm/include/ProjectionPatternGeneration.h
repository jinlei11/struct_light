#pragma once
//
//#include <iostream>
//#include <deque>
//#include <math.h>
//#include <opencv2/opencv.hpp>
//
//using namespace std;
////�ı����һλ�ַ�ֵ
//#define change_bit(x) x=((x)=='0'?'1':'0')
//#define pi CV_PI
//
////Projector Models
//enum DLPModel {
//	DLP4710 = 0,
//	DLP6500 = 1,
//	DLP3010 = 2,
//	DLP4500 = 3,
//};
//
//
//class BuildGrayPattern {
//public:
//	// Gray code bits width of the image and height of the image 
//	BuildGrayPattern(int bit, int width, int height);
//	//Constructor, pass in projector model and number of images generated
//	BuildGrayPattern(int bit, enum DLPModel DLP, int direction);
//
//	//Generate Gray code values
//	deque<string> BuildGrayCode();
//
//	//retrieve the opposite of what one wants
//	void qufan(char* a, char* b, int n);
//
//	//Generate Gray Code Images
//	vector<cv::Mat> BuildGrayP();
//
//	//Generate grayscale values of 0 or 255 based on Gray code values.
//	int jugevalue(char a);
//
//	//Storing the generated Gray code
//	std::deque<std::string>m_GrayArr;
//
//	//For storing Gray Code images
//	vector<cv::Mat>M_GrayPicture;
//
//	//Produces  all-black image    Existing in Projection_Pattern file
//	void buildBlack(cv::Size Image_size);
//	//Generate full-white images   Existing in Projection_Pattern file
//	void buildWhite(cv::Size Image_size);
//
//	//Produces  all-black image    Existing in the Specify path file
//	void buildBlack(const std::string& path);
//	//Generate full-white images   Existing in the Specify path file
//	void buildWhite(const std::string& path);
//
//	//Produces an all-black image with the specified gray value, less than 20
//	void buildBlack(const std::string& path, int value);
//	//Generates an all-white image with the specified gray value, greater than 240
//	void buildWhite(const std::string& path, int value);
//
//
//	//Generate and save images
//	void BuildAndSave_GrayPattren(const std::string path);
//
//	//Saving Gray Code Images
//	void SaveGray(const std::string path);
//
//
//
//
//	//�ع�����
//	//���ɸ������м��ַ�����
//	//
//	
//
//
//
//private:
//	//Number of bits used to hold the Gray code image
//	int m_Bit;
//	//Size, width and height of the generated image
//	int m_Width;
//	int m_Height;
//	//Determine two scenarios
//	bool flag;
//	int b;
//
//	int M_direction;
//};
//
//class CreateTW {
//
//public:
//	//constructer
//	CreateTW(float f, int bit, int width, int height, int direction);
//	//constructer
//	CreateTW(float f, int bit, enum DLPModel DLP, int direction);
//
//	//Generate phase-shifted stripe patterns
//	vector<cv::Mat> creatP();
//
//	//Saving phase-shifted streak images
//	void SavePhase(std::string& path);
//
//	//For storing phase-shifted stripe patterns
//	vector<cv::Mat>M_PhasePicture;
//
//	//Generate and save phase shift patterns
//	void BuildAndSave_phase_shift(const std::string& path);
//
//private:
//	//Frequency of the phase-shifted image
//	float f;
//	//image size
//	int m_bit;
//	int m_width;
//	int m_height;
//	//Number of phase shift patterns generated
//	int N = 4;
//	int M_direction;
//};




#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <bitset>
#include "FileManipulation.h"

#define G_GrayCodeMaxNum (1 << 10)
#define G_GrayCodeMaxBit 20    
#define pi CV_PI


enum DLPModel {
	DLP4710 = 0,
	DLP6500 = 1,
	DLP3010 = 2,
	DLP4500 = 3,
};

//--------------------------------------������������-------------------------------------------------//
/*******************************��һ�ַ���**********************************/
//���㲽������Ҫ����nλ���ݵĸ�����ֵ�����һֵΪ 00000 ��n����
//��һ�����ı����ұߵ�λԪֵ��
//�ڶ������ı������һ��Ϊ1��λԪ�����λԪ��
//�����������Ĳ��ظ���һ���͵ڶ�����ֱ�����еĸ����������ϣ����仰˵���Ѿ�����(2 ^ n) - 1 ������


/********************************�ڶ��ַ���******************************************/
//�Ƽ�ʹ�����ַ���
//���ݶ����ƺ͸�����֮���ת����ϵ
//��һ�������ɶ������루�����ʵ�ܼ򵥣�������ڴ洢�ľ��Ƕ����Ƶ���ʽ�������м�λ���ͼ���� 2^n �����Ϳ����ˣ�
//�ڶ��������ݶ������룬ͨ��ת����ϵת���ɸ�����
// ����������һλ��ͬ������õ�һλ��ֵ���Ҳ�����λ������������ͬΪ�㣬��ͬΪһ��
//���������������ɸ�������ֵ


/**********************************�����ַ���***************************************/
//�����������ɸ����뷨
//��ֻ��һλ��ʱ�򣬸�����ֵҪô��0��Ҫô��1��
//������λʱ,�����������ǣ�0
//                         1
// ���0��1����ֻ��һλʱ�ĸ�����ֲ�
//�����õ�             1
//                         0
//Ȼ�����ϰ���������0���°���������1
//               00
//               01
//               11
//               10
//������λʱ���Ҳ���λ�ĸ������������λʱ����ֵ
//               00
//               01
//               11
//               10
//�����о���-----------
//               10
//               11
//               01
//               00
// Ȼ�����ϰ���������0���°���������1  �Դ�����

/************************************�����ַ���********************************************/
//�����������ɸ����뷽�����ԳƷ���





//--------------------------------------������������-------------------------------------------------//
//��һ�ַ���
//��һ�������һ��������ֵ����ߵ�һλ���ֲ��䣬��Ϊ�����Ķ�������ĵ�һλ
//�ڶ������Ը��������ߵڶ�λ��ʼ��������һλ�������������õ��ڶ�λ��������ֵ
//���������ظ�����������ֱ����������һλ����������



//����Ϊ����ϰһ���鹹�����Լ�����͸������Ĳ���������������ֲ�ͬ���������Ϊ�������������д�����麯��
class GrayCodeOperation
{
public:
	GrayCodeOperation();

	virtual ~GrayCodeOperation();

	virtual std::vector<int> Generate_value() = 0;

	virtual std::vector<cv::Mat> Generate_pics() = 0;

	virtual std::vector<cv::Mat> Generate_pics_reverse() = 0;

	virtual std::vector<cv::Mat> GenerateAndSave_pics(std::string& file_path) = 0;

	virtual std::vector<cv::Mat> GenerateAndSave_pics_reverse(std::string& file_path) = 0;

	virtual std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) = 0;

	virtual std::vector<cv::Mat> GenerateAndSave_pics_block(std::string& file_path, int bit, int rowBlock, int colBlock) = 0;

	virtual cv::Mat Decode_Pics(std::string& file_path, int PicNum) = 0;

	virtual cv::Mat Decode_PicsPositiveAndNegative(std::string& file_path, int PicNum) = 0;

	virtual cv::Mat GeneratePatternWithValue(const std::string& path, int value) = 0;



protected:

	inline void ShowImages(const std::vector<cv::Mat>& images, const std::string& windowName = "GraycodePics") {
		cv::namedWindow(windowName, cv::WINDOW_NORMAL);
		for (size_t i = 0; i < images.size(); ++i) {
			if (images[i].empty()) {
				std::cerr << "Warning: Image at index " << i << " is empty, skipping." << std::endl;
				continue;
			}
			cv::imshow(windowName, images[i]);
			cv::waitKey(0);
		}
		cv::destroyWindow(windowName);
	}


};



class GrayMethod_1 : public GrayCodeOperation {

public:
	/*
	 * @ brief  ���෽��һ���캯��
	 * @ prarm  bit      ���ɸ������λ�� 
	 * @ prarm  width    ����ͼƬ�Ŀ��
	 * @ param  height   ����ͼƬ�ĸ߶�
	 * @ return 
	*/
	GrayMethod_1(int bit, int width, int height);


	/*
	 * @ brief   ���ɸ�����ֵ
	 * @ return  ͨ��int���͵����������洢���ɵĸ�����ֵ
	*/
	std::vector<int> Generate_value() override;

	/*
	 * @ brief   ���ɸ�����ͼ��
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> Generate_pics() override;

	/*
	 * @ brief   ���ɸ�����ͼ�񣬰��ջ��ֿ�ķ�ʽ���ɣ�ÿһ�еĸ�����ֵ������һ���ġ�
	 * @ param   ���ɸ������λ�����͹��캯���еĸ�����λ����һ����˼��
	 * @ param   ���ɸ�����������п���
	 * @ param   ���ɸ�����������п���
 	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) override;


private:
	int M_bit;//������λ��
	int M_width;//ͼ����
	int M_height;//ͼ��߶�

	char M_grayCodes[G_GrayCodeMaxNum][G_GrayCodeMaxBit];    //Gray's code sequence
};


//�Ƽ�ʹ�õķ���
class GrayMethod_2 : public GrayCodeOperation {

public:
	/*
	 * @ brief  ���෽�������캯��
	 * @ prarm  bit      ���ɸ������λ��
	 * @ prarm  width    ����ͼƬ�Ŀ��
	 * @ param  height   ����ͼƬ�ĸ߶�
	 * @ return
	*/
	GrayMethod_2(int bit, enum DLPModel DLP, int dircetion);

	/*
	 * @ brief  ���෽���� ��������
	 * @ return
	*/
	~GrayMethod_2();

	/*
	 * @ brief   ���ɸ�����ֵ
	 * @ return  ͨ��int���͵����������洢���ɵĸ�����ֵ
	*/
	std::vector<int> Generate_value() override;


	/*
	 * @ brief   ����ָ���Ҷ�ֵ������ͼ��
	 * @ param  file_path   ����ͼ���·��
	 * @ param  value       ָ���ĻҶ�ֵ
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	cv::Mat GeneratePatternWithValue(const std::string& path, int value);

	/*
	 * @ brief   ����һϵ�о��в�ͬ���ȻҶ�ֵ��ͼ������
	 * @ param   filepath    ͼƬ���·��
	 * @ param   step        ����ͼ��ʱ�Ҷ�ֵ�Ĳ���
	*/
	void generationDiffGrayValuePattern(const std::string& filepath, int step);


	/*
	 * @ brief   ���ɸ�����ͼ��
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> Generate_pics() override;

	/*
	 * @ brief   ���ɸ�����ͼ�񡪡�����ͼ��
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> Generate_pics_reverse() override;


	/*
	 * @ brief   ���ɲ����������ͼ��
	 * @ param  file_path   ���������ͼ���·��
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> GenerateAndSave_pics(std::string& file_path) override;

	/*
	 * @ brief   ���ɲ����������ͼ�񡪡�����ͼ��
	 * @ param  file_path   ���������ͼ���·��
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> GenerateAndSave_pics_reverse(std::string& file_path) override;

	/*
	 * @ brief   ���ɸ�����ͼ�񣬰��ջ��ֿ�ķ�ʽ���ɣ�ÿһ�еĸ�����ֵ������һ���ġ�
	 * @ param  bit        ���ɸ������λ�����͹��캯���еĸ�����λ����һ����˼��
	 * @ param  rowBlock   ���ɸ�����������п���
	 * @ param  colBlock   ���ɸ�����������п���
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) override;

	/*
	 * @ brief   ���ɸ�����ͼ�񣬰��ջ��ֿ�ķ�ʽ���ɣ�ÿһ�еĸ�����ֵ������һ���ġ�
	 * @ param  file_path  ���ͼ���·��
	 * @ param  bit        ���ɸ������λ�����͹��캯���еĸ�����λ����һ����˼��
	 * @ param  rowBlock   ���ɸ�����������п���
	 * @ param  colBlock   ���ɸ�����������п���
	 * @ return  ͨ��Mat���͵����������洢���ɵĸ�����ͼ��
	*/
	std::vector<cv::Mat> GenerateAndSave_pics_block(std::string& file_path, int bit, int rowBlock, int colBlock) override;

	/*
	 * @ brief   ���������ͼ����Ҫ�Ƕ�ֵ�����ͼ��
	 * @ param   ������ͼ���ŵ�·��
	 * @ param   ������ͼ�������
	 * @ return  ͨ��Mat���͵��������洢������ĸ�����ͼ��
	*/
	cv::Mat Decode_Pics(std::string& file_path, int PicNum) override;

	/*
	 * @ brief   ���������ͼ��������ͼ��,�����ͼ����ж�ֵ������
	 * @ param   ������ͼ���ŵ�·��
	 * @ param   ������ͼ�������
	 * @ return  ͨ��Mat���͵��������洢������ĸ�����ͼ��
	*/
	cv::Mat Decode_PicsPositiveAndNegative(std::string& file_path, int PicNum) override;


private:
	//�ļ�������Ա����
	FileManipulation M_ImageOperation;

	int M_bit;//������λ��
	int M_Width;//ͼ����
	int M_Height;//ͼ��߶�
	int M_Direction;//ͼ��ķ���
};



class GrayMethod_3 : public GrayCodeOperation {

public:

	GrayMethod_3();

	std::vector<int> Generate_value() override;

	std::vector<cv::Mat> Generate_pics() override;

};



class PhasePatternOperation {

public:
	/*
	 * @ brief  ����ͼ��������Ĺ��캯��
	 * @ prarm  f         ���Ƶ�ͼ����Ƶ��
	 * @ prarm  bit       ������ͼ���λ��
	 * @ param  DLP       DLPͶӰ�ǵ��ͺţ���Ҫ������������ͼƬ�ĸ߶ȺͿ��
	 * @ param  direction ��������ͼ���ķ���ˮƽ������ֱ
	 * @ param  num       ��������ͼ���Ĳ�����Ĭ��Ϊ�Ĳ�����
 	 * @ return
	*/
	PhasePatternOperation(float f, enum DLPModel DLP, int direction = 0, int num = 4);
	
	~PhasePatternOperation();

	//Generate phase-shifted stripe patterns
	std::vector<cv::Mat> CreatPhasePatterns();


	/*
	 * @ brief  ���ɲ�������������ͼ��
	 * @ param  filepath     ������������ͼ����·��
	 * @ return ���ش��������������ͼ����ͼ������
	*/
	std::vector<cv::Mat> BuildAndSave_PhaseShiftingPattern(const std::string& filepath);

private:

	inline void ShowPhaseImages(const std::vector<cv::Mat>& images, const std::string& windowName = "PhaseShiftingPics") {
		cv::namedWindow(windowName, cv::WINDOW_NORMAL);
		for (size_t i = 0; i < images.size(); ++i) {
			if (images[i].empty()) {
				std::cerr << "Warning: Image at index " << i << " is empty, skipping." << std::endl;
				continue;
			}
			cv::imshow(windowName, images[i]);
			cv::waitKey(0);
		}
		cv::destroyWindow(windowName);
	}


	//�ļ�������Ա����
	FileManipulation M_ImageOperation;

	//For storing phase-shifted stripe patterns
	std::vector<cv::Mat>M_PhasePicture;

	int M_PhaseShiftingNum;//����ͼ��Ĳ���
	int M_GrayCodeBit;//������ͼ���λ��
	float M_F;//����ͼ���Ƶ��
	int M_Direction;//����ͼ��ķ���
	int M_Width;//����ͼ��Ŀ��
	int M_Height;//����ͼ��ĸ߶�


};
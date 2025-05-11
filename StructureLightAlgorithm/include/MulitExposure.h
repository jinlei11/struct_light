#pragma once
//����ļ��н����ֲ�ͬ�����عⷽ��

#include "ALLHeads.h"
#include <opencv2/opencv.hpp>


//��������عⷽ������
class HDR_Reconstruction {


public:

	HDR_Reconstruction();

	~HDR_Reconstruction();

	///*
	// * @ brief   �����ع�ʱ��
	// * @ return  �������ɵ��ع�ʱ������
	//*/
	//virtual std::vector<double> GenerateExposure() = 0;

	///*
	// * @ brief   ��������Ӧ����ͼ������
	// * @ return  �������ɵ�����ͼ������
	//*/
	//virtual std::vector<cv::Mat> GenerateAdaptiveStripePattern() = 0;

	///*
	// * @ brief   ���ɻҶ�ֱ��ͼ
	// * @ return  ���ػҶ�ֱ��ͼͼ��
	//*/
	//virtual cv::Mat GenerationGrayHistogram() = 0;


private:


};


//�ο����ף�������,¬��ʤ,����ï.�߶�̬��Χ��������Ӧ����ͶӰ��������[J].��ѧѧ��,2021,41(19):145-154.
/*
	��������һ�ֻ��ڵݹ������Ӧ����ͶӰ������������ʹ�õ��ǵ�Ŀ����
*/
class HDR_AdaptiveStripePattern : public HDR_Reconstruction {

public:

	HDR_AdaptiveStripePattern();


	~HDR_AdaptiveStripePattern();


	//cv::Mat GenerationGrayHistogram() override;







private:




};



//�ο����ף���ά,���پ�,������,��.��������Ӧ���Ƶĸ߷��������ά���β�������[J].��ѧѧ��,2020,40(05):119-127.
/*
	һ�ֻ���ͼ���ںϺͲ�ֵԤ�������Ӧ����ͶӰ����
*/

class HDR_AdaptiveStripePattern_FW : public HDR_Reconstruction {

public:

	HDR_AdaptiveStripePattern_FW();


	~HDR_AdaptiveStripePattern_FW();


	//cv::Mat GenerationGrayHistogram() override;




private:
 



public:



	/*
	 * @ brief   ������Ĥͼ�����У�����ͼ�����У���ǰͼ��ĻҶ�ֵ������ͼ�������е����ֵ������������ֵ������ֵ
	 *			 ����ǰ��Ĥͼ������ֵ��Ϊ0�������ǰͼ��Ҷ�ֵ�������ֵ�������ֵδ������ֵ������Ĥͼ������ֵ��Ϊ1
	 * @ param   DiffGrayValuePatterns     �������Ĳ�ͬ�Ҷ�ֵ�ľ��ȻҶ�ͼ��
	 * @ param   MaskImages                ���ɵ���Ĥͼ������         
	*/
	void getMaskImages(std::vector<cv::Mat>& DiffGrayValuePatterns, std::vector<cv::Mat>& MaskImages);

	/*
	 * @ brief   ��ȡ��Ч���ȻҶ�ͼ�����У����Ҷ�ֵ������ֵ������ֵ��Ϊ0���������ֲ���
	 * @ param   DiffGrayValuePatterns      �������Ĳ�ͬ�Ҷ�ֵ�ľ��ȻҶ�ͼ��
	 * @ param   GraySequences              ���ɵ���Ч�Ҷ�ͼ������
	*/
	void getEffectiveGraySequences(std::vector<cv::Mat>& DiffGrayValuePatterns, std::vector<cv::Mat>& GraySequences);

	/*
	 * @ brief   ������Ĥͼ�����Ч���ȻҶ�ͼ����������ͶӰͼ�񣺽���������ͼ����������ۼ�
	 * @ param   MaskImages             ���ɵ���Ĥͼ������
	 * @ param   GraySequences          ���ɵ���Ч�Ҷ�ͼ������
	*/
	void generateProjectedImages(std::vector<cv::Mat>& MaskImages, std::vector<cv::Mat>& GraySequences);


	//�����ͼ�񲢲������ͶӰ�Ҷ�ֵ�����ò�ֵԤ����ٲ��ҷ�������ͶӰ�Ҷ�ֵ


	//�������-ͶӰ�Ǳ궨


	//�ܽ᣺��ʵ���վ���������һ����һ���Ҷ�ֵ���������У������ֵ����255���Ǽ��������ĳһ��ֵ


};


//�ο����ף�[1]����,����Ӣ,������,��.����߷������Ľṹ�����β�������[J].�й���ѧ(��Ӣ��),2025,18(01):42-52.
/*
	���������Ӧ���ߣ����㷴Ӧ������淴��ǿ�ȵķ��ն�ͼ��ͨ��ģ��C��ֵ���෽��������Ӧ�ָͬ���ն�����
	����ø��������ĵķ��նȣ����������Ӧ����Ԥ�ⲻͬ�������������ع�ʱ��
*/
class HDR_CameraResponseCurve : public HDR_Reconstruction {

public:

	/*
	 * @ brief   ���������Ӧ���ߵĸ߶�̬��Χ�ؽ������캯��
	 * @ param   filepath               ��Ų�ͬ�ع�ʱ��ͼ����ļ�·��
	 * @ param   Exposures              �ع�ʱ������
	*/
	HDR_CameraResponseCurve(std::string& filepath,
							std::string& resultpath,
							std::vector<float>& Exposures,
							int width,
							int height);


	~HDR_CameraResponseCurve();
	 
public:




	/*
	 * @ brief   ��ȡ�������Ӧ���ߣ����벻ͬ�ع��µ�ͼ������
	 * @ param   DiffExposureImages     ��ͬ�ع��µ�ͼ������
	 * @ param   Exposures              �ع�ʱ������
	 * @ param   ResponseCurveImg       ���������Ӧ�������飬�����նȵĶ����������Ӧ�ĻҶ�ֵ֮��Ĺ�ϵ
	 * PS �����μ����������Ӧ���߲ο���https://blog.csdn.net/wc781708249/article/details/78490976
	 * https://blog.csdn.net/u013049912/article/details/85157402
	*/
	bool getCameraResponseCurve(std::vector<cv::Mat>& DiffExposureImages, 
								std::vector<float>& Exposures,
								cv::Mat& ResponseCurveImg);


	/*
	 * @ brief   ������ն�ͼ�񣬴��뵥��ͼ��,����������ܻ�����
	 * @ param   Image           �����ع������㵽��ͼ��
	 * @ param   Exposure        ���ε��ع�ʱ��        
	 * @ param   gain            ͼ������
	*/
	bool CalculateIrradianceImage(cv::Mat& Image, float Exposurem, float gain);




private:
	//����ع�ͼ�����е��ļ�
	std::string M_ExposureImagesFilepath;

	//��������Ӧͼ����
	std::string M_CamResponseResultpath;

	//�ع�ʱ������
	std::vector<float> M_Exposures;

	int M_Height;
	int M_Width;

};
#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include "FileManipulation.h"

class DecodePattern {

public:
	/*
	 * @ brief  ͶӰͼ�������๹�캯��
	 * @ prarm  GrayBit             ������ͼ����λ��
	 * @ prarm  PhaseShiftingNum    ����ͼ���Ĳ���
	 * @ param  PhaseFrequency      ����ͼ����Ƶ��
	 * @ param  filepath            ͶӰͼ�����浽·��
	 * @ return
	*/
	DecodePattern(int GrayBit,int PhaseShiftingNum,int PhaseFrequency, const std::string& filepath);

	~DecodePattern();

	cv::Mat DecodeGrayPattern();

	cv::Mat DecodePhasePattern();

	cv::Mat DecodeProjectionPattern();

private:
	int M_GrayBit;//������ͼ����λ��
	int M_PhaseShiftingNum;//����ͼ���Ĳ���
	int M_PhaseFrequency;//����ͼ����Ƶ��
	std::string M_Filepath;//ͶӰͼ�����浽·��

};
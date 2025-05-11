#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include "FileManipulation.h"

class DecodePattern {

public:
	/*
	 * @ brief  投影图案解相类构造函数
	 * @ prarm  GrayBit             格雷码图案的位数
	 * @ prarm  PhaseShiftingNum    相移图案的步数
	 * @ param  PhaseFrequency      相移图案的频率
	 * @ param  filepath            投影图案保存到路径
	 * @ return
	*/
	DecodePattern(int GrayBit,int PhaseShiftingNum,int PhaseFrequency, const std::string& filepath);

	~DecodePattern();

	cv::Mat DecodeGrayPattern();

	cv::Mat DecodePhasePattern();

	cv::Mat DecodeProjectionPattern();

private:
	int M_GrayBit;//格雷码图案的位数
	int M_PhaseShiftingNum;//相移图案的步数
	int M_PhaseFrequency;//相移图案的频率
	std::string M_Filepath;//投影图案保存到路径

};
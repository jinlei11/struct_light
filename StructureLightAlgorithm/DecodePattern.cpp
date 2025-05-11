#include"DecodePattern.h"

DecodePattern::DecodePattern(int GrayBit, 
							 int PhaseShiftingNum, 
							 int PhaseFrequency, 
							 const std::string& filepath) 
	:M_GrayBit(GrayBit),M_PhaseShiftingNum(PhaseShiftingNum),M_PhaseFrequency(PhaseFrequency),M_Filepath(filepath)
{

};

DecodePattern::~DecodePattern() {

};

cv::Mat DecodePattern::DecodeGrayPattern() {
	cv::Mat DecodePattern;


	return	DecodePattern;
};

cv::Mat DecodePattern::DecodePhasePattern() {
	cv::Mat DecodePattern;


	return	DecodePattern;
};

cv::Mat DecodePattern::DecodeProjectionPattern() {
	cv::Mat DecodePattern;


	return	DecodePattern;
};
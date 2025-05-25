#include"ImageReprocess.h"



PhaseNoiseFilter::PhaseNoiseFilter(double minArea)
    : areaThreshold(minArea) {}


// ֱ����ԭͼ���޸ģ�С������Ϊ0��
void PhaseNoiseFilter::filterInPlace(cv::Mat& phaseMap) {
    CV_Assert(phaseMap.type() == CV_32FC1);

    // ֻҪ��Ϊ0����Ϊ�Ǻ�ѡ����
    cv::Mat nonZeroMask = phaseMap != 0.0f;

    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(nonZeroMask, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // ����һ����ʱ��Ĥ���ڱ���������
    cv::Mat cleanMask = cv::Mat::zeros(phaseMap.size(), CV_8U);
    for (size_t i = 0; i < contours.size(); ++i) {
        double area = cv::contourArea(contours[i]);
        if (area >= areaThreshold) {
            cv::drawContours(cleanMask, contours, static_cast<int>(i), cv::Scalar(255), cv::FILLED);
        }
    }

    // С��������Ĥ���� 0��ֱ����Ϊ 0
    phaseMap.setTo(0.0f, cleanMask == 0);
}



cv::Mat PhaseNoiseFilter::filter(const cv::Mat& phaseMap) {
    CV_Assert(phaseMap.type() == CV_32FC1);

    // ������������Ĥ��ֻҪ����0������Ϊ�Ǻ�ѡ����
    cv::Mat nonZeroMask = phaseMap != 0.0f;

    // ����ͨ����
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(nonZeroMask, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // ��������Ĥ
    cv::Mat cleanMask = cv::Mat::zeros(phaseMap.size(), CV_8U);
    for (size_t i = 0; i < contours.size(); ++i) {
        double area = cv::contourArea(contours[i]);
        if (area >= areaThreshold) {
            cv::drawContours(cleanMask, contours, static_cast<int>(i), cv::Scalar(255), cv::FILLED);
        }
    }

    // Ӧ����Ĥ
    cv::Mat result = cv::Mat::zeros(phaseMap.size(), phaseMap.type());
    phaseMap.copyTo(result, cleanMask);
    return result;
}

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




void PhaseNoiseFilter::filterPhaseByGradient(const cv::Mat& phaseMap,
                                             float gradThreshold,
                                             cv::Mat& mask,
                                             cv::Mat& filteredPhaseMap){
    CV_Assert(phaseMap.type() == CV_32FC1);

    cv::Mat gradX, gradY, gradMag;

    // ���� x �� y ������ݶ�
    cv::Sobel(phaseMap, gradX, CV_32F, 1, 0, 3);
    cv::Sobel(phaseMap, gradY, CV_32F, 0, 1, 3);

    // �ݶȷ�ֵ
    cv::magnitude(gradX, gradY, gradMag);

    // ��ֵɸѡ
    cv::threshold(gradMag, mask, gradThreshold, 255.0, cv::THRESH_BINARY_INV);
    mask.convertTo(mask, CV_8UC1);

    // ��ѡ����̬ѧ����
    cv::morphologyEx(mask, mask, cv::MORPH_OPEN, cv::Mat(), cv::Point(-1, -1), 1);
    cv::morphologyEx(mask, mask, cv::MORPH_DILATE, cv::Mat(), cv::Point(-1, -1), 1);

    // Ӧ����Ĥ
    filteredPhaseMap = cv::Mat::zeros(phaseMap.size(), phaseMap.type());
    phaseMap.copyTo(filteredPhaseMap, mask);
}
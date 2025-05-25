#pragma once
#include <opencv2/opencv.hpp>

class PhaseNoiseFilter {
public:
    // ���캯����minArea �Ǳ�������С��ͨ�������
    PhaseNoiseFilter(double minArea = 300.0);


    void filterInPlace(cv::Mat& phaseMap);

    // �˳�С��ͨ��������������������������
    cv::Mat filter(const cv::Mat& phaseMap);

private:
    double areaThreshold;  // �����ֵ
};



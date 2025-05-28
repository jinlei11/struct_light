#pragma once
#include <opencv2/opencv.hpp>

class PhaseNoiseFilter {
public:
    // 构造函数：minArea 是保留的最小连通区域面积
    PhaseNoiseFilter(double minArea = 300.0);


    void filterInPlace(cv::Mat& phaseMap);

    // 滤除小连通区域噪声，保留大区域（主区域）
    cv::Mat filter(const cv::Mat& phaseMap);


    //背景相位滤除，根据相位梯度滤除背景
    void filterPhaseByGradient(const cv::Mat& phaseMap,
                               float gradThreshold,
                               cv::Mat& mask,
                               cv::Mat& filteredPhaseMap);


private:
    double areaThreshold;  // 面积阈值
};



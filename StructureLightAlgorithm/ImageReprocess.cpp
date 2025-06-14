#include"ImageReprocess.h"



PhaseNoiseFilter::PhaseNoiseFilter(double minArea)
    : areaThreshold(minArea) {}


// 直接在原图上修改（小区域设为0）
void PhaseNoiseFilter::filterInPlace(cv::Mat& phaseMap) {
    CV_Assert(phaseMap.type() == CV_32FC1);

    // 只要不为0都认为是候选像素
    cv::Mat nonZeroMask = phaseMap != 0.0f;

    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(nonZeroMask, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // 创建一个临时掩膜用于保留大区域
    cv::Mat cleanMask = cv::Mat::zeros(phaseMap.size(), CV_8U);
    for (size_t i = 0; i < contours.size(); ++i) {
        double area = cv::contourArea(contours[i]);
        if (area >= areaThreshold) {
            cv::drawContours(cleanMask, contours, static_cast<int>(i), cv::Scalar(255), cv::FILLED);
        }
    }

    // 小区域在掩膜上是 0，直接置为 0
    phaseMap.setTo(0.0f, cleanMask == 0);
}



cv::Mat PhaseNoiseFilter::filter(const cv::Mat& phaseMap) {
    CV_Assert(phaseMap.type() == CV_32FC1);

    // 先制作非零掩膜（只要不是0，就认为是候选区域）
    cv::Mat nonZeroMask = phaseMap != 0.0f;

    // 找连通区域
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(nonZeroMask, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // 创建新掩膜
    cv::Mat cleanMask = cv::Mat::zeros(phaseMap.size(), CV_8U);
    for (size_t i = 0; i < contours.size(); ++i) {
        double area = cv::contourArea(contours[i]);
        if (area >= areaThreshold) {
            cv::drawContours(cleanMask, contours, static_cast<int>(i), cv::Scalar(255), cv::FILLED);
        }
    }

    // 应用掩膜
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

    // 计算 x 和 y 方向的梯度
    cv::Sobel(phaseMap, gradX, CV_32F, 1, 0, 3);
    cv::Sobel(phaseMap, gradY, CV_32F, 0, 1, 3);

    // 梯度幅值
    cv::magnitude(gradX, gradY, gradMag);

    // 阈值筛选
    cv::threshold(gradMag, mask, gradThreshold, 255.0, cv::THRESH_BINARY_INV);
    mask.convertTo(mask, CV_8UC1);

    // 可选：形态学操作
    cv::morphologyEx(mask, mask, cv::MORPH_OPEN, cv::Mat(), cv::Point(-1, -1), 1);
    cv::morphologyEx(mask, mask, cv::MORPH_DILATE, cv::Mat(), cv::Point(-1, -1), 1);

    // 应用掩膜
    filteredPhaseMap = cv::Mat::zeros(phaseMap.size(), phaseMap.type());
    phaseMap.copyTo(filteredPhaseMap, mask);
}
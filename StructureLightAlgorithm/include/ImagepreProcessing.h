#pragma once
#include <opencv2/opencv.hpp>
#include <iostream>

class ImageProcessor {
public:
    ImageProcessor(const std::string& path) : imagePath(path) {}

    bool loadImage();

    //提取刀具ROI区域
    std::vector<cv::Point> process();

    void showResult();

    void saveResult(const std::string& outputPath);

private:
    std::string imagePath;
    cv::Mat image, gray, denoised, adaptiveThresh, morph;
};

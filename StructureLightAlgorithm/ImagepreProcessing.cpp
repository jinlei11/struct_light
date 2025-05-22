#include"ImagepreProcessing.h"


bool ImageProcessor::loadImage() {
    image = cv::imread(imagePath);
    if (image.empty()) {
        std::cerr << "无法加载图像: " << imagePath << std::endl;
        return false;
    }
    return true;
}


std::vector<cv::Point> ImageProcessor::process() {

    // 转灰度
    cv::cvtColor(image, gray, cv::COLOR_BGR2GRAY);

    // 中值滤波
    cv::medianBlur(gray, denoised, 5);

    // 自适应阈值
    cv::adaptiveThreshold(denoised, adaptiveThresh, 255,
        cv::ADAPTIVE_THRESH_GAUSSIAN_C,
        cv::THRESH_BINARY_INV, 21, 10);

    // 闭操作
    cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(5, 5));
    cv::morphologyEx(adaptiveThresh, morph, cv::MORPH_CLOSE, kernel);

    // 查找轮廓
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(morph, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    int height = gray.rows;

    // 过滤轮廓
    std::vector<std::vector<cv::Point>> filteredContours;
    for (const auto& cnt : contours) {
        double area = cv::contourArea(cnt);
        if (area < 1000) continue;

        cv::Moments M = cv::moments(cnt);
        if (M.m00 == 0) continue;

        int cy = static_cast<int>(M.m01 / M.m00);
        if (cy < height * 0.85)
            filteredContours.push_back(cnt);
    }

    // 绘制最小外接矩形
    std::vector<cv::Point> box;
    if (!filteredContours.empty()) {
        std::vector<cv::Point> allPoints;
        for (const auto& cnt : filteredContours) {
            allPoints.insert(allPoints.end(), cnt.begin(), cnt.end());
        }

        cv::RotatedRect rect = cv::minAreaRect(allPoints);
        cv::Point2f boxPoints[4];
        rect.points(boxPoints);

        for (int i = 0; i < 4; ++i) {
            box.push_back(boxPoints[i]);
            std::cout << "X: " << boxPoints[i].x << "Y: " << boxPoints[i].y << std::endl;
        }

        cv::polylines(image, box, true, cv::Scalar(0, 165, 255), 4); // 橙色框
    }

    return box;
}



void ImageProcessor::showResult() {
    if (!image.empty()) {
        cv::namedWindow("ROI", cv::WINDOW_NORMAL);
        cv::imshow("ROI", image);
        cv::waitKey(0);
    }
}


void ImageProcessor::saveResult(const std::string& outputPath) {
    if (!image.empty()) {
        cv::imwrite(outputPath, image);
    }
}
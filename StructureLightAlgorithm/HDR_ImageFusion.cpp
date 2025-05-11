#include"include/HDR_ImageFusion.h"


//double OptimalinterclassVarianceMethod(cv::Mat& Image, int GrayValue) {
//
//	float TotalPixelNum = static_cast<float>(Image.rows * Image.cols);//总像素数
//
//    //统计灰度直方图
//    cv::Mat Hist;
//    int histSize = 256;    // 0-255级灰度
//    float range[] = { 0, 256 };
//    const float* histRange = { range };
//    cv::calcHist(&Image,         // 输入图像
//                 1,              // 图像数量
//                 0,              // 通道索引（灰度图为0）
//                 cv::Mat(),      // 掩码（空表示全图）
//                 Hist,           // 输出直方图
//                 1,              // 直方图维度
//                 &histSize,      // bin数量
//                 &histRange);    // 值范围
//
//    //低于阈值灰度级概率
//    double Pc0 = 0;
//    for (int value = 0; value <= GrayValue; value++) {
//        Pc0 += (Hist.at<float>(value, 0)/ TotalPixelNum);
//    }
//    //高于阈值灰度级概率
//    double Pc1 = 0;
//    for (int value = GrayValue + 1; value <= 255; value++) {
//        Pc1 += (Hist.at<float>(value, 0) / TotalPixelNum);
//    }
//
//    //低于阈值灰度级的平均灰度值
//    double mu0 = 0;
//    for (int value = 0; value <= GrayValue; value++) {
//        mu0 += ((Hist.at<float>(value, 0) / TotalPixelNum) * value / Pc0);
//    }
//
//    //高于阈值灰度级的平均灰度值
//    double mu1 = 0;
//    for (int value = GrayValue + 1; value <= 255; value++) {
//        mu1 += ((Hist.at<float>(value, 0) / TotalPixelNum) * value / Pc1);
//    }
//
//    //整幅图像的平均灰度
//    double mu = 0;
//    for (int value = 0; value <= 255; value++) {
//        mu += value * (Hist.at<float>(value, 0) / TotalPixelNum);
//    }
//
//
//    double sigma = Pc0 * (mu0 - mu) * (mu0 - mu) + Pc1 * (mu1 - mu) * (mu1 - mu);
//
//    return sigma;
//
//};



double OptimalinterclassVarianceMethod(cv::Mat& Image, int GrayValue) {
    CV_Assert(Image.type() == CV_8UC1); // 确保是灰度图

    const float TotalPixelNum = Image.rows * Image.cols;
    if (TotalPixelNum == 0) return 0;

    // 计算直方图（同原代码）
    cv::Mat Hist;
    const int histSize = 256;
    const float range[] = { 0, 256 };
    const float* histRange = { range };
    cv::calcHist(&Image, 1, 0, cv::Mat(), Hist, 1, &histSize, &histRange);

    // 合并计算所有参数
    double Pc0 = 0, Pc1 = 0, mu0_part = 0, mu1_part = 0, mu = 0;

    for (int value = 0; value < 256; ++value) {
        const double prob = Hist.at<float>(value) / TotalPixelNum;
        mu += value * prob; // 全局均值

        if (value <= GrayValue) {
            Pc0 += prob;
            mu0_part += value * prob;
        }
        else {
            Pc1 += prob;
            mu1_part += value * prob;
        }
    }

    // 处理除零异常
    const double mu0 = (Pc0 > 1e-7) ? (mu0_part / Pc0) : 0;
    const double mu1 = (Pc1 > 1e-7) ? (mu1_part / Pc1) : 0;

    return Pc0 * (mu0 - mu) * (mu0 - mu) + Pc1 * (mu1 - mu) * (mu1 - mu);
}
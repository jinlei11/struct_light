#include"include/HDR_ImageFusion.h"


//double OptimalinterclassVarianceMethod(cv::Mat& Image, int GrayValue) {
//
//	float TotalPixelNum = static_cast<float>(Image.rows * Image.cols);//��������
//
//    //ͳ�ƻҶ�ֱ��ͼ
//    cv::Mat Hist;
//    int histSize = 256;    // 0-255���Ҷ�
//    float range[] = { 0, 256 };
//    const float* histRange = { range };
//    cv::calcHist(&Image,         // ����ͼ��
//                 1,              // ͼ������
//                 0,              // ͨ���������Ҷ�ͼΪ0��
//                 cv::Mat(),      // ���루�ձ�ʾȫͼ��
//                 Hist,           // ���ֱ��ͼ
//                 1,              // ֱ��ͼά��
//                 &histSize,      // bin����
//                 &histRange);    // ֵ��Χ
//
//    //������ֵ�Ҷȼ�����
//    double Pc0 = 0;
//    for (int value = 0; value <= GrayValue; value++) {
//        Pc0 += (Hist.at<float>(value, 0)/ TotalPixelNum);
//    }
//    //������ֵ�Ҷȼ�����
//    double Pc1 = 0;
//    for (int value = GrayValue + 1; value <= 255; value++) {
//        Pc1 += (Hist.at<float>(value, 0) / TotalPixelNum);
//    }
//
//    //������ֵ�Ҷȼ���ƽ���Ҷ�ֵ
//    double mu0 = 0;
//    for (int value = 0; value <= GrayValue; value++) {
//        mu0 += ((Hist.at<float>(value, 0) / TotalPixelNum) * value / Pc0);
//    }
//
//    //������ֵ�Ҷȼ���ƽ���Ҷ�ֵ
//    double mu1 = 0;
//    for (int value = GrayValue + 1; value <= 255; value++) {
//        mu1 += ((Hist.at<float>(value, 0) / TotalPixelNum) * value / Pc1);
//    }
//
//    //����ͼ���ƽ���Ҷ�
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
    CV_Assert(Image.type() == CV_8UC1); // ȷ���ǻҶ�ͼ

    const float TotalPixelNum = Image.rows * Image.cols;
    if (TotalPixelNum == 0) return 0;

    // ����ֱ��ͼ��ͬԭ���룩
    cv::Mat Hist;
    const int histSize = 256;
    const float range[] = { 0, 256 };
    const float* histRange = { range };
    cv::calcHist(&Image, 1, 0, cv::Mat(), Hist, 1, &histSize, &histRange);

    // �ϲ��������в���
    double Pc0 = 0, Pc1 = 0, mu0_part = 0, mu1_part = 0, mu = 0;

    for (int value = 0; value < 256; ++value) {
        const double prob = Hist.at<float>(value) / TotalPixelNum;
        mu += value * prob; // ȫ�־�ֵ

        if (value <= GrayValue) {
            Pc0 += prob;
            mu0_part += value * prob;
        }
        else {
            Pc1 += prob;
            mu1_part += value * prob;
        }
    }

    // ��������쳣
    const double mu0 = (Pc0 > 1e-7) ? (mu0_part / Pc0) : 0;
    const double mu1 = (Pc1 > 1e-7) ? (mu1_part / Pc1) : 0;

    return Pc0 * (mu0 - mu) * (mu0 - mu) + Pc1 * (mu1 - mu) * (mu1 - mu);
}
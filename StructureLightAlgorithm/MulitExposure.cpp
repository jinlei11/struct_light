#include"MulitExposure.h"


/***************************************    基类     ***************************************/
HDR_Reconstruction::HDR_Reconstruction() {

};

HDR_Reconstruction::~HDR_Reconstruction() {

};


/****************************************************    自适应条纹投影测量法      *******************************************/

//  汪锦航,卢荣胜,刘端茂.高动态范围表面自适应条纹投影测量方法[J].光学学报,2021,41(19):145-154.
HDR_AdaptiveStripePattern::HDR_AdaptiveStripePattern() {

};


HDR_AdaptiveStripePattern::~HDR_AdaptiveStripePattern() {


};


//cv::Mat HDR_AdaptiveStripePattern::GenerationGrayHistogram() {
//
//	cv::Mat Hist;
//
//	return Hist;
//}










// 冯维,汤少靖,赵晓冬,等.基于自适应条纹的高反光表面三维面形测量方法[J].光学学报,2020,40(05):119-127.

HDR_AdaptiveStripePattern_FW::HDR_AdaptiveStripePattern_FW() {

};


HDR_AdaptiveStripePattern_FW::~HDR_AdaptiveStripePattern_FW() {


};



//cv::Mat HDR_AdaptiveStripePattern_FW::GenerationGrayHistogram() {
//
//	cv::Mat Hist;
//
//	return Hist;
//}







/****************************************************    自适应条纹投影测量法      *******************************************/




/****************************************************    基于相机响应曲线测量法      *******************************************/


HDR_CameraResponseCurve::HDR_CameraResponseCurve(std::string& filepath, 
                                                 std::string& resultpath,
                                                 std::vector<float>& Exposures,
                                                 int width,
                                                 int height)
    :M_ExposureImagesFilepath(filepath), M_Exposures(Exposures),M_Width(width),M_Height(height),
    M_CamResponseResultpath(resultpath)
{
    FileManipulation GetImages;
    std::vector<cv::Mat> DiffExposuresImgs;

    DiffExposuresImgs = GetImages.ReadPics(filepath);

    cv::Mat ResponseCurveImg = cv::Mat::zeros(height, width, CV_32FC1);
    getCameraResponseCurve(DiffExposuresImgs, Exposures, ResponseCurveImg);
};


HDR_CameraResponseCurve::~HDR_CameraResponseCurve() {



};


bool HDR_CameraResponseCurve::getCameraResponseCurve(std::vector<cv::Mat>& DiffExposureImages,
                            std::vector<float>& Exposures,
                            cv::Mat& ResponseCurveImg) {
    cv::Mat times_mat(Exposures.size(), 1, CV_32F);
    for (size_t i = 0; i < Exposures.size(); ++i) {
        times_mat.at<float>(static_cast<int>(i), 0) = Exposures[i];
    }

    // 计算相机响应曲线（使用Debevec算法）
    cv::Ptr<cv::CalibrateDebevec> calibrate = cv::createCalibrateDebevec();
    try {
        // 使用 CalibrateDebevec 进行响应曲线估计
        calibrate->process(DiffExposureImages, ResponseCurveImg, times_mat);
    }
    catch (const cv::Exception& e) {
        std::cerr << "处理响应曲线时发生错误: " << e.what() << std::endl;
        return false;
    }

    // 验证响应曲线结构（应为256x1）
    if (ResponseCurveImg.empty() || ResponseCurveImg.rows != 256 || ResponseCurveImg.cols != 1) {
        std::cerr << "响应曲线计算失败或格式错误！" << std::endl;
        return false;
    }

    //绘制相机响应曲线曲线
    int plot_width = 850, plot_height = 650;
    cv::Mat plot = cv::Mat(plot_height, plot_width, CV_8UC3, cv::Scalar(255, 255, 255)); // 白色背景
    std::vector<float> logE_values;

    double min_logE = ResponseCurveImg.at<float>(0, 0);
    double max_logE = ResponseCurveImg.at<float>(255, 0);

    // 绘制曲线（横轴：logE，纵轴：像素值）
    std::vector<cv::Point> points;
    for (int pixel_value = 0; pixel_value < 256; pixel_value++) {
        float logE = ResponseCurveImg.at<float>(pixel_value, 0);
        // 计算坐标
        float x = (logE - min_logE) / (max_logE - min_logE) * (plot_width - 100) + 50;
        float y = plot_height - static_cast<float>(pixel_value) / 255 * (plot_height - 100) - 50;
        points.push_back(cv::Point(x, y));
    }
    // 绘制折线（黑色曲线）
    polylines(plot, points, false, cv::Scalar(0, 0, 0), 2);

    //----------------------------  绘制坐标轴和标签  ------------------------------
    // X轴（logE）
    line(plot, cv::Point(50, plot_height - 50), cv::Point(plot_width - 50, plot_height - 50), cv::Scalar(0, 0, 0), 2);
    // 标签（logE范围，间隔为 (max_logE - min_logE)/5）
    for (int i = 0; i <= 5; i++) {
        float logE = min_logE + i * (max_logE - min_logE) / 5;
        float x = i * (plot_width - 100) / 5 + 50;
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << logE;
        putText(plot, ss.str(), cv::Point(x, plot_height - 30), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    }
    putText(plot, "Log Exposure (logE)", cv::Point(plot_width - 200, plot_height - 10), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 0, 0), 1);

    // Y轴（像素值0~255）
    line(plot, cv::Point(50, 50), cv::Point(50, plot_height - 50), cv::Scalar(0, 0, 0), 2);
    for (int i = 0; i <= 255; i += 50) {
        float y = plot_height - 50 - static_cast<float>(i) / 255 * (plot_height - 85);
        putText(plot, std::to_string(i), cv::Point(10, y), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1);
    }
    putText(plot, "Pixel Value (Z)", cv::Point(50, 30), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 0, 0), 1);


    // 显示并保存图像
    imshow("Gray Camera Response Curve", plot);
    imwrite(M_CamResponseResultpath + "result.png", plot);
    cv::waitKey(0);
    return true;

};




bool HDR_CameraResponseCurve::CalculateIrradianceImage(cv::Mat& Image, float Exposure, float gain) {
    cv::Mat image_float;
    Image.convertTo(image_float, CV_32FC1);

    float D = 25.0f;        // 暗电流
    float k = 0.1f;         // 校准系数 (W·s/m² per DN)

    cv::Mat irradiance(image_float.size(), CV_32FC1);
    for (int y = 0; y < image_float.rows; y++) {
        for (int x = 0; x < image_float.cols; x++) {
            // 获取当前像素值（浮点型）
            float V = image_float.at<float>(y, x);

            // 应用辐照度公式: E = (V - D) / (k * exposure * gain)
            float E = (V - D) / (k * Exposure * gain);

            // 存储结果
            irradiance.at<float>(y, x) = E;
        }
    }
    double minVal, maxVal;
    minMaxLoc(irradiance, &minVal, &maxVal);
    cv::Mat irradiance_normalized;
    irradiance.convertTo(irradiance_normalized, CV_8UC1, 255.0 / (maxVal - minVal), -minVal * 255.0 / (maxVal - minVal));

    // 保存结果
    imwrite("irradiance_output.png", irradiance_normalized);
    return true;
};



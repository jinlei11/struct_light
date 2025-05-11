// StructureLightAlgorithms.cpp: 定义应用程序的入口点。
//

#include "StructureLightAlgorithms.h"
#include"ProjectionPatternGeneration.h"

//using namespace std;
//
//int main()
//{
//
//    // 图像路径（根据你的目录结构调整）
//    vector<string> image_files = {
//        "../../../../CamResponseCurveImgs/1.bmp",
//        "../../../../CamResponseCurveImgs/2.bmp",
//        "../../../../CamResponseCurveImgs/3.bmp",
//        "../../../../CamResponseCurveImgs/4.bmp",
//        "../../../../CamResponseCurveImgs/5.bmp",
//        "../../../../CamResponseCurveImgs/6.bmp"
//    };
//
//    // 曝光时间（单位：秒）
//    vector<float> exposure_times = { 0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.06f };
//
//    // 确保图像数量和曝光时间数量一致
//    if (image_files.size() != exposure_times.size()) {
//        cerr << "图像数量与曝光时间数量不一致！" << endl;
//        return -1;
//    }
//
//    // 读取图像并保持为灰度图（保持单通道）
//    vector<cv::Mat> images;
//    for (const auto& path : image_files) {
//        cv::Mat gray = cv::imread(path, cv::IMREAD_GRAYSCALE);
//        if (gray.empty()) {
//            cerr << "无法读取图像: " << path << endl;
//            return -1;
//        }
//
//        // 将灰度图像转换为三通道图像
//        cv::Mat color_img;
//        cv::cvtColor(gray, color_img, cv::COLOR_GRAY2BGR);
//        images.push_back(color_img);
//    }
//
//    // 创建曝光时间矩阵（确保是单列矩阵，且类型为 CV_32F）
//    cv::Mat times_mat(exposure_times.size(), 1, CV_32F);
//    for (size_t i = 0; i < exposure_times.size(); ++i)
//        times_mat.at<float>(static_cast<int>(i), 0) = exposure_times[i];
//
//    // 打印曝光时间矩阵，确保它正确
//    cout << "曝光时间矩阵：" << endl;
//    for (int i = 0; i < times_mat.rows; ++i) {
//        cout << times_mat.at<float>(i, 0) << " ";
//    }
//    cout << endl;
//
//    // 创建 CRF 估计器
//    cv::Ptr<cv::CalibrateDebevec> calibrate = cv::createCalibrateDebevec();
//
//    // 估计响应曲线
//    cv::Mat response;
//    try {
//        // 使用 CalibrateDebevec 进行响应曲线估计
//        calibrate->process(images, response, times_mat);
//    }
//    catch (const cv::Exception& e) {
//        cerr << "处理响应曲线时发生错误: " << e.what() << endl;
//        return -1;
//    }
//
//    // 检查输出尺寸
//    if (response.empty()) {
//        cerr << "响应曲线计算失败！" << endl;
//        return -1;
//    }
//
//    // 打印响应曲线前几个值
//    cout << "响应曲线（前几个灰度值）:" << endl;
//    for (int i = 0; i < 256; i += 25)
//        cout << "灰度 " << i << ": " << response.at<float>(i, 0) << endl;
//
//    // 可视化响应曲线（灰度通道）
//    int width = 500, height = 300; // 增加宽度和高度，以便显示坐标系
//    cv::Mat crf_plot(height, width, CV_8UC3, cv::Scalar(255, 255, 255)); // 白底
//
//    double minVal, maxVal;
//    cv::minMaxLoc(response.col(0), &minVal, &maxVal);
//
//    // 绘制响应曲线（直接使用响应值绘制）
//    for (int i = 0; i < 256; ++i) {  // 遍历灰度值范围
//        float gray_value = static_cast<float>(i);
//        float v = response.at<float>(i, 0);  // 响应值（直接用响应值）
//
//        // 归一化响应值
//        float norm_v = (v - minVal) / (maxVal - minVal);
//
//        // Y坐标表示归一化后的响应值
//        int y = height - cvRound(norm_v * height);
//
//        // X坐标根据灰度值映射
//        int x = cvRound((gray_value / 255.0) * width);
//
//        // 在图像上画点
//        cv::circle(crf_plot, cv::Point(x, y), 2, cv::Scalar(0, 0, 0), -1);  // 画点
//    }
//
//    // 连接这些点以绘制曲线
//    for (int i = 1; i < 256; ++i) {  // 遍历灰度值
//        float gray_value1 = static_cast<float>(i - 1);
//        float gray_value2 = static_cast<float>(i);
//        float v1 = response.at<float>(i - 1, 0);
//        float v2 = response.at<float>(i, 0);
//
//        float norm_v1 = (v1 - minVal) / (maxVal - minVal);
//        float norm_v2 = (v2 - minVal) / (maxVal - minVal);
//
//        int y1 = height - cvRound(norm_v1 * height);
//        int y2 = height - cvRound(norm_v2 * height);
//
//        int x1 = cvRound((gray_value1 / 255.0) * width);
//        int x2 = cvRound((gray_value2 / 255.0) * width);
//
//        // 连接相邻的点绘制曲线
//        cv::line(crf_plot, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(0, 0, 0), 2);  // 黑色曲线
//    }
//
//    // 绘制坐标系（X轴和Y轴）
//    cv::line(crf_plot, cv::Point(0, height - 1), cv::Point(width - 1, height - 1), cv::Scalar(0, 0, 0), 2); // X轴
//    cv::line(crf_plot, cv::Point(0, 0), cv::Point(0, height - 1), cv::Scalar(0, 0, 0), 2); // Y轴
//
//    // 添加坐标轴标记（X轴为灰度值）
//    for (int i = 0; i <= 255; i += 50) {  // 每隔50个灰度值添加标记（0, 50, 100,...255）
//        int x = cvRound((i / 255.0) * width);
//        cv::putText(crf_plot, to_string(i), cv::Point(x, height - 5), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1, 8);
//    }
//
//    // 添加Y轴标记（响应值）
//    for (int i = 0; i <= 255; i += 50) { // 每隔50个响应值添加标记
//        int y = height - cvRound((i / 255.0) * height);
//        cv::putText(crf_plot, to_string(i), cv::Point(5, y), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0), 1, 8);
//    }
//
//    // 保存和显示
//    cv::imwrite("crf_curve_full.jpg", crf_plot); // 保存完整的响应曲线图
//    cv::imshow("Full Gray Camera Response Curve", crf_plot); // 显示完整的响应曲线图
//    cv::waitKey(0);
//
//    return 0;
//
//}
//
//
//
//    
//
//
//
//
//
//  
//
//
//
//
//

  















 


#pragma once
#include <iostream>
#include <opencv2/opencv.hpp>
#include <fstream>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <numeric>
#include "FileManipulation.h"
#include "Calibration.h"




using namespace std;
using namespace cv;






class OpticalPlaneCalibration {
public:

	OpticalPlaneCalibration(const std::string& BoardPath, const std::string& OpticalPlanePath, enum CalibrationPlate boardClass);


public:
	// 计算光平面标定图像的r t 
	void getRT(const Mat& cameraMatrix, const Mat& distCoeffs,
		vector<Mat>& rotationMat, vector<Mat>& translationMat,
		bool show);

	void Get_Camera_Data(const std::string& cameraParamsPath, Mat& instrinsic_Mat, Mat& distCoeff);


private:
	//Number of corner points of a checkerboard grid or circular calibration board (Width, Height)
	cv::Size M_board_size;

	//Checkerboard or circular calibration board with square dimensions or circle spacing
	cv::Size M_squareSize;

	//Distinguishing calibration plate types
	int M_CalibPlateType;

	//保存标定板图片的路径
	std::string M_PlaneBoardFilePath;
	//保存激光线图片的路径
	std::string M_PlaneLineFilePath;


};





// 获取激光条纹中心线
class Centerline {
public:
	// 预处理  将原图 src 进行处理，得到处理后的图像 imgHSVMask
	void picture_blur(Mat& src, Mat& imgHSVMask, vector<Rect>& rois) {
		Mat img;
		Mat imgGauss, imgRrode, imgDilate;

		threshold(src, img, 0, 255, THRESH_OTSU);

		GaussianBlur(img, imgGauss, Size(3, 3), 0);

		Mat element = getStructuringElement(MORPH_RECT, Size(3, 3));

		morphologyEx(imgGauss, imgHSVMask, MORPH_OPEN, element);
		//threshold(imgDilate, imgHSVMask, 50, 255, THRESH_TOZERO);

		// 自动框选roi
		int sSize = 2000;	// 面积阈值100	
		imgHSVMask = autoroi(imgHSVMask, sSize, rois);	// 好像相机标定的激光图不太适用

		// 保存ROI区域的灰度值
		//saveMultipleRowsToTxtInROIs(imgHSVMask, rois);
		// 
		// 手动框选roi
		//imgHSVMask = handroi(imgHSVMask);
	}

	// 提取激光条纹中心线  修改：应该将结果作为返回值传回，而不是用成员变量存储
	std::vector<cv::Point2d> ExtractCenters(Mat src, Mat& imgHSVMask, vector<Rect>& rois) {
		int w, brightValue;

		// 评估光条宽度w
		estimateWidth(imgHSVMask, rois, w, brightValue);

		vector<Point2d> centers = Centers(imgHSVMask, rois, w);
		//vector<Point2d> centers = extractLaserCentersGravity(imgHSVMask, rois, w);
		//cout << "Centers: " << centers << endl;

		Mat huidu_img0;
		src.copyTo(huidu_img0);
		drawPointsOnImage(centers, huidu_img0);

		// 定义半径k
		double k = 5.0;
		//以 半径k 对中心点遍历校正
		for (const auto& center_point : centers) {
			Point2d smoothed_point = SmoothedPoint(centers, center_point, k);
			smooth_points.push_back(smoothed_point);
		}

		//// 创建一个与 src 同尺寸的白色图像
		//Mat white_img = Mat::zeros(g_midImage.size(), CV_8UC3); // 创建一个全黑图像
		//white_img.setTo(Scalar(255, 255, 255)); // 将图像设置为全白
		//drawPointsOnImage(smooth_points, white_img, (0, 0, 0));

		Mat huidu_img1;
		src.copyTo(huidu_img1);
		//cvtColor(huidu_img1, huidu_img1, COLOR_GRAY2BGR);
		drawPointsOnImage(smooth_points, huidu_img1);
		return smooth_points;
	}

	vector<Point2d> Center_points; //激光中心线
	vector<Point2d> smooth_points; //K邻域校正后中心点

private:
	Mat src;//原图
	Mat distCoeff;//相机畸变系数
	Mat instrinsic_Mat;//相机内参
	Mat g_dstImage, g_midImage;

	// ------------------------在图像上绘制点集------------------------
	void drawPointsOnImage(vector<Point2d>& points, Mat& image, Scalar color = Scalar(0, 0, 255)) {

		for (const auto& point : points) {
			// 将小数坐标转换为图像中的像素坐标（保持小数精度）
			Point2d pixelPoint(point.x, point.y);

			// 缩放小数坐标到图像坐标，以绘制近似半径为0.1像素的实心点
			Point scaledPoint(cvRound(pixelPoint.x), cvRound(pixelPoint.y));

			// 绘制实心圆
			circle(image, scaledPoint, 0, color, 1);  // 半径为0的实心圆

		}
	}


	// ------------------------ROI------------------------
	// 自动框选roi 	 修改：传入引用即可
	Mat autoroi(Mat img, int size, vector<Rect>& rois) {
		Mat cimg;
		img.copyTo(cimg);
		threshold(cimg, cimg, 0, 255, THRESH_BINARY);
		//cvtColor(cimg, cimg, COLOR_BGR2GRAY);

		// 查找轮廓（链码，采用 CV_CHAIN_APPROX_NONE 表示不压缩边缘点）
		vector<vector<Point>> contours;
		vector<Vec4i> hierarchy;
		findContours(cimg, contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_NONE);

		//Mat contourImg = Mat::zeros(cimg.size(), CV_8UC3);
		//drawContours(contourImg, contours, -1, Scalar(0, 255, 0), 1);
		//imshow("Contours", contourImg);
		//waitKey(0);

		// 创建黑底图像用于绘制ROI区域
		Mat out_Image = Mat::zeros(img.size(), img.type());

		int padding = 2; // 控制四周扩展像素

		for (size_t i = 0; i < contours.size(); i++) {
			// 计算轮廓面积
			double area = contourArea(contours[i]);

			if (area > size) { // 面积阈值筛选
				Rect rect = boundingRect(contours[i]);

				// 扩展ROI四周
				Rect roi(rect.x - padding, rect.y - padding,
					rect.width + 2 * padding, rect.height + 2 * padding);

				// 限制在图像范围内
				roi &= Rect(0, 0, img.cols, img.rows);

				rois.push_back(roi);  // 保存 ROI 块

				// 拷贝原图像到输出图
				img(roi).copyTo(out_Image(roi));

				// 画扩展后的矩形框
				//rectangle(out_Image, roi, Scalar(255), 1);
			}
		}
		return out_Image;
	}

	// 手动框选roi
	struct ROIData {
		Mat temp_Image;
		Mat roi_Image;
		Mat roi_result_Image;
		Point pre_pt;
		Point cur_pt;
	};

	static void on_mouse(int event, int x, int y, int flags, void* userdata)
	{
		ROIData* data = static_cast<ROIData*>(userdata); // 通过userdata访问结构体
		Mat tmp, dst;

		if (event == EVENT_LBUTTONDOWN)
		{
			data->pre_pt = Point(x, y);  // 记录鼠标按下位置
		}
		else if (event == EVENT_MOUSEMOVE && flags) // 鼠标移动时，实时绘制矩形框
		{
			data->temp_Image.copyTo(tmp);
			data->cur_pt = Point(x, y);
			rectangle(tmp, data->pre_pt, data->cur_pt, Scalar(255, 255, 255), 1, 8, 0);
			imshow("img", tmp);  // 实时显示ROI矩形
		}
		else if (event == EVENT_LBUTTONUP)  // 松开鼠标左键，完成ROI框选
		{
			data->temp_Image.copyTo(data->roi_Image);
			rectangle(data->roi_Image, data->pre_pt, data->cur_pt, Scalar(255, 255, 255), 1, 8, 0);
			imshow("img", data->roi_Image);

			int width = abs(data->pre_pt.x - data->cur_pt.x);
			int height = abs(data->pre_pt.y - data->cur_pt.y);
			dst = data->temp_Image(Rect(min(data->cur_pt.x, data->pre_pt.x), min(data->cur_pt.y, data->pre_pt.y), width, height));

			Mat out_Image(data->temp_Image.size(), data->temp_Image.type(), Scalar(0)); // 创建全黑图
			Rect roi_rect = Rect(min(data->cur_pt.x, data->pre_pt.x), min(data->cur_pt.y, data->pre_pt.y), width, height);
			dst.copyTo(out_Image(roi_rect));
			out_Image.copyTo(data->roi_result_Image);
		}
	}

	Mat handroi(Mat img)
	{
		ROIData data;  // 创建结构体以存储回调数据
		img.copyTo(data.temp_Image);  // 将输入图像复制到结构体
		namedWindow("img", WINDOW_NORMAL);
		imshow("img", img);

		setMouseCallback("img", on_mouse, &data);  // 传递结构体指针
		waitKey(0);
		destroyWindow("img");

		return data.roi_result_Image;  // 返回结果图像
	}


	// ----------------------查看像素分布-------------------------
	// 将图像中某一行的灰度值保存为 txt 文件
	void saveRowToTxt(const Mat& gray, int y) {
		CV_Assert(gray.channels() == 1 && y >= 0 && y < gray.rows);

		// 构建文件名
		string filename = "第" + to_string(y) + "行.txt";

		ofstream fout(filename);
		if (!fout.is_open()) {
			cerr << "无法打开文件 " << filename << " 进行写入！" << endl;
			return;
		}

		const uchar* row = gray.ptr<uchar>(y);
		for (int x = 0; x < gray.cols; ++x) {
			fout << static_cast<int>(row[x]) << endl;  // 一行一个灰度值
		}

		fout.close();
	}

	// 将roi区域每行像素进行保存
	void saveMultipleRowsToTxtInROIs(const  Mat& gray, const vector< Rect>& rois, const string& filename = "多行灰度值.txt") {
		CV_Assert(gray.channels() == 1);

		ofstream fout(filename);
		if (!fout.is_open()) {
			cerr << "无法打开文件 " << filename << " 进行写入！" << endl;
			return;
		}

		for (const auto& roi : rois) {
			int startY = max(0, roi.y);
			int endY = min(gray.rows - 1, roi.y + roi.height - 1);
			int startX = max(0, roi.x);
			int endX = min(gray.cols - 1, roi.x + roi.width - 1);

			for (int y = startY; y <= endY; ++y) {
				const uchar* row = gray.ptr<uchar>(y);
				for (int x = startX; x <= endX; ++x) {
					fout << static_cast<int>(row[x]);
					if (x != endX)
						fout << " ";
				}
				fout << endl; // 换行
			}
		}

		fout.close();
		cout << "已保存所选 ROI 中的多行灰度值到文件：" << filename << endl;
	}

	// 区域 大于0
	void saveMultipleRowsToTxtInROI(const  Mat& gray, const vector< Rect>& rois, const string& filename = "多行灰度值.txt") {
		CV_Assert(gray.channels() == 1);

		ofstream fout(filename);
		if (!fout.is_open()) {
			cerr << "无法打开文件 " << filename << " 进行写入！" << endl;
			return;
		}

		for (const auto& roi : rois) {
			int startY = max(0, roi.y);
			int endY = min(gray.rows - 1, roi.y + roi.height - 1);
			int startX = max(0, roi.x);
			int endX = min(gray.cols - 1, roi.x + roi.width - 1);

			for (int y = startY; y <= endY; ++y) {
				const uchar* row = gray.ptr<uchar>(y);
				bool first = true;

				for (int x = startX; x <= endX; ++x) {
					int val = static_cast<int>(row[x]);
					if (val > 0) {
						if (!first) fout << " ";
						fout << val;
						first = false;
					}
				}

				fout << endl; // 每行一个 ROI 子区域
			}
		}

		fout.close();
		cout << "已保存所有 ROI 区域中非零灰度值到文件：" << filename << endl;
	}


	// ----------------------自适应提取中心线------------------------
	// 光条宽度w
	void estimateWidth(const Mat& gray, const vector<Rect>& rois, int& w, int& brightValue) {
		vector<uchar> nonZeroValues;
		//存入ROI区域内所有的非0灰度值
		for (const auto& roi : rois) {
			for (int y = roi.y; y < roi.y + roi.height; ++y) {
				const uchar* row = gray.ptr<uchar>(y);
				for (int x = roi.x; x < roi.x + roi.width; ++x) {
					uchar val = row[x];
					if (val > 0) nonZeroValues.push_back(val);
				}
			}
		}

		if (nonZeroValues.empty()) {
			cout << "无非零像素，无法估计亮度阈值与宽度。" << endl;
			w = 0;
			brightValue = 255;
			return;
		}

		sort(nonZeroValues.begin(), nonZeroValues.end(), greater<uchar>());
		int index = max(1, static_cast<int>(nonZeroValues.size() * 0.2)) - 1;
		brightValue = nonZeroValues[index];
		cout << "亮度阈值（前20%的最小亮度）= " << brightValue << endl;

		// 使用此 brightValue 计算宽度
		int totalPixels = 0;
		int rowCount = 0;

		//对所有 ROI 中的每一行，统计该行中大于等于 brightValue 的像素个数，累加求平均，最终得到条纹宽度 w
		for (const auto& roi : rois) {
			for (int y = roi.y; y < roi.y + roi.height; ++y) {
				const uchar* row = gray.ptr<uchar>(y);
				int count = 0;
				for (int x = roi.x; x < roi.x + roi.width; ++x) {
					if (row[x] >= brightValue)
						++count;
				}
				if (count > 0) {
					++rowCount;
					totalPixels += count;
				}
			}
		}

		w = (rowCount > 0) ? (totalPixels / rowCount) : 0;
		cout << "估计的激光条纹宽度 w = " << w << " 像素" << endl;
	}

	// 初始中心点选取
	int initialCenter(const Mat& gray, int y, const Rect& roi) {
		const uchar* row = gray.ptr<uchar>(y);

		int maxVal = 0;
		vector<int> maxPositions;

		for (int x = roi.x; x < roi.x + roi.width; ++x) {
			int val = row[x];
			if (val > maxVal) {
				maxVal = val;
				maxPositions.clear();
				maxPositions.push_back(x);
			}
			else if (val == maxVal) {
				maxPositions.push_back(x);
			}
		}

		if (maxPositions.empty()) return -1;

		// 取最大值位置平均
		int sumX = accumulate(maxPositions.begin(), maxPositions.end(), 0);
		return static_cast<int>(sumX / maxPositions.size());
	}

	// 高斯拟合：y = A * exp(-(x - mu)^2 / (2 * sigma^2))，提取中心位置 mu
	double gaussFit(const vector<uchar>& profile) {
		vector<Point2f> points;

		for (int i = 0; i < profile.size(); ++i) {
			if (profile[i] > 30) // 可调阈值
				points.emplace_back(i, profile[i]);
		}

		if (points.size() < 5) return -1;

		double A = *max_element(profile.begin(), profile.end());
		double mu = distance(profile.begin(), max_element(profile.begin(), profile.end()));
		double sigma = 1.5; // 固定或预设值更稳定

		// 简化 Gauss-Newton：手工迭代求解 mu
		for (int iter = 0; iter < 10; ++iter) {
			double num_mu = 0, denom = 0;
			for (auto& pt : points) {
				double x = pt.x;
				double y = pt.y;
				double exp_term = exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
				double weight = A * exp_term;
				num_mu += weight * x * (y - weight);
				denom += weight * (y - weight);
			}
			if (abs(denom) < 1e-6) break;
			mu += num_mu / denom;
		}
		return mu;
	}

	// ===== 灰度重心 =====
	Point2f computeGravityCenter(const vector<uchar>& profile, int offsetX) {
		double sumI = 0, sumIX = 0;
		for (int i = 0; i < profile.size(); ++i) {
			sumI += profile[i];
			sumIX += i * profile[i];
		}
		if (sumI <= 0) return Point2f(-1, -1);
		return Point2f(sumIX / sumI + offsetX, 0);
	}

	// ===== 高斯拟合中心 =====
	double gaussFullFitCenter(const vector<uchar>& profile) {
		if (profile.size() < 5) return -1;

		double A = *max_element(profile.begin(), profile.end());
		if (A < 30) return -1; // 灰度太低忽略

		// 仅用高亮部分
		vector<Point2f> points;
		for (int i = 0; i < profile.size(); ++i) {
			if (profile[i] > 0.6 * A)
				points.emplace_back(i, profile[i]);
		}

		if (points.size() < 5) return -1;

		double mu = distance(profile.begin(), max_element(profile.begin(), profile.end()));
		double sigma = 1.5;

		for (int iter = 0; iter < 10; ++iter) {
			double num_mu = 0, denom = 0;
			for (auto& pt : points) {
				double x = pt.x, y = pt.y;
				double exp_term = exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
				double weight = A * exp_term;
				num_mu += weight * x * (y - weight);
				denom += weight * (y - weight);
			}
			if (abs(denom) < 1e-6) break;
			mu += num_mu / denom;
		}

		return mu;
	}

	// ===== 中心线提取主函数 =====
	vector<Point2d> extractLaserCentersGravity(const Mat& gray, const vector<Rect>& rois, int w) {
		vector<Point2d> centers;

		for (const Rect& roi : rois) {
			for (int row = roi.y; row < roi.y + roi.height; ++row) {
				// 1. 找初始中心点
				int x0 = initialCenter(gray, row, roi); // 用户已有实现
				if (x0 < 0) continue;

				// 2. 设定左右 2w 范围
				int left = max(roi.x, x0 - 2 * w);
				int right = min(roi.x + roi.width - 1, x0 + 2 * w);
				if (right <= left + 5) continue;

				// 3. 提取 profile
				vector<uchar> profile;
				for (int x = left; x <= right; ++x)
					profile.push_back(gray.at<uchar>(row, x));

				// 4. 平滑
				Mat profMat(profile, true);        // 转成行向量 Mat
				profMat = profMat.reshape(1, 1);   // 显式 reshape 成一行
				GaussianBlur(profMat, profMat, Size(5, 1), 1.0);  // 高斯滤波
				profile.assign(profMat.begin<uchar>(), profMat.end<uchar>());  // 赋值回 vector


				// 5. 灰度重心估计初值
				Point2f gravity = computeGravityCenter(profile, left);
				int center_x = static_cast<int>(gravity.x - left);

				// 6. 局部截取用于高斯拟合
				int sub_left = max(0, center_x - w);
				int sub_right = min((int)profile.size() - 1, center_x + w);
				// 此处的sub_left  sub_right  数据存在问题，越界了。需要检查修改。

				vector<uchar> sub_profile(profile.begin() + sub_left, profile.begin() + sub_right + 1);

				// 7. 拟合与融合
				double mu = gaussFullFitCenter(sub_profile);
				if (mu >= 0)
					centers.emplace_back(mu + left + sub_left, row);
				else
					centers.emplace_back(gravity.x, row);
			}
		}

		return centers;
	}

	// 计算灰度重心
	Point2f gravityCenter(const vector<uchar>& profile, int x_start) {
		double weightedSum = 0;
		double sum = 0;
		for (int i = 0; i < profile.size(); ++i) {
			weightedSum += (x_start + i) * profile[i];
			sum += profile[i];
		}
		if (sum <= 0) return Point2f(-1, -1);
		return Point2f(weightedSum / sum, 0);
	}

	// 亚像素化y坐标 (未使用)
	double refineYSubpixel(const  Mat& gray, int x, int y) {
		if (y <= 0 || y >= gray.rows - 1 || x <= 0 || x >= gray.cols - 1)
			return static_cast<double>(y);  // 边界不做处理

		// 获取上下三行中该x处灰度值
		double I1 = gray.at<uchar>(y - 1, x);
		double I2 = gray.at<uchar>(y, x);
		double I3 = gray.at<uchar>(y + 1, x);

		double denom = (I1 - 2 * I2 + I3);
		if (abs(denom) < 1e-3) return static_cast<double>(y);

		// 抛物线极值点 y_sub = y - 0.5 * (I3 - I1) / (I1 - 2*I2 + I3)
		double delta = 0.5 * (I1 - I3) / denom;
		return y + delta;
	}

	// 提取中心线
	vector<Point2d> Centers(const Mat& gray, const vector<Rect>& rois, int w) {
		vector<Point2d> centers;

		for (const Rect& roi : rois) {
			for (int row = roi.y; row < roi.y + roi.height; ++row) {

				// 1. 找初始中心点
				int x0 = initialCenter(gray, row, roi);
				if (x0 < 0) continue;

				// 2. 确定中心点左右 2w 范围
				int left = max(roi.x, x0 - 2 * w);
				int right = min(roi.x + roi.width - 1, x0 + 2 * w);

				if (right <= left + 5) continue; // 保证区域宽度足够

				// 3. 提取刚刚确定的 4w 区域的灰度值
				vector<uchar> profile;
				for (int x = left; x <= right; ++x) {
					profile.push_back(gray.at<uchar>(row, x));
				}

				// 4. 在区域内高斯滤波 并 平滑处理降低噪声影响
				Mat profMat(profile);
				profMat = profMat.reshape(1, 1);
				GaussianBlur(profMat, profMat, Size(5, 1), 1.0);
				profile.assign(profMat.begin<uchar>(), profMat.end<uchar>());

				// 5. 灰度重心
				Point2f pt = gravityCenter(profile, left);
				if (pt.x >= 0) {
					/*int intX = static_cast<int>(round(pt.x));
					double refinedY = refineYSubpixel(gray, intX, row);
					centers.emplace_back(pt.x, refinedY);*/
					centers.emplace_back(pt.x, row);
				}
			}
		}
		return centers;
	}


	//------------------------获取中心线------------------------
	void huidu(Mat PointImage, Mat src, vector<Point2d>& centers)
	{
		// 读取PointImage中的中心点
		vector<Point2d> points;
		point(PointImage, points);
		vector<double> kcal;

		for (int i = 0; i < points.size() - 1; i++) {

			double pfCosSita = 0, pfSinSita = 0;
			CalcNormVec(Point2d(points[i - 1].x, points[i - 1].y), Point2d(points[i].x, points[i].y), Point2d(points[i + 1].x, points[i + 1].y), pfCosSita, pfSinSita);
			// 三点求法向量函数CalcNormVec
			//-------------------灰度中心法gdd---------------------//
			double sum = 0, sum_sumx = 0, sum_sumy = 0;
			for (int j = 0; j < 2; j++)		// 求该点以及法向量上两个点，一共三个点的灰度重心，得到亚像素级别的中心点
			{
				if (j == 0)
				{
					double cj = points[i].x;
					double ci = points[i].y;
					sum = ijpixel(cj, ci, src);//灰度值总和    此时ci、cj是整数，ijpixel得到的就是（x，y）处的灰度值
					sum_sumx = ijpixel(cj, ci, src) * cj;//（x*灰度值）的总和
					sum_sumy = ijpixel(cj, ci, src) * ci;//（y*灰度值）的总和
				}
				else
				{
					//cout << "pfCosSita = " << pfCosSita << endl;
					double x_cor = points[i].x + j * pfCosSita;	//cor点(在j=0(第一个)点的法向量后)的x坐标：x + cos
					double y_cor = points[i].y + j * pfSinSita;
					double x_cor1 = points[i].x - j * pfCosSita;	//cor_1点的x坐标
					double y_cor1 = points[i].y - j * pfSinSita;//此时x_cor等是小数，ijpixel得到的就是（x_cor，y_cor）等处的灰度值
					sum = sum + ijpixel(x_cor, y_cor, src) + ijpixel(x_cor1, y_cor1, src);
					sum_sumx = sum_sumx + ijpixel(x_cor, y_cor, src) * x_cor + ijpixel(x_cor1, y_cor1, src) * x_cor1;
					sum_sumy = sum_sumy + ijpixel(x_cor, y_cor, src) * y_cor + ijpixel(x_cor1, y_cor1, src) * y_cor1;
				}
			}
			Point2d p;
			p.x = sum_sumx / sum;
			p.y = sum_sumy / sum;
			//cout << "i = " << i << endl;
			//cout << "sumx = " << sum_sumx << "sumy = " << sum_sumy << "sum = " << sum << endl;
			//cout << "p.x = " << p.x << "p.y = " << p.y << endl;
			centers.push_back(p);

		}
	}

	//第i，j个点像素值
	double ijpixel(double& x, double& y, Mat& m)
	{
		int x_0 = int(x);
		int x_1 = int(x + 1);
		int y_0 = int(y);
		int y_1 = int(y + 1);
		int px_0y_0 = int(m.at<uchar>(y_0, x_0));
		int px_0y_1 = int(m.at<uchar>(y_1, x_0));
		int px_1y_0 = int(m.at<uchar>(y_0, x_1));
		int px_1y_1 = int(m.at<uchar>(y_1, x_1));
		double x_y0 = px_0y_0 + (x - double(x_0)) * (px_1y_0 - px_0y_0);
		double x_y1 = px_0y_1 + (x - double(x_0)) * (px_1y_1 - px_0y_1);
		double x_y = x_y0 + (y - double(y_0)) * (x_y1 - x_y0);
		return x_y;
	}
	//计算法向量
	void CalcNormVec(Point2d ptA, Point2d ptB, Point2d ptC, double& pfCosSita, double& pfSinSita)
	{
		double fVec1_x, fVec1_y, fVec2_x, fVec2_y;
		if (ptA.x == 0 && ptA.y == 0)
		{
			ptA.x = ptC.x;
			ptA.y = ptC.y;
			//先用B点坐标减A点坐标。
			fVec1_x = -(ptB.x - ptA.x);
			fVec1_y = -(ptB.y - ptA.y);
		}
		else
		{
			//先用B点坐标减A点坐标。
			fVec1_x = ptB.x - ptA.x;
			fVec1_y = ptB.y - ptA.y;
		}

		if (ptC.x == 0 && ptC.y == 0)
		{
			ptC.x = ptA.x;
			ptC.y = ptA.y;
			//再用C点坐标减B点坐标。
			fVec2_x = (ptB.x - ptC.x);
			fVec2_y = (ptB.y - ptC.y);
		}
		else
		{
			//再用C点坐标减B点坐标。
			fVec2_x = ptC.x - ptB.x;
			fVec2_y = ptC.y - ptB.y;
		}

		//单位化。
		double fMod = sqrt(fVec1_x * fVec1_x + fVec1_y * fVec1_y);
		fVec1_x /= fMod;
		fVec1_y /= fMod;
		//计算垂线。
		double fPerpendicularVec1_x = -fVec1_y;
		double fPerpendicularVec1_y = fVec1_x;


		//单位化。
		fMod = sqrt(fVec2_x * fVec2_x + fVec2_y * fVec2_y);
		fVec2_x /= fMod;
		fVec2_y /= fMod;
		//计算垂线。
		double fPerpendicularVec2_x = -fVec2_y;
		double fPerpendicularVec2_y = fVec2_x;

		//求和。
		double fSumX = fPerpendicularVec1_x + fPerpendicularVec2_x;
		double fSumY = fPerpendicularVec1_y + fPerpendicularVec2_y;
		//单位化。
		fMod = sqrt(fSumX * fSumX + fSumY * fSumY);
		double fCosSita = fSumX / fMod;
		double fSinSita = fSumY / fMod;
		pfCosSita = fCosSita;
		pfSinSita = fSinSita;
	}










	//------------------------校正中心点------------------------
	Point2d SmoothedPoint(const vector<Point2d>& center_points, Point2d a, double k) {
		double numerator_x = 0.0;
		double numerator_y = 0.0;
		double denominator = 0.0;

		for (const auto& p : center_points) {
			double distance = norm(a - p);
			if (distance < k && distance != 0) { // 只考虑半径k以内的点且不包括当前中心点
				double delta = 1.0 / distance; // 计算权重Δi
				numerator_x += delta * p.x;    // 使用p.x计算加权总和
				numerator_y += delta * p.y;    // 使用p.y计算加权总和
				denominator += delta;          // 计算权重总和
			}
		}

		if (denominator != 0) {
			double corrected_x = numerator_x / denominator; // 计算加权平均x坐标
			double corrected_y = numerator_y / denominator; // 计算加权平均y坐标
			return Point2d(corrected_x, corrected_y);
		}
		else {
			return a; // 如果半径内没有点，则返回原始点
		}
	}



	//中心线坐标
	void point(Mat& inputimg, vector<Point2d>& pt)
	{
		for (int i = 0; i < inputimg.cols; i++) {
			for (int j = 0; j < inputimg.rows; j++) {
				if (inputimg.at<uchar>(j, i) >= 95) {
					Point2d curr = Point2d(i, j);
					pt.push_back(curr);
				}
			}
		}
	}

};


//获取激光线与标点板交点，并转为相机坐标系
// 第一步：Set_param 传入标定板内参、平移向量、标定板旋转向量、标定板图片、2D中心线数据、宽---与激光线垂直、高---与激光线平行，完成初始化。
// 第二步：Get_plane 获取标定板平面方程。
// 第三步：Get_Intersection 获取激光线与标定板直线的交点。
// 第四步：LaserPoints_cam 将交点数据转为相机坐标系。
// 第五步：访问 laser_cam 获取相机坐标系下点。
class Turn_to_3D {
public:
	void Set_param(Mat instrinsic_Mat, Mat distCoeff, Mat translationMat, Mat rotationMat, Mat image,
		vector<Point2d> laser_pix, enum CalibrationPlate boardClass){
		this->instrinsic_Mat = instrinsic_Mat;
		this->distCoeff = distCoeff;
		this->translationMat = translationMat;
		this->rotationMat = rotationMat;
		this->image = image;
		this->laser_pix = laser_pix;

		switch (boardClass)
		{
		case Plate_GP050:
			this->chessboard_width = 11;
			this->chessboard_height = 8;
			break;
		case Plate_doubleCircle:
			this->chessboard_width = 7;
			this->chessboard_height = 5;
			break;
		case Plate_bigCircle:
			this->chessboard_width = 13;
			this->chessboard_height = 10;
			break;
		case Plate_Big_chess:
			this->chessboard_width = 9;
			this->chessboard_height = 6;
			break;
		case Plate_smallCircle:
			this->chessboard_width = 12;
			this->chessboard_height = 9;
			break;
		case Plate_chess811:
			this->chessboard_width = 8;
			this->chessboard_height = 11;
			break;
		default:
			break;
		}
	}
	/**
	* @brief 标定板外参转为标定板平面方程 Ax+By+Cz+D=0
	*/
	void Get_plane() {
		//将相机外参数向量整合成一个4×4矩阵,external_Mat
		Mat raMat, tempMat, externalMat;
		Mat bottomline = (Mat_<double>(1, 4) << 0, 0, 0, 1);
		//罗德里格斯公式
		Rodrigues(rotationMat, raMat);
		hconcat(raMat, translationMat, tempMat);//水平拼接
		vconcat(tempMat, bottomline, externalMat);//垂直拼接

		//求外参数矩阵的逆，并求出标定板所在平面的平面方程，坐标系为相机坐标系
		Mat externalMat_invert = externalMat.inv();	//求逆
		Mat ex_param = (Mat_<double>(1, 4) << 0, 0, 1, 0);
		Mat planeEquation = ex_param * externalMat_invert;	//平面方程=外参*

		plane = planeEquation;
	}
	/**
	* @brief 获取激光线与标定板直线的交点
	*/
	//在图像上绘制点集
	void drawPointsOnImage(vector<Point2d>& points, Mat& image, Scalar color = Scalar(0, 0, 255)) {  // 默认红色
		for (const auto& point : points) {
			// 将小数坐标转换为图像中的像素坐标（保持小数精度）
			Point2d pixelPoint(point.x, point.y);

			// 缩放小数坐标到图像坐标，以绘制近似半径为0.1像素的实心点
			Point scaledPoint(cvRound(pixelPoint.x), cvRound(pixelPoint.y));

			// 绘制实心圆
			circle(image, scaledPoint, 0, color, -1);  // 半径为0的实心圆
		}
	}


	void Get_Intersection(){
		Size boardSize;
		boardSize.width = chessboard_width;//宽---与激光线垂直
		boardSize.height = chessboard_height;//高---与激光线平行

		vector<double> laserline;	//激光线
		// 中心点去畸变
		undistortPoints(laser_pix, undistorted_centers, instrinsic_Mat, distCoeff);
		// 归一化点转回像素坐标系
		for (auto& pt : undistorted_centers) {
			pt.x = pt.x * instrinsic_Mat.at<double>(0, 0) + instrinsic_Mat.at<double>(0, 2);
			pt.y = pt.y * instrinsic_Mat.at<double>(1, 1) + instrinsic_Mat.at<double>(1, 2);
		}

		laserline = line_fitting(undistorted_centers);

		vector<Point2d> allp;

		vector<Point2d> pointbuf;
		bool found;

		found = findChessboardCornersSB(image, boardSize, pointbuf, CALIB_CB_EXHAUSTIVE | CALIB_CB_ACCURACY);
		//cout << "found = " << found << endl;
		//drawPointsOnImage(pointbuf, img2);

		if (found) {
			vector<Point2d> ipoint;
			vector<vector<double>> myline;
			for (size_t i = 0; i < boardSize.height; i++) // 遍历棋盘格的每一行
			{
				vector<Point2d> tps;
				for (size_t j = 0; j < boardSize.width; j++) // 遍历当前行的每个角点
				{
					Point2d tp;
					tp = pointbuf[j + i * boardSize.width]; // 获取当前角点坐标
					tps.push_back(tp); // 将角点加入当前行的点集合
				}

				vector<double> tl;
				tl = line_fitting(tps); // 对当前行的角点集合进行直线拟合，得到该行角点的直线方程
				myline.push_back(tl); // 将拟合得到的直线参数存入 myline

				Point2d ips = intersection(laserline, tl); // 计算当前行直线与激光线的交点
				allp.push_back(ips); // 将交点加入到交点集合 allp 中
			}

			laser_chess_p = allp;

			Mat huidu_img2;
			image.copyTo(huidu_img2);
			drawPointsOnImage(laser_chess_p, huidu_img2);

			ofstream ofs("交点.txt", ios::trunc);
			for (const auto& pt : laser_chess_p) {
				ofs << pt.x << "\t" << pt.y << endl;
			}
			ofs.close();
		}
		else {
			cout << "Chessboard corners not found." << endl;
		}
	}


	//** @brief 像素坐标系转相机坐标系

	vector<Point3d> LaserPoints_cam(){
		//平面方程的方向向量
		double a1 = plane.at<double>(0, 0);
		double b1 = plane.at<double>(0, 1);
		double c1 = plane.at<double>(0, 2);
		double d1 = plane.at<double>(0, 3);
		cout << "平面方程的方向向量：" << a1 << "\t" << b1 << "\t" << c1 << "\t" << d1 << endl;

		Mat instrinsic_Mat_invert = instrinsic_Mat.inv();	//内参矩阵求逆

		vector<Point3d> camPoints;
		for (int i = 0; i < laser_chess_p.size(); i++){
			//激光条纹上的点在像素坐标系下为pixelPoint，（u,v,1）
			Mat pixelPoint = (Mat_<double>(3, 1) << laser_chess_p[i].x, laser_chess_p[i].y, 1);
			//(u,v,1)乘相机的内参数矩阵的逆，就是激光条纹特征点在归一化平面的坐标，位于相机坐标系
			Mat imagePoint = instrinsic_Mat_invert * pixelPoint;

			double x1 = imagePoint.at<double>(0, 0);
			double y1 = imagePoint.at<double>(1, 0);
			double z1 = 1.0;  // 归一化平面上 z 坐标设为 1

			 //计算射线和平面交点
			double denominator = a1 * x1 + b1 * y1 + c1 * z1;
			if (denominator == 0) {
				cout << "warning：射线与平面平行，无法计算交点" << endl;
				continue;
			}

			// 计算交点的比例因子 lambda
			double lambda = -d1 / denominator;

			// 得到相机坐标系中的三维坐标
			double xc = lambda * x1;
			double yc = lambda * y1;
			double zc = lambda;

			Point3d p(xc, yc, zc);
			camPoints.push_back(p);

		}
		return camPoints;
	}

private:
	Mat instrinsic_Mat;
	Mat distCoeff;
	Mat translationMat;
	Mat rotationMat;
	Mat image;//标定板图片
	Mat plane;
	vector<Point2d> laser_pix;//像素坐标系下的激光中心线
	vector<Point2d> laser_chess_p;
	int chessboard_width, chessboard_height;
	vector<Point2d> undistorted_centers;	// 中心线去畸变后

	vector<double> line_fitting(vector<Point2d> input){
		Mat dst = Mat(2, 2, CV_64FC1, Scalar(0));
		//最小二乘法
		int n = input.size();
		for (int i = 0; i < n; i++)
		{
			dst.at<double>(0, 0) += pow(input[i].x, 2);
			dst.at<double>(0, 1) += input[i].x * input[i].y;
			dst.at<double>(1, 0) += input[i].x;
			dst.at<double>(1, 1) += input[i].y;
		}
		double k = (dst.at<double>(0, 1) - dst.at<double>(1, 0) * dst.at<double>(1, 1) / n) / (dst.at<double>(0, 0) - dst.at<double>(1, 0) * dst.at<double>(1, 0) / n);

		double x = dst.at<double>(1, 0) / n;
		double y = dst.at<double>(1, 1) / n;
		double b = y - k * x;

		vector<double> line(2);
		
		line[0] = k;
		line[1] = b;

		return line;
	}
	Point2d intersection(vector<double> laserline, vector<double> chessline)
	{
		double x = (chessline[1] - laserline[1]) / (laserline[0] - chessline[0]);
		double y = laserline[0] * x + laserline[1];

		return Point2d(x, y);
	}
};


void OptPlaneCalibration(const std::string& PlaneBoardFilename, 
						 const std::string& PlaneLineFilename,
						 const std::string& camParamPath,
						 const std::string& PointsPath,
						 const std::string& PlanePath,
	                     enum CalibrationPlate BoardClass);
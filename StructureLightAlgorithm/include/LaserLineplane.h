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
	// �����ƽ��궨ͼ���r t 
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

	//����궨��ͼƬ��·��
	std::string M_PlaneBoardFilePath;
	//���漤����ͼƬ��·��
	std::string M_PlaneLineFilePath;


};





// ��ȡ��������������
class Centerline {
public:
	// Ԥ����  ��ԭͼ src ���д����õ�������ͼ�� imgHSVMask
	void picture_blur(Mat& src, Mat& imgHSVMask, vector<Rect>& rois) {
		Mat img;
		Mat imgGauss, imgRrode, imgDilate;

		threshold(src, img, 0, 255, THRESH_OTSU);

		GaussianBlur(img, imgGauss, Size(3, 3), 0);

		Mat element = getStructuringElement(MORPH_RECT, Size(3, 3));

		morphologyEx(imgGauss, imgHSVMask, MORPH_OPEN, element);
		//threshold(imgDilate, imgHSVMask, 50, 255, THRESH_TOZERO);

		// �Զ���ѡroi
		int sSize = 2000;	// �����ֵ100	
		imgHSVMask = autoroi(imgHSVMask, sSize, rois);	// ��������궨�ļ���ͼ��̫����

		// ����ROI����ĻҶ�ֵ
		//saveMultipleRowsToTxtInROIs(imgHSVMask, rois);
		// 
		// �ֶ���ѡroi
		//imgHSVMask = handroi(imgHSVMask);
	}

	// ��ȡ��������������  �޸ģ�Ӧ�ý������Ϊ����ֵ���أ��������ó�Ա�����洢
	std::vector<cv::Point2d> ExtractCenters(Mat src, Mat& imgHSVMask, vector<Rect>& rois) {
		int w, brightValue;

		// �����������w
		estimateWidth(imgHSVMask, rois, w, brightValue);

		vector<Point2d> centers = Centers(imgHSVMask, rois, w);
		//vector<Point2d> centers = extractLaserCentersGravity(imgHSVMask, rois, w);
		//cout << "Centers: " << centers << endl;

		Mat huidu_img0;
		src.copyTo(huidu_img0);
		drawPointsOnImage(centers, huidu_img0);

		// ����뾶k
		double k = 5.0;
		//�� �뾶k �����ĵ����У��
		for (const auto& center_point : centers) {
			Point2d smoothed_point = SmoothedPoint(centers, center_point, k);
			smooth_points.push_back(smoothed_point);
		}

		//// ����һ���� src ͬ�ߴ�İ�ɫͼ��
		//Mat white_img = Mat::zeros(g_midImage.size(), CV_8UC3); // ����һ��ȫ��ͼ��
		//white_img.setTo(Scalar(255, 255, 255)); // ��ͼ������Ϊȫ��
		//drawPointsOnImage(smooth_points, white_img, (0, 0, 0));

		Mat huidu_img1;
		src.copyTo(huidu_img1);
		//cvtColor(huidu_img1, huidu_img1, COLOR_GRAY2BGR);
		drawPointsOnImage(smooth_points, huidu_img1);
		return smooth_points;
	}

	vector<Point2d> Center_points; //����������
	vector<Point2d> smooth_points; //K����У�������ĵ�

private:
	Mat src;//ԭͼ
	Mat distCoeff;//�������ϵ��
	Mat instrinsic_Mat;//����ڲ�
	Mat g_dstImage, g_midImage;

	// ------------------------��ͼ���ϻ��Ƶ㼯------------------------
	void drawPointsOnImage(vector<Point2d>& points, Mat& image, Scalar color = Scalar(0, 0, 255)) {

		for (const auto& point : points) {
			// ��С������ת��Ϊͼ���е��������꣨����С�����ȣ�
			Point2d pixelPoint(point.x, point.y);

			// ����С�����굽ͼ�����꣬�Ի��ƽ��ư뾶Ϊ0.1���ص�ʵ�ĵ�
			Point scaledPoint(cvRound(pixelPoint.x), cvRound(pixelPoint.y));

			// ����ʵ��Բ
			circle(image, scaledPoint, 0, color, 1);  // �뾶Ϊ0��ʵ��Բ

		}
	}


	// ------------------------ROI------------------------
	// �Զ���ѡroi 	 �޸ģ��������ü���
	Mat autoroi(Mat img, int size, vector<Rect>& rois) {
		Mat cimg;
		img.copyTo(cimg);
		threshold(cimg, cimg, 0, 255, THRESH_BINARY);
		//cvtColor(cimg, cimg, COLOR_BGR2GRAY);

		// �������������룬���� CV_CHAIN_APPROX_NONE ��ʾ��ѹ����Ե�㣩
		vector<vector<Point>> contours;
		vector<Vec4i> hierarchy;
		findContours(cimg, contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_NONE);

		//Mat contourImg = Mat::zeros(cimg.size(), CV_8UC3);
		//drawContours(contourImg, contours, -1, Scalar(0, 255, 0), 1);
		//imshow("Contours", contourImg);
		//waitKey(0);

		// �����ڵ�ͼ�����ڻ���ROI����
		Mat out_Image = Mat::zeros(img.size(), img.type());

		int padding = 2; // ����������չ����

		for (size_t i = 0; i < contours.size(); i++) {
			// �����������
			double area = contourArea(contours[i]);

			if (area > size) { // �����ֵɸѡ
				Rect rect = boundingRect(contours[i]);

				// ��չROI����
				Rect roi(rect.x - padding, rect.y - padding,
					rect.width + 2 * padding, rect.height + 2 * padding);

				// ������ͼ��Χ��
				roi &= Rect(0, 0, img.cols, img.rows);

				rois.push_back(roi);  // ���� ROI ��

				// ����ԭͼ�����ͼ
				img(roi).copyTo(out_Image(roi));

				// ����չ��ľ��ο�
				//rectangle(out_Image, roi, Scalar(255), 1);
			}
		}
		return out_Image;
	}

	// �ֶ���ѡroi
	struct ROIData {
		Mat temp_Image;
		Mat roi_Image;
		Mat roi_result_Image;
		Point pre_pt;
		Point cur_pt;
	};

	static void on_mouse(int event, int x, int y, int flags, void* userdata)
	{
		ROIData* data = static_cast<ROIData*>(userdata); // ͨ��userdata���ʽṹ��
		Mat tmp, dst;

		if (event == EVENT_LBUTTONDOWN)
		{
			data->pre_pt = Point(x, y);  // ��¼��갴��λ��
		}
		else if (event == EVENT_MOUSEMOVE && flags) // ����ƶ�ʱ��ʵʱ���ƾ��ο�
		{
			data->temp_Image.copyTo(tmp);
			data->cur_pt = Point(x, y);
			rectangle(tmp, data->pre_pt, data->cur_pt, Scalar(255, 255, 255), 1, 8, 0);
			imshow("img", tmp);  // ʵʱ��ʾROI����
		}
		else if (event == EVENT_LBUTTONUP)  // �ɿ������������ROI��ѡ
		{
			data->temp_Image.copyTo(data->roi_Image);
			rectangle(data->roi_Image, data->pre_pt, data->cur_pt, Scalar(255, 255, 255), 1, 8, 0);
			imshow("img", data->roi_Image);

			int width = abs(data->pre_pt.x - data->cur_pt.x);
			int height = abs(data->pre_pt.y - data->cur_pt.y);
			dst = data->temp_Image(Rect(min(data->cur_pt.x, data->pre_pt.x), min(data->cur_pt.y, data->pre_pt.y), width, height));

			Mat out_Image(data->temp_Image.size(), data->temp_Image.type(), Scalar(0)); // ����ȫ��ͼ
			Rect roi_rect = Rect(min(data->cur_pt.x, data->pre_pt.x), min(data->cur_pt.y, data->pre_pt.y), width, height);
			dst.copyTo(out_Image(roi_rect));
			out_Image.copyTo(data->roi_result_Image);
		}
	}

	Mat handroi(Mat img)
	{
		ROIData data;  // �����ṹ���Դ洢�ص�����
		img.copyTo(data.temp_Image);  // ������ͼ���Ƶ��ṹ��
		namedWindow("img", WINDOW_NORMAL);
		imshow("img", img);

		setMouseCallback("img", on_mouse, &data);  // ���ݽṹ��ָ��
		waitKey(0);
		destroyWindow("img");

		return data.roi_result_Image;  // ���ؽ��ͼ��
	}


	// ----------------------�鿴���طֲ�-------------------------
	// ��ͼ����ĳһ�еĻҶ�ֵ����Ϊ txt �ļ�
	void saveRowToTxt(const Mat& gray, int y) {
		CV_Assert(gray.channels() == 1 && y >= 0 && y < gray.rows);

		// �����ļ���
		string filename = "��" + to_string(y) + "��.txt";

		ofstream fout(filename);
		if (!fout.is_open()) {
			cerr << "�޷����ļ� " << filename << " ����д�룡" << endl;
			return;
		}

		const uchar* row = gray.ptr<uchar>(y);
		for (int x = 0; x < gray.cols; ++x) {
			fout << static_cast<int>(row[x]) << endl;  // һ��һ���Ҷ�ֵ
		}

		fout.close();
	}

	// ��roi����ÿ�����ؽ��б���
	void saveMultipleRowsToTxtInROIs(const  Mat& gray, const vector< Rect>& rois, const string& filename = "���лҶ�ֵ.txt") {
		CV_Assert(gray.channels() == 1);

		ofstream fout(filename);
		if (!fout.is_open()) {
			cerr << "�޷����ļ� " << filename << " ����д�룡" << endl;
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
				fout << endl; // ����
			}
		}

		fout.close();
		cout << "�ѱ�����ѡ ROI �еĶ��лҶ�ֵ���ļ���" << filename << endl;
	}

	// ���� ����0
	void saveMultipleRowsToTxtInROI(const  Mat& gray, const vector< Rect>& rois, const string& filename = "���лҶ�ֵ.txt") {
		CV_Assert(gray.channels() == 1);

		ofstream fout(filename);
		if (!fout.is_open()) {
			cerr << "�޷����ļ� " << filename << " ����д�룡" << endl;
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

				fout << endl; // ÿ��һ�� ROI ������
			}
		}

		fout.close();
		cout << "�ѱ������� ROI �����з���Ҷ�ֵ���ļ���" << filename << endl;
	}


	// ----------------------����Ӧ��ȡ������------------------------
	// �������w
	void estimateWidth(const Mat& gray, const vector<Rect>& rois, int& w, int& brightValue) {
		vector<uchar> nonZeroValues;
		//����ROI���������еķ�0�Ҷ�ֵ
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
			cout << "�޷������أ��޷�����������ֵ���ȡ�" << endl;
			w = 0;
			brightValue = 255;
			return;
		}

		sort(nonZeroValues.begin(), nonZeroValues.end(), greater<uchar>());
		int index = max(1, static_cast<int>(nonZeroValues.size() * 0.2)) - 1;
		brightValue = nonZeroValues[index];
		cout << "������ֵ��ǰ20%����С���ȣ�= " << brightValue << endl;

		// ʹ�ô� brightValue ������
		int totalPixels = 0;
		int rowCount = 0;

		//������ ROI �е�ÿһ�У�ͳ�Ƹ����д��ڵ��� brightValue �����ظ������ۼ���ƽ�������յõ����ƿ�� w
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
		cout << "���Ƶļ������ƿ�� w = " << w << " ����" << endl;
	}

	// ��ʼ���ĵ�ѡȡ
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

		// ȡ���ֵλ��ƽ��
		int sumX = accumulate(maxPositions.begin(), maxPositions.end(), 0);
		return static_cast<int>(sumX / maxPositions.size());
	}

	// ��˹��ϣ�y = A * exp(-(x - mu)^2 / (2 * sigma^2))����ȡ����λ�� mu
	double gaussFit(const vector<uchar>& profile) {
		vector<Point2f> points;

		for (int i = 0; i < profile.size(); ++i) {
			if (profile[i] > 30) // �ɵ���ֵ
				points.emplace_back(i, profile[i]);
		}

		if (points.size() < 5) return -1;

		double A = *max_element(profile.begin(), profile.end());
		double mu = distance(profile.begin(), max_element(profile.begin(), profile.end()));
		double sigma = 1.5; // �̶���Ԥ��ֵ���ȶ�

		// �� Gauss-Newton���ֹ�������� mu
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

	// ===== �Ҷ����� =====
	Point2f computeGravityCenter(const vector<uchar>& profile, int offsetX) {
		double sumI = 0, sumIX = 0;
		for (int i = 0; i < profile.size(); ++i) {
			sumI += profile[i];
			sumIX += i * profile[i];
		}
		if (sumI <= 0) return Point2f(-1, -1);
		return Point2f(sumIX / sumI + offsetX, 0);
	}

	// ===== ��˹������� =====
	double gaussFullFitCenter(const vector<uchar>& profile) {
		if (profile.size() < 5) return -1;

		double A = *max_element(profile.begin(), profile.end());
		if (A < 30) return -1; // �Ҷ�̫�ͺ���

		// ���ø�������
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

	// ===== ��������ȡ������ =====
	vector<Point2d> extractLaserCentersGravity(const Mat& gray, const vector<Rect>& rois, int w) {
		vector<Point2d> centers;

		for (const Rect& roi : rois) {
			for (int row = roi.y; row < roi.y + roi.height; ++row) {
				// 1. �ҳ�ʼ���ĵ�
				int x0 = initialCenter(gray, row, roi); // �û�����ʵ��
				if (x0 < 0) continue;

				// 2. �趨���� 2w ��Χ
				int left = max(roi.x, x0 - 2 * w);
				int right = min(roi.x + roi.width - 1, x0 + 2 * w);
				if (right <= left + 5) continue;

				// 3. ��ȡ profile
				vector<uchar> profile;
				for (int x = left; x <= right; ++x)
					profile.push_back(gray.at<uchar>(row, x));

				// 4. ƽ��
				Mat profMat(profile, true);        // ת�������� Mat
				profMat = profMat.reshape(1, 1);   // ��ʽ reshape ��һ��
				GaussianBlur(profMat, profMat, Size(5, 1), 1.0);  // ��˹�˲�
				profile.assign(profMat.begin<uchar>(), profMat.end<uchar>());  // ��ֵ�� vector


				// 5. �Ҷ����Ĺ��Ƴ�ֵ
				Point2f gravity = computeGravityCenter(profile, left);
				int center_x = static_cast<int>(gravity.x - left);

				// 6. �ֲ���ȡ���ڸ�˹���
				int sub_left = max(0, center_x - w);
				int sub_right = min((int)profile.size() - 1, center_x + w);
				// �˴���sub_left  sub_right  ���ݴ������⣬Խ���ˡ���Ҫ����޸ġ�

				vector<uchar> sub_profile(profile.begin() + sub_left, profile.begin() + sub_right + 1);

				// 7. ������ں�
				double mu = gaussFullFitCenter(sub_profile);
				if (mu >= 0)
					centers.emplace_back(mu + left + sub_left, row);
				else
					centers.emplace_back(gravity.x, row);
			}
		}

		return centers;
	}

	// ����Ҷ�����
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

	// �����ػ�y���� (δʹ��)
	double refineYSubpixel(const  Mat& gray, int x, int y) {
		if (y <= 0 || y >= gray.rows - 1 || x <= 0 || x >= gray.cols - 1)
			return static_cast<double>(y);  // �߽粻������

		// ��ȡ���������и�x���Ҷ�ֵ
		double I1 = gray.at<uchar>(y - 1, x);
		double I2 = gray.at<uchar>(y, x);
		double I3 = gray.at<uchar>(y + 1, x);

		double denom = (I1 - 2 * I2 + I3);
		if (abs(denom) < 1e-3) return static_cast<double>(y);

		// �����߼�ֵ�� y_sub = y - 0.5 * (I3 - I1) / (I1 - 2*I2 + I3)
		double delta = 0.5 * (I1 - I3) / denom;
		return y + delta;
	}

	// ��ȡ������
	vector<Point2d> Centers(const Mat& gray, const vector<Rect>& rois, int w) {
		vector<Point2d> centers;

		for (const Rect& roi : rois) {
			for (int row = roi.y; row < roi.y + roi.height; ++row) {

				// 1. �ҳ�ʼ���ĵ�
				int x0 = initialCenter(gray, row, roi);
				if (x0 < 0) continue;

				// 2. ȷ�����ĵ����� 2w ��Χ
				int left = max(roi.x, x0 - 2 * w);
				int right = min(roi.x + roi.width - 1, x0 + 2 * w);

				if (right <= left + 5) continue; // ��֤�������㹻

				// 3. ��ȡ�ո�ȷ���� 4w ����ĻҶ�ֵ
				vector<uchar> profile;
				for (int x = left; x <= right; ++x) {
					profile.push_back(gray.at<uchar>(row, x));
				}

				// 4. �������ڸ�˹�˲� �� ƽ������������Ӱ��
				Mat profMat(profile);
				profMat = profMat.reshape(1, 1);
				GaussianBlur(profMat, profMat, Size(5, 1), 1.0);
				profile.assign(profMat.begin<uchar>(), profMat.end<uchar>());

				// 5. �Ҷ�����
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


	//------------------------��ȡ������------------------------
	void huidu(Mat PointImage, Mat src, vector<Point2d>& centers)
	{
		// ��ȡPointImage�е����ĵ�
		vector<Point2d> points;
		point(PointImage, points);
		vector<double> kcal;

		for (int i = 0; i < points.size() - 1; i++) {

			double pfCosSita = 0, pfSinSita = 0;
			CalcNormVec(Point2d(points[i - 1].x, points[i - 1].y), Point2d(points[i].x, points[i].y), Point2d(points[i + 1].x, points[i + 1].y), pfCosSita, pfSinSita);
			// ��������������CalcNormVec
			//-------------------�Ҷ����ķ�gdd---------------------//
			double sum = 0, sum_sumx = 0, sum_sumy = 0;
			for (int j = 0; j < 2; j++)		// ��õ��Լ��������������㣬һ��������ĻҶ����ģ��õ������ؼ�������ĵ�
			{
				if (j == 0)
				{
					double cj = points[i].x;
					double ci = points[i].y;
					sum = ijpixel(cj, ci, src);//�Ҷ�ֵ�ܺ�    ��ʱci��cj��������ijpixel�õ��ľ��ǣ�x��y�����ĻҶ�ֵ
					sum_sumx = ijpixel(cj, ci, src) * cj;//��x*�Ҷ�ֵ�����ܺ�
					sum_sumy = ijpixel(cj, ci, src) * ci;//��y*�Ҷ�ֵ�����ܺ�
				}
				else
				{
					//cout << "pfCosSita = " << pfCosSita << endl;
					double x_cor = points[i].x + j * pfCosSita;	//cor��(��j=0(��һ��)��ķ�������)��x���꣺x + cos
					double y_cor = points[i].y + j * pfSinSita;
					double x_cor1 = points[i].x - j * pfCosSita;	//cor_1���x����
					double y_cor1 = points[i].y - j * pfSinSita;//��ʱx_cor����С����ijpixel�õ��ľ��ǣ�x_cor��y_cor���ȴ��ĻҶ�ֵ
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

	//��i��j��������ֵ
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
	//���㷨����
	void CalcNormVec(Point2d ptA, Point2d ptB, Point2d ptC, double& pfCosSita, double& pfSinSita)
	{
		double fVec1_x, fVec1_y, fVec2_x, fVec2_y;
		if (ptA.x == 0 && ptA.y == 0)
		{
			ptA.x = ptC.x;
			ptA.y = ptC.y;
			//����B�������A�����ꡣ
			fVec1_x = -(ptB.x - ptA.x);
			fVec1_y = -(ptB.y - ptA.y);
		}
		else
		{
			//����B�������A�����ꡣ
			fVec1_x = ptB.x - ptA.x;
			fVec1_y = ptB.y - ptA.y;
		}

		if (ptC.x == 0 && ptC.y == 0)
		{
			ptC.x = ptA.x;
			ptC.y = ptA.y;
			//����C�������B�����ꡣ
			fVec2_x = (ptB.x - ptC.x);
			fVec2_y = (ptB.y - ptC.y);
		}
		else
		{
			//����C�������B�����ꡣ
			fVec2_x = ptC.x - ptB.x;
			fVec2_y = ptC.y - ptB.y;
		}

		//��λ����
		double fMod = sqrt(fVec1_x * fVec1_x + fVec1_y * fVec1_y);
		fVec1_x /= fMod;
		fVec1_y /= fMod;
		//���㴹�ߡ�
		double fPerpendicularVec1_x = -fVec1_y;
		double fPerpendicularVec1_y = fVec1_x;


		//��λ����
		fMod = sqrt(fVec2_x * fVec2_x + fVec2_y * fVec2_y);
		fVec2_x /= fMod;
		fVec2_y /= fMod;
		//���㴹�ߡ�
		double fPerpendicularVec2_x = -fVec2_y;
		double fPerpendicularVec2_y = fVec2_x;

		//��͡�
		double fSumX = fPerpendicularVec1_x + fPerpendicularVec2_x;
		double fSumY = fPerpendicularVec1_y + fPerpendicularVec2_y;
		//��λ����
		fMod = sqrt(fSumX * fSumX + fSumY * fSumY);
		double fCosSita = fSumX / fMod;
		double fSinSita = fSumY / fMod;
		pfCosSita = fCosSita;
		pfSinSita = fSinSita;
	}










	//------------------------У�����ĵ�------------------------
	Point2d SmoothedPoint(const vector<Point2d>& center_points, Point2d a, double k) {
		double numerator_x = 0.0;
		double numerator_y = 0.0;
		double denominator = 0.0;

		for (const auto& p : center_points) {
			double distance = norm(a - p);
			if (distance < k && distance != 0) { // ֻ���ǰ뾶k���ڵĵ��Ҳ�������ǰ���ĵ�
				double delta = 1.0 / distance; // ����Ȩ�ئ�i
				numerator_x += delta * p.x;    // ʹ��p.x�����Ȩ�ܺ�
				numerator_y += delta * p.y;    // ʹ��p.y�����Ȩ�ܺ�
				denominator += delta;          // ����Ȩ���ܺ�
			}
		}

		if (denominator != 0) {
			double corrected_x = numerator_x / denominator; // �����Ȩƽ��x����
			double corrected_y = numerator_y / denominator; // �����Ȩƽ��y����
			return Point2d(corrected_x, corrected_y);
		}
		else {
			return a; // ����뾶��û�е㣬�򷵻�ԭʼ��
		}
	}



	//����������
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


//��ȡ����������彻�㣬��תΪ�������ϵ
// ��һ����Set_param ����궨���ڲΡ�ƽ���������궨����ת�������궨��ͼƬ��2D���������ݡ���---�뼤���ߴ�ֱ����---�뼤����ƽ�У���ɳ�ʼ����
// �ڶ�����Get_plane ��ȡ�궨��ƽ�淽�̡�
// ��������Get_Intersection ��ȡ��������궨��ֱ�ߵĽ��㡣
// ���Ĳ���LaserPoints_cam ����������תΪ�������ϵ��
// ���岽������ laser_cam ��ȡ�������ϵ�µ㡣
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
	* @brief �궨�����תΪ�궨��ƽ�淽�� Ax+By+Cz+D=0
	*/
	void Get_plane() {
		//�����������������ϳ�һ��4��4����,external_Mat
		Mat raMat, tempMat, externalMat;
		Mat bottomline = (Mat_<double>(1, 4) << 0, 0, 0, 1);
		//�޵����˹��ʽ
		Rodrigues(rotationMat, raMat);
		hconcat(raMat, translationMat, tempMat);//ˮƽƴ��
		vconcat(tempMat, bottomline, externalMat);//��ֱƴ��

		//�������������棬������궨������ƽ���ƽ�淽�̣�����ϵΪ�������ϵ
		Mat externalMat_invert = externalMat.inv();	//����
		Mat ex_param = (Mat_<double>(1, 4) << 0, 0, 1, 0);
		Mat planeEquation = ex_param * externalMat_invert;	//ƽ�淽��=���*

		plane = planeEquation;
	}
	/**
	* @brief ��ȡ��������궨��ֱ�ߵĽ���
	*/
	//��ͼ���ϻ��Ƶ㼯
	void drawPointsOnImage(vector<Point2d>& points, Mat& image, Scalar color = Scalar(0, 0, 255)) {  // Ĭ�Ϻ�ɫ
		for (const auto& point : points) {
			// ��С������ת��Ϊͼ���е��������꣨����С�����ȣ�
			Point2d pixelPoint(point.x, point.y);

			// ����С�����굽ͼ�����꣬�Ի��ƽ��ư뾶Ϊ0.1���ص�ʵ�ĵ�
			Point scaledPoint(cvRound(pixelPoint.x), cvRound(pixelPoint.y));

			// ����ʵ��Բ
			circle(image, scaledPoint, 0, color, -1);  // �뾶Ϊ0��ʵ��Բ
		}
	}


	void Get_Intersection(){
		Size boardSize;
		boardSize.width = chessboard_width;//��---�뼤���ߴ�ֱ
		boardSize.height = chessboard_height;//��---�뼤����ƽ��

		vector<double> laserline;	//������
		// ���ĵ�ȥ����
		undistortPoints(laser_pix, undistorted_centers, instrinsic_Mat, distCoeff);
		// ��һ����ת����������ϵ
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
			for (size_t i = 0; i < boardSize.height; i++) // �������̸��ÿһ��
			{
				vector<Point2d> tps;
				for (size_t j = 0; j < boardSize.width; j++) // ������ǰ�е�ÿ���ǵ�
				{
					Point2d tp;
					tp = pointbuf[j + i * boardSize.width]; // ��ȡ��ǰ�ǵ�����
					tps.push_back(tp); // ���ǵ���뵱ǰ�еĵ㼯��
				}

				vector<double> tl;
				tl = line_fitting(tps); // �Ե�ǰ�еĽǵ㼯�Ͻ���ֱ����ϣ��õ����нǵ��ֱ�߷���
				myline.push_back(tl); // ����ϵõ���ֱ�߲������� myline

				Point2d ips = intersection(laserline, tl); // ���㵱ǰ��ֱ���뼤���ߵĽ���
				allp.push_back(ips); // ��������뵽���㼯�� allp ��
			}

			laser_chess_p = allp;

			Mat huidu_img2;
			image.copyTo(huidu_img2);
			drawPointsOnImage(laser_chess_p, huidu_img2);

			ofstream ofs("����.txt", ios::trunc);
			for (const auto& pt : laser_chess_p) {
				ofs << pt.x << "\t" << pt.y << endl;
			}
			ofs.close();
		}
		else {
			cout << "Chessboard corners not found." << endl;
		}
	}


	//** @brief ��������ϵת�������ϵ

	vector<Point3d> LaserPoints_cam(){
		//ƽ�淽�̵ķ�������
		double a1 = plane.at<double>(0, 0);
		double b1 = plane.at<double>(0, 1);
		double c1 = plane.at<double>(0, 2);
		double d1 = plane.at<double>(0, 3);
		cout << "ƽ�淽�̵ķ���������" << a1 << "\t" << b1 << "\t" << c1 << "\t" << d1 << endl;

		Mat instrinsic_Mat_invert = instrinsic_Mat.inv();	//�ڲξ�������

		vector<Point3d> camPoints;
		for (int i = 0; i < laser_chess_p.size(); i++){
			//���������ϵĵ�����������ϵ��ΪpixelPoint����u,v,1��
			Mat pixelPoint = (Mat_<double>(3, 1) << laser_chess_p[i].x, laser_chess_p[i].y, 1);
			//(u,v,1)��������ڲ���������棬���Ǽ��������������ڹ�һ��ƽ������꣬λ���������ϵ
			Mat imagePoint = instrinsic_Mat_invert * pixelPoint;

			double x1 = imagePoint.at<double>(0, 0);
			double y1 = imagePoint.at<double>(1, 0);
			double z1 = 1.0;  // ��һ��ƽ���� z ������Ϊ 1

			 //�������ߺ�ƽ�潻��
			double denominator = a1 * x1 + b1 * y1 + c1 * z1;
			if (denominator == 0) {
				cout << "warning��������ƽ��ƽ�У��޷����㽻��" << endl;
				continue;
			}

			// ���㽻��ı������� lambda
			double lambda = -d1 / denominator;

			// �õ��������ϵ�е���ά����
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
	Mat image;//�궨��ͼƬ
	Mat plane;
	vector<Point2d> laser_pix;//��������ϵ�µļ���������
	vector<Point2d> laser_chess_p;
	int chessboard_width, chessboard_height;
	vector<Point2d> undistorted_centers;	// ������ȥ�����

	vector<double> line_fitting(vector<Point2d> input){
		Mat dst = Mat(2, 2, CV_64FC1, Scalar(0));
		//��С���˷�
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
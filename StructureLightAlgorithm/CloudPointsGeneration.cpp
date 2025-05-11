#include"CloudPointsGeneration.h"

CloudPointsGeneration::CloudPointsGeneration(const std::string& CloudPointsSavePath)
	:M_CloudPointsSavePath(CloudPointsSavePath)
{
	
};


cv::Mat CloudPointsGeneration::CalculateDisparity(cv::Mat& UnwrappedLeft,
												  cv::Mat& UnwrappedRight,
												  cv::Mat& disparity,
												  int maxDisparity) {
	disparity = cv::Mat::zeros(UnwrappedLeft.size(), CV_32FC1);

	for (int row = 0; row < UnwrappedLeft.rows; row++) {
		float* leftRowPtr = UnwrappedLeft.ptr<float>(row);
		float* rightRowPtr = UnwrappedRight.ptr<float>(row);

		for (int col = 0; col < UnwrappedLeft.cols; col++) {
			float leftValue = leftRowPtr[col];
			float rightValue = rightRowPtr[col];
			//这一步有一些问题，这样回导致一些么问题的视野被滤除掉
			//if (leftValue == 0 || rightValue == 0) {
			//	// Remove background
			//	disparity.ptr<float>(row)[col] = 0.0;
			//	continue;
			//}

			if (leftValue == 0 && rightValue == 0) {
				// Remove background
				disparity.ptr<float>(row)[col] = 0.0;
				continue;
			}

			int flag = col; // Initialize flag with the current column
			float minDiff = std::numeric_limits<float>::max();

			// Define the search range based on maxDisparity
			int minSearchCol = std::max(0, col - maxDisparity);
			int maxSearchCol = std::min(UnwrappedRight.cols - 1, col + maxDisparity);

			for (int j = minSearchCol; j <= maxSearchCol; j++) {
				float diff = std::abs(leftValue - rightRowPtr[j]);

				if (diff < minDiff) {
					minDiff = diff;
					flag = j;
				}
			}
			//Control parallax values not to be too large
			if (minDiff > 0.5) {
				disparity.ptr<float>(row)[col] = 0;
			}
			else if (flag == col || flag == 0) { // Check if flag is at the border
				disparity.ptr<float>(row)[col] = static_cast<float>(col - flag);
			}
			else {
				// Ensure flag-1 is within the valid range
				float rightValuePrev = (flag > 0) ? rightRowPtr[flag - 1] : rightRowPtr[flag];
				double k = (leftValue - rightValuePrev) / (rightRowPtr[flag] - rightValuePrev);
				disparity.ptr<float>(row)[col] = static_cast<float>(col - (flag - 1 + k));
			}
		}
	}
	return disparity;
};


cv::Mat CloudPointsGeneration::CalculateDepthMap(cv::Mat& Disparity, 
												 cv::Mat& Q,	
												 float mindep, 
												 float maxdep,
												 float minDisparity,
												 float maxDisparity) {
	//Coordinates of the intersection of the camera's main optical axis 
	//and the image plane (principal point)
	const double cx = -Q.at<double>(0, 3);
	const double cy = -Q.at<double>(1, 3);
	//Focal length after reprojection
	const double f = Q.at<double>(2, 3);
	//Distance between projection centers of two cameras
	const double Tx = -1 / Q.at<double>(3, 2);
	//Difference between the x-coordinate of the main point of the left camera 
	//and the x-coordinate of the main point of the right camera	
	const double cxlr = Q.at<double>(3, 3) * Tx;
	//The z-coordinate of the point cloud
	double z;

	cv::Mat Depth = cv::Mat::zeros(Disparity.size(), Disparity.type());
	//Depth calculation based on parallax
	for (int row = 0; row < Disparity.rows; row++) {
		for (int col = 0; col < Disparity.cols; col++) {
			//过滤视差图中的错误信息
			if (Disparity.ptr<float>(row)[col] == 0 || Disparity.ptr<float>(row)[col] == NAN ||
				Disparity.ptr<float>(row)[col] > maxDisparity || Disparity.ptr<float>(row)[col] < minDisparity)
			{
				continue;
			}
			float d = Disparity.ptr<float>(row)[col];
			z = f * Tx / (-d + cxlr);//depth = Z/W = f * T / -disparity + (Cx - Cx')
			if (z > maxdep || z < mindep) {
				continue;
			}
			if (isnan(z)) {
				Depth.ptr<float>(row)[col] = 0;
				continue;
			}
			Depth.ptr<float>(row)[col] = z;
		}
	}
	//Do some filtering on the image
	cv::Mat th;
	cv::Mat s;
	cv::threshold(Depth, th, 1, 255, cv::THRESH_BINARY);
	th.convertTo(th, CV_8UC1, 255);
	cv::Mat kerenl = cv::getStructuringElement(0, cv::Size(3, 3));
	cv::erode(th, th, kerenl, cv::Point(-1, -1), 3);
	Depth.copyTo(s, th);

	return Depth;
};

//优化版：提升运算速度
cv::Mat CloudPointsGeneration::CalculateDepthMap_(cv::Mat& Disparity,
												  cv::Mat& Q,
												  cv::Mat& DepthMap,
												  float mindep,
												  float maxdep,
												  float minDisparity,
												  float maxDisparity) {
	//摄像机主光轴与图像平面（主点）的交点坐标
	const double cx = -Q.at<double>(0, 3);
	const double cy = -Q.at<double>(1, 3);
	//重影后的焦距
	const double f = Q.at<double>(2, 3);
	//两台摄像机投影中心之间的距离
	const double Tx = -1 / Q.at<double>(3, 2);
	//左摄像机主点 x 坐标与右摄像机主点 x 坐标之差
	const double cxlr = Q.at<double>(3, 3) * Tx;
	//点云的 Z 坐标
	double z;

	cv::Mat Depth = cv::Mat::zeros(Disparity.size(), Disparity.type());
	//Depth calculation based on parallax
	for (int row = 0; row < Disparity.rows; row++) {
		for (int col = 0; col < Disparity.cols; col++) {
			//过滤视差图中的错误信息
			if (Disparity.ptr<float>(row)[col] == 0 || Disparity.ptr<float>(row)[col] == NAN ||
				Disparity.ptr<float>(row)[col] > maxDisparity || Disparity.ptr<float>(row)[col] < minDisparity)
			{
				continue;
			}
			float d = Disparity.ptr<float>(row)[col];
			z = f * Tx / (-d + cxlr);//depth = Z/W = f * T / -disparity + (Cx - Cx')
			if (z > maxdep || z < mindep) {
				continue;
			}
			if (isnan(z)) {
				Depth.ptr<float>(row)[col] = 0;
				continue;
			}
			Depth.ptr<float>(row)[col] = z;
		}
	}
	//Do some filtering on the image
	cv::Mat th;
	cv::Mat s;
	cv::threshold(Depth, th, 1, 255, cv::THRESH_BINARY);
	th.convertTo(th, CV_8UC1, 255);
	cv::Mat kerenl = cv::getStructuringElement(0, cv::Size(3, 3));
	cv::erode(th, th, kerenl, cv::Point(-1, -1), 3);
	Depth.copyTo(s, th);

	return Depth;





	//// 预计算固定值
	//const double cx = -Q.at<double>(0, 3);
	//const double cy = -Q.at<double>(1, 3);
	//const double f = Q.at<double>(2, 3);
	//const double Tx = -1 / Q.at<double>(3, 2);
	//const double cxlr = Q.at<double>(3, 3) * Tx;

	//DepthMap = cv::Mat::zeros(Disparity.size(), Disparity.type());

	//// 循环处理视差图
	//for (int row = 0; row < Disparity.rows; row++) {
	//	for (int col = 0; col < Disparity.cols; col++) {
	//		float d = Disparity.at<float>(row, col);

	//		if (d == 0 || std::isnan(d) || d < minDisparity || d > maxDisparity) {
	//			continue;
	//		}

	//		if (-d + cxlr == 0) {
	//			continue;
	//		}

	//		double z = f * Tx / (-d + cxlr);

	//		if (z > maxdep || z < mindep || std::isnan(z)) {
	//			Depth.at<float>(row, col) = 0;
	//			continue;
	//		}

	//		Depth.at<float>(row, col) = static_cast<float>(z);
	//	}
	//}

	//// 形态学处理
	//cv::Mat th;
	//cv::threshold(Depth, th, 1, 255, cv::THRESH_BINARY);
	//th.convertTo(th, CV_8UC1, 255);
	//cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3));
	//cv::erode(th, th, kernel, cv::Point(-1, -1), 3);
	//Depth.copyTo(Depth, th);

	//return Depth;
};





void CloudPointsGeneration::SaveCloudPointsToTxt(cv::Mat& DepthMap,
												 cv::Mat& Q, 
												 std::string& CloudPointsPath) {
	float ratios = 1;
	//Coordinates of the intersection of the camera's main optical axis 
	//and the image plane (principal point)
	const double cx = -Q.at<double>(0, 3);
	const double cy = -Q.at<double>(1, 3);
	const double q00 = -Q.at<double>(0, 0);
	const double q11 = -Q.at<double>(1, 1);
	//Focal length after reprojection
	const double f = Q.at<double>(2, 3);
	//Distance between projection centers of two cameras
	const double Tx = -1 / Q.at<double>(3, 2);
	//Difference between the x-coordinate of the main point of the left camera 
	//and the x-coordinate of the main point of the right camera
	const double cxlr = Q.at<double>(3, 3) * Tx;

	std::fstream outfile;
	outfile.open(CloudPointsPath, std::ios::trunc | std::ios::out);
	if (!outfile) {
		std::cout << "file can not open" << std::endl;
		return;
	}
	for (int row = 0; row < DepthMap.rows; row++) {
		for (int col = 0; col < DepthMap.cols; col++) {
			if ((DepthMap.at<float>(row, col) == 0) || (DepthMap.at<float>(row, col) == NAN)) {
				//continue;
			}
			else {
				float z = DepthMap.at<float>(row, col);
				uint32_t rgb;
				std::string str;
				//Generate 3D coordinates
				double x = (z / f * (col - cx) * 1.00) * ratios;
				double y = (z / f * (row - cy) * 1.00) * ratios;
				double z1 = z * ratios;
				outfile << x << " " << y << " " << z1 << std::endl;
			}
		}
	}
	outfile.close();
};



CloudPointsGeneration::~CloudPointsGeneration() {

};
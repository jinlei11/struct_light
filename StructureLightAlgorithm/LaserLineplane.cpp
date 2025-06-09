// 光平面标定，修改了中心线的提取方法。2025.05.15

#include"include/LaserLineplane.h"
#include <Eigen/Dense>


using namespace std;
using namespace cv;



OpticalPlaneCalibration::OpticalPlaneCalibration(const std::string& BoardPath, const std::string& OpticalPlanePath, enum CalibrationPlate boardClass)
	:M_PlaneBoardFilePath(BoardPath), M_PlaneLineFilePath(OpticalPlanePath)
{
	switch (boardClass)
	{
	case Plate_GP050:
		this->M_board_size = cv::Size(11, 8);
		this->M_squareSize = cv::Size(3, 3);
		this->M_CalibPlateType = 0;
		break;
	case Plate_doubleCircle:
		this->M_board_size = cv::Size(7, 5);
		this->M_squareSize = cv::Size(9, 9);
		this->M_CalibPlateType = 1;
		break;
	case Plate_bigCircle:
		this->M_board_size = cv::Size(13, 10);
		this->M_squareSize = cv::Size(15, 15);
		this->M_CalibPlateType = 2;
		break;
	case Plate_Big_chess:
		this->M_board_size = cv::Size(9, 6);
		this->M_squareSize = cv::Size(50, 50);
		this->M_CalibPlateType = 3;
	case Plate_smallCircle:
		this->M_board_size = cv::Size(12, 9);
		this->M_squareSize = cv::Size(3, 3);
		this->M_CalibPlateType = 4;

	case Plate_chess811:
		this->M_board_size = cv::Size(8, 11);
		this->M_squareSize = cv::Size(3, 3);
		this->M_CalibPlateType = 5;

	default:
		break;
	}

};




// 计算光平面标定图像的r t 
void OpticalPlaneCalibration::getRT(const Mat& cameraMatrix, const Mat& distCoeffs,
	       std::vector<cv::Mat>& rotationMat, std::vector<cv::Mat>& translationMat, 
		   bool show = false) {

	// 构建棋盘格世界坐标 objectPoints
	vector<Point3f> objPoints;
	for (int i = 0; i < M_board_size.height; ++i) {
		for (int j = 0; j < M_board_size.width; ++j) {
			objPoints.emplace_back(j * M_squareSize.width, i * M_squareSize.width, 0.0f);
		}
	}

	// 打开文件并清空旧内容
	std::ofstream outFile("extrinsic_parameters.txt", std::ios::out | std::ios::trunc);
	if (!outFile.is_open()) {
		std::cerr << "无法打开文件 extrinsic_parameters.txt 进行写入。" << std::endl;
		return;
	}


	//获取图像
	FileManipulation ImageGet;
	std::vector<cv::Mat> BoardPics = ImageGet.ReadPics(M_PlaneBoardFilePath);
	//根据标定板类型，进行分类处理
	if (M_CalibPlateType == 1 || M_CalibPlateType == 4) {
		for (auto& img : BoardPics) {
			vector<Point2f> imgPoints;

			cv::SimpleBlobDetector::Params params;
			cv::Ptr<cv::FeatureDetector> blob = cv::SimpleBlobDetector::create(params);
			params.maxArea = 1000;
			params.minArea = 100;
			params.minDistBetweenBlobs = 10;

			bool found = cv::findCirclesGrid(img, M_board_size, imgPoints, cv::CALIB_CB_SYMMETRIC_GRID, blob);
			if (show) {
				cv::Mat imgCopy;
				if (img.channels() == 1) {
					cv::cvtColor(img, imgCopy, cv::COLOR_GRAY2BGR);  // 将灰度图转为彩色图
				}
				else {
					imgCopy = img.clone();
				}
				cv::drawChessboardCorners(imgCopy, M_board_size, imgPoints, found);
				cv::namedWindow("corner", cv::WINDOW_NORMAL);
				cv::imshow("corner", imgCopy);
				cv::waitKey(500);
			}

			// 求解 R/T
			Mat rvec, tvec;
			solvePnP(objPoints, imgPoints, cameraMatrix, distCoeffs, rvec, tvec);

			// 写入 TXT 文件
			outFile << "Image " << ":" << std::endl;
			outFile << "Rotation Vector (Rodrigues):" << std::endl;
			for (int j = 0; j < 3; ++j) {
				outFile << rvec.at<double>(j, 0) << " ";
			}
			outFile << std::endl;

			outFile << "Translation Vector:" << std::endl;
			for (int j = 0; j < 3; ++j) {
				outFile << tvec.at<double>(j, 0) << " ";
			}
			outFile << std::endl << std::endl;
		}
	}
	else {
		//提取角点
		for (auto& img : BoardPics) {
			vector<Point2f> imgPoints;
			bool found = findChessboardCornersSB(img, M_board_size, imgPoints);
			// 精细化角点
			Mat gray;
			//cvtColor(img, gray, COLOR_BGR2GRAY);
			//cornerSubPix(gray, imgPoints, Size(11, 11), Size(-1, -1),
			//	         TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 30, 0.01));

			if (show) {
				cv::Mat imgCopy;
				if (img.channels() == 1) {
					cv::cvtColor(img, imgCopy, cv::COLOR_GRAY2BGR);  // 将灰度图转为彩色图
				}
				else {
					imgCopy = img.clone();
				}
				cv::drawChessboardCorners(imgCopy, M_board_size, imgPoints, found);
				cv::namedWindow("corner", cv::WINDOW_NORMAL);
				cv::imshow("corner", imgCopy);
				cv::waitKey(500);
			}

			// 求解 R/T
			Mat rvec, tvec;
			solvePnP(objPoints, imgPoints, cameraMatrix, distCoeffs, rvec, tvec);

			rotationMat.push_back(rvec);
			translationMat.push_back(tvec);

			// 写入 TXT 文件
			outFile << "Image " << ":" << std::endl;
			outFile << "Rotation Vector (Rodrigues):" << std::endl;
			for (int j = 0; j < 3; ++j) {
				outFile << rvec.at<double>(j, 0) << " ";
			}
			outFile << std::endl;

			outFile << "Translation Vector:" << std::endl;
			for (int j = 0; j < 3; ++j) {
				outFile << tvec.at<double>(j, 0) << " ";
			}
			outFile << std::endl << std::endl;

		}
	}
	outFile.close(); // 关闭文件
	destroyAllWindows();
	cout << "所有图像的外参提取完成，共处理 " << rotationMat.size() << " 张图像。" << endl;
}



/**
* @brief 通过读取out_camera_data.yml获取内参和外参
* @param instrinsic_Mat	存放相机内参的矩阵
* @param distCoeff		存放相机畸变的矩阵
* @param translationMat	存放外参平移的容器
* @param rotationMat	存放外参旋转的容器
*
*/
void OpticalPlaneCalibration::Get_Camera_Data(const std::string& cameraParamsPath, 
											  Mat& instrinsic_Mat, 
											  Mat& distCoeff) {
	//读取相机标定的结果
	FileStorage camPara(cameraParamsPath, FileStorage::READ);
	if (!camPara.isOpened())
	{
		cout << "--------------相机参数打开失败-----------" << endl;
		exit(1);
	}
	camPara["camera_matrix"] >> instrinsic_Mat;
	camPara["distortion_coefficients"] >> distCoeff;
	camPara.release();
}


/**
* @brief 通过最小二乘法，通过像素坐标系下激光平面上的点来拟合激光平面方程
* @param input	所有拟合平面的点
*
* @return 平面方程Ax+By+Cz+D=0
*/
cv::Mat plane_fitting(const vector<Point3f>& input) {
	cv::Mat dst = Mat(3, 3, CV_32F, Scalar(0));  // 3x3 矩阵
	cv::Mat out = Mat(3, 1, CV_32F, Scalar(0));

	// 计算 dst 和 out 矩阵
	for (int i = 0; i < input.size(); i++) {
		dst.at<float>(0, 0) += pow(input[i].x, 2);
		dst.at<float>(0, 1) += input[i].x * input[i].y;
		dst.at<float>(0, 2) += input[i].x;
		dst.at<float>(1, 0) += input[i].x * input[i].y;
		dst.at<float>(1, 1) += pow(input[i].y, 2);
		dst.at<float>(1, 2) += input[i].y;
		dst.at<float>(2, 0) += input[i].x;
		dst.at<float>(2, 1) += input[i].y;
		dst.at<float>(2, 2) = input.size();

		out.at<float>(0, 0) += input[i].x * input[i].z;
		out.at<float>(1, 0) += input[i].y * input[i].z;
		out.at<float>(2, 0) += input[i].z;
	}

	// 判断矩阵是否为奇异
	if (abs(determinant(dst)) < 0.001) {
		cout << "矩阵奇异，无法解算" << endl;
		return Mat(); // 返回空矩阵
	}

	// 通过 solve 函数求解，避免求逆
	Mat output;
	solve(dst, out, output);

	// 提取平面方程系数
	double a = output.at<float>(0, 0);
	double b = output.at<float>(1, 0);
	double c = -1.0;  // 强制 c = -1
	double d = output.at<float>(2, 0); // d 值表示常数项

	cout << "拟合平面方程: " << a << " * x + " << b << " * y + " << c << " * z + " << d << " = 0" << endl;
	return output;
}


// 使用拟合平面的方程生成一系列平面点
vector<Point3f> generatePlanePoints(double a, double b, double c, double d, float range_min, float range_max, double step) {
	vector<Point3f> planePoints;
	for (float x = range_min; x <= range_max; x += step) {
		for (float y = range_min; y <= range_max; y += step) {
			float z = a * x + b * y + d;  // 计算z坐标
			planePoints.push_back(Point3f(x, y, z));
		}
	}
	return planePoints;
}


/**
* @brief 计算平面拟合残差
* @param threeCoordinate	所有拟合平面的点
* @param a					平面方程Ax+By+Cz+D=0 中的A
* @param b					平面方程Ax+By+Cz+D=0 中的B
* @param c					平面方程Ax+By+Cz+D=0 中的C
* @param d					平面方程Ax+By+Cz+D=0 中的D
*
* @return 拟合残差
*/
double Residual(const vector<Point3f>& threeCoordinate, double a, double b, double c, double d) {
	double distance = 0.0;
	int num = threeCoordinate.size();

	for (const auto& point : threeCoordinate) {
		double x = point.x, y = point.y, z = point.z;
		distance += (a * x + b * y + c * z + d) * (a * x + b * y + c * z + d) / (a * a + b * b + c * c);
	}
	cout << "系数： " << a << "\t" << b << "\t" << c << "\t" << d << endl;
	double derror = sqrt(distance / num);  // 计算均方根误差
	cout << "RMS：" << derror << endl;
	return derror;
}




void OptPlaneCalibration(const std::string& PlaneBoardFilename, const std::string& PlaneLineFilename,
						 const std::string& camParamPath, const std::string& PointsPath,
						 const std::string& PlanePath, enum CalibrationPlate BoardClass) {

	//-----------------------------------获取内参和外参-----------------------------------
	//获取基础参数信息
	Mat instrinsic_Mat, distCoeff;
	vector<Mat>translationMat, rotationMat;
	OpticalPlaneCalibration getBaseParam(PlaneBoardFilename, PlaneLineFilename, BoardClass);
	getBaseParam.Get_Camera_Data(camParamPath, instrinsic_Mat, distCoeff);
	getBaseParam.getRT(instrinsic_Mat, distCoeff, rotationMat, translationMat);


	//------------------------------------------------------------------------------------
	vector<Point3f> laserpoint_cam_all;//存放所有中心点的集合
	//读取图片
	FileManipulation ImageGet;
	std::vector<cv::Mat>BoardPics;
	std::vector<cv::Mat>OptPlanePics;
	ImageGet.ReadPics(PlaneBoardFilename, BoardPics);
	ImageGet.ReadPics(PlaneLineFilename, OptPlanePics);


	for (size_t i = 0; i < BoardPics.size(); i++){
		cv::Mat AfterDis, Board_AfterDis;
		//-----------------------------------获取激光中心线-----------------------------------
		cv::Mat imgHSVMask;
		std::vector<cv::Rect> rois;

		Centerline centerline;
		centerline.picture_blur(OptPlanePics[i], imgHSVMask, rois);
		std::vector<cv::Point2d> center_points = centerline.ExtractCenters(OptPlanePics[i], imgHSVMask, rois);


		//------------------------------------转相机坐标系------------------------------------
		Turn_to_3D turn3d;
		turn3d.Set_param(instrinsic_Mat, distCoeff, translationMat[i], rotationMat[i], BoardPics[i], center_points, BoardClass);
		turn3d.Get_plane();	//获取标定板平面方程
		turn3d.Get_Intersection();	//获取激光线与标定板交点
		vector<Point3d> centers_cam = turn3d.LaserPoints_cam();	//交点坐标由像素坐标系转换到相机坐标系

		//将读到的交点存放至laserpoint_cam_all
		laserpoint_cam_all.insert(laserpoint_cam_all.end(), centers_cam.begin(), centers_cam.end());
		cout << "第" << to_string(i + 1) << "张保存成功" << endl;
	}
	//------------------------------------保存点------------------------------------
	ofstream ofs;
	ofs.open(PointsPath, ios::trunc);
	for (int i = 0; i < laserpoint_cam_all.size(); i++){
		ofs << laserpoint_cam_all[i].x << "\t" << laserpoint_cam_all[i].y << "\t" << laserpoint_cam_all[i].z << endl;
	}
	//------------------------------------------------------------------------------

	//-----------------------------------拟合平面-----------------------------------
	cv::Mat plane0102 = plane_fitting(laserpoint_cam_all);
	// 平面方程 ax + by + cz + d = 0

	// 提取拟合平面方程的系数
	double a = plane0102.at<float>(0, 0);
	double b = plane0102.at<float>(1, 0);
	double c = -1.0;
	double d = plane0102.at<float>(2, 0);

	// 计算均方误差
	double derror = Residual(laserpoint_cam_all, a, b, c, d);

	// 生成拟合平面上的点
	vector<Point3f> planePoints = generatePlanePoints(a, b, c, d, -50, 50, 0.1);
	ofstream ofs_plane(PlanePath, ios::trunc);
	for (int i = 0; i < planePoints.size(); i++) {
		ofs_plane << planePoints[i].x << "\t" << planePoints[i].y << "\t" << planePoints[i].z << endl;
	}
	ofs_plane.close();
}
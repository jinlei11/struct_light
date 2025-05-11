#include"Calibration.h"


CameraCalibration::CameraCalibration(enum CalibrationPlate Plate, 
									 const std::string& filepath_left,
									 const std::string& filepath_Right,
									 const std::string& ResultPath,
									 bool isDrawResult = true)
	:M_CalibrationPicsPath_Left(filepath_left),M_CalibrationPicsPath_Right(filepath_Right),
	M_CalibrationResultPath(ResultPath),M_IsDrawResult(isDrawResult)
{
	switch (Plate)
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

	default:
		break;
	}
};

CameraCalibration::CameraCalibration(const cv::Size& board_size,
									 const cv::Size& squareSize, 
									 const std::string& filepath_left,
									 const std::string& filepath_Right,
									 const std::string& ResultPath,
									 int CalibrationPlateMode,
									 bool isDrawResult = true)
	:M_board_size(board_size), M_squareSize(squareSize),M_CalibrationPicsPath_Left(filepath_left),
	M_CalibrationPicsPath_Right(filepath_Right),M_CalibrationResultPath(ResultPath),
	M_IsDrawResult(isDrawResult), M_CalibPlateType(CalibrationPlateMode)
{

};




void CameraCalibration::ProcessImagePointsCircle(const std::vector<cv::Mat>& imgs,
										   std::vector<std::vector<cv::Point2f>>& imgPoints,
										   const cv::Size& boardSize,
										   const cv::Ptr<cv::FeatureDetector>& blob,
										   bool showCorners){
	for (const auto& img : imgs) {
		std::vector<cv::Point2f> points;
		bool found = cv::findCirclesGrid(img, boardSize, points, cv::CALIB_CB_SYMMETRIC_GRID, blob);
		if (showCorners) {
			cv::Mat imgCopy = img.clone();
			cv::drawChessboardCorners(imgCopy, boardSize, points, found);
			cv::namedWindow("corner", cv::WINDOW_NORMAL);
			cv::imshow("corner", imgCopy);
			cv::waitKey(500);
		}
		std::cout << "[Circle Grid] Found " << points.size() << " points." << std::endl;
		imgPoints.push_back(points);
	}
}

void CameraCalibration::GetImagePointsCircle(bool IsDrawResult) {
	// 1. 设置 blob 参数（需在创建 detector 前设置）
	cv::SimpleBlobDetector::Params params;
	cv::Ptr<cv::FeatureDetector> blob = cv::SimpleBlobDetector::create(params);
	params.maxArea = 1000;
	params.minArea = 100;
	params.minDistBetweenBlobs = 10;

	// 3. 提取左右图像的圆心点
	ProcessImagePointsCircle(M_imgLs, M_imgsPointsL, M_board_size, blob, IsDrawResult);
	ProcessImagePointsCircle(M_imgRs, M_imgsPointsR, M_board_size, blob, IsDrawResult);
};



void CameraCalibration::ProcessImagePointsChess(const std::vector<cv::Mat>& imgs,
												std::vector<std::vector<cv::Point2f>>& imgPoints,
												const cv::Size& boardSize,
												bool showCorners) {
	for (const auto& img : imgs) {
		std::vector<cv::Point2f> points;
		bool found = cv::findChessboardCornersSB(img, boardSize, points);
		if (showCorners) {
			cv::Mat imgCopy = img.clone();
			cv::drawChessboardCorners(imgCopy, boardSize, points, found);
			cv::namedWindow("corner", cv::WINDOW_NORMAL);
			cv::imshow("corner", imgCopy);
			cv::waitKey(500);
		}
		std::cout << "[Chess Grid] Found " << points.size() << " points." << std::endl;
		imgPoints.push_back(points);
	}
}

void CameraCalibration::GetImagePointsChess(bool IsDrawResult) {
	ProcessImagePointsChess(M_imgLs, M_imgsPointsL, M_board_size, IsDrawResult);
	ProcessImagePointsChess(M_imgRs, M_imgsPointsR, M_board_size, IsDrawResult);
};


void CameraCalibration::GetObjectPoints() {
	// 创建一张标定板的三维坐标点（z=0）
	std::vector<cv::Point3f> objectTemplate;
	for (int i = 0; i < M_board_size.height; ++i) {
		for (int j = 0; j < M_board_size.width; ++j) {
			objectTemplate.emplace_back(
				j * M_squareSize.width,
				i * M_squareSize.height,
				0.0f
			);
		}
	}
	// 将模板复制到所有图像帧对应的位置
	M_objectPoints.assign(M_imgsPointsL.size(), objectTemplate);
}



cv::Mat CameraCalibration::MyCameraCalibrationAndSave() {

	//Initialize image size
	cv::Size Imagesize;
	Imagesize.width = M_imgLs[0].cols;
	Imagesize.height = M_imgLs[0].rows;

	//Define Mat matrix for storing internal reference matrix and distortion coefficients
	cv::Mat MatrixL = cv::Mat(3, 3, CV_32FC1, cv::Scalar::all(0.f));
	cv::Mat MatrixR = cv::Mat(3, 3, CV_32FC1, cv::Scalar::all(0.f));
	cv::Mat distL = cv::Mat(1, 5, CV_32FC1, cv::Scalar::all(0.f));
	cv::Mat distR = cv::Mat(1, 5, CV_32FC1, cv::Scalar::all(0.f));

	//Define matrices to hold rotation transformation matrices and translation transformation matrices
	std::vector<cv::Mat>rvecsL;
	std::vector<cv::Mat>tvecsL;
	std::vector<cv::Mat>rvecsR;
	std::vector<cv::Mat>tvecsR;

	//Camera Calibration，error demonstrates the reprojection error 
	auto errorL = calibrateCamera(M_objectPoints, M_imgsPointsL, Imagesize, MatrixL, distL, rvecsL, tvecsL, 0);
	auto errorR = calibrateCamera(M_objectPoints, M_imgsPointsR, Imagesize, MatrixR, distR, rvecsR, tvecsR, 0);

	std::cout << "calibration_errorL: " << errorL << std::endl;
	std::cout << "calibration_errorR: " << errorR << std::endl;

	////calcualte re-projection error----------------------------------------------------------------------
	double errL = 0.f;
	double err_sumL = 0;
	std::vector<cv::Point2f>pointL;
	std::cout << "每幅图像的标定误差" << std::endl;

	for (int i = 0; i < M_imgLs.size(); i++) {

		cv::projectPoints(M_objectPoints[i], rvecsL[i], tvecsL[i], MatrixL, distL, pointL);

		std::vector<cv::Point2f>tempImagepoint = M_imgsPointsL[i];
		cv::Mat ImagepointsFL = cv::Mat(1, tempImagepoint.size(), CV_32FC2);
		cv::Mat ImagepointsAL = cv::Mat(1, pointL.size(), CV_32FC2);

		for (int j = 0; j < tempImagepoint.size(); j++) {
			ImagepointsFL.at<cv::Vec2f>(0, j) = cv::Vec2f(tempImagepoint[j].x, tempImagepoint[j].y);
			ImagepointsAL.at<cv::Vec2f>(0, j) = cv::Vec2f(pointL[j].x, pointL[j].y);
		}
		errL = cv::norm(ImagepointsFL, ImagepointsAL, cv::NORM_L2);
		//err_sumL += errL /= 88;
		//float result = sqrt(errL * errL / 130 / 130);
		float result = errL / 130;
		std::cout << "左相机第" << i + 1 << "幅图像的平均误差：" << result << std::endl;
	}


	double errR = 0.f;
	double err_sumR = 0;
	std::vector<cv::Point2f>pointR;
	for (int i = 0; i < M_imgRs.size(); i++) {

		cv::projectPoints(M_objectPoints[i], rvecsR[i], tvecsR[i], MatrixR, distR, pointR);

		std::vector<cv::Point2f>tempImagepoint = M_imgsPointsR[i];
		cv::Mat ImagepointsFR = cv::Mat(1, tempImagepoint.size(), CV_32FC2);
		cv::Mat ImagepointsAR = cv::Mat(1, pointR.size(), CV_32FC2);

		for (int j = 0; j < tempImagepoint.size(); j++) {
			ImagepointsFR.at<cv::Vec2f>(0, j) = cv::Vec2f(tempImagepoint[j].x, tempImagepoint[j].y);
			ImagepointsAR.at<cv::Vec2f>(0, j) = cv::Vec2f(pointR[j].x, pointR[j].y);
		}
		errR = cv::norm(ImagepointsFR, ImagepointsAR, cv::NORM_L2);
		//err_sumR += errR /= 88;
		//float result = sqrt(errR * errR / 130 / 130);
		float result = errR / 130;
		std::cout << "右相机第" << i + 1 << "幅图像的平均误差：" << result << std::endl;
	}



	//After calibrating a single camera, the rotation translation transformation matrix between the two cameras is computed;
	//R is the rotation matrix between the cameras
	//T is the translation matrix between cameras
	//E is the eigenmatrix between cameras (rotation + translation)
	//F is the fundamental matrix between cameras Eigenmatrix + internal reference between two cameras;
	cv::Mat R, T, E, F;

	//Perform binocular calibration
	auto error2 = stereoCalibrate(M_objectPoints, M_imgsPointsL, M_imgsPointsR, MatrixL, distL, MatrixR, distR, Imagesize, R, T, E, F);

	//Compute the corrective transformation matrix
	cv::Mat R1, R2, P1, P2, Q;
	//R1 矫正旋转矩阵,将第一个相机坐标系下未矫正的点变换到第一个相机矫正坐标系下，即 R_{左矫正坐标系}{左未矫正坐标系}。
	//R2 矫正旋转矩阵,将第二个相机坐标系下未矫正的点变换到第二个相机矫正坐标系下，即 R_{右矫正坐标系}{右未矫正坐标系}。
	//P1 3x4左相机投影矩阵。将左矫正坐标系下的点投影到左矫正坐标系图像平面坐标系。
	//P2 3x4右相机投影矩阵。将左矫正坐标系下的点投影到右矫正坐标系图像平面坐标系。
	//Q  4x4的视差深度映射矩阵。

	cv::stereoRectify(MatrixL, distL, MatrixR, distR, Imagesize, R, T, R1, R2, P1, P2, Q, 0);

	//cv::Mat M_map11, M_map12, M_map21, M_map22;

	//initUndistortRectifyMap(MatrixL, distL, R1, P1, Imagesize, CV_32FC1, M_map11, M_map12);
	//initUndistortRectifyMap(MatrixR, distR, R2, P2, Imagesize, CV_32FC1, M_map21, M_map22);

	//cv::Mat rectifyImg(M_imgLs[0].rows, M_imgLs[0].cols * 2, CV_8UC1, cv::Scalar(0.f));
	//remap(M_imgLs[0], rectifyImg(cv::Rect(0, 0 ,rectifyImg.cols / 2, rectifyImg.rows)), M_map11, M_map12, cv::INTER_LINEAR);
	//remap(M_imgRs[0], rectifyImg(cv::Rect(rectifyImg.cols / 2, 0, rectifyImg.cols / 2, rectifyImg.rows)), M_map21, M_map22, cv::INTER_LINEAR);
	//
	//cv::Mat temp;
	//cv::cvtColor(rectifyImg, temp, cv::COLOR_GRAY2BGR);
	//for (size_t i = 0; i < rectifyImg.rows; i += 100) {
	//	cv::line(temp, cv::Point(0, i), cv::Point(rectifyImg.cols - 1, i), cv::Scalar(255, 0, 0));
	//}

	//Save calibration parameters
	cv::FileStorage FStore(M_CalibrationResultPath, cv::FileStorage::WRITE);
	FStore << "MatrixL" << MatrixL;
	FStore << "distL" << distL;
	FStore << "MatrixR" << MatrixR;
	FStore << "distR" << distR;
	FStore << "Q" << Q;
	FStore << "R1" << R1;
	FStore << "R2" << R2;
	FStore << "P1" << P1;
	FStore << "P2" << P2;
	FStore.release();
	std::cout << "保存完毕" << std::endl;




	//map11 输出的X坐标重映射参数
	//map12 输出的Y坐标重映射参数
	//map21 输出的X坐标重映射参数
	//map22 输出的Y坐标重映射参数
	//initUndistortRectifyMap(MatrixL, distL, R1, P1, Imagesize, CV_16SC2, M_map11, M_map12);	
	//initUndistortRectifyMap(MatrixR, distR, R2, P2, Imagesize, CV_16SC2, M_map21, M_map22);



	////Save calibration results  TXT
	//std::ofstream Fout("C:\\Users\\86187\\Pictures\\Calibration\\Calibration.txt");
	//Fout << "左相机内参矩阵：" << endl;
	//Fout << MatrixL << endl << endl;
	//Fout << "左相机畸变系数矩阵：" << endl;
	//Fout << distL << endl << endl;
	//Fout << "右相机内参矩阵：" << endl;
	//Fout << MatrixR << endl << endl;
	//Fout << "右相机畸变系数矩阵：" << endl;
	//Fout << distR << endl << endl;
	//Fout << "相机之间的旋转矩阵R：" << endl;
	//Fout << R << endl << endl;
	//Fout << "相机之间的平移矩阵T:" << endl;
	//Fout << T << endl;
	//Fout << "相机之间的本征矩阵（旋转+平移）E:" << endl;
	//Fout << E << endl << endl;
	//Fout << "相机之间的基本矩阵  本征矩阵+两个相机之间的内参F:" << endl;
	//Fout << F << endl;
	//Fout << "进行校正变换之后-------------------------------------------------------" << endl;
	//Fout << "4x4的视差深度映射矩阵:" << endl;
	//Fout << Q << endl << endl;
	//Fout << "矫正旋转矩阵,将第一个相机坐标系下未矫正的点变换到第一个相机矫正坐标系下，即 R_{左矫正坐标系}{左未矫正坐标系}:R1" << endl;
	//Fout << R1 << endl << endl;
	//Fout << "矫正旋转矩阵,将第二个相机坐标系下未矫正的点变换到第二个相机矫正坐标系下，即 R_{右矫正坐标系}{右未矫正坐标系}:R2" << endl;
	//Fout << R2 << endl << endl;
	//Fout << "3x4左相机投影矩阵。将左矫正坐标系下的点投影到左矫正坐标系图像平面坐标系:P1" << endl;
	//Fout << P1 << endl << endl;
	//Fout << "3x4右相机投影矩阵。将左矫正坐标系下的点投影到右矫正坐标系图像平面坐标系:P2" << endl;
	//Fout << P2 << endl << endl;



	//进行图像的校正
	//for (int i = 0; i < M_imgLs.size(); i++)
	//{
	//	//进行校正映射
	//	cv::Mat img1r, img2r;
	//	remap(M_imgLs[i], img1r, M_map11, M_map12, cv::INTER_LINEAR);
	//	remap(M_imgRs[i], img2r, M_map21, M_map22, cv::INTER_LINEAR);
	//	//用于显示校正结果
	//	//拼接图像
	//	cv::Mat result;
	//	hconcat(img1r, img2r, result);
	//	//绘制直线，用于比较同一个内角点y轴是否一致
	//	cv::namedWindow("result", cv::WINDOW_NORMAL);
	//	line(result, cv::Point(-1, M_imgsPointsL[i][0].y), cv::Point(result.cols, M_imgsPointsL[i][0].y), cv::Scalar(0, 0, 255), 2);
	//	imshow("result", result);
	//	cv::waitKey(0);
	//}
	return Q;
};


cv::Mat CameraCalibration::M_MyCameraCalibrationAndSave() {
	//读取图像
	M_GetImage.ReadPics(M_CalibrationPicsPath_Left, M_imgLs);
	M_GetImage.ReadPics(M_CalibrationPicsPath_Right, M_imgRs);
	//计算图像坐标
	if (M_CalibPlateType == 1 || M_CalibPlateType == 2 || M_CalibPlateType == 4) {
		GetImagePointsCircle(true);
	}
	else {
		GetImagePointsChess(true);
	}
	//计算世界坐标
	GetObjectPoints();

	//标定
	cv::Size imageSize = M_imgLs[0].size();

	cv::Mat MatrixL, MatrixR, distL, distR;
	std::vector<cv::Mat> rvecsL, tvecsL, rvecsR, tvecsR;

	double errorL = calibrateCamera(M_objectPoints, M_imgsPointsL, imageSize, MatrixL, distL, rvecsL, tvecsL);
	double errorR = calibrateCamera(M_objectPoints, M_imgsPointsR, imageSize, MatrixR, distR, rvecsR, tvecsR);

	std::cout << "calibration_errorL: " << errorL << std::endl;
	std::cout << "calibration_errorR: " << errorR << std::endl;

	// 计算重投影误差的 Lambda 简化逻辑
	auto computeReprojectionError = [](const std::vector<std::vector<cv::Point3f>>& objectPoints,
		const std::vector<std::vector<cv::Point2f>>& imagePoints,
		const std::vector<cv::Mat>& rvecs,
		const std::vector<cv::Mat>& tvecs,
		const cv::Mat& cameraMatrix,
		const cv::Mat& distCoeffs,
		const std::string& cameraName) {

			for (size_t i = 0; i < objectPoints.size(); ++i) {
				std::vector<cv::Point2f> projectedPoints;
				cv::projectPoints(objectPoints[i], rvecs[i], tvecs[i], cameraMatrix, distCoeffs, projectedPoints);

				double err = cv::norm(imagePoints[i], projectedPoints, cv::NORM_L2);
				double avgErr = err / imagePoints[i].size();

				std::cout << cameraName << "第" << i + 1 << "幅图像的平均误差：" << avgErr << std::endl;
			}
	};

	std::cout << "每幅图像的标定误差" << std::endl;
	computeReprojectionError(M_objectPoints, M_imgsPointsL, rvecsL, tvecsL, MatrixL, distL, "左相机");
	computeReprojectionError(M_objectPoints, M_imgsPointsR, rvecsR, tvecsR, MatrixR, distR, "右相机");

	// 双目标定
	cv::Mat R, T, E, F;
	stereoCalibrate(M_objectPoints, M_imgsPointsL, M_imgsPointsR, MatrixL, distL, MatrixR, distR, imageSize, R, T, E, F);

	// 立体矫正
	cv::Mat R1, R2, P1, P2, Q;
	stereoRectify(MatrixL, distL, MatrixR, distR, imageSize, R, T, R1, R2, P1, P2, Q);

	// 保存校准参数
	cv::FileStorage fs(M_CalibrationResultPath, cv::FileStorage::WRITE);
	fs << "MatrixL" << MatrixL << "distL" << distL
		<< "MatrixR" << MatrixR << "distR" << distR
		<< "Q" << Q << "R1" << R1 << "R2" << R2
		<< "P1" << P1 << "P2" << P2;
	fs.release();

	std::cout << "保存完毕" << std::endl;
	return Q;
}

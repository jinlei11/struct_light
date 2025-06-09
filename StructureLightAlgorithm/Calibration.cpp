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
		 
	case Plate_chess811:
		this->M_board_size = cv::Size(8, 11);
		this->M_squareSize = cv::Size(3, 3);
		this->M_CalibPlateType = 5;

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
			cv::Mat imgCopy;
			if (img.channels() == 1) {
				cv::cvtColor(img, imgCopy, cv::COLOR_GRAY2BGR);  // 将灰度图转为彩色图
			}
			else {
				imgCopy = img.clone();
			}
			cv::drawChessboardCorners(imgCopy, boardSize, points, found);
			cv::namedWindow("corner", cv::WINDOW_NORMAL);
			cv::imshow("corner", imgCopy);
			cv::waitKey(500);
		}
		std::cout << "[Circle Grid] Found " << points.size() << " points." << std::endl;
		imgPoints.push_back(points);
	}
	cv::destroyAllWindows();
}


void CameraCalibration::ProcessImagePointsCircleInvert(const std::vector<cv::Mat>& imgs,
														std::vector<std::vector<cv::Point2f>>& imgPoints,
														const cv::Size& boardSize,
														const cv::Ptr<cv::FeatureDetector>& blob,
														bool showCorners) {
	for (const auto& img : imgs) {
		std::vector<cv::Point2f> points;
		bool found = cv::findCirclesGrid(~img, boardSize, points, cv::CALIB_CB_SYMMETRIC_GRID, blob);
		if (showCorners) {
			cv::Mat imgCopy;
			if (img.channels() == 1) {
				cv::cvtColor(img, imgCopy, cv::COLOR_GRAY2BGR);  // 将灰度图转为彩色图
			}
			else {
				imgCopy = img.clone();
			}
			cv::drawChessboardCorners(imgCopy, boardSize, points, found);
			cv::namedWindow("corner", cv::WINDOW_NORMAL);
			cv::imshow("corner", imgCopy);
			cv::waitKey(500);
		}
		std::cout << "[Circle Grid] Found " << points.size() << " points." << std::endl;
		imgPoints.push_back(points);
	}
	cv::destroyAllWindows();
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



void CameraCalibration::GetImagePointsCircleInvert(bool IsDrawResult) {
	// 1. 设置 blob 参数（需在创建 detector 前设置）
	cv::SimpleBlobDetector::Params params;
	cv::Ptr<cv::FeatureDetector> blob = cv::SimpleBlobDetector::create(params);
	params.maxArea = 1000;
	params.minArea = 100;
	params.minDistBetweenBlobs = 10;

	// 3. 提取左右图像的圆心点
	ProcessImagePointsCircleInvert(M_imgLs, M_imgsPointsL, M_board_size, blob, IsDrawResult);
	ProcessImagePointsCircleInvert(M_imgRs, M_imgsPointsR, M_board_size, blob, IsDrawResult);
};



void CameraCalibration::ProcessImagePointsChess(const std::vector<cv::Mat>& imgs,
												std::vector<std::vector<cv::Point2f>>& imgPoints,
												const cv::Size& boardSize,
												bool showCorners) {
	for (const auto& img : imgs) {
		std::vector<cv::Point2f> points;
		bool found = cv::findChessboardCornersSB(img, boardSize, points);
		if (showCorners) {
			cv::Mat imgCopy;
			if (img.channels() == 1) {
				cv::cvtColor(img, imgCopy, cv::COLOR_GRAY2BGR);  // 将灰度图转为彩色图
			}
			else {
				imgCopy = img.clone();
			}
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




////修改板
//cv::Mat CameraCalibration::M_MyCameraCalibrationAndSave() {
//	//读取图像
//	M_GetImage.ReadPics(M_CalibrationPicsPath_Left, M_imgLs);
//	M_GetImage.ReadPics(M_CalibrationPicsPath_Right, M_imgRs);
//
//	// 检查第一张图像的属性
//	std::cout << "Left image size: " << M_imgLs[0].size() << ", channels: " << M_imgLs[0].channels() << ", type: " << M_imgLs[0].type() << std::endl;
//	std::cout << "Right image size: " << M_imgRs[0].size() << ", channels: " << M_imgRs[0].channels() << ", type: " << M_imgRs[0].type() << std::endl;
//
//	//计算图像坐标
//	if (M_CalibPlateType == 1 || M_CalibPlateType == 4) {
//		GetImagePointsCircle(true);
//	}
//	else if (M_CalibPlateType == 2) {
//		GetImagePointsCircleInvert(true);
//	}
//	else {
//		GetImagePointsChess(true);
//	}
//
//	// 检查特征点数量
//	if (M_imgsPointsL.size() != M_imgsPointsR.size() || M_imgsPointsL.empty()) {
//		std::cerr << "Error: Invalid image points detected!" << std::endl;
//		return cv::Mat();
//	}
//
//	//计算世界坐标
//	GetObjectPoints();
//
//	//标定
//	cv::Size imageSize = M_imgLs[0].size();
//
//	cv::Mat MatrixL, MatrixR, distL, distR;
//	std::vector<cv::Mat> rvecsL, tvecsL, rvecsR, tvecsR;
//
//	double errorL = calibrateCamera(M_objectPoints, M_imgsPointsL, imageSize, MatrixL, distL, rvecsL, tvecsL,0);
//	double errorR = calibrateCamera(M_objectPoints, M_imgsPointsR, imageSize, MatrixR, distR, rvecsR, tvecsR,0);
//
//	std::cout << "calibration_errorL: " << errorL << std::endl;
//	std::cout << "calibration_errorR: " << errorR << std::endl;
//
//	// 检查标定误差
//	if (errorL > 1.0 || errorR > 1.0) {
//		std::cout << "Warning: High calibration error detected!" << std::endl;
//	}
//
//	// 计算重投影误差的 Lambda 简化逻辑
//	auto computeReprojectionError = [](const std::vector<std::vector<cv::Point3f>>& objectPoints,
//		const std::vector<std::vector<cv::Point2f>>& imagePoints,
//		const std::vector<cv::Mat>& rvecs,
//		const std::vector<cv::Mat>& tvecs,
//		const cv::Mat& cameraMatrix,
//		const cv::Mat& distCoeffs,
//		const std::string& cameraName) {
//
//			for (size_t i = 0; i < objectPoints.size(); ++i) {
//				std::vector<cv::Point2f> projectedPoints;
//				cv::projectPoints(objectPoints[i], rvecs[i], tvecs[i], cameraMatrix, distCoeffs, projectedPoints);
//
//				double err = cv::norm(imagePoints[i], projectedPoints, cv::NORM_L2);
//				double avgErr = err / imagePoints[i].size();
//
//				std::cout << cameraName << "第" << i + 1 << "幅图像的平均误差：" << avgErr << std::endl;
//			}
//	};
//
//	std::cout << "每幅图像的标定误差" << std::endl;
//	computeReprojectionError(M_objectPoints, M_imgsPointsL, rvecsL, tvecsL, MatrixL, distL, "左相机");
//	computeReprojectionError(M_objectPoints, M_imgsPointsR, rvecsR, tvecsR, MatrixR, distR, "右相机");
//
//	// 双目标定 - 使用更合适的标志
//	cv::Mat R, T, E, F;
//	//int stereoFlags = cv::CALIB_FIX_INTRINSIC | cv::CALIB_USE_INTRINSIC_GUESS | cv::CALIB_SAME_FOCAL_LENGTH;
//
//	double stereo_error = stereoCalibrate(M_objectPoints, M_imgsPointsL, M_imgsPointsR,
//											MatrixL, distL, MatrixR, distR, imageSize,
//											R, T, E, F);
//	std::cout << "Stereo calibration error: " << stereo_error << std::endl;
//
//	// 检查双目标定结果
//	std::cout << "Translation vector T: " << T << std::endl;
//	std::cout << "Rotation matrix R: " << R << std::endl;
//
//	// 立体矫正 - 使用alpha参数控制视野
//	cv::Mat R1, R2, P1, P2, Q;
//	double alpha = 0.0; // 0表示裁剪掉无效像素，1表示保留所有像素
//	cv::stereoRectify(MatrixL, distL, MatrixR, distR, imageSize, R, T,
//					  R1, R2, P1, P2, Q);
//
//
//	// 保存校准参数
//	cv::FileStorage fs(M_CalibrationResultPath, cv::FileStorage::WRITE);
//	fs << "MatrixL" << MatrixL << "distL" << distL
//		<< "MatrixR" << MatrixR << "distR" << distR
//		<< "Q" << Q << "R1" << R1 << "R2" << R2
//		<< "P1" << P1 << "P2" << P2
//		<< "R" << R << "T" << T << "E" << E << "F" << F;
//	fs.release();
//
//	// 计算重映射参数 - 使用CV_32FC1类型更精确
//	cv::Mat M_map11, M_map12, M_map21, M_map22;
//	initUndistortRectifyMap(MatrixL, distL, R1, P1, imageSize, CV_32FC1, M_map11, M_map12);
//	initUndistortRectifyMap(MatrixR, distR, R2, P2, imageSize, CV_32FC1, M_map21, M_map22);
//
//	// 检查重映射参数
//	std::cout << "Map sizes - Left: " << M_map11.size() << ", Right: " << M_map21.size() << std::endl;
//
//	for (int i = 0; i < M_imgLs.size(); i++) {
//		std::cout << "Processing image pair " << i + 1 << "/" << M_imgLs.size() << std::endl;
//
//		// 检查原始图像
//		if (M_imgLs[i].empty() || M_imgRs[i].empty()) {
//			std::cerr << "Warning: Empty image at index " << i << std::endl;
//			continue;
//		}
//
//		// 确保图像格式一致
//		cv::Mat leftImg = M_imgLs[i];
//		cv::Mat rightImg = M_imgRs[i];
//
//		// 如果是彩色图像，确保都是3通道
//		if (leftImg.channels() == 1 && rightImg.channels() == 3) {
//			cv::cvtColor(leftImg, leftImg, cv::COLOR_GRAY2BGR);
//		}
//		else if (leftImg.channels() == 3 && rightImg.channels() == 1) {
//			cv::cvtColor(rightImg, rightImg, cv::COLOR_GRAY2BGR);
//		}
//
//		//进行校正映射
//		cv::Mat imgL, imgR;
//		remap(leftImg, imgL, M_map11, M_map12, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0));
//		remap(rightImg, imgR, M_map21, M_map22, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0));
//
//		// 检查矫正后的图像
//		if (imgL.empty() || imgR.empty()) {
//			std::cerr << "Error: Failed to rectify image pair " << i << std::endl;
//			continue;
//		}
//
//		std::cout << "Rectified left image size: " << imgL.size() << ", type: " << imgL.type() << std::endl;
//		std::cout << "Rectified right image size: " << imgR.size() << ", type: " << imgR.type() << std::endl;
//
//		// 确保两个图像具有相同的尺寸和类型
//		if (imgL.size() != imgR.size() || imgL.type() != imgR.type()) {
//			std::cerr << "Error: Left and right images have different sizes or types after rectification!" << std::endl;
//			std::cerr << "Left: " << imgL.size() << ", type: " << imgL.type() << std::endl;
//			std::cerr << "Right: " << imgR.size() << ", type: " << imgR.type() << std::endl;
//			continue;
//		}
//
//		// 检查图像内容是否有效
//		cv::Scalar meanL = cv::mean(imgL);
//		cv::Scalar meanR = cv::mean(imgR);
//		std::cout << "Mean pixel values - Left: " << meanL << ", Right: " << meanR << std::endl;
//
//		// 拼接图像
//		cv::Mat result;
//		try {
//			hconcat(imgL, imgR, result);
//
//			// 检查拼接结果
//			if (result.empty()) {
//				std::cerr << "Error: Failed to concatenate images" << std::endl;
//				continue;
//			}
//
//			std::cout << "Result image size: " << result.size() << ", type: " << result.type() << std::endl;
//
//			// 绘制极线 - 使用多条线来验证矫正效果
//			int step = result.rows / 20; // 绘制20条极线
//			for (int y = step; y < result.rows; y += step) {
//				line(result, cv::Point(0, y), cv::Point(result.cols, y), cv::Scalar(0, 255, 0), 1);
//			}
//
//			// 在图像中央绘制一条红色的极线用于重点参考
//			int centerY = result.rows / 2;
//			line(result, cv::Point(0, centerY), cv::Point(result.cols, centerY), cv::Scalar(0, 0, 255), 2);
//
//			// 调整显示窗口大小
//			cv::namedWindow("Stereo Rectification Result", cv::WINDOW_NORMAL);
//			cv::resizeWindow("Stereo Rectification Result", 1200, 400);
//			imshow("Stereo Rectification Result", result);
//
//			std::cout << "Press any key to continue, 'q' to quit, 's' to save current result..." << std::endl;
//			char key = cv::waitKey(0);
//			if (key == 'q' || key == 'Q') {
//				break;
//			}
//			else if (key == 's' || key == 'S') {
//				std::string filename = "rectified_result_" + std::to_string(i) + ".jpg";
//				cv::imwrite(filename, result);
//				std::cout << "Saved: " << filename << std::endl;
//			}
//
//		}
//		catch (const cv::Exception& e) {
//			std::cerr << "OpenCV Exception: " << e.what() << std::endl;
//		}
//	}
//
//	cv::destroyAllWindows();
//	std::cout << "双目标定完成" << std::endl;
//	return Q;
//}




/* ****************************************   标定函数   ************************************************/
cv::Mat CameraCalibration::M_MyCameraCalibrationAndSave() {
	//读图像
	ReadCalibrationImages();    
	//提取特征点
	if (!ExtractFeaturePoints()) return cv::Mat(); 
	// 相机标定
	CalibrateSingleCameras();         
	// 重投影误差计算
	PrintReprojectionError(M_objectPoints, M_imgsPointsL, M_rvecsL, M_tvecsL, M_MatrixL, M_DistL, "左相机");
	PrintReprojectionError(M_objectPoints, M_imgsPointsR, M_rvecsR, M_tvecsR, M_MatrixR, M_DistR, "右相机");
	//双目标定
	StereoCalibration();        
	//矫正与拼接显示
	StereoRectificationAndRemap();  
	std::cout << "双目标定完成" << std::endl;
	return M_Q;
}



void CameraCalibration::ReadCalibrationImages() {
	M_GetImage.ReadPics(M_CalibrationPicsPath_Left, M_imgLs);
	M_GetImage.ReadPics(M_CalibrationPicsPath_Right, M_imgRs);
	std::cout << "读取图像完成，尺寸检查:" << M_imgLs[0].size() << " | " << M_imgRs[0].size() << std::endl;
}


bool CameraCalibration::ExtractFeaturePoints() {
	if (M_CalibPlateType == 1 || M_CalibPlateType == 4)
		GetImagePointsCircle(true);
	else if (M_CalibPlateType == 2)
		GetImagePointsCircleInvert(true);
	else
		GetImagePointsChess(true);

	if (M_imgsPointsL.size() != M_imgsPointsR.size() || M_imgsPointsL.empty()) {
		std::cerr << "特征点提取失败！" << std::endl;
		return false;
	}

	GetObjectPoints();
	return true;
}


void CameraCalibration::CalibrateSingleCameras() {
	cv::Size imageSize = M_imgLs[0].size();
	auto m_errorL = calibrateCamera(M_objectPoints, M_imgsPointsL, imageSize, M_MatrixL, M_DistL, M_rvecsL, M_tvecsL, 0);
	auto m_errorR = calibrateCamera(M_objectPoints, M_imgsPointsR, imageSize, M_MatrixR, M_DistR, M_rvecsR, M_tvecsR, 0);
	std::cout << "左相机误差：" << m_errorL << " | 右相机误差：" << m_errorR << std::endl;
}


void CameraCalibration::PrintReprojectionError(const std::vector<std::vector<cv::Point3f>>& objPts,
												const std::vector<std::vector<cv::Point2f>>& imgPts,
												const std::vector<cv::Mat>& rvecs,
												const std::vector<cv::Mat>& tvecs,
												const cv::Mat& camMat, const cv::Mat& distCoeffs,
												const std::string& camName) {
	for (size_t i = 0; i < objPts.size(); ++i) {
		std::vector<cv::Point2f> projPts;
		cv::projectPoints(objPts[i], rvecs[i], tvecs[i], camMat, distCoeffs, projPts);
		double avgErr = cv::norm(imgPts[i], projPts, cv::NORM_L2) / imgPts[i].size();
		std::cout << camName << " 第" << i + 1 << "图像误差: " << avgErr << std::endl;
	}
}

void CameraCalibration::StereoCalibration() {
	cv::Size imageSize = M_imgLs[0].size();
	auto stereoError = cv::stereoCalibrate(
		M_objectPoints, M_imgsPointsL, M_imgsPointsR,
		M_MatrixL, M_DistL, M_MatrixR, M_DistR, imageSize,
		M_R, M_T, M_E, M_F);
	std::cout << "Stereo calibration error: " << stereoError << std::endl;

	cv::stereoRectify(M_MatrixL, M_DistL, M_MatrixR, M_DistR, imageSize, M_R, M_T,
		M_R1, M_R2, M_P1, M_P2, M_Q, 0);

	cv::FileStorage fs(M_CalibrationResultPath, cv::FileStorage::WRITE);
	fs << "MatrixL" << M_MatrixL << "distL" << M_DistL
		<< "MatrixR" << M_MatrixR << "distR" << M_DistR
		<< "Q" << M_Q << "R1" << M_R1 << "R2" << M_R2
		<< "P1" << M_P1 << "P2" << M_P2
		<< "R" << M_R << "T" << M_T << "E" << M_E << "F" << M_F;
	fs.release();
}




void CameraCalibration::StereoRectificationAndRemap() {
	cv::Size imageSize = M_imgLs[0].size();
	cv::initUndistortRectifyMap(M_MatrixL, M_DistL, M_R1, M_P1, imageSize, CV_32FC1, M_map11, M_map12);
	cv::initUndistortRectifyMap(M_MatrixR, M_DistR, M_R2, M_P2, imageSize, CV_32FC1, M_map21, M_map22);

	for (size_t i = 0; i < M_imgLs.size(); ++i) {
		cv::Mat imgL_rect, imgR_rect;
		cv::remap(M_imgLs[i], imgL_rect, M_map11, M_map12, cv::INTER_LINEAR);
		cv::remap(M_imgRs[i], imgR_rect, M_map21, M_map22, cv::INTER_LINEAR);

		// 拼接显示
		cv::Mat result;
		hconcat(imgL_rect, imgR_rect, result);
		std::cout << result.rows<< std::endl;
		DrawEpipolarLines(result);

		//// 评估极线误差
		//EpipolarErrorStats stats = evaluateEpipolarError(imgL_rect, imgR_rect);
		//std::cout << "Image Pair " << i << " - Matches: " << stats.num_matches
		//	<< ", Mean: " << stats.mean_error
		//	<< ", StdDev: " << stats.std_error
		//	<< ", Max: " << stats.max_error << std::endl;
	}
}


void CameraCalibration::DrawEpipolarLines(cv::Mat& img) {
	if (img.channels() == 1) {
		cv::cvtColor(img, img, cv::COLOR_GRAY2BGR);
	}
	int h = img.rows;
	int step = h / 20;
	std::cout << "step: " << step << std::endl;
	for (int y = step; y < h; y += step) {
		cv::line(img, { 0, y }, { img.cols, y }, { 0, 255, 0 }, 3); // 整体图画线
	}

	cv::namedWindow("Stereo Rectification Result", cv::WINDOW_NORMAL);
	cv::imshow("Stereo Rectification Result", img);
	char key = cv::waitKey(0);
}



void CameraCalibration::ShowAndSaveResult(const cv::Mat& result, int index) {
	cv::namedWindow("Stereo Rectification Result", cv::WINDOW_NORMAL);
	cv::resizeWindow("Stereo Rectification Result", 1200, 400);
	cv::imshow("Stereo Rectification Result", result);
	char key = cv::waitKey(0);
	if (key == 's' || key == 'S') {
		std::string name = "rectified_result_" + std::to_string(index) + ".jpg";
		cv::imwrite(name, result);
		std::cout << "保存图像：" << name << std::endl;
	}
	if (key == 'q' || key == 'Q') exit(0);
}




/* ****************************************   标定函数   ************************************************/


EpipolarErrorStats CameraCalibration::evaluateEpipolarError(const cv::Mat& img_left_rect, const cv::Mat& img_right_rect) {
	
	cv::Ptr<cv::SIFT> detector = cv::SIFT::create();
	
	std::vector<cv::KeyPoint> keypoints_left, keypoints_right;
	cv::Mat descriptors_left, descriptors_right;

	detector->detectAndCompute(img_left_rect, cv::noArray(), keypoints_left, descriptors_left);
	detector->detectAndCompute(img_right_rect, cv::noArray(), keypoints_right, descriptors_right);

	cv::BFMatcher matcher;
	std::vector<std::vector<cv::DMatch>> knn_matches;
	matcher.knnMatch(descriptors_left, descriptors_right, knn_matches, 2);

	std::vector<cv::DMatch> good_matches;
	const float ratio_thresh = 0.75f;
	for (size_t i = 0; i < knn_matches.size(); i++) {
		if (knn_matches[i].size() < 2) continue;
		if (knn_matches[i][0].distance < ratio_thresh * knn_matches[i][1].distance) {
			good_matches.push_back(knn_matches[i][0]);
		}
	}

	std::vector<double> epipolar_errors;
	for (const auto& match : good_matches) {
		cv::Point2f pt_left = keypoints_left[match.queryIdx].pt;
		cv::Point2f pt_right = keypoints_right[match.trainIdx].pt;
		epipolar_errors.push_back(std::abs(pt_left.y - pt_right.y));
	}

	EpipolarErrorStats stats;
	stats.num_matches = epipolar_errors.size();

	if (epipolar_errors.empty()) {
		stats.mean_error = stats.std_error = stats.max_error = 0.0;
		return stats;
	}

	double sum = 0.0, max_err = 0.0;
	for (double err : epipolar_errors) {
		sum += err;
		max_err = std::max(max_err, err);
	}
	stats.mean_error = sum / epipolar_errors.size();

	double variance = 0.0;
	for (double err : epipolar_errors) {
		variance += (err - stats.mean_error) * (err - stats.mean_error);
	}
	stats.std_error = std::sqrt(variance / epipolar_errors.size());
	stats.max_error = max_err;

	return stats;
}


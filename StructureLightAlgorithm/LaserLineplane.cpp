// ��ƽ��궨���޸��������ߵ���ȡ������2025.05.15

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




// �����ƽ��궨ͼ���r t 
void OpticalPlaneCalibration::getRT(const Mat& cameraMatrix, const Mat& distCoeffs,
	       std::vector<cv::Mat>& rotationMat, std::vector<cv::Mat>& translationMat, 
		   bool show = false) {

	// �������̸��������� objectPoints
	vector<Point3f> objPoints;
	for (int i = 0; i < M_board_size.height; ++i) {
		for (int j = 0; j < M_board_size.width; ++j) {
			objPoints.emplace_back(j * M_squareSize.width, i * M_squareSize.width, 0.0f);
		}
	}

	// ���ļ�����վ�����
	std::ofstream outFile("extrinsic_parameters.txt", std::ios::out | std::ios::trunc);
	if (!outFile.is_open()) {
		std::cerr << "�޷����ļ� extrinsic_parameters.txt ����д�롣" << std::endl;
		return;
	}


	//��ȡͼ��
	FileManipulation ImageGet;
	std::vector<cv::Mat> BoardPics = ImageGet.ReadPics(M_PlaneBoardFilePath);
	//���ݱ궨�����ͣ����з��ദ��
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
					cv::cvtColor(img, imgCopy, cv::COLOR_GRAY2BGR);  // ���Ҷ�ͼתΪ��ɫͼ
				}
				else {
					imgCopy = img.clone();
				}
				cv::drawChessboardCorners(imgCopy, M_board_size, imgPoints, found);
				cv::namedWindow("corner", cv::WINDOW_NORMAL);
				cv::imshow("corner", imgCopy);
				cv::waitKey(500);
			}

			// ��� R/T
			Mat rvec, tvec;
			solvePnP(objPoints, imgPoints, cameraMatrix, distCoeffs, rvec, tvec);

			// д�� TXT �ļ�
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
		//��ȡ�ǵ�
		for (auto& img : BoardPics) {
			vector<Point2f> imgPoints;
			bool found = findChessboardCornersSB(img, M_board_size, imgPoints);
			// ��ϸ���ǵ�
			Mat gray;
			//cvtColor(img, gray, COLOR_BGR2GRAY);
			//cornerSubPix(gray, imgPoints, Size(11, 11), Size(-1, -1),
			//	         TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 30, 0.01));

			if (show) {
				cv::Mat imgCopy;
				if (img.channels() == 1) {
					cv::cvtColor(img, imgCopy, cv::COLOR_GRAY2BGR);  // ���Ҷ�ͼתΪ��ɫͼ
				}
				else {
					imgCopy = img.clone();
				}
				cv::drawChessboardCorners(imgCopy, M_board_size, imgPoints, found);
				cv::namedWindow("corner", cv::WINDOW_NORMAL);
				cv::imshow("corner", imgCopy);
				cv::waitKey(500);
			}

			// ��� R/T
			Mat rvec, tvec;
			solvePnP(objPoints, imgPoints, cameraMatrix, distCoeffs, rvec, tvec);

			rotationMat.push_back(rvec);
			translationMat.push_back(tvec);

			// д�� TXT �ļ�
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
	outFile.close(); // �ر��ļ�
	destroyAllWindows();
	cout << "����ͼ��������ȡ��ɣ������� " << rotationMat.size() << " ��ͼ��" << endl;
}



/**
* @brief ͨ����ȡout_camera_data.yml��ȡ�ڲκ����
* @param instrinsic_Mat	�������ڲεľ���
* @param distCoeff		����������ľ���
* @param translationMat	������ƽ�Ƶ�����
* @param rotationMat	��������ת������
*
*/
void OpticalPlaneCalibration::Get_Camera_Data(const std::string& cameraParamsPath, 
											  Mat& instrinsic_Mat, 
											  Mat& distCoeff) {
	//��ȡ����궨�Ľ��
	FileStorage camPara(cameraParamsPath, FileStorage::READ);
	if (!camPara.isOpened())
	{
		cout << "--------------���������ʧ��-----------" << endl;
		exit(1);
	}
	camPara["camera_matrix"] >> instrinsic_Mat;
	camPara["distortion_coefficients"] >> distCoeff;
	camPara.release();
}


/**
* @brief ͨ����С���˷���ͨ����������ϵ�¼���ƽ���ϵĵ�����ϼ���ƽ�淽��
* @param input	�������ƽ��ĵ�
*
* @return ƽ�淽��Ax+By+Cz+D=0
*/
cv::Mat plane_fitting(const vector<Point3f>& input) {
	cv::Mat dst = Mat(3, 3, CV_32F, Scalar(0));  // 3x3 ����
	cv::Mat out = Mat(3, 1, CV_32F, Scalar(0));

	// ���� dst �� out ����
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

	// �жϾ����Ƿ�Ϊ����
	if (abs(determinant(dst)) < 0.001) {
		cout << "�������죬�޷�����" << endl;
		return Mat(); // ���ؿվ���
	}

	// ͨ�� solve ������⣬��������
	Mat output;
	solve(dst, out, output);

	// ��ȡƽ�淽��ϵ��
	double a = output.at<float>(0, 0);
	double b = output.at<float>(1, 0);
	double c = -1.0;  // ǿ�� c = -1
	double d = output.at<float>(2, 0); // d ֵ��ʾ������

	cout << "���ƽ�淽��: " << a << " * x + " << b << " * y + " << c << " * z + " << d << " = 0" << endl;
	return output;
}


// ʹ�����ƽ��ķ�������һϵ��ƽ���
vector<Point3f> generatePlanePoints(double a, double b, double c, double d, float range_min, float range_max, double step) {
	vector<Point3f> planePoints;
	for (float x = range_min; x <= range_max; x += step) {
		for (float y = range_min; y <= range_max; y += step) {
			float z = a * x + b * y + d;  // ����z����
			planePoints.push_back(Point3f(x, y, z));
		}
	}
	return planePoints;
}


/**
* @brief ����ƽ����ϲв�
* @param threeCoordinate	�������ƽ��ĵ�
* @param a					ƽ�淽��Ax+By+Cz+D=0 �е�A
* @param b					ƽ�淽��Ax+By+Cz+D=0 �е�B
* @param c					ƽ�淽��Ax+By+Cz+D=0 �е�C
* @param d					ƽ�淽��Ax+By+Cz+D=0 �е�D
*
* @return ��ϲв�
*/
double Residual(const vector<Point3f>& threeCoordinate, double a, double b, double c, double d) {
	double distance = 0.0;
	int num = threeCoordinate.size();

	for (const auto& point : threeCoordinate) {
		double x = point.x, y = point.y, z = point.z;
		distance += (a * x + b * y + c * z + d) * (a * x + b * y + c * z + d) / (a * a + b * b + c * c);
	}
	cout << "ϵ���� " << a << "\t" << b << "\t" << c << "\t" << d << endl;
	double derror = sqrt(distance / num);  // ������������
	cout << "RMS��" << derror << endl;
	return derror;
}




void OptPlaneCalibration(const std::string& PlaneBoardFilename, const std::string& PlaneLineFilename,
						 const std::string& camParamPath, const std::string& PointsPath,
						 const std::string& PlanePath, enum CalibrationPlate BoardClass) {

	//-----------------------------------��ȡ�ڲκ����-----------------------------------
	//��ȡ����������Ϣ
	Mat instrinsic_Mat, distCoeff;
	vector<Mat>translationMat, rotationMat;
	OpticalPlaneCalibration getBaseParam(PlaneBoardFilename, PlaneLineFilename, BoardClass);
	getBaseParam.Get_Camera_Data(camParamPath, instrinsic_Mat, distCoeff);
	getBaseParam.getRT(instrinsic_Mat, distCoeff, rotationMat, translationMat);


	//------------------------------------------------------------------------------------
	vector<Point3f> laserpoint_cam_all;//����������ĵ�ļ���
	//��ȡͼƬ
	FileManipulation ImageGet;
	std::vector<cv::Mat>BoardPics;
	std::vector<cv::Mat>OptPlanePics;
	ImageGet.ReadPics(PlaneBoardFilename, BoardPics);
	ImageGet.ReadPics(PlaneLineFilename, OptPlanePics);


	for (size_t i = 0; i < BoardPics.size(); i++){
		cv::Mat AfterDis, Board_AfterDis;
		//-----------------------------------��ȡ����������-----------------------------------
		cv::Mat imgHSVMask;
		std::vector<cv::Rect> rois;

		Centerline centerline;
		centerline.picture_blur(OptPlanePics[i], imgHSVMask, rois);
		std::vector<cv::Point2d> center_points = centerline.ExtractCenters(OptPlanePics[i], imgHSVMask, rois);


		//------------------------------------ת�������ϵ------------------------------------
		Turn_to_3D turn3d;
		turn3d.Set_param(instrinsic_Mat, distCoeff, translationMat[i], rotationMat[i], BoardPics[i], center_points, BoardClass);
		turn3d.Get_plane();	//��ȡ�궨��ƽ�淽��
		turn3d.Get_Intersection();	//��ȡ��������궨�彻��
		vector<Point3d> centers_cam = turn3d.LaserPoints_cam();	//������������������ϵת�����������ϵ

		//�������Ľ�������laserpoint_cam_all
		laserpoint_cam_all.insert(laserpoint_cam_all.end(), centers_cam.begin(), centers_cam.end());
		cout << "��" << to_string(i + 1) << "�ű���ɹ�" << endl;
	}
	//------------------------------------�����------------------------------------
	ofstream ofs;
	ofs.open(PointsPath, ios::trunc);
	for (int i = 0; i < laserpoint_cam_all.size(); i++){
		ofs << laserpoint_cam_all[i].x << "\t" << laserpoint_cam_all[i].y << "\t" << laserpoint_cam_all[i].z << endl;
	}
	//------------------------------------------------------------------------------

	//-----------------------------------���ƽ��-----------------------------------
	cv::Mat plane0102 = plane_fitting(laserpoint_cam_all);
	// ƽ�淽�� ax + by + cz + d = 0

	// ��ȡ���ƽ�淽�̵�ϵ��
	double a = plane0102.at<float>(0, 0);
	double b = plane0102.at<float>(1, 0);
	double c = -1.0;
	double d = plane0102.at<float>(2, 0);

	// ����������
	double derror = Residual(laserpoint_cam_all, a, b, c, d);

	// �������ƽ���ϵĵ�
	vector<Point3f> planePoints = generatePlanePoints(a, b, c, d, -50, 50, 0.1);
	ofstream ofs_plane(PlanePath, ios::trunc);
	for (int i = 0; i < planePoints.size(); i++) {
		ofs_plane << planePoints[i].x << "\t" << planePoints[i].y << "\t" << planePoints[i].z << endl;
	}
	ofs_plane.close();
}
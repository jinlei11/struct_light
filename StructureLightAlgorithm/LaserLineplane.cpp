// ��ƽ��궨���޸��������ߵ���ȡ������2025.05.15

#include"include/LaserLineplane.h"
#include <Eigen/Dense>


using namespace std;
using namespace cv;


// �����ƽ��궨ͼ���r t 
void getRT(int num, const string& imagePath, const Size& boardSize, float squareSize,
	const Mat& cameraMatrix, const Mat& distCoeffs, vector<Mat>& rotationMat, vector<Mat>& translationMat, bool show = false) {

	// �������̸��������� objectPoints
	vector<Point3f> objPoints;
	for (int i = 0; i < boardSize.height; ++i) {
		for (int j = 0; j < boardSize.width; ++j) {
			objPoints.emplace_back(j * squareSize, i * squareSize, 0.0f);
		}
	}

	rotationMat.clear();
	translationMat.clear();

	// ���ļ�����վ�����
	std::ofstream outFile("extrinsic_parameters.txt", std::ios::out | std::ios::trunc);
	if (!outFile.is_open()) {
		std::cerr << "�޷����ļ� extrinsic_parameters.txt ����д�롣" << std::endl;
		return;
	}

	for (int i = 1; i <= num; ++i) {
		string filename = imagePath + to_string(i) + ".bmp";
		Mat img = imread(filename);
		if (img.empty()) {
			cout << "�޷���ȡͼ��: " << filename << endl;
			continue;
		}

		vector<Point2f> imgPoints;
		bool found = findChessboardCornersSB(img, boardSize, imgPoints);
		if (!found) {
			cout << "δ�ҵ��ǵ�: " << filename << endl;
			continue;
		}

		// ��ϸ���ǵ�
		Mat gray;
		cvtColor(img, gray, COLOR_BGR2GRAY);
		cornerSubPix(gray, imgPoints, Size(11, 11), Size(-1, -1),
			TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 30, 0.01));

		// ��� R/T
		Mat rvec, tvec;
		solvePnP(objPoints, imgPoints, cameraMatrix, distCoeffs, rvec, tvec);

		rotationMat.push_back(rvec);
		translationMat.push_back(tvec);

		// д�� TXT �ļ�
		outFile << "Image " << i << ":" << std::endl;
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

		if (show) {
			drawChessboardCorners(img, boardSize, imgPoints, found);
			imshow("corners", img);
			waitKey(100);
		}

		//cout << "ͼ�� " << i << " ������ɡ�" << endl;
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
void Get_Camera_Data(Mat& instrinsic_Mat, Mat& distCoeff) {
	//��ȡ����궨�Ľ��
	//string camPath = "D:\\ʵ��������\\24.09.25�ؽ�\\�ڲ�\\out_camera_data.yml";
	string camPath = "C:\\Users\\W\\Desktop\\0530\\camera\\out_camera_data.yml";
	//�������Լ��ṩ�ĸ�Ŀ¼�µ��������
	//string camPath = "out_camera_data���.yml";
	FileStorage camPara(camPath, FileStorage::READ);
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
Mat plane_fitting(const vector<Point3f>& input) {
	Mat dst = Mat(3, 3, CV_32F, Scalar(0));  // 3x3 ����
	Mat out = Mat(3, 1, CV_32F, Scalar(0));

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




//int main(int argc, char* argv){
//	Eigen::MatrixXd m(2, 2);
//	m << 0, 1, 2, 3;
//	Eigen::Vector3d x;
//
//	//-----------------------------------��ȡ�ڲκ����-----------------------------------
//	Mat instrinsic_Mat, distCoeff;
//	vector<Mat>translationMat, rotationMat;
//	int num = 11;
//
//	string PlaneBoardFilename = "C:\\Users\\W\\Desktop\\0530\\plane\\boards\\";
//	string PlaneLineFilename = "C:\\Users\\W\\Desktop\\0530\\plane\\lines\\";
//	Size boardSize(8, 11);		// �궨��ߴ�
//	float squareSize = 3.0f;     // �궨�巽��ߴ磨���ף�
//
//	Get_Camera_Data(instrinsic_Mat, distCoeff);
//
//	getRT(num, PlaneBoardFilename, boardSize, squareSize, instrinsic_Mat, distCoeff, rotationMat, translationMat);
//
//
//	//------------------------------------------------------------------------------------
//
//	vector<Point3f> laserpoint_cam_all;//����������ĵ�ļ���
//
//	for (size_t i = 0; i < num; i++){
//		//��ȡ����궨ͼƬ
//		Mat image = imread(PlaneBoardFilename + to_string(i + 1) + ".bmp");
//		Mat src = imread(PlaneLineFilename + to_string(i + 1) + ".bmp");
//
//		Mat AfterDis, Board_AfterDis;
//		//-----------------------------------��ȡ����������-----------------------------------
//		Mat imgHSVMask;
//		int startY, endY;
//		vector<Rect> rois;
//
//		Centerline centerline;
//		centerline.picture_blur(src, imgHSVMask, rois);
//		centerline.ExtractCenters(src, imgHSVMask, rois);
//
//		vector<Point2d> center_points;
//		center_points = centerline.smooth_points;
//		//------------------------------------------------------------------------------------
//
//		//------------------------------------ת�������ϵ------------------------------------
//		Turn_to_3D turn3d;
//		turn3d.Set_param(instrinsic_Mat, distCoeff, translationMat[i], rotationMat[i], image, center_points, 8, 11);
//		turn3d.Get_plane();	//��ȡ�궨��ƽ�淽��
//		turn3d.Get_Intersection(i + 1);	//��ȡ��������궨�彻��
//
//		turn3d.LaserPoints_cam();	//������������������ϵת�����������ϵ
//
//		vector<Point3d> centers_cam;
//		centers_cam = turn3d.laser_cam;	//��ȡ�������ϵ�µĵ�
//
//
//		//�������Ľ�������laserpoint_cam_all
//		laserpoint_cam_all.insert(laserpoint_cam_all.end(), centers_cam.begin(), centers_cam.end());
//		cout << "��" << to_string(i + 1) << "�ű���ɹ�" << endl;
//	}
//	//------------------------------------�����------------------------------------
//	ofstream ofs;
//	ofs.open("Intersections.txt", ios::trunc);
//	for (int i = 0; i < laserpoint_cam_all.size(); i++){
//		ofs << laserpoint_cam_all[i].x << "\t" << laserpoint_cam_all[i].y << "\t" << laserpoint_cam_all[i].z << endl;
//	}
//	//------------------------------------------------------------------------------
//
//	//-----------------------------------���ƽ��-----------------------------------
//	Mat plane0102 = plane_fitting(laserpoint_cam_all);
//	// ƽ�淽�� ax + by + cz + d = 0
//
//	// ��ȡ���ƽ�淽�̵�ϵ��
//	double a = plane0102.at<float>(0, 0);
//	double b = plane0102.at<float>(1, 0);
//	double c = -1.0;
//	double d = plane0102.at<float>(2, 0);
//
//
//	// ����������
//	double derror = Residual(laserpoint_cam_all, a, b, c, d);
//
//	// �������ƽ���ϵĵ�
//	vector<Point3f> planePoints = generatePlanePoints(a, b, c, d, -50, 50, 0.1);
//
//	ofstream ofs_plane("PLANE.txt", ios::trunc);
//
//	for (int i = 0; i < planePoints.size(); i++) {
//		ofs_plane << planePoints[i].x << "\t" << planePoints[i].y << "\t" << planePoints[i].z << endl;
//	}
//	ofs_plane.close();
//
//	waitKey(0);	system("pause");
//	return 0;
//}
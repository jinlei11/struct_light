// 光平面标定，修改了中心线的提取方法。2025.05.15

#include"include/LaserLineplane.h"
#include <Eigen/Dense>


using namespace std;
using namespace cv;


// 计算光平面标定图像的r t 
void getRT(int num, const string& imagePath, const Size& boardSize, float squareSize,
	const Mat& cameraMatrix, const Mat& distCoeffs, vector<Mat>& rotationMat, vector<Mat>& translationMat, bool show = false) {

	// 构建棋盘格世界坐标 objectPoints
	vector<Point3f> objPoints;
	for (int i = 0; i < boardSize.height; ++i) {
		for (int j = 0; j < boardSize.width; ++j) {
			objPoints.emplace_back(j * squareSize, i * squareSize, 0.0f);
		}
	}

	rotationMat.clear();
	translationMat.clear();

	// 打开文件并清空旧内容
	std::ofstream outFile("extrinsic_parameters.txt", std::ios::out | std::ios::trunc);
	if (!outFile.is_open()) {
		std::cerr << "无法打开文件 extrinsic_parameters.txt 进行写入。" << std::endl;
		return;
	}

	for (int i = 1; i <= num; ++i) {
		string filename = imagePath + to_string(i) + ".bmp";
		Mat img = imread(filename);
		if (img.empty()) {
			cout << "无法读取图像: " << filename << endl;
			continue;
		}

		vector<Point2f> imgPoints;
		bool found = findChessboardCornersSB(img, boardSize, imgPoints);
		if (!found) {
			cout << "未找到角点: " << filename << endl;
			continue;
		}

		// 精细化角点
		Mat gray;
		cvtColor(img, gray, COLOR_BGR2GRAY);
		cornerSubPix(gray, imgPoints, Size(11, 11), Size(-1, -1),
			TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 30, 0.01));

		// 求解 R/T
		Mat rvec, tvec;
		solvePnP(objPoints, imgPoints, cameraMatrix, distCoeffs, rvec, tvec);

		rotationMat.push_back(rvec);
		translationMat.push_back(tvec);

		// 写入 TXT 文件
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

		//cout << "图像 " << i << " 处理完成。" << endl;
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
void Get_Camera_Data(Mat& instrinsic_Mat, Mat& distCoeff) {
	//读取相机标定的结果
	//string camPath = "D:\\实验室资料\\24.09.25重建\\内参\\out_camera_data.yml";
	string camPath = "C:\\Users\\W\\Desktop\\0530\\camera\\out_camera_data.yml";
	//或者是自己提供的根目录下的相机参数
	//string camPath = "out_camera_data相机.yml";
	FileStorage camPara(camPath, FileStorage::READ);
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
Mat plane_fitting(const vector<Point3f>& input) {
	Mat dst = Mat(3, 3, CV_32F, Scalar(0));  // 3x3 矩阵
	Mat out = Mat(3, 1, CV_32F, Scalar(0));

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




//int main(int argc, char* argv){
//	Eigen::MatrixXd m(2, 2);
//	m << 0, 1, 2, 3;
//	Eigen::Vector3d x;
//
//	//-----------------------------------获取内参和外参-----------------------------------
//	Mat instrinsic_Mat, distCoeff;
//	vector<Mat>translationMat, rotationMat;
//	int num = 11;
//
//	string PlaneBoardFilename = "C:\\Users\\W\\Desktop\\0530\\plane\\boards\\";
//	string PlaneLineFilename = "C:\\Users\\W\\Desktop\\0530\\plane\\lines\\";
//	Size boardSize(8, 11);		// 标定板尺寸
//	float squareSize = 3.0f;     // 标定板方格尺寸（毫米）
//
//	Get_Camera_Data(instrinsic_Mat, distCoeff);
//
//	getRT(num, PlaneBoardFilename, boardSize, squareSize, instrinsic_Mat, distCoeff, rotationMat, translationMat);
//
//
//	//------------------------------------------------------------------------------------
//
//	vector<Point3f> laserpoint_cam_all;//存放所有中心点的集合
//
//	for (size_t i = 0; i < num; i++){
//		//读取相机标定图片
//		Mat image = imread(PlaneBoardFilename + to_string(i + 1) + ".bmp");
//		Mat src = imread(PlaneLineFilename + to_string(i + 1) + ".bmp");
//
//		Mat AfterDis, Board_AfterDis;
//		//-----------------------------------获取激光中心线-----------------------------------
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
//		//------------------------------------转相机坐标系------------------------------------
//		Turn_to_3D turn3d;
//		turn3d.Set_param(instrinsic_Mat, distCoeff, translationMat[i], rotationMat[i], image, center_points, 8, 11);
//		turn3d.Get_plane();	//获取标定板平面方程
//		turn3d.Get_Intersection(i + 1);	//获取激光线与标定板交点
//
//		turn3d.LaserPoints_cam();	//交点坐标由像素坐标系转换到相机坐标系
//
//		vector<Point3d> centers_cam;
//		centers_cam = turn3d.laser_cam;	//获取相机坐标系下的点
//
//
//		//将读到的交点存放至laserpoint_cam_all
//		laserpoint_cam_all.insert(laserpoint_cam_all.end(), centers_cam.begin(), centers_cam.end());
//		cout << "第" << to_string(i + 1) << "张保存成功" << endl;
//	}
//	//------------------------------------保存点------------------------------------
//	ofstream ofs;
//	ofs.open("Intersections.txt", ios::trunc);
//	for (int i = 0; i < laserpoint_cam_all.size(); i++){
//		ofs << laserpoint_cam_all[i].x << "\t" << laserpoint_cam_all[i].y << "\t" << laserpoint_cam_all[i].z << endl;
//	}
//	//------------------------------------------------------------------------------
//
//	//-----------------------------------拟合平面-----------------------------------
//	Mat plane0102 = plane_fitting(laserpoint_cam_all);
//	// 平面方程 ax + by + cz + d = 0
//
//	// 提取拟合平面方程的系数
//	double a = plane0102.at<float>(0, 0);
//	double b = plane0102.at<float>(1, 0);
//	double c = -1.0;
//	double d = plane0102.at<float>(2, 0);
//
//
//	// 计算均方误差
//	double derror = Residual(laserpoint_cam_all, a, b, c, d);
//
//	// 生成拟合平面上的点
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
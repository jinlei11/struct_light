#include"ALLHeads.h"



std::string currentPath = "C:\\mycode\\StructureLightAlgorithms\\";
std::string CalibrationResultPath = currentPath + "CalibrationImage\\CalirationInfo.yml";
std::string Calibration_LeftCamImgsPath = currentPath + "CalibrationImage\\left\\";
std::string Calibration_RightCamImgsPath = currentPath + "CalibrationImage\\right\\";

std::string CamResponseCurveImgsPath = currentPath + "CamResponseCurveImgs\\";
std::string CamResponseCurveReuslt = currentPath + "CamResponseResult\\";

std::string RecoveryPahsePath = currentPath + "RecoveryPahse\\";

std::string CapturePictures_Left = currentPath + "CapturePictures\\A7500CU35_DC38534AAK00017\\";
std::string CapturePictures_Right = currentPath + "CapturePictures\\A7500MU35_DF13188AAK00019\\";
std::string CloudPointsSavedPath = currentPath + "CloudPoints\\result.txt";


std::string ProjectionPattern = currentPath + "ProjectionPattern\\";

//������ͼ����
int GrayImgsNum = 7;
//����ͼ����
int PhaseImgsNum = 4;
//�ع����
int exposureNums = 6;


//����ͼ������
int auxiliaryImageNum = 1;

int main() {

	/********************  ���ɸ�����ͼ��   *************************/
	//GrayCodeOperation* TEST = new GrayMethod_2(6, DLP4500,0);
	///*TEST->GenerateAndSave_pics(ProjectionPattern);*/
	//TEST->GeneratePatternWithValue(ProjectionPattern, 255);
	/********************  ���ɸ�����ͼ��   *************************/


	/********************   ��������ͼ��   *************************/
	//std::string fielpath = R"(../../../../ProjectionPattern/DLP3010/)";
	//PhasePatternOperation Phase(32.0, DLP4710);
	//Phase.BuildAndSave_PhaseShiftingPattern(fielpath);
	/********************   ��������ͼ��   *************************/


	/***********************    ����궨   *********************************/
	//CameraCalibration calibration(Plate_smallCircle, 
	//								Calibration_LeftCamImgsPath, 
	//								Calibration_RightCamImgsPath, 
	//								CalibrationResultPath, 
	//								true);
	//calibration.M_MyCameraCalibrationAndSave();
	/***********************    ����궨   *********************************/


	/***********************    ͼ��Ԥ����   *********************************/
	ImageProcessor processor(CapturePictures_Right + "12.bmp");
	if (processor.loadImage()) {
		processor.process();
		processor.showResult();
		processor.saveResult(CapturePictures_Right + "processed_output.png");
	}

	/***********************    ͼ��Ԥ����   *********************************/







	/***********************    ��ά�ؽ�   *********************************/
	////�����ع�
	//SolvingPacel Reconstruction_SingleExposure(CapturePictures_Left,
	//										   CapturePictures_Right, 
	//										   CalibrationResultPath,
	//										   GrayImgsNum,
	//										   PhaseImgsNum,
	//										   A7500,
	//										   1);
	 
	 
	//bool firstExectue = false;
	////����ع�
	//SolvingPacel Reconstruction_MulitExposure( CapturePictures_Left,
	//										   CapturePictures_Right,
	//										   CalibrationResultPath,
	//										   GrayImgsNum,
	//										   PhaseImgsNum,
	//										   A7500,
	//										   exposureNums,
	//										   auxiliaryImageNum,
	//										   true);
	//cv::Mat CorrectionLeft = cv::Mat::zeros(2048, 2448, CV_32FC1);
	//cv::Mat CorrectionRight = cv::Mat::zeros(2048, 2448, CV_32FC1);
	//cv::Mat Q = Reconstruction_MulitExposure.BinocularPhaseRecovery(CorrectionLeft, CorrectionRight,254);


	//CloudPointsGeneration Generate3DPoints(CloudPointsSavedPath); 
	//cv::Mat Disparity = cv::Mat::zeros(2048, 2448, CV_32FC1); ;
	//Generate3DPoints.CalculateDisparity(CorrectionLeft, CorrectionRight, Disparity,1500);
	//cv::Mat DepthMap = Generate3DPoints.CalculateDepthMap(Disparity, Q, -10000, 20000, -2000, 1000);
	//Generate3DPoints.SaveCloudPointsToTxt(DepthMap, Q, CloudPointsSavedPath);






	
	/***********************    ��ά�ؽ�   *********************************/


	/***********************    ��ά�ؽ�:���������Ӧ����   *********************************/
	//std::vector<float> Exposures = { 0.02f, 0.04f, 0.05f, 0.07f, 0.09f, 0.12f ,0.15f, 0.20f, 0.25f ,0.30f,0.40f,0.50f,0.70f };

	//HDR_CameraResponseCurve test(CamResponseCurveImgsPath, CamResponseCurveReuslt, Exposures, 2448, 2048);

	//cv::Mat image = cv::imread(CamResponseCurveImgsPath + "10.bmp");
	//test.CalculateIrradianceImage(image, 0.30f, 1.0f);


	/***********************    ��ά�ؽ�:���������Ӧ����   *********************************/


	return 0;
}
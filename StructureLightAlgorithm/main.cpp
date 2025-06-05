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

std::string ROIPathLeft = currentPath + "ROI\\left\\";
std::string ROIPathRight = currentPath + "ROI\\right\\";

std::string ProjectionPattern = currentPath + "ProjectionPattern\\";

//格雷码图像数
int GrayImgsNum = 7;
//相移图像数
int PhaseImgsNum = 4;
//曝光次数
int exposureNums = 4;


//辅助图像数量
int auxiliaryImageNum = 1;

int main() {

	/********************  生成格雷码图案   *************************/
	//GrayCodeOperation* TEST = new GrayMethod_2(6, DLP4500,0);
	///*TEST->GenerateAndSave_pics(ProjectionPattern);*/
	//TEST->GeneratePatternWithValue(ProjectionPattern, 255);
	/********************  生成格雷码图案   *************************/


	/********************   生成相移图案   *************************/
	//std::string fielpath = R"(../../../../ProjectionPattern/DLP3010/)";
	//PhasePatternOperation Phase(32.0, DLP4710);
	//Phase.BuildAndSave_PhaseShiftingPattern(fielpath);
	/********************   生成相移图案   *************************/


	/***********************    相机标定   *********************************/
	CameraCalibration calibration(Plate_doubleCircle,
									Calibration_LeftCamImgsPath, 
									Calibration_RightCamImgsPath, 
									CalibrationResultPath, 
									true);
	calibration.M_MyCameraCalibrationAndSave();
	/***********************    相机标定   *********************************/


	/***********************    图像预处理   *********************************/
	//ImageProcessor processor(ROIPathLeft + "11.bmp");
	//if (processor.loadImage()) {
	//	processor.process();
	//	processor.showResult();
	//	processor.saveResult(ROIPathLeft + "ROI.bmp");
	//}

	//ImageProcessor processor_right(ROIPathRight + "12.bmp");
	//if (processor_right.loadImage()) {
	//	processor_right.process();
	//	processor_right.showResult();
	//	processor_right.saveResult(ROIPathRight + "ROI.bmp");
	//}
	/***********************    图像预处理   *********************************/





	/***********************    三维重建   *********************************/	 
	//bool firstExectue = false;
	////多次曝光
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
	//std::vector<cv::Mat> StripesLevel_Left;
	//std::vector<cv::Mat> StripesLevel_Right;
	//cv::Mat Q = Reconstruction_MulitExposure.BinocularPhaseRecovery(CorrectionLeft, CorrectionRight,
	//																StripesLevel_Left, StripesLevel_Right,
	//																254);
	////图像滤波
	//PhaseNoiseFilter Filter;
	//cv::Mat mask_left, mask_right;
	//cv::Mat FilterCorrectionLeft, FilterCorrectionRight;
	////基于相位梯度滤波
	//Filter.filterPhaseByGradient(CorrectionLeft, 3.f, mask_left, FilterCorrectionLeft);
	//Filter.filterPhaseByGradient(CorrectionRight, 3.f, mask_right, FilterCorrectionRight);
	////连通域滤波
	//Filter.filterInPlace(CorrectionLeft);
	//Filter.filterInPlace(CorrectionRight);


	////点云计算
	//CloudPointsGeneration Generate3DPoints(CloudPointsSavedPath); 
	//cv::Mat Disparity = cv::Mat::zeros(2048, 2448, CV_32FC1); ;
	//Generate3DPoints.CalculateDisparity(FilterCorrectionLeft, FilterCorrectionRight, Disparity,1500);
	//cv::Mat DepthMap = Generate3DPoints.CalculateDepthMap(Disparity, Q, -10000, 20000, -2000, 1000);
	//Generate3DPoints.SaveCloudPointsToTxt(DepthMap, Q, CloudPointsSavedPath);


	/***********************    三维重建   *********************************/


	/***********************    三维重建:基于相机响应曲线   *********************************/
	//std::vector<float> Exposures = { 0.02f, 0.04f, 0.05f, 0.07f, 0.09f, 0.12f ,0.15f, 0.20f, 0.25f ,0.30f,0.40f,0.50f,0.70f };

	//HDR_CameraResponseCurve test(CamResponseCurveImgsPath, CamResponseCurveReuslt, Exposures, 2448, 2048);

	//cv::Mat image = cv::imread(CamResponseCurveImgsPath + "10.bmp");
	//test.CalculateIrradianceImage(image, 0.30f, 1.0f);


	/***********************    三维重建:基于相机响应曲线   *********************************/


	return 0;
}
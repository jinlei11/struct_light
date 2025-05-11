#include "SolvingPacel.h"



SolvingPacel::SolvingPacel(const std::string& ImagepathLeft,
						   const std::string& ImagepathRight,
	                       const std::string& CalibrationResultPath,
						   int GrayImgsNum,
						   int PhaseImgsNum,
						   int width,
						   int height,
						   int AncillaryPatternNum,
						   bool isRename)
	:M_Width(width),M_Height(height),M_GrayImgsNum(GrayImgsNum),M_PhaseImgsNum(PhaseImgsNum), 
	M_ImagesPathLeft(ImagepathLeft), M_ImagesPathRight(ImagepathRight),
	M_AncillaryPatternNum(AncillaryPatternNum), M_CalibrationResultPath(CalibrationResultPath),
	M_isRename(isRename)
{


}


SolvingPacel::SolvingPacel(const std::string& ImagepathLeft,
						   const std::string& ImagepathRight,
	                       const std::string& CalibrationResultPath,
						   int GrayImgsNum,
						   int PhaseImgsNum,
						   enum CameraSolutioin Form,
						   int AncillaryPatternNum,
						   bool isRename)
	:M_GrayImgsNum(GrayImgsNum), M_PhaseImgsNum(PhaseImgsNum),
	M_ImagesPathLeft(ImagepathLeft), M_ImagesPathRight(ImagepathRight),
	M_AncillaryPatternNum(AncillaryPatternNum), M_CalibrationResultPath(CalibrationResultPath),
	M_isRename(isRename)
{
	switch (Form)
	{
	case A5131:
		this->M_Width = 1280;
		this->M_Height = 1024;
		break;
	case A5201:
		this->M_Width = 1920;
		this->M_Height = 1200;
		break;
	case A7300:
		this->M_Width = 2048;
		this->M_Height = 1536;
		break;
	case A7500:
		this->M_Width = 2448;
		this->M_Height = 2048;
		break;
	default:
		break;
	}
}


SolvingPacel::SolvingPacel(const std::string& ImagepathLeft,
						   const std::string& ImagepathRight,
	                       const std::string& CalibrationResultPath,
						   int GrayImgsNum,
						   int PhaseImgsNum,
						   int width,
					       int height,
						   int ExposureNums,
						   int AncillaryPatternNum,
						   bool isRename)
	:M_GrayImgsNum(GrayImgsNum), M_PhaseImgsNum(PhaseImgsNum),
	M_ImagesPathLeft(ImagepathLeft), M_ImagesPathRight(ImagepathRight),
	M_Width(width), M_Height(height), M_AncillaryPatternNum(AncillaryPatternNum),
	M_ExposureNums(ExposureNums), M_CalibrationResultPath(CalibrationResultPath),
	M_isRename(isRename)
{
	M_IsMulExposure = true;
}


SolvingPacel::SolvingPacel(const std::string& ImagepathLeft,
						   const std::string& ImagepathRight,
	                       const std::string& CalibrationResultPath,
						   int NumOfGray,
						   int NumOfPhase,
						   enum CameraSolutioin Form,
						   int ExposureNums,
						   int AncillaryPatternNum,
						   bool isRename)
	:M_AncillaryPatternNum(AncillaryPatternNum), M_ExposureNums(ExposureNums),
	M_GrayImgsNum(NumOfGray), M_PhaseImgsNum(NumOfPhase),
	M_CalibrationResultPath(CalibrationResultPath),
	M_ImagesPathLeft(ImagepathLeft), M_ImagesPathRight(ImagepathRight),
	M_isRename(isRename)
{
	M_IsMulExposure = true;
	switch (Form)
	{
	case A5131:
		this->M_Width = 1280;
		this->M_Height = 1024;
		break;
	case A5201:
		this->M_Width = 1920;
		this->M_Height = 1200;
		break;
	case A7300:
		this->M_Width = 2048;
		this->M_Height = 1536;
		break;
	case A7500:
		this->M_Width = 2448;
		this->M_Height = 2048;
		break;
	default:
		break;
	}
}



void SolvingPacel::SingleDePhase(const std::vector<cv::Mat>& ProjectionPatterns, 
								 cv::Mat& Phase, 
								 int regulatory = 5) {
	//����������ʱ��������Ƿ�Ϊ�棬�������Ϊ�٣�����ֹ����ִ�в���ӡ������Ϣ��ʹ��OpenCV�еĺ�
	CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum));
	//���������λ,�ڶ�ȡͼ��ʱ�Ѿ������rows��cols��ֵ
	cv::Mat wrappedphase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	//ÿ�����Ƶ�ֵ
	const float shiftVal = static_cast<float>(CV_2PI) / M_PhaseImgsNum;

	cv::Mat testPhase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range)
		{
			//�����������ͼ��ÿһ�е���ʼ��ַ�������뵽imgsPtrs������
			std::vector<const float*> imgsPtrs(M_PhaseImgsNum);
			for (int i = range.start; i < range.end; ++i)
			{
				for (int k = 0; k < M_PhaseImgsNum; ++k)
				{
					imgsPtrs[k] = ProjectionPatterns[M_GrayImgsNum +k].ptr<float>(i);
				}
				//��Ű�����λÿһ�е���ʼ��ַ
				auto wrappedphasePtr0 = wrappedphase.ptr<float>(i);
				//�����б���
				for (int j = 0; j < M_Width; ++j)
				{
					float molecules = 0.f, denominator = 0.f;
					for (int k = 0; k < M_PhaseImgsNum; ++k)
					{
						molecules += imgsPtrs[k][j] * sin(k * shiftVal);
						denominator += imgsPtrs[k][j] * cos(k * shiftVal);
					}
					(sqrt((molecules * molecules) + (denominator * denominator)) * 2 / M_PhaseImgsNum) > regulatory ? wrappedphasePtr0[j] = -atan2(molecules, denominator) : wrappedphasePtr0[j] = NAN;
				}
			}
		});

	cv::Mat A = cv::Mat::zeros(ProjectionPatterns[0].size(), CV_32FC1);
	for (int i = 0; i < M_PhaseImgsNum; ++i)
	{
		A += ProjectionPatterns[M_GrayImgsNum +i] / M_PhaseImgsNum;
	}

	//�����
	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range)
		{
			//��Ÿ�����ͼ���ÿһ�е���ʼ��ַ����ʼ��СΪ�߸�
			std::vector<const float*> imgsPtrs1(M_GrayImgsNum);

			for (int i = range.start; i < range.end; ++i)
			{
				for (int j = 0; j < M_GrayImgsNum; ++j)
				{
					imgsPtrs1[j] = ProjectionPatterns[j].ptr<float>(i);
				}
				//���������λ����Ҫ�õ�����ǿ��A��������λwrappedphase�����ǻ�ȡ���еĵ�ַ
				auto PhasePtr = testPhase.ptr<float>(i);
				auto confidencePtr = A.ptr<float>(i);
				auto wrappedPhasePtr = wrappedphase.ptr<float>(i);

				for (int j = 0; j < M_Width; ++j)
				{
					int K1 = 0, tempVal = 0;
					for (int k = 0; k < M_GrayImgsNum - 1; ++k)
					{
						//������ͼ��ֵ���ڻ���ǿ��A��ֵ���õ�0����1���õ������Ƶ�ֵ
						//tempVal���뷵�صĽ����������������ͬΪ0����ͬΪ1��
						tempVal ^= imgsPtrs1[k][j] > confidencePtr[j];
						// K1�õ��ľ��Ƕ�����ת��Ϊʮ���Ƶ���������Ĳ������ݶ����Ƶ�ֵת��Ϊʮ����
						//���յõ���K1����ʮ���Ƶ���
						K1 = (K1 << 1) + tempVal;
					}
					//�ж����Ż��������룬���ǽ�������Ĳ������������ɵĹ�ʽ�е㲻һ��
					tempVal ^= imgsPtrs1[M_GrayImgsNum - 1][j] > confidencePtr[j];
					const int K2 = ((K1 << 1) + tempVal + 1) / 2;

					//���յõ���K1��K2��ֵ����ʮ���Ƶ��������Խ�����λչ���ļ���
					if (wrappedPhasePtr[j] <= CV_PI / 2)
					{
						PhasePtr[j] = wrappedPhasePtr[j] + 2 * CV_PI * K2;
					}
					else if (wrappedPhasePtr[j] >= CV_PI * 3 / 2)
					{
						PhasePtr[j] = wrappedPhasePtr[j] + 2 * CV_PI * (K2 - 1);
					}
					else
					{
						PhasePtr[j] = wrappedPhasePtr[j] + 2 * CV_PI * K1;
					}
				}
			}
		});
	Phase = testPhase.clone();
	return;
};



void SolvingPacel::MultiDePhase(const std::vector<cv::Mat>& ProjectionPatterns,
								 std::vector<cv::Mat>& SovelPacelPhases,
								 int ExposureNums, 
								 int regulatory = 5) {
	//����������ʱ��������Ƿ�Ϊ�棬�������Ϊ�٣�����ֹ����ִ�в���ӡ������Ϣ��ʹ��OpenCV�еĺ�
	if (M_AncillaryPatternNum == 0) {
		//δʹ�ø���ͼ��
		CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum)* ExposureNums);
		int OnceProjectionPatternNum = M_GrayImgsNum + M_PhaseImgsNum;
		for (int ExposureNum = 0; ExposureNum < ExposureNums; ExposureNum++) {
			std::vector<cv::Mat>OnceExposure(ProjectionPatterns.begin() + OnceProjectionPatternNum * ExposureNum, ProjectionPatterns.begin() + OnceProjectionPatternNum * (ExposureNum + 1) - 1);
			SingleDePhase(OnceExposure, SovelPacelPhases[ExposureNum]);
		}
	}
	else if (M_AncillaryPatternNum == 1) {
		//ʹ��һ�Ÿ���ͼ��
		CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum + M_AncillaryPatternNum) * ExposureNums);
		int OnceProjectionPatternNum = M_GrayImgsNum + M_PhaseImgsNum + M_AncillaryPatternNum;
		for (int ExposureNum = 0; ExposureNum < ExposureNums; ExposureNum++) {
			std::vector<cv::Mat>OnceExposure(ProjectionPatterns.begin() + OnceProjectionPatternNum * ExposureNum, ProjectionPatterns.begin() + OnceProjectionPatternNum * (ExposureNum + 1) - M_AncillaryPatternNum);//ʹ�����ַ�������ʱ��������಻�����Ҳ�
			SingleDePhase(OnceExposure, SovelPacelPhases[ExposureNum]);
			//�洢����ͼ��
			M_AncillaryPatterns.push_back(ProjectionPatterns[OnceProjectionPatternNum * (ExposureNum+1) - 1]);
		}
	}
	return;
};


//�Ľ��汾
void SolvingPacel::PhasesIntegration(std::vector<cv::Mat>& SovelPacelPhases,
									 int ExposureNums,
									 cv::Mat& IntegrationPhase,
								   	 int GrayThresholds,
									 CameraPosition position) {
	std::vector<cv::Mat> AncillaryPatterns;
	if (position == Left) {
		AncillaryPatterns.assign(M_AncillaryPatterns.begin(), M_AncillaryPatterns.begin() + ExposureNums);
	}
	else {
		AncillaryPatterns.assign(M_AncillaryPatterns.begin()+ ExposureNums, M_AncillaryPatterns.end());
	}
	//����Ա����M_AncillaryPatterns�Ƿ���и���ͼ��
	CV_Assert(AncillaryPatterns.size() >= static_cast<size_t>(ExposureNums));
	//����ƽ����λ���õļ������󣬴�Ų�ͬͼ������ͬ����λ��Ӧ����Ĵ���
	cv::Mat PhasesIntegrationMask = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat IntegrationMask = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);

	//�����ع����������ѭ��
	for (int n = 0; n < ExposureNums; n++) {
		std::cout << n << std::endl;
		for (int i = 0; i < M_Height; i++) {
			for (int j = 0; j < M_Width; j++) {
				if (isnan(SovelPacelPhases[n].at<float>(i, j))) {
					continue;
				}
				else if (AncillaryPatterns[n].at<float>(i, j) < GrayThresholds &&
					AncillaryPatterns[n].at<float>(i, j) > 0.5 &&
					IntegrationMask.at<float>(i, j) == 0) {
					IntegrationPhase.at<float>(i, j) = SovelPacelPhases[n].at<float>(i, j);
					IntegrationMask.at<float>(i, j) = 1;
				}
			}
		}
	}
};



//��ʼ�汾
//void SolvingPacel::PhasesIntegration(std::vector<cv::Mat>& SovelPacelPhases,
//									 int ExposureNums,
//									 cv::Mat& IntegrationPhase,
//									 int GrayThresholds) {
//	//����Ա����M_AncillaryPatterns�Ƿ���и���ͼ��
//	CV_Assert(M_AncillaryPatterns.size() >= static_cast<size_t>(ExposureNums));
//	//����ƽ����λ���õļ������󣬴�Ų�ͬͼ������ͬ����λ��Ӧ����Ĵ���
//	cv::Mat PhasesIntegrationMask = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
//
//	//�����ع����������ѭ��
//	for (int n = 0; n < ExposureNums; n++) {
//		std::cout << n << std::endl;
//		for (int i = 0; i < M_Height; i++) {
//			for (int j = 0; j < M_Width; j++) {
//				if (isnan(SovelPacelPhases[n].at<float>(i, j))){
//					continue;
//				}
//				else if(M_AncillaryPatterns[n].at<float>(i, j) < GrayThresholds && M_AncillaryPatterns[n].at<float>(i, j) > 0.5){
//					IntegrationPhase.at<float>(i, j) += SovelPacelPhases[n].at<float>(i, j);
//					PhasesIntegrationMask.at<float>(i, j) += 1;
//				}
//			}
//		}
//
//	}
//	IntegrationPhase /= PhasesIntegrationMask;
//};




cv::Mat SolvingPacel::readCameraInfoAndPolarrectification(const std::string& filename,
											cv::Mat& unwrappedLeft,
											cv::Mat& unwrappedRight,
											cv::Mat& jiaozhengL,
											cv::Mat& jiaozhengR,
											bool isShowPolarrectification) {
	cv::FileStorage FSget(filename, cv::FileStorage::READ);
	cv::Mat MatrixL;
	cv::Mat MatrixR;
	cv::Mat distL;
	cv::Mat distR;
	cv::Mat R1;
	cv::Mat R2;
	cv::Mat P1;
	cv::Mat P2;
	cv::Mat Q;

	FSget["MatrixL"] >> MatrixL;
	FSget["MatrixR"] >> MatrixR;
	FSget["distL"] >> distL;
	FSget["distR"] >> distR;
	FSget["R1"] >> R1;
	FSget["R2"] >> R2;
	FSget["P1"] >> P1;
	FSget["P2"] >> P2;
	FSget["Q"] >> Q;

	//����һ��ת������֪����û����
	R1.convertTo(R1, CV_32FC1);
	R2.convertTo(R2, CV_32FC1);
	P1.convertTo(P1, CV_32FC1);
	P2.convertTo(P2, CV_32FC1);

	FSget.release();

	cv::Mat M_map11, M_map12, M_map21, M_map22;
	initUndistortRectifyMap(MatrixL, distL, R1, P1, unwrappedLeft.size(), CV_32FC1, M_map11, M_map12);
	initUndistortRectifyMap(MatrixR, distR, R2, P2, unwrappedRight.size(), CV_32FC1, M_map21, M_map22);

	remap(unwrappedLeft, jiaozhengL, M_map11, M_map12, cv::INTER_LINEAR, cv::BORDER_CONSTANT);
	remap(unwrappedRight, jiaozhengR, M_map21, M_map22, cv::INTER_LINEAR, cv::BORDER_CONSTANT);

	if(isShowPolarrectification){
		cv::Mat horizontallyConcatenated;
		cv::hconcat(jiaozhengL, jiaozhengR, horizontallyConcatenated);

		cv::Mat temp;
		cv::cvtColor(horizontallyConcatenated, temp, cv::COLOR_GRAY2BGR);

		for (size_t i = 0; i < horizontallyConcatenated.rows; i += 100) {
			cv::line(temp, cv::Point(0, i), cv::Point(horizontallyConcatenated.cols - 1, i), cv::Scalar(255, 0, 0));
		}
		cv::namedWindow("jiaozheng", 2);
		cv::imshow("jiaozheng", temp);
		cv::waitKey(0);
	}
	return Q;
};




cv::Mat SolvingPacel::BinocularPhaseRecovery(cv::Mat& correction_left, 
											 cv::Mat& correction_right, 
										 	 int GrayThreshold,
											 bool isSaved) {
	//��ȡ����ͼƬ�ļ�
	M_GetImage.ReadPics_openmp(M_ImagesPathLeft, M_Patterns_Left, CV_32FC1, M_isRename);
	M_GetImage.ReadPics_openmp(M_ImagesPathRight, M_Patterns_Right, CV_32FC1 , M_isRename);
	//��ʼ�����еĵ��Կ��Բ��ó�ʼ���������еı���Ҫ��ʼ��
	cv::Mat Dephase_Left = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat Dephase_Right = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	/*std::vector<cv::Mat> Dephases_Left, Dephases_Right;*/

	std::vector<cv::Mat> Dephases_Left(
		M_ExposureNums,
		cv::Mat::zeros(M_Height, M_Width, CV_32FC1)  // ����Ԫ�صĳ�ʼֵ
	);
	std::vector<cv::Mat> Dephases_Right(
		M_ExposureNums, 
		cv::Mat::zeros(M_Height, M_Width, CV_32FC1)  // ����Ԫ�صĳ�ʼֵ
	);


	//��λ����
	if (M_IsMulExposure) {
		//�����ع�
		MultiDePhase(M_Patterns_Left, Dephases_Left, M_ExposureNums);
		PhasesIntegration(Dephases_Left, M_ExposureNums, Dephase_Left, GrayThreshold, Left);

		MultiDePhase(M_Patterns_Right, Dephases_Right, M_ExposureNums);
		PhasesIntegration(Dephases_Right, M_ExposureNums, Dephase_Right, GrayThreshold, Right);
	}
	else {
		//�����ع�
		SingleDePhase(M_Patterns_Left, Dephase_Left);
		SingleDePhase(M_Patterns_Right, Dephase_Right);
	}
	//���߽���
	cv::Mat Q = readCameraInfoAndPolarrectification(M_CalibrationResultPath,
													Dephase_Left, Dephase_Right,
													correction_left, correction_right);
	if (isSaved) {
		cv::imwrite("C:\\mycode\\download_source\\StructureLightAlgorithm\\RecoveryPahse\\correction_left.bmp", correction_left);
		cv::imwrite("C:\\mycode\\download_source\\StructureLightAlgorithm\\RecoveryPahse\\correction_right.bmp", correction_right);
		cv::imwrite("C:\\mycode\\download_source\\StructureLightAlgorithm\\RecoveryPahse\\Q.bmp", Q);
	}
	return Q;
};


//��ʵ����λ����ʱ��û��Ҫһֱ��ȡ������ͼ��ֻ��Ҫ��ȡ����ͼ�񼴿ɣ�����������㷨
bool SolvingPacel::BinocularPhaseRecovery(cv::Mat& correction_left, cv::Mat& correction_right, cv::Mat& Q){
	return true;
}



//û���ã���֧�ֱ��渡������ʽ
void SolvingPacel::readRecoveryPhaseImageAndQ(std::string& path, cv::Mat& correction_left, cv::Mat& correction_right, cv::Mat& Q) {
	correction_left = cv::imread(path + "correction_left.bmp", cv::IMREAD_UNCHANGED);
	correction_right = cv::imread(path + "correction_right.bmp", cv::IMREAD_UNCHANGED);
	Q = cv::imread(path + "Q.bmp", cv::IMREAD_UNCHANGED);
};
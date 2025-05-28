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


//初始版本
//void SolvingPacel::SingleDePhase(const std::vector<cv::Mat>& ProjectionPatterns, 
//								 cv::Mat& Phase, 
//								 int regulatory = 5) {
//	//用于在运行时检查条件是否为真，如果条件为假，则终止程序执行并打印错误信息。使用OpenCV中的宏
//	CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum));
//	//定义包裹相位,在读取图像时已经获得了rows和cols的值
//	cv::Mat wrappedphase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
//	//每次相移的值
//	const float shiftVal = static_cast<float>(CV_2PI) / M_PhaseImgsNum;
//
//	cv::Mat testPhase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
//	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range)
//		{
//			//存放所有相移图像每一行的起始地址，并存入到imgsPtrs容器中
//			std::vector<const float*> imgsPtrs(M_PhaseImgsNum);
//			for (int i = range.start; i < range.end; ++i)
//			{
//				for (int k = 0; k < M_PhaseImgsNum; ++k)
//				{
//					imgsPtrs[k] = ProjectionPatterns[M_GrayImgsNum +k].ptr<float>(i);
//				}
//				//存放包裹相位每一行的起始地址
//				auto wrappedphasePtr0 = wrappedphase.ptr<float>(i);
//				//进行列遍历
//				for (int j = 0; j < M_Width; ++j)
//				{
//					float molecules = 0.f, denominator = 0.f;
//					for (int k = 0; k < M_PhaseImgsNum; ++k)
//					{
//						molecules += imgsPtrs[k][j] * sin(k * shiftVal);
//						denominator += imgsPtrs[k][j] * cos(k * shiftVal);
//					}
//					(sqrt((molecules * molecules) + (denominator * denominator)) * 2 / M_PhaseImgsNum) > regulatory ? wrappedphasePtr0[j] = -atan2(molecules, denominator) : wrappedphasePtr0[j] = NAN;
//				}
//			}
//		});
//
//	cv::Mat A = cv::Mat::zeros(ProjectionPatterns[0].size(), CV_32FC1);
//	for (int i = 0; i < M_PhaseImgsNum; ++i)
//	{
//		A += ProjectionPatterns[M_GrayImgsNum +i] / M_PhaseImgsNum;
//	}
//
//	//解包裹
//	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range)
//		{
//			//存放格雷码图像的每一行的起始地址，起始大小为七个
//			std::vector<const float*> imgsPtrs1(M_GrayImgsNum);
//
//			for (int i = range.start; i < range.end; ++i)
//			{
//				for (int j = 0; j < M_GrayImgsNum; ++j)
//				{
//					imgsPtrs1[j] = ProjectionPatterns[j].ptr<float>(i);
//				}
//				//计算绝对相位、需要用到环境强度A，包裹相位wrappedphase，还是获取首行的地址
//				auto PhasePtr = testPhase.ptr<float>(i);
//				auto confidencePtr = A.ptr<float>(i);
//				auto wrappedPhasePtr = wrappedphase.ptr<float>(i);
//
//				for (int j = 0; j < M_Width; ++j)
//				{
//					int K1 = 0, tempVal = 0;
//					for (int k = 0; k < M_GrayImgsNum - 1; ++k)
//					{
//						//格雷码图像值大于环境强度A的值，得到0或者1，得到二进制的值
//						//tempVal在与返回的结果进行异或操作，相同为0，不同为1。
//						tempVal ^= imgsPtrs1[k][j] > confidencePtr[j];
//						// K1得到的就是二进制转换为十进制的数，下面的操作根据二进制的值转换为十进制
//						//最终得到的K1就是十进制的数
//						K1 = (K1 << 1) + tempVal;
//					}
//					//判断那张互补格雷码，还是进行上面的操作，但是生成的公式有点不一样
//					tempVal ^= imgsPtrs1[M_GrayImgsNum - 1][j] > confidencePtr[j];
//					const int K2 = ((K1 << 1) + tempVal + 1) / 2;
//
//					//最终得到的K1和K2的值都是十进制的数，可以进行相位展开的计算
//					if (wrappedPhasePtr[j] <= CV_PI / 2)
//					{
//						PhasePtr[j] = wrappedPhasePtr[j] + 2 * CV_PI * K2;
//					}
//					else if (wrappedPhasePtr[j] >= CV_PI * 3 / 2)
//					{
//						PhasePtr[j] = wrappedPhasePtr[j] + 2 * CV_PI * (K2 - 1);
//					}
//					else
//					{
//						PhasePtr[j] = wrappedPhasePtr[j] + 2 * CV_PI * K1;
//					}
//				}
//			}
//		});
//	Phase = testPhase.clone();
//	return;
//};




//初始版本 - 添加调制度存储，太长
//void SolvingPacel::SingleDePhase(const std::vector<cv::Mat>& ProjectionPatterns,
//								 cv::Mat& Phase,
//								 cv::Mat& FringeOrder,
//								 int regulatory = 5) {
//	// 使用OpenCV中的宏进行运行时检查
//	CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum));
//	CV_Assert(!ProjectionPatterns.empty() && !ProjectionPatterns[0].empty());
//
//	// 预分配内存，避免重复分配
//	cv::Mat wrappedPhase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);//包裹相位
//	cv::Mat absolutePhase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);//绝对相位
//	cv::Mat averageIntensity = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);//背景强度
//	cv::Mat modulation = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);    // 新增：调制度矩阵
//	cv::Mat fringeOrder = cv::Mat::zeros(M_Height, M_Width, CV_32SC1);  // 存储条纹级次K
//
//	// 每次相移的值：根据 M_PhaseImgsNum判断是几步相移
//	const float shiftVal = static_cast<float>(CV_2PI) / M_PhaseImgsNum;
//	const float invPhaseImgsNum = 1.0f / M_PhaseImgsNum;//相移步数的倒数
//	const float regulatoryThreshold = static_cast<float>(regulatory);
//
//	// 预计算三角函数值，避免重复计算
//	std::vector<float> sinVals(M_PhaseImgsNum), cosVals(M_PhaseImgsNum);
//	for (int k = 0; k < M_PhaseImgsNum; ++k) {
//		const float angle = k * shiftVal;
//		sinVals[k] = sin(angle);
//		cosVals[k] = cos(angle);
//	}
//
//	// 第一步：计算包裹相位、平均强度和调制度
//	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range) {
//		// 预分配指针数组，避免在内层循环中重复分配
//		std::vector<const float*> phaseImgPtrs(M_PhaseImgsNum);
//
//		for (int i = range.start; i < range.end; ++i) {
//			// 获取所有相移图像的行指针
//			for (int k = 0; k < M_PhaseImgsNum; ++k) {
//				phaseImgPtrs[k] = ProjectionPatterns[M_GrayImgsNum + k].ptr<float>(i);
//			}
//
//			float* wrappedPhasePtr = wrappedPhase.ptr<float>(i);
//			float* avgIntensityPtr = averageIntensity.ptr<float>(i);
//			float* modulationPtr = modulation.ptr<float>(i);  // 新增：调制度指针
//
//			// 列遍历
//			for (int j = 0; j < M_Width; ++j) {
//				float numerator = 0.0f, denominator = 0.0f, sumIntensity = 0.0f;
//
//				// 同时计算相位和平均强度
//				for (int k = 0; k < M_PhaseImgsNum; ++k) {
//					const float intensity = phaseImgPtrs[k][j];
//					numerator += intensity * sinVals[k];
//					denominator += intensity * cosVals[k];
//					sumIntensity += intensity;
//				}
//
//				// 计算调制度并存储
//				const float modulationValue = sqrt(numerator * numerator + denominator * denominator) * 2 * invPhaseImgsNum;
//				modulationPtr[j] = modulationValue;  // 存储调制度值
//				avgIntensityPtr[j] = sumIntensity * invPhaseImgsNum;
//
//				// 判断是否有效并计算相位
//				if (modulationValue > regulatoryThreshold) {
//					wrappedPhasePtr[j] = -atan2(numerator, denominator);
//				}
//				else {
//					wrappedPhasePtr[j] = NAN;
//				}
//			}
//		}
//		});
//
//	// 第二步：相位解包裹并记录条纹级次
//	int minK = INT_MAX, maxK = INT_MIN;
//	std::mutex minMaxMutex;  // 用于线程安全更新最值
//
//	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range) {
//		// 预分配格雷码图像指针数组
//		std::vector<const float*> grayImgPtrs(M_GrayImgsNum);
//		int localMinK = INT_MAX, localMaxK = INT_MIN;
//
//		for (int i = range.start; i < range.end; ++i) {
//			// 获取所有格雷码图像的行指针
//			for (int j = 0; j < M_GrayImgsNum; ++j) {
//				grayImgPtrs[j] = ProjectionPatterns[j].ptr<float>(i);
//			}
//
//			float* absolutePhasePtr = absolutePhase.ptr<float>(i);
//			int* fringeOrderPtr = fringeOrder.ptr<int>(i);
//			const float* avgIntensityPtr = averageIntensity.ptr<float>(i);
//			const float* wrappedPhasePtr = wrappedPhase.ptr<float>(i);
//
//			for (int j = 0; j < M_Width; ++j) {
//				// 跳过无效像素
//				if (std::isnan(wrappedPhasePtr[j])) {
//					absolutePhasePtr[j] = NAN;
//					fringeOrderPtr[j] = -1;  // 无效像素标记为-1
//					continue;
//				}
//
//				int K1 = 0, tempVal = 0;
//				const float threshold = avgIntensityPtr[j];
//
//				// 解码格雷码获取K1，初始值为零，因为与零异或保持不变
//				for (int k = 0; k < M_GrayImgsNum - 1; ++k) {
//					tempVal ^= (grayImgPtrs[k][j] > threshold) ? 1 : 0;
//					K1 = (K1 << 1) + tempVal;
//				}
//
//				// 处理互补格雷码获取K2
//				tempVal ^= (grayImgPtrs[M_GrayImgsNum - 1][j] > threshold) ? 1 : 0;
//				const int K2 = ((K1 << 1) + tempVal + 1) >> 1; // 使用位移代替除法
//
//				// 相位展开并记录使用的K值
//				const float wrappedVal = wrappedPhasePtr[j];
//				int finalK;
//				if (wrappedVal <= CV_PI * 0.5f) {
//					absolutePhasePtr[j] = wrappedVal + CV_2PI * K2;
//					finalK = K2;
//				}
//				else if (wrappedVal >= CV_PI * 1.5f) {
//					absolutePhasePtr[j] = wrappedVal + CV_2PI * (K2 - 1);
//					finalK = K2 - 1;
//				}
//				else {
//					absolutePhasePtr[j] = wrappedVal + CV_2PI * K1;
//					finalK = K1;
//				}
//
//				fringeOrderPtr[j] = finalK;
//				localMinK = std::min(localMinK, finalK);
//				localMaxK = std::max(localMaxK, finalK);
//			}
//		}
//		// 线程安全地更新全局最值
//		std::lock_guard<std::mutex> lock(minMaxMutex);
//		if (localMinK != INT_MAX) minK = std::min(minK, localMinK);
//		if (localMaxK != INT_MIN) maxK = std::max(maxK, localMaxK);
//		});
//
//	//条纹级次
//	FringeOrder = std::move(fringeOrder);
//	//输出绝对相位
//	Phase = std::move(absolutePhase);
//
//	//检测相位误差
//	cv::Mat phaseError = PhaseErrorAnalyzer::calculateTheoreticalPhaseError(modulation, averageIntensity, 1.0f);
//	// 根据相位误差进行质量分级
//	cv::Mat qualityMask = cv::Mat::zeros(M_Height, M_Width, CV_8UC1);
//	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range) {
//		for (int i = range.start; i < range.end; ++i) {
//			const float* errorPtr = phaseError.ptr<float>(i);
//			uchar* qualityPtr = qualityMask.ptr<uchar>(i);
//
//			for (int j = 0; j < M_Width; ++j) {
//				if (errorPtr[j] < 0.1) {
//					qualityPtr[j] = 255; // 高质量
//				}
//				else if (errorPtr[j] < 0.5) {
//					qualityPtr[j] = 128; // 中等质量
//				}
//				else {
//					qualityPtr[j] = 0;   // 低质量
//				}
//			}
//		}
//	});
//}



// 性能优化版本：减少函数调用，保持数据局部性
void SolvingPacel::SingleDePhase(const std::vector<cv::Mat>& ProjectionPatterns,
								 cv::Mat& Phase, 
								 cv::Mat& FringeOrder, 
								 int regulatory = 5) {
	// 输入验证
	CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum));
	CV_Assert(!ProjectionPatterns.empty() && !ProjectionPatterns[0].empty());

	// 初始化矩阵
	cv::Mat wrappedPhase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat absolutePhase = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat averageIntensity = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat modulation = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat fringeOrder = cv::Mat::zeros(M_Height, M_Width, CV_32SC1);

	const float regulatoryThreshold = static_cast<float>(regulatory);
	const TrigValues trig(M_PhaseImgsNum);

	// 合并两个主要步骤到一个并行循环中，保持数据局部性
	cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range) {
		// 预分配指针数组
		std::vector<const float*> phaseImgPtrs(M_PhaseImgsNum);
		std::vector<const float*> grayImgPtrs(M_GrayImgsNum);

		for (int i = range.start; i < range.end; ++i) {
			// 获取所有图像的行指针
			for (int k = 0; k < M_PhaseImgsNum; ++k) {
				phaseImgPtrs[k] = ProjectionPatterns[M_GrayImgsNum + k].ptr<float>(i);
			}
			for (int j = 0; j < M_GrayImgsNum; ++j) {
				grayImgPtrs[j] = ProjectionPatterns[j].ptr<float>(i);
			}

			// 获取输出矩阵指针
			float* wrappedPtr = wrappedPhase.ptr<float>(i);
			float* absPhasePtr = absolutePhase.ptr<float>(i);
			float* avgPtr = averageIntensity.ptr<float>(i);
			float* modPtr = modulation.ptr<float>(i);
			int* fringeOrderPtr = fringeOrder.ptr<int>(i);

			for (int j = 0; j < M_Width; ++j) {
				// 步骤1：计算包裹相位、调制度和平均强度（内联处理）
				float numerator = 0.0f, denominator = 0.0f, sumIntensity = 0.0f;
				for (int k = 0; k < M_PhaseImgsNum; ++k) {
					const float intensity = phaseImgPtrs[k][j];
					numerator += intensity * trig.sinVals[k];
					denominator += intensity * trig.cosVals[k];
					sumIntensity += intensity;
				}

				const float modulationValue = sqrt(numerator * numerator + denominator * denominator) * 2 * trig.invPhaseImgsNum;
				modPtr[j] = modulationValue;
				avgPtr[j] = sumIntensity * trig.invPhaseImgsNum;

				if (modulationValue > regulatoryThreshold) {
					wrappedPtr[j] = -atan2(numerator, denominator);

					// 步骤2：相位展开（内联处理）
					const float threshold = avgPtr[j];
					auto [K1, K2] = decodeGrayCodeInline(grayImgPtrs, j, threshold);

					auto [unwrappedPhase, finalK] = unwrapPhaseInline(wrappedPtr[j], K1, K2);
					absPhasePtr[j] = unwrappedPhase;
					fringeOrderPtr[j] = finalK;
				}
				else {
					wrappedPtr[j] = NAN;
					absPhasePtr[j] = NAN;
					fringeOrderPtr[j] = -1;
				}
			}
		}
		});

	// 输出结果
	Phase = std::move(absolutePhase);
	FringeOrder = std::move(fringeOrder);

	std::vector<cv::Mat> Phases(ProjectionPatterns.begin() + M_GrayImgsNum - 1, ProjectionPatterns.end() - M_AncillaryPatternNum);
	std::cout << Phases.size() << std::endl;
	cv::Mat PhaseError = PhaseErrorAnalyzer::calculateComprehensivePhaseError(Phases,
																			  wrappedPhase,
																			  modulation,
																			  averageIntensity,
																			  M_PhaseImgsNum);
}



void SolvingPacel::MultiDePhase(const std::vector<cv::Mat>& ProjectionPatterns,
								 std::vector<cv::Mat>& SovelPacelPhases,
								 std::vector<cv::Mat>& StripesLevel,
								 int ExposureNums, 
								 int regulatory = 5) {
	cv::Mat FringeOrderGray = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);

	if (M_AncillaryPatternNum == 0) {
		//未使用辅助图像
		CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum)* ExposureNums);
		int OnceProjectionPatternNum = M_GrayImgsNum + M_PhaseImgsNum;
		for (int ExposureNum = 0; ExposureNum < ExposureNums; ExposureNum++) {
			std::vector<cv::Mat>OnceExposure(ProjectionPatterns.begin() + OnceProjectionPatternNum * ExposureNum, ProjectionPatterns.begin() + OnceProjectionPatternNum * (ExposureNum + 1) - 1);
			SingleDePhase(OnceExposure, SovelPacelPhases[ExposureNum], FringeOrderGray);
			StripesLevel.push_back(FringeOrderGray);
		}
	}
	else if (M_AncillaryPatternNum == 1) {
		//使用一张辅助图像
		CV_Assert(ProjectionPatterns.size() >= static_cast<size_t>(M_GrayImgsNum + M_PhaseImgsNum + M_AncillaryPatternNum) * ExposureNums);
		int OnceProjectionPatternNum = M_GrayImgsNum + M_PhaseImgsNum + M_AncillaryPatternNum;
		for (int ExposureNum = 0; ExposureNum < ExposureNums; ExposureNum++) {
			std::vector<cv::Mat>OnceExposure(ProjectionPatterns.begin() + OnceProjectionPatternNum * ExposureNum, ProjectionPatterns.begin() + OnceProjectionPatternNum * (ExposureNum + 1) - M_AncillaryPatternNum);//使用这种方法构造时，包含左侧不包含右侧
			SingleDePhase(OnceExposure, SovelPacelPhases[ExposureNum], FringeOrderGray);
			StripesLevel.push_back(FringeOrderGray);
			//存储辅助图案
			M_AncillaryPatterns.push_back(ProjectionPatterns[OnceProjectionPatternNum * (ExposureNum+1) - 1]);
		}
	}
	return;
};


//改进版本
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
	//检测成员变量M_AncillaryPatterns是否存有辅助图像
	CV_Assert(AncillaryPatterns.size() >= static_cast<size_t>(ExposureNums));
	//计算平均相位所用的计数矩阵，存放不同图像中相同像素位置应计算的次数
	cv::Mat PhasesIntegrationMask = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat IntegrationMask = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);

	//根据曝光次数来进行循环
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



//初始版本
//void SolvingPacel::PhasesIntegration(std::vector<cv::Mat>& SovelPacelPhases,
//									 int ExposureNums,
//									 cv::Mat& IntegrationPhase,
//									 int GrayThresholds) {
//	//检测成员变量M_AncillaryPatterns是否存有辅助图像
//	CV_Assert(M_AncillaryPatterns.size() >= static_cast<size_t>(ExposureNums));
//	//计算平均相位所用的计数矩阵，存放不同图像中相同像素位置应计算的次数
//	cv::Mat PhasesIntegrationMask = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
//
//	//根据曝光次数来进行循环
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

	//新增一步转换，不知道有没有用
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
											 std::vector<cv::Mat>& StripesLevel_Left,
											 std::vector<cv::Mat>& StripesLevel_Right,
										 	 int GrayThreshold,
											 bool isSaved) {
	//获取所有图片文件
	M_GetImage.ReadPics_openmp(M_ImagesPathLeft, M_Patterns_Left, CV_32FC1, M_isRename);
	M_GetImage.ReadPics_openmp(M_ImagesPathRight, M_Patterns_Right, CV_32FC1 , M_isRename);
	//初始化，有的电脑可以不用初始化，但是有的必须要初始化
	cv::Mat Dephase_Left = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat Dephase_Right = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);
	cv::Mat FringeOrderGray = cv::Mat::zeros(M_Height, M_Width, CV_32FC1);


	std::vector<cv::Mat> Dephases_Left(
		M_ExposureNums,
		cv::Mat::zeros(M_Height, M_Width, CV_32FC1)  // 所有元素的初始值
	);
	std::vector<cv::Mat> Dephases_Right(
		M_ExposureNums, 
		cv::Mat::zeros(M_Height, M_Width, CV_32FC1)  // 所有元素的初始值
	);


	//相位解算
	if (M_IsMulExposure) {
		//多重曝光
		MultiDePhase(M_Patterns_Left, Dephases_Left, StripesLevel_Left, M_ExposureNums);
		PhasesIntegration(Dephases_Left, M_ExposureNums, Dephase_Left, GrayThreshold, Left);

		MultiDePhase(M_Patterns_Right, Dephases_Right, StripesLevel_Right, M_ExposureNums);
		PhasesIntegration(Dephases_Right, M_ExposureNums, Dephase_Right, GrayThreshold, Right);
	}
	else {
		//单次曝光
		SingleDePhase(M_Patterns_Left, Dephase_Left, FringeOrderGray);
		SingleDePhase(M_Patterns_Right, Dephase_Right, FringeOrderGray);
	}
	//极线矫正
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


//其实在相位解相时，没必要一直获取格雷码图像，只需要获取相移图像即可，完善上面的算法
bool SolvingPacel::BinocularPhaseRecovery(cv::Mat& correction_left, cv::Mat& correction_right, cv::Mat& Q){
	return true;
}



//没有用，不支持保存浮点数格式
void SolvingPacel::readRecoveryPhaseImageAndQ(std::string& path, cv::Mat& correction_left, cv::Mat& correction_right, cv::Mat& Q) {
	correction_left = cv::imread(path + "correction_left.bmp", cv::IMREAD_UNCHANGED);
	correction_right = cv::imread(path + "correction_right.bmp", cv::IMREAD_UNCHANGED);
	Q = cv::imread(path + "Q.bmp", cv::IMREAD_UNCHANGED);
};
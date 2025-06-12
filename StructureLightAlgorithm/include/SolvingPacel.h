#pragma once
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <math.h>
#include <filesystem>
#include <string>
#include <algorithm>
#include <functional>
#include "FileManipulation.h"
#include "PhaseErrorAnalyze.h"


#define pi CV_PI

enum CameraSolutioin {
	A5131 = 0,
	A5201 = 1,
	A7300 = 2,
	A7500 = 3,
};


enum CameraPosition {
	Left = 0,
	Right = 1,
};



enum SolveMethod {

	DPWC = 0,    //��Ƶ��
	HBGLM = 1,   //����������


};





/*
 * @ brief  ����Ӧ��ֵ���㷽��ö��
 * @ OTSU                Otsu����Ӧ��ֵ��������䷽�����
 * @ LOCAL_MEAN          �ֲ���ֵ��ֵ��ʹ���������ؾ�ֵ
 * @ ADAPTIVE_GAUSSIAN   ��˹��Ȩ����Ӧ��ֵ��ʹ�ø�˹Ȩ�ؼ���������ֵ
 * @ HYBRID              ��Ϸ��������ݾֲ���������Ӧ����Ȩ��
 */
enum class ThresholdMethod {
	OTSU,
	LOCAL_MEAN,
	ADAPTIVE_GAUSSIAN,
	HYBRID
};


class SolvingPacel {
public:
	/*
	 * @ brief  ��λ��Ϣ�����๹�캯��������ڵ����ع�
	 * @ param  ImagepathLeft             �������������ͶӰͼ����·��
	 * @ param  ImagepathRight            �������������ͶӰͼ����·��
	 * @ param  CalibrationResultPath     �������궨���
	 * @ param  GrayImgsNum               ������ͼ�������
	 * @ param  PhaseImgsNum              ����ͼ������
	 * @ param  width                     ͼ��Ŀ��
	 * @ param  height                    ͼ��ĸ߶�
	 * @ param  AncillaryPatternNum       ͶӰ�ĸ���ͼƬ������
	 * @ param  isRename                  �Ƿ���Ҫ������ͼƬ
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
				 const std::string& CalibrationResultPath,
				 int GrayImgsNum,
				 int PhaseImgsNum,
				 int width,
				 int height,
				 int AncillaryPatternNum,
				 bool isRename);

	/*
	 * @ brief  ��λ��Ϣ�����๹�캯��������ڵ����ع�
	 * @ param  imagepath            ��������ͶӰͼ����·��
	 * @ param  GrayImgsNum          ������ͼ�������
	 * @ param  PhaseImgsNum         ����ͼ������
	 * @ param  Form                 ������ͺţ��ɴ�������ͼ��Ŀ�Ⱥ͸߶�
	 * @ param  AncillaryPatternNum  ͶӰ�ĸ���ͼƬ������
 	 * @ param  isRename             �Ƿ���Ҫ������ͼƬ
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
		         const std::string& CalibrationResultPath,
				 int NumOfGray,
				 int NumOfPhase,
				 enum CameraSolutioin Form,
				 int AncillaryPatternNum,
				 bool isRename);

	/*
	 * @ brief  ��λ��Ϣ�����๹�캯������������ع�
	 * @ param  ImagepathLeft             �������������ͶӰͼ����·��
	 * @ param  ImagepathRight            �������������ͶӰͼ����·��
	 * @ param  CalibrationResultPath     �������궨���
	 * @ param  GrayImgsNum               ������ͼ�������
	 * @ param  PhaseImgsNum              ����ͼ������
	 * @ param  width                     ͼ��Ŀ��
	 * @ param  height                    ͼ��ĸ߶�
	 * @ param  ExposureNums              �����ع���ع����
	 * @ param  AncillaryPatternNum       ͶӰ�ĸ���ͼƬ������
	 * @ param  isRename                  �Ƿ���Ҫ������ͼƬ
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
		         const std::string& CalibrationResultPath,
				 int GrayImgsNum,
				 int PhaseImgsNum,
				 int width,
				 int height,
				 int ExposureNums,
				 int AncillaryPatternNum,
		         bool isRename);

	/*
	 * @ brief  ��λ��Ϣ�����๹�캯������������ع�
	 * @ param  ImagepathLeft             �������������ͶӰͼ����·��
	 * @ param  ImagepathRight            �������������ͶӰͼ����·��
	 * @ param  CalibrationResultPath     �������궨���
	 * @ param  GrayImgsNum               ������ͼ�������
	 * @ param  PhaseImgsNum              ����ͼ������
	 * @ param  Form                      ������ͺţ��ɴ�������ͼ��Ŀ�Ⱥ͸߶�
	 * @ param  ExposureNums              �����ع���ع����
	 * @ param  AncillaryPatternNum       ͶӰ�ĸ���ͼƬ������
	 * @ param  isRename                  �Ƿ���Ҫ������ͼƬ
	*/
	SolvingPacel(const std::string& ImagepathLeft,
				 const std::string& ImagepathRight,
		         const std::string& CalibrationResultPath,
				 int NumOfGray,
				 int NumOfPhase,
				 enum CameraSolutioin Form,
				 int ExposureNums,
				 int AncillaryPatternNum,
				 bool isRename);



	/*
	 * @ brief  �Ե���ͶӰͼ��������λ���㣬������λչ�����������Ż��汾��
	 * @ param  ProjectionPatterns  ͶӰͼ�����飬����������ͼ����λͼ
	 * @ param  Phase               ����ľ�����λ����
	 * @ param  FringeOrder         ��������Ƽ��ξ���
	 * @ param  regulatory          ���ƶ���ֵ��Ĭ��Ϊ5
	 * @ param  thresholdMethod     ����Ӧ��ֵ���㷽����Ĭ��Ϊ��Ϸ���
	 */
	void SingleDePhase(const std::vector<cv::Mat>& ProjectionPatterns,
									 cv::Mat& Phase,
									 cv::Mat& FringeOrder,
									 int regulatory,
									 ThresholdMethod thresholdMethod);
	 

	/*
	 * @ brief  �Զ���ͶӰͼ��������λ����,�ú���ֱ�ӵ��õ��ν��㺯����ÿ���ع��ͶӰͼ��ֱ�ӽ��໥������
	 * @ param  ProjectionPatterns   ����ͶӰͼ��
	 * @ param  SovelPacelPhases     ������Ľ������λ����
	 * @ param  StripesLevel         ���Ƽ���
	 * @ param  ExposureNums         �����ع���ع����
	 * @ param  regulatory           ���ƶ�B��ֵ
	*/
	void MultiDePhase(const std::vector<cv::Mat>& ProjectionPatterns, 
					  std::vector<cv::Mat>& SovelPacelPhases, 
					  std::vector<cv::Mat>& StripesLevel,
					  int ExposureNums, 
					  int regulatory);


	/*
	 * @ brief  �Զ����ع�Ľ������λ�����ں�
	 * @ param  SovelPacelPhases     ������Ľ������λ����
	 * @ param  ExposureNums         �����ع���ع����
	 * @ param  IntegrationPhase     �ں�֮�����λ
	 * @ param  GrayThresholds       �Ҷ���ֵ
	 * @ param  position             ��������������������
	*/
	void PhasesIntegration(std::vector<cv::Mat>& SovelPacelPhases, 
						   int ExposureNums,
						   cv::Mat& IntegrationPhase,
						   int GrayThresholds,
		                   CameraPosition position);


	/*
	 * @ brief  ���߽�������
	 * @ param  filename                  �������궨�������ļ�·��
	 * @ param  unwrappedLeft             �������������Ľ������λ
	 * @ param  unwrappedRight            �������������Ľ������λ
	 * @ param  jiaozhengL                ���߽������������������λ
	 * @ param  jiaozhengR                ���߽������������������λ
	 * @ param  isShowPolarrectification  �Ƿ���ʾ���߽����Ľ��
	 * @ return  ˫Ŀ�궨��Q����
	*/
	cv::Mat readCameraInfoAndPolarrectification(const std::string& filename,
												cv::Mat& unwrappedLeft,
												cv::Mat& unwrappedRight,
												cv::Mat& jiaozhengL,
												cv::Mat& jiaozhengR,
												bool isShowPolarrectification = false);

	/*
	 * @ brief  ˫Ŀ�ṹ����λ�ָ�
	 * @ param  correction_left      ���߽�����������չ����λ
	 * @ param  correction_right     ���߽�����������չ����λ
	 * @ param  StripesLevel_Left    ��������Ƽ���
	 * @ param  StripesLevel_Right   ��������Ƽ���
	 * @ param  GrayThreshold        �Ҷ���ֵ
	 * @ param  isSaved              �Ƿ񱣴漫��У�����ͼ��
	 * @ return  ˫Ŀ�궨��Q����
	*/
	cv::Mat BinocularPhaseRecovery(cv::Mat& correction_left, 
								   cv::Mat& correction_right, 
								   std::vector<cv::Mat>& StripesLevel_Left,
								   std::vector<cv::Mat>& StripesLevel_Right,
								   int GrayThreshold, bool isSaved = false);


	/*
	 * @ brief  ˫Ŀ�ṹ����λ�ָ�,������������ظ����������ͼ��ֻ����һ���ع��µĸ�����ͼ�񼴿�
	 * @ param  correction_left      ���߽�����������չ����λ
	 * @ param  correction_right     ���߽�����������չ����λ
	 * @ param  Q                     ˫Ŀ�궨��Q����
	*/
	bool BinocularPhaseRecovery(cv::Mat& correction_left, cv::Mat& correction_right ,cv::Mat& Q);



	/*
	 * @ brief  ��ȡ�Ѿ������������������߽������ͼ���Լ�Q����
	 * @ param  path                 ���ͼƬ��·��
	 * @ param  correction_left      ���߽�����������չ����λ
	 * @ param  correction_right     ���߽�����������չ����λ
	 * @ param  Q                     ˫Ŀ�궨��Q����
	*/
	void readRecoveryPhaseImageAndQ(std::string& path, cv::Mat& correction_left, cv::Mat& correction_right, cv::Mat& Q);





	//Binarise the image Change the pixel value to 0 or 255
	std::vector<cv::Mat> changeValue(std::vector<cv::Mat>& grays,
									 cv::Mat& threshold,
									 int width,
									 int height);


	//��float����ת��Ϊdouble������˵
	cv::Mat Get_wrappedphase_double();
	cv::Mat getunwtapped_double();
	cv::Mat RemoveBackground_double();


private:
	/*
	 * @ brief  ���º���������λչ������
	 * @ struct    TrigValues                         Ԥ�������Ǻ���ֵ�Ľṹ��
	 * @ function  generateQualityMask                ������������
	 * @ function  decodeGrayCodeInline               ����������ȡ���Ƽ���
	 * @ function  unwrapPhaseInline                  ��λչ��
	*/

	// Ԥ�������Ǻ���ֵ�Ľṹ��
	struct TrigValues {
		std::vector<float> sinVals, cosVals;
		float shiftVal, invPhaseImgsNum;

		TrigValues(int phaseImgsNum) {
			shiftVal = static_cast<float>(CV_2PI) / phaseImgsNum;
			invPhaseImgsNum = 1.0f / phaseImgsNum;
			sinVals.resize(phaseImgsNum);
			cosVals.resize(phaseImgsNum);

			for (int k = 0; k < phaseImgsNum; ++k) {
				const float angle = k * shiftVal;
				sinVals[k] = sin(angle);
				cosVals[k] = cos(angle);
			}
		}
	};

	// ������������
	cv::Mat generateQualityMask(const cv::Mat& modulation, const cv::Mat& averageIntensity) {
		cv::Mat phaseError = PhaseErrorAnalyzer::calculateTheoreticalPhaseError(modulation, averageIntensity, 1.0f);
		cv::Mat qualityMask = cv::Mat::zeros(M_Height, M_Width, CV_8UC1);

		cv::parallel_for_(cv::Range(0, M_Height), [&](const cv::Range& range) {
			for (int i = range.start; i < range.end; ++i) {
				const float* errorPtr = phaseError.ptr<float>(i);
				uchar* qualityPtr = qualityMask.ptr<uchar>(i);

				for (int j = 0; j < M_Width; ++j) {
					qualityPtr[j] = (errorPtr[j] < 0.1) ? 255 :
						(errorPtr[j] < 0.5) ? 128 : 0;
				}
			}
			});

		return qualityMask;
	}

	// �����������������ٺ������ÿ���
	inline std::pair<int, int> decodeGrayCodeInline(const std::vector<const float*>& grayImgPtrs, 
													int j,
													float threshold) {
		int K1 = 0, tempVal = 0;
		// ����������ȡK1����ʼֵΪ�㣬��Ϊ������򱣳ֲ���
		for (int k = 0; k < M_GrayImgsNum - 1; ++k) {
			tempVal ^= (grayImgPtrs[k][j] > threshold) ? 1 : 0;
			K1 = (K1 << 1) + tempVal;
		}
		// �������������ȡK2
		tempVal ^= (grayImgPtrs[M_GrayImgsNum - 1][j] > threshold) ? 1 : 0;
		const int K2 = ((K1 << 1) + tempVal + 1) >> 1; // ʹ��λ�ƴ������
		return { K1, K2 };
	}

	inline std::pair<float, int> unwrapPhaseInline(float wrappedPhase, int K1, int K2) {
		if (wrappedPhase <= CV_PI * 0.5f) {
			return { wrappedPhase + CV_2PI * K2, K2 };
		}
		else if (wrappedPhase >= CV_PI * 1.5f) {
			return { wrappedPhase + CV_2PI * (K2 - 1), K2 - 1 };
		}
		else {
			return { wrappedPhase + CV_2PI * K1, K1 };
		}
	}


	/*
	 * @ brief     ���º������ڸ�����Ķ�ֵ��
	 * @ struct    calculateOtsuThreshold                 ����Otsu����Ӧ��ֵ����Ե���������Χ����
	 * @ function  calculateLocalAdaptiveThreshold        ����ֲ�����Ӧ��ֵ
	 * @ function  calculateGaussianAdaptiveThreshold     �����˹��Ȩ����Ӧ��ֵ
	 * @ function  calculateHybridThreshold               ����������Ӧ��ֵ����
	 * @ function  decodeGrayCodeInlineOptimized          �Ż��ĸ�������뺯����ʹ������Ӧ��ֵ
	*/


	/*
	 * @ brief  ����Otsu����Ӧ��ֵ����Ե���������Χ����
	 * @ param  grayImgPtrs     ������ͼ����ָ������
	 * @ param  centerJ         �������ص�������
	 * @ param  width           ͼ����
	 * @ param  windowSize      ���㴰�ڴ�С��Ĭ��Ϊ5
	 * @ return float           ����õ���Otsu��ֵ
	 */
	inline float calculateOtsuThreshold(const std::vector<const float*>& grayImgPtrs,
										int centerJ, 
										int width, 
										int windowSize = 5) {
		const int halfWindow = windowSize / 2;
		const int startJ = std::max(0, centerJ - halfWindow);
		const int endJ = std::min(width - 1, centerJ + halfWindow);

		// �ռ����������и�����ͼ�������ֵ
		std::vector<float> pixels;
		pixels.reserve((endJ - startJ + 1) * grayImgPtrs.size());

		for (size_t imgIdx = 0; imgIdx < grayImgPtrs.size(); ++imgIdx) {
			for (int j = startJ; j <= endJ; ++j) {
				pixels.push_back(grayImgPtrs[imgIdx][j]);
			}
		}

		if (pixels.empty()) return 128.0f;

		// �򻯵�Otsu�㷨
		std::sort(pixels.begin(), pixels.end());

		float bestThreshold = pixels[pixels.size() / 2]; // Ĭ����λ��
		float maxVariance = 0.0f;

		const size_t totalPixels = pixels.size();

		// ������ֵ��Χ�����������ֵ
		for (size_t t = 1; t < totalPixels - 1; ++t) {
			const float threshold = pixels[t];

			// ����ǰ���ͱ�����Ȩ�غ;�ֵ
			float w0 = static_cast<float>(t) / totalPixels;
			float w1 = 1.0f - w0;

			if (w0 == 0 || w1 == 0) continue;

			float sum0 = 0, sum1 = 0;
			for (size_t i = 0; i < t; ++i) sum0 += pixels[i];
			for (size_t i = t; i < totalPixels; ++i) sum1 += pixels[i];

			float mean0 = sum0 / t;
			float mean1 = sum1 / (totalPixels - t);

			// ��䷽��
			float variance = w0 * w1 * (mean0 - mean1) * (mean0 - mean1);

			if (variance > maxVariance) {
				maxVariance = variance;
				bestThreshold = threshold;
			}
		}
		return bestThreshold;
	}


	/*
	 * @ brief  ����ֲ�����Ӧ��ֵ
	 * @ param  grayImgPtrs     ������ͼ����ָ������
	 * @ param  centerJ         �������ص�������
	 * @ param  width           ͼ����
	 * @ param  windowSize      ���㴰�ڴ�С��Ĭ��Ϊ7
	 * @ return float           ����õ��ľֲ���ֵ��ֵ
	 */
	inline float calculateLocalAdaptiveThreshold(const std::vector<const float*>& grayImgPtrs,
												 int centerJ, 
												 int width, 
												 int windowSize = 7) {
		const int halfWindow = windowSize / 2;
		const int startJ = std::max(0, centerJ - halfWindow);
		const int endJ = std::min(width - 1, centerJ + halfWindow);

		float sum = 0.0f;
		int count = 0;

		// ���㴰�������и�����ͼ���ƽ��ֵ
		for (size_t imgIdx = 0; imgIdx < grayImgPtrs.size(); ++imgIdx) {
			for (int j = startJ; j <= endJ; ++j) {
				sum += grayImgPtrs[imgIdx][j];
				count++;
			}
		}
		return count > 0 ? sum / count : 128.0f;
	}


	/*
	 * @ brief  �����˹��Ȩ����Ӧ��ֵ
	 * @ param  grayImgPtrs     ������ͼ����ָ������
	 * @ param  centerJ         �������ص�������
	 * @ param  width           ͼ����
	 * @ param  windowSize      ���㴰�ڴ�С��Ĭ��Ϊ7
	 * @ return float           ����õ��ĸ�˹��Ȩ��ֵ
	 */
	inline float calculateGaussianAdaptiveThreshold(const std::vector<const float*>& grayImgPtrs,
													int centerJ, 
													int width, 
													int windowSize = 7) {
		const int halfWindow = windowSize / 2;
		const int startJ = std::max(0, centerJ - halfWindow);
		const int endJ = std::min(width - 1, centerJ + halfWindow);

		// Ԥ�����˹Ȩ��
		std::vector<float> gaussWeights(windowSize);
		const float sigma = windowSize / 3.0f;
		float weightSum = 0.0f;

		for (int i = 0; i < windowSize; ++i) {
			float dist = (i - halfWindow) * (i - halfWindow);
			gaussWeights[i] = exp(-dist / (2 * sigma * sigma));
			weightSum += gaussWeights[i];
		}

		// ��һ��Ȩ��
		for (float& w : gaussWeights) {
			w /= weightSum;
		}

		float weightedSum = 0.0f;

		// �����Ȩƽ��ֵ
		for (size_t imgIdx = 0; imgIdx < grayImgPtrs.size(); ++imgIdx) {
			for (int j = startJ; j <= endJ; ++j) {
				int weightIdx = j - startJ;
				if (weightIdx < static_cast<int>(gaussWeights.size())) {
					weightedSum += grayImgPtrs[imgIdx][j] * gaussWeights[weightIdx];
				}
			}
		}

		return weightedSum;
	}


	/*
	 * @ brief  ����������Ӧ��ֵ����
	 * @ param  grayImgPtrs     ������ͼ����ָ������
	 * @ param  centerJ         �������ص�������
	 * @ param  width           ͼ����
	 * @ param  avgIntensity    ����ƽ��ǿ��ֵ
	 * @ return float           ����õ��Ļ����ֵ
	 */
	inline float calculateHybridThreshold(const std::vector<const float*>& grayImgPtrs,
										  int centerJ, 
										  int width, 
										  float avgIntensity) {
		// ��Ͼֲ�����Ӧ��ȫ��ƽ��ֵ
		float localThreshold = calculateLocalAdaptiveThreshold(grayImgPtrs, centerJ, width);
		float globalThreshold = avgIntensity;

		// ���ݾֲ��������Ȩ��
		float localVariance = 0.0f;
		int count = 0;
		const int windowSize = 5;
		const int halfWindow = windowSize / 2;
		const int startJ = std::max(0, centerJ - halfWindow);
		const int endJ = std::min(width - 1, centerJ + halfWindow);

		for (size_t imgIdx = 0; imgIdx < grayImgPtrs.size(); ++imgIdx) {
			for (int j = startJ; j <= endJ; ++j) {
				float diff = grayImgPtrs[imgIdx][j] - localThreshold;
				localVariance += diff * diff;
				count++;
			}
		}

		localVariance = count > 0 ? sqrt(localVariance / count) : 0.0f;

		// �����ʱ�������ֲ���ֵ������Сʱ������ȫ����ֵ
		float localWeight = std::min(1.0f, localVariance / 50.0f); // 50�Ǿ���ֵ���ɵ���
		return localWeight * localThreshold + (1.0f - localWeight) * globalThreshold;
	}


	/*
	 * @ brief  �Ż��ĸ�������뺯����ʹ������Ӧ��ֵ
	 * @ param  grayImgPtrs     ������ͼ����ָ������
	 * @ param  j               ��ǰ���ص�������
	 * @ param  width           ͼ����
	 * @ param  avgIntensity    ����ƽ��ǿ��ֵ
	 * @ param  method          ����Ӧ��ֵ���㷽����Ĭ��Ϊ��Ϸ���
	 * @ return std::pair<int, int>  ���ؽ���õ���K1��K2ֵ
	 */
	inline std::pair<int, int> decodeGrayCodeInlineOptimized(const std::vector<const float*>& grayImgPtrs,
															 int j, int width,
															 float avgIntensity,
															 ThresholdMethod method = ThresholdMethod::HYBRID) {
		// ����ѡ��ķ�����������Ӧ��ֵ
		float threshold;
		switch (method) {
		case ThresholdMethod::OTSU:
			threshold = calculateOtsuThreshold(grayImgPtrs, j, width);
			break;
		case ThresholdMethod::LOCAL_MEAN:
			threshold = calculateLocalAdaptiveThreshold(grayImgPtrs, j, width);
			break;
		case ThresholdMethod::ADAPTIVE_GAUSSIAN:
			threshold = calculateGaussianAdaptiveThreshold(grayImgPtrs, j, width);
			break;
		case ThresholdMethod::HYBRID:
		default:
			threshold = calculateHybridThreshold(grayImgPtrs, j, width, avgIntensity);
			break;
		}

		int K1 = 0, tempVal = 0;
		// ����������ȡK1
		for (int k = 0; k < static_cast<int>(grayImgPtrs.size()) - 1; ++k) {
			tempVal ^= (grayImgPtrs[k][j] > threshold) ? 1 : 0;
			K1 = (K1 << 1) + tempVal;
		}

		// �������������ȡK2
		tempVal ^= (grayImgPtrs[grayImgPtrs.size() - 1][j] > threshold) ? 1 : 0;
		const int K2 = ((K1 << 1) + tempVal + 1) >> 1;
		return { K1, K2 };
	}





private:
	//�Ƿ�ʹ�ö����ع��־λ
	bool M_IsMulExposure = false;
	//�����ع����
	int M_ExposureNums;

	//ͶӰ�ĸ���ͼƬ������
	int M_AncillaryPatternNum;

	//�����ͼƬ���·��
	std::string M_ImagesPathLeft;
	//�����ͼƬ���·��
	std::string M_ImagesPathRight;
	//����궨�������·��
	std::string M_CalibrationResultPath;

	//������ͼ�������
	int M_GrayImgsNum;
	//����ͼ�������
	int M_PhaseImgsNum;
	//ͼ��ĳߴ�
	int M_Width;
	int M_Height;

	//��������ͼ������
	std::vector<cv::Mat> M_Patterns_Left;
	//��������ͼ������
	std::vector<cv::Mat> M_Patterns_Right;


	std::vector<cv::Mat> M_AncillaryPatterns;


	//�������ͼ��
	std::vector<cv::Mat> M_Allp;
	//��Ŷ�ֵ�����ͼ��
	std::vector<cv::Mat> M_BinaryP;
	//��Ű�����λ
	cv::Mat M_pha;
	cv::Mat M_wrapped;
	//���K1����K2����
	cv::Mat M_K1;
	cv::Mat M_K2;
	//��Ž����ͼ��
	cv::Mat M_unwrapped;
	//���Ŷ�
	cv::Mat M_confidence;

	//�ļ�����
	FileManipulation M_GetImage;
	//�Ƿ���Ҫ������ͼƬ
	bool M_isRename;



};

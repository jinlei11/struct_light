#pragma once

// ��λ�����㷽������
#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>

class PhaseErrorAnalyzer {

public:
    // ����1�����ڵ��ƶȵ�������λ������
    static cv::Mat calculateTheoreticalPhaseError(const cv::Mat& modulation,
                                                  const cv::Mat& averageIntensity,
                                                  float noiseStd = 1.0f) {
        cv::Mat phaseError = cv::Mat::zeros(modulation.size(), CV_32FC1);

        cv::parallel_for_(cv::Range(0, modulation.rows), [&](const cv::Range& range) {
            for (int i = range.start; i < range.end; ++i) {
                const float* modPtr = modulation.ptr<float>(i);
                const float* avgPtr = averageIntensity.ptr<float>(i);
                float* errorPtr = phaseError.ptr<float>(i);

                for (int j = 0; j < modulation.cols; ++j) {
                    if (modPtr[j] > 0 && avgPtr[j] > 0) {
                        // ���۹�ʽ: ��_�� �� ��_I / (M * I_avg)
                        errorPtr[j] = noiseStd / (modPtr[j] * avgPtr[j]);
                    }
                    else {
                        errorPtr[j] = std::numeric_limits<float>::max();
                    }
                }
            }
         });

        return phaseError;
    }

    // ����2�����ڲв��������λ������
    static cv::Mat calculateResidualPhaseError(const std::vector<cv::Mat>& phaseImages,
                                               const cv::Mat& wrappedPhase,
                                               const cv::Mat& averageIntensity,
                                               int phaseSteps) {
        cv::Mat phaseError = cv::Mat::zeros(wrappedPhase.size(), CV_32FC1);
        const float shiftVal = static_cast<float>(CV_2PI) / phaseSteps;

        cv::parallel_for_(cv::Range(0, wrappedPhase.rows), [&](const cv::Range& range) {
            std::vector<const float*> phasePtrs(phaseSteps);

            for (int i = range.start; i < range.end; ++i) {
                // ��ȡ��������ͼ�����ָ��
                for (int k = 0; k < phaseSteps; ++k) {
                    phasePtrs[k] = phaseImages[k].ptr<float>(i);
                }

                const float* phasePtr = wrappedPhase.ptr<float>(i);
                const float* avgPtr = averageIntensity.ptr<float>(i);
                float* errorPtr = phaseError.ptr<float>(i);

                for (int j = 0; j < wrappedPhase.cols; ++j) {
                    if (std::isnan(phasePtr[j]) || avgPtr[j] <= 0) {
                        errorPtr[j] = std::numeric_limits<float>::max();
                        continue;
                    }

                    // �����ع����
                    float sumSquaredResidual = 0.0f;
                    const float phase = phasePtr[j];
                    const float avgIntensity = avgPtr[j];

                    for (int k = 0; k < phaseSteps; ++k) {
                        // ���۹�ǿֵ
                        float theoreticalIntensity = avgIntensity * (1.0f + cos(phase + k * shiftVal));
                        // ʵ�ʹ�ǿֵ
                        float actualIntensity = phasePtrs[k][j];
                        // �в�
                        float residual = actualIntensity - theoreticalIntensity;
                        sumSquaredResidual += residual * residual;
                    }

                    // RMS�в���Ϊ��λ���Ĺ���
                    errorPtr[j] = sqrt(sumSquaredResidual / phaseSteps);
                }
            }
            });

        return phaseError;
    }

    // ����3��������λһ���Ե�������
    static cv::Mat calculatePhaseConsistencyError(const cv::Mat& wrappedPhase,
                                                  const cv::Mat& modulation,
                                                  int windowSize = 5) {
        cv::Mat phaseError = cv::Mat::zeros(wrappedPhase.size(), CV_32FC1);
        const int halfWindow = windowSize / 2;

        cv::parallel_for_(cv::Range(halfWindow, wrappedPhase.rows - halfWindow),
            [&](const cv::Range& range) {
                for (int i = range.start; i < range.end; ++i) {
                    float* errorPtr = phaseError.ptr<float>(i);

                    for (int j = halfWindow; j < wrappedPhase.cols - halfWindow; ++j) {
                        const float centerPhase = wrappedPhase.at<float>(i, j);
                        const float centerMod = modulation.at<float>(i, j);

                        if (std::isnan(centerPhase) || centerMod <= 0) {
                            errorPtr[j] = std::numeric_limits<float>::max();
                            continue;
                        }

                        // ������������λ�ı�׼��
                        std::vector<float> neighborPhases;
                        float sumWeight = 0.0f, weightedSum = 0.0f, weightedSumSq = 0.0f;

                        for (int di = -halfWindow; di <= halfWindow; ++di) {
                            for (int dj = -halfWindow; dj <= halfWindow; ++dj) {
                                if (di == 0 && dj == 0) continue;

                                const float neighborPhase = wrappedPhase.at<float>(i + di, j + dj);
                                const float neighborMod = modulation.at<float>(i + di, j + dj);

                                if (!std::isnan(neighborPhase) && neighborMod > 0) {
                                    // ������λ����
                                    float phaseDiff = neighborPhase - centerPhase;
                                    if (phaseDiff > CV_PI) phaseDiff -= CV_2PI;
                                    if (phaseDiff < -CV_PI) phaseDiff += CV_2PI;

                                    float adjustedPhase = centerPhase + phaseDiff;
                                    float weight = neighborMod; // �õ��ƶ���ΪȨ��

                                    weightedSum += adjustedPhase * weight;
                                    weightedSumSq += adjustedPhase * adjustedPhase * weight;
                                    sumWeight += weight;
                                }
                            }
                        }

                        if (sumWeight > 0) {
                            float mean = weightedSum / sumWeight;
                            float variance = (weightedSumSq / sumWeight) - (mean * mean);
                            errorPtr[j] = sqrt(std::max(0.0f, variance));
                        }
                        else {
                            errorPtr[j] = std::numeric_limits<float>::max();
                        }
                    }
                }
            });

        return phaseError;
    }

    // ����4�����ڸ߽�г������λ�����
    static cv::Mat calculateHarmonicPhaseError(const std::vector<cv::Mat>& phaseImages,
        const cv::Mat& wrappedPhase,
        int phaseSteps) {
        cv::Mat phaseError = cv::Mat::zeros(wrappedPhase.size(), CV_32FC1);
        const float shiftVal = static_cast<float>(CV_2PI) / phaseSteps;

        cv::parallel_for_(cv::Range(0, wrappedPhase.rows), [&](const cv::Range& range) {
            std::vector<const float*> phasePtrs(phaseSteps);

            for (int i = range.start; i < range.end; ++i) {
                for (int k = 0; k < phaseSteps; ++k) {
                    phasePtrs[k] = phaseImages[k].ptr<float>(i);
                }

                const float* phasePtr = wrappedPhase.ptr<float>(i);
                float* errorPtr = phaseError.ptr<float>(i);

                for (int j = 0; j < wrappedPhase.cols; ++j) {
                    if (std::isnan(phasePtr[j])) {
                        errorPtr[j] = std::numeric_limits<float>::max();
                        continue;
                    }

                    // �������г������
                    float sum2ndSin = 0.0f, sum2ndCos = 0.0f;
                    float sum1stSin = 0.0f, sum1stCos = 0.0f;

                    for (int k = 0; k < phaseSteps; ++k) {
                        const float intensity = phasePtrs[k][j];
                        const float angle = k * shiftVal;

                        // һ��г��
                        sum1stSin += intensity * sin(angle);
                        sum1stCos += intensity * cos(angle);

                        // ����г��
                        sum2ndSin += intensity * sin(2 * angle);
                        sum2ndCos += intensity * cos(2 * angle);
                    }

                    // һ�׺Ͷ���г���ķ�ֵ
                    float firstHarmonic = sqrt(sum1stSin * sum1stSin + sum1stCos * sum1stCos);
                    float secondHarmonic = sqrt(sum2ndSin * sum2ndSin + sum2ndCos * sum2ndCos);

                    // ����г����һ��г���ı�ֵ��Ϊ���ָ��
                    if (firstHarmonic > 0) {
                        errorPtr[j] = secondHarmonic / firstHarmonic;
                    }
                    else {
                        errorPtr[j] = std::numeric_limits<float>::max();
                    }
                }
            }
            });

        return phaseError;
    }

    // �ۺ���λ�������
    static cv::Mat calculateComprehensivePhaseError(const std::vector<cv::Mat>& phaseImages,
                                                    const cv::Mat& wrappedPhase,
                                                    const cv::Mat& modulation,
                                                    const cv::Mat& averageIntensity,
                                                    int phaseSteps,
                                                    float noiseStd = 1.0f) {
        // ����������
        cv::Mat theoreticalError = calculateTheoreticalPhaseError(modulation, averageIntensity, noiseStd);
        cv::Mat residualError = calculateResidualPhaseError(phaseImages, wrappedPhase, averageIntensity, phaseSteps);
        cv::Mat consistencyError = calculatePhaseConsistencyError(wrappedPhase, modulation);
        cv::Mat harmonicError = calculateHarmonicPhaseError(phaseImages, wrappedPhase, phaseSteps);

        // �ۺ�����
        cv::Mat comprehensiveError = cv::Mat::zeros(wrappedPhase.size(), CV_32FC1);

        cv::parallel_for_(cv::Range(0, wrappedPhase.rows), [&](const cv::Range& range) {
            for (int i = range.start; i < range.end; ++i) {
                float* compPtr = comprehensiveError.ptr<float>(i);
                const float* theoPtr = theoreticalError.ptr<float>(i);
                const float* residPtr = residualError.ptr<float>(i);
                const float* consPtr = consistencyError.ptr<float>(i);
                const float* harmPtr = harmonicError.ptr<float>(i);

                for (int j = 0; j < wrappedPhase.cols; ++j) {
                    // ��Ȩ��ϸ���������
                    float weightedError = 0.4f * theoPtr[j] +
                        0.3f * residPtr[j] +
                        0.2f * consPtr[j] +
                        0.1f * harmPtr[j];
                    compPtr[j] = weightedError;
                }
            }
            });

        return comprehensiveError;
    }
};
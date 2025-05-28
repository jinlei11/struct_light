#pragma once

// 相位误差计算方法集合
#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>

class PhaseErrorAnalyzer {

public:
    // 方法1：基于调制度的理论相位误差计算
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
                        // 理论公式: σ_φ ≈ σ_I / (M * I_avg)
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

    // 方法2：基于残差分析的相位误差估计
    static cv::Mat calculateResidualPhaseError(const std::vector<cv::Mat>& phaseImages,
                                               const cv::Mat& wrappedPhase,
                                               const cv::Mat& averageIntensity,
                                               int phaseSteps) {
        cv::Mat phaseError = cv::Mat::zeros(wrappedPhase.size(), CV_32FC1);
        const float shiftVal = static_cast<float>(CV_2PI) / phaseSteps;

        cv::parallel_for_(cv::Range(0, wrappedPhase.rows), [&](const cv::Range& range) {
            std::vector<const float*> phasePtrs(phaseSteps);

            for (int i = range.start; i < range.end; ++i) {
                // 获取所有相移图像的行指针
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

                    // 计算重构误差
                    float sumSquaredResidual = 0.0f;
                    const float phase = phasePtr[j];
                    const float avgIntensity = avgPtr[j];

                    for (int k = 0; k < phaseSteps; ++k) {
                        // 理论光强值
                        float theoreticalIntensity = avgIntensity * (1.0f + cos(phase + k * shiftVal));
                        // 实际光强值
                        float actualIntensity = phasePtrs[k][j];
                        // 残差
                        float residual = actualIntensity - theoreticalIntensity;
                        sumSquaredResidual += residual * residual;
                    }

                    // RMS残差作为相位误差的估计
                    errorPtr[j] = sqrt(sumSquaredResidual / phaseSteps);
                }
            }
            });

        return phaseError;
    }

    // 方法3：基于相位一致性的误差估计
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

                        // 计算邻域内相位的标准差
                        std::vector<float> neighborPhases;
                        float sumWeight = 0.0f, weightedSum = 0.0f, weightedSumSq = 0.0f;

                        for (int di = -halfWindow; di <= halfWindow; ++di) {
                            for (int dj = -halfWindow; dj <= halfWindow; ++dj) {
                                if (di == 0 && dj == 0) continue;

                                const float neighborPhase = wrappedPhase.at<float>(i + di, j + dj);
                                const float neighborMod = modulation.at<float>(i + di, j + dj);

                                if (!std::isnan(neighborPhase) && neighborMod > 0) {
                                    // 处理相位包裹
                                    float phaseDiff = neighborPhase - centerPhase;
                                    if (phaseDiff > CV_PI) phaseDiff -= CV_2PI;
                                    if (phaseDiff < -CV_PI) phaseDiff += CV_2PI;

                                    float adjustedPhase = centerPhase + phaseDiff;
                                    float weight = neighborMod; // 用调制度作为权重

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

    // 方法4：基于高阶谐波的相位误差检测
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

                    // 计算二阶谐波分量
                    float sum2ndSin = 0.0f, sum2ndCos = 0.0f;
                    float sum1stSin = 0.0f, sum1stCos = 0.0f;

                    for (int k = 0; k < phaseSteps; ++k) {
                        const float intensity = phasePtrs[k][j];
                        const float angle = k * shiftVal;

                        // 一阶谐波
                        sum1stSin += intensity * sin(angle);
                        sum1stCos += intensity * cos(angle);

                        // 二阶谐波
                        sum2ndSin += intensity * sin(2 * angle);
                        sum2ndCos += intensity * cos(2 * angle);
                    }

                    // 一阶和二阶谐波的幅值
                    float firstHarmonic = sqrt(sum1stSin * sum1stSin + sum1stCos * sum1stCos);
                    float secondHarmonic = sqrt(sum2ndSin * sum2ndSin + sum2ndCos * sum2ndCos);

                    // 二阶谐波与一阶谐波的比值作为误差指标
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

    // 综合相位误差评估
    static cv::Mat calculateComprehensivePhaseError(const std::vector<cv::Mat>& phaseImages,
                                                    const cv::Mat& wrappedPhase,
                                                    const cv::Mat& modulation,
                                                    const cv::Mat& averageIntensity,
                                                    int phaseSteps,
                                                    float noiseStd = 1.0f) {
        // 计算各种误差
        cv::Mat theoreticalError = calculateTheoreticalPhaseError(modulation, averageIntensity, noiseStd);
        cv::Mat residualError = calculateResidualPhaseError(phaseImages, wrappedPhase, averageIntensity, phaseSteps);
        cv::Mat consistencyError = calculatePhaseConsistencyError(wrappedPhase, modulation);
        cv::Mat harmonicError = calculateHarmonicPhaseError(phaseImages, wrappedPhase, phaseSteps);

        // 综合评估
        cv::Mat comprehensiveError = cv::Mat::zeros(wrappedPhase.size(), CV_32FC1);

        cv::parallel_for_(cv::Range(0, wrappedPhase.rows), [&](const cv::Range& range) {
            for (int i = range.start; i < range.end; ++i) {
                float* compPtr = comprehensiveError.ptr<float>(i);
                const float* theoPtr = theoreticalError.ptr<float>(i);
                const float* residPtr = residualError.ptr<float>(i);
                const float* consPtr = consistencyError.ptr<float>(i);
                const float* harmPtr = harmonicError.ptr<float>(i);

                for (int j = 0; j < wrappedPhase.cols; ++j) {
                    // 加权组合各种误差估计
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
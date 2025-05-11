#pragma once
#include<opencv2/opencv.hpp>
#include<iostream>




//最大类间方法差，进行灰度值阈值选取
double OptimalinterclassVarianceMethod(cv::Mat& Image, int GrayValue);
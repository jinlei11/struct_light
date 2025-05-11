#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <filesystem>
#include <algorithm>


class FileManipulation {

public:

	FileManipulation();
	~FileManipulation();

	/*
	 * @ brief  读取文件夹中的图像，无需传入图片数量，自动读取文件夹下所有图片
	 * @ param  filepath         投影图像保存路径
	 * @ return  读取到的图片，使用容器存储
	*/
	std::vector<cv::Mat> ReadPics(const std::string& filepath);

	/*
	 * @ brief  读取文件夹中的图像，需要指定图片的数量
	 * @ param  filepath         投影图像保存路径
	 * @ param  ImageNum         指定读取的图像数量
	 * @ param  isRename         是否需要重命名
	 * @ return  读取到的图片，使用容器存储
	*/
	std::vector<cv::Mat> ReadPics(const std::string& filepath, int ImageNum, bool isRename = false);

	/*
	 * @ brief  读取文件夹中的图像
	 * @ param  filepath          投影图像保存路径
	 * @ param  Originate         指定读取的图像的起始序号
	 * @ param  end               指定读取的图像的结束序号
	 * @ return  读取到的图片，使用容器存储
	*/
	std::vector<cv::Mat> ReadPics(const std::string& filepath, int Originate, int end);


	/*
	 * @ brief  读取文件夹中的图像
	 * @ param  filepath          投影图像保存路径
	 * @ param  Imgs_even         存放图像序列中的偶数图像
	 * @ param  Imgs_odd          存放图像序列中的奇数图像
	 * @ return  是否读取成功的结果
	*/
	bool ReadPicsToTwoSequence(const std::string& filepath, 
							   std::vector<cv::Mat>& Imgs_even,
							   std::vector<cv::Mat>& Imgs_odd);


	/*
	 * @ brief  读取文件夹中的图像
	 * @ param  filepath         投影图像保存路径
	 * @ param  Imgs             存放图像的容器
	 * @ return  是否读取成功的结果
	*/
	bool ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs);

	/*
	 * @ brief  读取文件夹中的图像，该函数可转换图片格式
	 * @ param  filepath          投影图像保存路径
	 * @ param  Imgs              存放图像的容器
	 * @ param  targetType        将读取到的图像转换为该类型：比如传入CV_32FC1
	 * @ param  isRename         是否需要重命名
	 * @ return  是否读取成功的结果
	*/
	bool ReadPics(const std::string& filepath,
				 std::vector<cv::Mat>& Imgs,
			  	 int targetType,
				 bool isRename = false);


	/*
	 * @ brief  读取文件夹中的图像，该函数可转换图片格式,使用并行库OpenPM
	 * @ param  filepath          投影图像保存路径
	 * @ param  Imgs              存放图像的容器
	 * @ param  targetType        将读取到的图像转换为该类型：比如传入CV_32FC1
	 * @ param  isRename         是否需要重命名
	 * @ return  是否读取成功的结果
	*/
	bool ReadPics_openmp(const std::string& filepath,
						 std::vector<cv::Mat>& Imgs,
						 int targetType,
						 bool isRename = false);





	/*
	 * @ brief  读取文件夹中的图像
	 * @ param  filepath          投影图像保存路径
	 * @ param  Imgs              存放图像的容器
	 * @ param  Originate         指定读取的图像的起始序号
	 * @ param  end               指定读取的图像的结束序号
	 * @ return  是否读取成功的结果
	*/
	bool ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs, int Originate, int end);


	/*
	 * @ brief  保存图像函数
	 * @ param  Imgs             存放图像的容器
	 * @ param  filepath         保存图片文件的路径
	 * @ return  是否保存成功的结果
	*/
	bool SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath);

	/*
	 * @ brief  保存图像函数
	 * @ param  Imgs             存放图像的容器
	 * @ param  filepath         保存图片文件的路径
	 * @ param  ImageNum         指定保存的图像数量
	 * @ return  是否保存成功的结果
	*/
	bool SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath, int ImageNum);

	/*
	 * @ brief  保存图像函数
	 * @ param  Imgs             存放图像的容器
	 * @ param  filepath         保存图片文件的路径
	 * @ param  Originate         指定读取的图像的起始序号
	 * @ param  end               指定读取的图像的结束序号
	 * @ return  是否保存成功的结果
	*/
	bool SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath, int Originate, int end);


private:

	//为文件夹中图片自动改名的函数，设置图片的序号为从零开始
	bool numericStringCompare(const std::string& s1, const std::string& s2);
	int renameImagesInFolder(const std::string& folderPath, int StartingValue = 0);


};
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
	 * @ brief  ��ȡ�ļ����е�ͼ�����贫��ͼƬ�������Զ���ȡ�ļ���������ͼƬ
	 * @ param  filepath         ͶӰͼ�񱣴�·��
	 * @ return  ��ȡ����ͼƬ��ʹ�������洢
	*/
	std::vector<cv::Mat> ReadPics(const std::string& filepath);

	/*
	 * @ brief  ��ȡ�ļ����е�ͼ����Ҫָ��ͼƬ������
	 * @ param  filepath         ͶӰͼ�񱣴�·��
	 * @ param  ImageNum         ָ����ȡ��ͼ������
	 * @ param  isRename         �Ƿ���Ҫ������
	 * @ return  ��ȡ����ͼƬ��ʹ�������洢
	*/
	std::vector<cv::Mat> ReadPics(const std::string& filepath, int ImageNum, bool isRename = false);

	/*
	 * @ brief  ��ȡ�ļ����е�ͼ��
	 * @ param  filepath          ͶӰͼ�񱣴�·��
	 * @ param  Originate         ָ����ȡ��ͼ�����ʼ���
	 * @ param  end               ָ����ȡ��ͼ��Ľ������
	 * @ return  ��ȡ����ͼƬ��ʹ�������洢
	*/
	std::vector<cv::Mat> ReadPics(const std::string& filepath, int Originate, int end);


	/*
	 * @ brief  ��ȡ�ļ����е�ͼ��
	 * @ param  filepath          ͶӰͼ�񱣴�·��
	 * @ param  Imgs_even         ���ͼ�������е�ż��ͼ��
	 * @ param  Imgs_odd          ���ͼ�������е�����ͼ��
	 * @ return  �Ƿ��ȡ�ɹ��Ľ��
	*/
	bool ReadPicsToTwoSequence(const std::string& filepath, 
							   std::vector<cv::Mat>& Imgs_even,
							   std::vector<cv::Mat>& Imgs_odd);


	/*
	 * @ brief  ��ȡ�ļ����е�ͼ��
	 * @ param  filepath         ͶӰͼ�񱣴�·��
	 * @ param  Imgs             ���ͼ�������
	 * @ return  �Ƿ��ȡ�ɹ��Ľ��
	*/
	bool ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs);

	/*
	 * @ brief  ��ȡ�ļ����е�ͼ�񣬸ú�����ת��ͼƬ��ʽ
	 * @ param  filepath          ͶӰͼ�񱣴�·��
	 * @ param  Imgs              ���ͼ�������
	 * @ param  targetType        ����ȡ����ͼ��ת��Ϊ�����ͣ����紫��CV_32FC1
	 * @ param  isRename         �Ƿ���Ҫ������
	 * @ return  �Ƿ��ȡ�ɹ��Ľ��
	*/
	bool ReadPics(const std::string& filepath,
				 std::vector<cv::Mat>& Imgs,
			  	 int targetType,
				 bool isRename = false);


	/*
	 * @ brief  ��ȡ�ļ����е�ͼ�񣬸ú�����ת��ͼƬ��ʽ,ʹ�ò��п�OpenPM
	 * @ param  filepath          ͶӰͼ�񱣴�·��
	 * @ param  Imgs              ���ͼ�������
	 * @ param  targetType        ����ȡ����ͼ��ת��Ϊ�����ͣ����紫��CV_32FC1
	 * @ param  isRename         �Ƿ���Ҫ������
	 * @ return  �Ƿ��ȡ�ɹ��Ľ��
	*/
	bool ReadPics_openmp(const std::string& filepath,
						 std::vector<cv::Mat>& Imgs,
						 int targetType,
						 bool isRename = false);





	/*
	 * @ brief  ��ȡ�ļ����е�ͼ��
	 * @ param  filepath          ͶӰͼ�񱣴�·��
	 * @ param  Imgs              ���ͼ�������
	 * @ param  Originate         ָ����ȡ��ͼ�����ʼ���
	 * @ param  end               ָ����ȡ��ͼ��Ľ������
	 * @ return  �Ƿ��ȡ�ɹ��Ľ��
	*/
	bool ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs, int Originate, int end);


	/*
	 * @ brief  ����ͼ����
	 * @ param  Imgs             ���ͼ�������
	 * @ param  filepath         ����ͼƬ�ļ���·��
	 * @ return  �Ƿ񱣴�ɹ��Ľ��
	*/
	bool SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath);

	/*
	 * @ brief  ����ͼ����
	 * @ param  Imgs             ���ͼ�������
	 * @ param  filepath         ����ͼƬ�ļ���·��
	 * @ param  ImageNum         ָ�������ͼ������
	 * @ return  �Ƿ񱣴�ɹ��Ľ��
	*/
	bool SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath, int ImageNum);

	/*
	 * @ brief  ����ͼ����
	 * @ param  Imgs             ���ͼ�������
	 * @ param  filepath         ����ͼƬ�ļ���·��
	 * @ param  Originate         ָ����ȡ��ͼ�����ʼ���
	 * @ param  end               ָ����ȡ��ͼ��Ľ������
	 * @ return  �Ƿ񱣴�ɹ��Ľ��
	*/
	bool SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath, int Originate, int end);


private:

	//Ϊ�ļ�����ͼƬ�Զ������ĺ���������ͼƬ�����Ϊ���㿪ʼ
	bool numericStringCompare(const std::string& s1, const std::string& s2);
	int renameImagesInFolder(const std::string& folderPath, int StartingValue = 0);


};
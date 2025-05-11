#include"FileManipulation.h"

namespace fs = std::filesystem;


FileManipulation::FileManipulation() {


};

FileManipulation::~FileManipulation() {

};


std::vector<cv::Mat> FileManipulation::ReadPics(const std::string& filepath) {
	int num = renameImagesInFolder(filepath);
	std::vector<cv::Mat> Imgs;
	for (int i = 0; i < num; i++) {
		cv::Mat image = cv::imread(filepath + std::to_string(i) + ".bmp", cv::IMREAD_GRAYSCALE);
		if (!image.data) {
			std::cout << "�����Ƿ�������ȷ��ͼ���ļ�" << std::endl;
		}
		else{
			Imgs.push_back(image);
		}
	}
	return Imgs;
};

std::vector<cv::Mat> FileManipulation::ReadPics(const std::string& filepath, int ImageNum, bool isRename) {
	if (isRename) {
		renameImagesInFolder(filepath);
	}
	std::vector<cv::Mat> Imgs;
	for (int i = 0; i < ImageNum; i++) {
		cv::Mat image = cv::imread(filepath + std::to_string(i) + ".bmp", cv::IMREAD_GRAYSCALE);
		if (!image.data) {
			std::cout << "�����Ƿ�������ȷ��ͼ���ļ�" << std::endl;
		}
		else {
			Imgs.push_back(image);
		}
	}
	return Imgs;
};

std::vector<cv::Mat> FileManipulation::ReadPics(const std::string& filepath, int Originate, int end) {
	renameImagesInFolder(filepath);
	std::vector<cv::Mat> Imgs;
	for (int i = Originate; i < end; i++) {
		cv::Mat image = cv::imread(filepath + std::to_string(i) + ".bmp", cv::IMREAD_GRAYSCALE);
		if (!image.data) {
			std::cout << "�����Ƿ�������ȷ��ͼ���ļ�" << std::endl;
		}
		else {
			Imgs.push_back(image);
		}
	}
	return Imgs;
};


bool FileManipulation::ReadPicsToTwoSequence(const std::string& filepath,
											 std::vector<cv::Mat>& Imgs_even,
											 std::vector<cv::Mat>& Imgs_odd) {
	int num = renameImagesInFolder(filepath);
	for (int i = 0; i < num; i++) {
		cv::Mat image = cv::imread(filepath + std::to_string(i) + ".bmp", cv::IMREAD_GRAYSCALE);
		if (!image.data) {
			std::cout << "�����Ƿ�������ȷ��ͼ���ļ�" << std::endl;
		}
		else {
			//����
			if ((i & 1) == 1) {
				Imgs_odd.push_back(image);
			}
			//ż��
			else {
				Imgs_even.push_back(image);
			}
		}
	}
	return true;
};



bool FileManipulation::ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs) {
	int num = renameImagesInFolder(filepath);
	for (int i = 0; i < num; i++) {
		cv::Mat image = cv::imread(filepath + std::to_string(i) + ".bmp", cv::IMREAD_GRAYSCALE);
		if (!image.data) {
			std::cout << "�����Ƿ�������ȷ��ͼ���ļ�" << std::endl;
		}
		else {
			Imgs.push_back(image);
		}
	}
	return true;
};




bool FileManipulation::ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs, int targetType, bool isRename) {
	int num = 0;
	//�Ƿ����������
	if (isRename) {
		num = renameImagesInFolder(filepath);
	}
	else {
		try {
			for (const auto& entry : fs::directory_iterator(filepath)) {
				if (entry.is_regular_file()) {
					++num;
				}
			}
		}
		catch (const fs::filesystem_error& e) {
			std::cerr << "Error: " << e.what() << std::endl;
		}
	}

	Imgs.clear();
	Imgs.reserve(num);

	for (int i = 0; i < num; i++) {
		std::string fullPath = filepath + std::to_string(i) + ".bmp";
		cv::Mat image = cv::imread(fullPath, cv::IMREAD_GRAYSCALE); // ��ȡ�Ҷ�ͼ��
		if (image.empty()) {
			std::cout << "�����Ƿ�������ȷ��ͼ���ļ���" << fullPath << std::endl;
			continue;
		}

		// ��ԭʼͼ��ͨ������Ŀ������ͨ������ͬ����Ҫת��ͨ����
		int srcChannels = image.channels();
		int dstChannels = CV_MAT_CN(targetType);

		if (srcChannels != dstChannels) {
			if (srcChannels == 1 && dstChannels == 3) {
				cv::cvtColor(image, image, cv::COLOR_GRAY2BGR);
			}
			else if (srcChannels == 3 && dstChannels == 1) {
				cv::cvtColor(image, image, cv::COLOR_BGR2GRAY);
			}
			else {
				std::cout << "��֧�ֵ�ͨ��ת����" << fullPath << std::endl;
				continue;
			}
		}
		// ����ת��������� CV_8U �� CV_32FC1��
		if (image.type() != targetType) {
			cv::Mat converted;
			int depth = CV_MAT_DEPTH(targetType);
			image.convertTo(converted, depth);
			Imgs.push_back(converted);
		}
		else {
			Imgs.push_back(image);
		}
	}
	return true;
}




//ʹ�ò��д����Ż���
bool FileManipulation::ReadPics_openmp(const std::string& filepath,
									   std::vector<cv::Mat>& Imgs,
									   int targetType,
									   bool isRename) {
	int num = 0;
	// ����ԭʼ�ļ��������߼�
	if (isRename) {
		num = renameImagesInFolder(filepath);
	}
	else {
		try {
			for (const auto& entry : fs::directory_iterator(filepath)) {
				if (entry.is_regular_file()) ++num;
			}
		}
		catch (const fs::filesystem_error& e) {
			std::cerr << "Error: " << e.what() << std::endl;
		}
	}

	// Ԥ���������ļ�·��
	std::vector<std::string> filePaths;
	filePaths.reserve(num);
	for (int i = 0; i < num; ++i) {
		filePaths.emplace_back(filepath + std::to_string(i) + ".bmp");
	}

	// �ؼ��޸ĵ㣺ֱ��Ԥ����Imgs���������
	Imgs.clear();
	Imgs.resize(num);  // Ԥ����num��λ��

	const int targetDepth = CV_MAT_DEPTH(targetType);
	const int targetChannels = CV_MAT_CN(targetType);
	const int loadFlags = (targetChannels == 3) ? cv::IMREAD_COLOR : cv::IMREAD_GRAYSCALE;

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < num; ++i) {
		cv::Mat image = cv::imread(filePaths[i], loadFlags);

		if (image.empty()) {
			#pragma omp critical
			std::cerr << "����ʧ��: " << filePaths[i] << std::endl;
			continue;  // Imgs[i]����Ĭ�Ϲ���Ŀ�Mat
		}

		// ͨ��/����ת��
		if (image.channels() != targetChannels) {
			cv::cvtColor(image, image,
				targetChannels == 1 ? cv::COLOR_BGR2GRAY : cv::COLOR_GRAY2BGR);
		}
		if (image.depth() != targetDepth) {
			image.convertTo(image, targetDepth);
		}

		Imgs[i] = std::move(image);  // ֱ�Ӵ���Imgs�Ķ�Ӧλ��
	}

	// �Ƴ�����ʧ�ܵ�ͼƬ����ѡ��
	Imgs.erase(
		std::remove_if(Imgs.begin(), Imgs.end(),
			[](const cv::Mat& m) { return m.empty(); }),
		Imgs.end()
	);
	return !Imgs.empty();
}





bool FileManipulation::ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs, int Originate, int end) {
	renameImagesInFolder(filepath);
	for (int i = Originate; i < end; i++) {
		cv::Mat image = cv::imread(filepath + std::to_string(i) + ".bmp", cv::IMREAD_GRAYSCALE);
		if (!image.data) {
			std::cout << "�����Ƿ�������ȷ��ͼ���ļ�" << std::endl;
		}
		else {
			Imgs.push_back(image);
		}
	}
	return true;
};


bool FileManipulation::SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath) {
	//�������õ��Ż�������ʹ�ö��̣߳����������ߡ�������ģ��
	for (int i = 0; i < Imgs.size(); i++) {
		std::string filename = filepath + std::to_string(i) + ".bmp";
		try {
			bool success = cv::imwrite(filename, Imgs[i]);
			if (!success) {
				std::cerr << "Failed to write image: " << filename << std::endl;
				return false;
			}
		}
		catch (const cv::Exception& e) {
			std::cerr << "Exception writing image " << filename << ": " << e.what() << std::endl;
			return false;
		}
	}
	return true;
};
bool FileManipulation::SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath, int ImageNum) {
	//�������õ��Ż�������ʹ�ö��̣߳����������ߡ�������ģ��
	for (int i = 0; i < ImageNum; i++) {
		std::string filename = filepath + std::to_string(i) + ".bmp";
		try {
			bool success = cv::imwrite(filename, Imgs[i]);
			if (!success) {
				std::cerr << "Failed to write image: " << filename << std::endl;
				return false;
			}
		}
		catch (const cv::Exception& e) {
			std::cerr << "Exception writing image " << filename << ": " << e.what() << std::endl;
			return false;
		}
	}
	return true;
};
bool FileManipulation::SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath, int Originate, int end) {
	//�������õ��Ż�������ʹ�ö��̣߳����������ߡ�������ģ��
	for (int i = Originate; i < end; i++) {
		std::string filename = filepath + std::to_string(i) + ".bmp";
		try {
			bool success = cv::imwrite(filename, Imgs[i]);
			if (!success) {
				std::cerr << "Failed to write image: " << filename << std::endl;
				return false;
			}
		}
		catch (const cv::Exception& e) {
			std::cerr << "Exception writing image " << filename << ": " << e.what() << std::endl;
			return false;
		}
	}
	return true;
};


//Make sure to read the image files in numerical order
bool FileManipulation::numericStringCompare(const std::string& s1, const std::string& s2) {
	return std::stoi(s1) < std::stoi(s2);
}
// Renaming functions under the C++ 17 standard
int FileManipulation::renameImagesInFolder(const std::string& folderPath , int StartingValue) {

	std::vector<std::string> filenames;
	for (const auto& entry : fs::directory_iterator(folderPath)) {
		if (entry.is_regular_file()) {
			filenames.push_back(entry.path().filename().string());
		}
	}
	//Sorting file names 
	//Since numericStringCompare is a member function of CameraCalibration.
	//Need to use member function pointers or lambda functions in std::sort calls
	std::sort(filenames.begin(), filenames.end(), [this](const std::string& s1, const std::string& s2) {
		return numericStringCompare(s1, s2);
		});

	for (const auto& filename : filenames) {
		try {
			fs::path currentPath = folderPath + "\\" + filename;
			std::string newFilename = std::to_string(StartingValue) + ".bmp";
			fs::path newPath = currentPath.parent_path() / newFilename;//�����ˡ�/����һ��·��ƴ�ӵĹ���
			StartingValue++;
			fs::rename(currentPath, newPath);
		}
		catch (const std::exception& e) {
			std::cerr << "Error: " << e.what() << std::endl;
		}
	}
	return filenames.size();
}
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
			std::cout << "请检查是否输入正确的图像文件" << std::endl;
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
			std::cout << "请检查是否输入正确的图像文件" << std::endl;
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
			std::cout << "请检查是否输入正确的图像文件" << std::endl;
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
			std::cout << "请检查是否输入正确的图像文件" << std::endl;
		}
		else {
			//奇数
			if ((i & 1) == 1) {
				Imgs_odd.push_back(image);
			}
			//偶数
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
			std::cout << "请检查是否输入正确的图像文件" << std::endl;
		}
		else {
			Imgs.push_back(image);
		}
	}
	return true;
};




bool FileManipulation::ReadPics(const std::string& filepath, std::vector<cv::Mat>& Imgs, int targetType, bool isRename) {
	int num = 0;
	//是否进行重命名
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
		cv::Mat image = cv::imread(fullPath, cv::IMREAD_GRAYSCALE); // 读取灰度图像
		if (image.empty()) {
			std::cout << "请检查是否输入正确的图像文件：" << fullPath << std::endl;
			continue;
		}

		// 若原始图像通道数与目标类型通道数不同，需要转换通道数
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
				std::cout << "不支持的通道转换：" << fullPath << std::endl;
				continue;
			}
		}
		// 类型转换（比如从 CV_8U 到 CV_32FC1）
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




//使用并行处理优化后
bool FileManipulation::ReadPics_openmp(const std::string& filepath,
									   std::vector<cv::Mat>& Imgs,
									   int targetType,
									   bool isRename) {
	int num = 0;
	// 保留原始文件名处理逻辑
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

	// 预生成所有文件路径
	std::vector<std::string> filePaths;
	filePaths.reserve(num);
	for (int i = 0; i < num; ++i) {
		filePaths.emplace_back(filepath + std::to_string(i) + ".bmp");
	}

	// 关键修改点：直接预分配Imgs并并行填充
	Imgs.clear();
	Imgs.resize(num);  // 预分配num个位置

	const int targetDepth = CV_MAT_DEPTH(targetType);
	const int targetChannels = CV_MAT_CN(targetType);
	const int loadFlags = (targetChannels == 3) ? cv::IMREAD_COLOR : cv::IMREAD_GRAYSCALE;

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < num; ++i) {
		cv::Mat image = cv::imread(filePaths[i], loadFlags);

		if (image.empty()) {
			#pragma omp critical
			std::cerr << "加载失败: " << filePaths[i] << std::endl;
			continue;  // Imgs[i]保持默认构造的空Mat
		}

		// 通道/类型转换
		if (image.channels() != targetChannels) {
			cv::cvtColor(image, image,
				targetChannels == 1 ? cv::COLOR_BGR2GRAY : cv::COLOR_GRAY2BGR);
		}
		if (image.depth() != targetDepth) {
			image.convertTo(image, targetDepth);
		}

		Imgs[i] = std::move(image);  // 直接存入Imgs的对应位置
	}

	// 移除加载失败的图片（可选）
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
			std::cout << "请检查是否输入正确的图像文件" << std::endl;
		}
		else {
			Imgs.push_back(image);
		}
	}
	return true;
};


bool FileManipulation::SavePics(std::vector<cv::Mat>& Imgs, const std::string& filepath) {
	//后续更好的优化方法，使用多线程，构建生产者、消费者模型
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
	//后续更好的优化方法，使用多线程，构建生产者、消费者模型
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
	//后续更好的优化方法，使用多线程，构建生产者、消费者模型
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
			fs::path newPath = currentPath.parent_path() / newFilename;//重载了‘/’是一个路径拼接的功能
			StartingValue++;
			fs::rename(currentPath, newPath);
		}
		catch (const std::exception& e) {
			std::cerr << "Error: " << e.what() << std::endl;
		}
	}
	return filenames.size();
}
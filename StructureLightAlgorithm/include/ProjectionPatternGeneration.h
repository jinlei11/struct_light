#pragma once
//
//#include <iostream>
//#include <deque>
//#include <math.h>
//#include <opencv2/opencv.hpp>
//
//using namespace std;
////改变最后一位字符值
//#define change_bit(x) x=((x)=='0'?'1':'0')
//#define pi CV_PI
//
////Projector Models
//enum DLPModel {
//	DLP4710 = 0,
//	DLP6500 = 1,
//	DLP3010 = 2,
//	DLP4500 = 3,
//};
//
//
//class BuildGrayPattern {
//public:
//	// Gray code bits width of the image and height of the image 
//	BuildGrayPattern(int bit, int width, int height);
//	//Constructor, pass in projector model and number of images generated
//	BuildGrayPattern(int bit, enum DLPModel DLP, int direction);
//
//	//Generate Gray code values
//	deque<string> BuildGrayCode();
//
//	//retrieve the opposite of what one wants
//	void qufan(char* a, char* b, int n);
//
//	//Generate Gray Code Images
//	vector<cv::Mat> BuildGrayP();
//
//	//Generate grayscale values of 0 or 255 based on Gray code values.
//	int jugevalue(char a);
//
//	//Storing the generated Gray code
//	std::deque<std::string>m_GrayArr;
//
//	//For storing Gray Code images
//	vector<cv::Mat>M_GrayPicture;
//
//	//Produces  all-black image    Existing in Projection_Pattern file
//	void buildBlack(cv::Size Image_size);
//	//Generate full-white images   Existing in Projection_Pattern file
//	void buildWhite(cv::Size Image_size);
//
//	//Produces  all-black image    Existing in the Specify path file
//	void buildBlack(const std::string& path);
//	//Generate full-white images   Existing in the Specify path file
//	void buildWhite(const std::string& path);
//
//	//Produces an all-black image with the specified gray value, less than 20
//	void buildBlack(const std::string& path, int value);
//	//Generates an all-white image with the specified gray value, greater than 240
//	void buildWhite(const std::string& path, int value);
//
//
//	//Generate and save images
//	void BuildAndSave_GrayPattren(const std::string path);
//
//	//Saving Gray Code Images
//	void SaveGray(const std::string path);
//
//
//
//
//	//重构程序
//	//生成格雷码有几种方法：
//	//
//	
//
//
//
//private:
//	//Number of bits used to hold the Gray code image
//	int m_Bit;
//	//Size, width and height of the generated image
//	int m_Width;
//	int m_Height;
//	//Determine two scenarios
//	bool flag;
//	int b;
//
//	int M_direction;
//};
//
//class CreateTW {
//
//public:
//	//constructer
//	CreateTW(float f, int bit, int width, int height, int direction);
//	//constructer
//	CreateTW(float f, int bit, enum DLPModel DLP, int direction);
//
//	//Generate phase-shifted stripe patterns
//	vector<cv::Mat> creatP();
//
//	//Saving phase-shifted streak images
//	void SavePhase(std::string& path);
//
//	//For storing phase-shifted stripe patterns
//	vector<cv::Mat>M_PhasePicture;
//
//	//Generate and save phase shift patterns
//	void BuildAndSave_phase_shift(const std::string& path);
//
//private:
//	//Frequency of the phase-shifted image
//	float f;
//	//image size
//	int m_bit;
//	int m_width;
//	int m_height;
//	//Number of phase shift patterns generated
//	int N = 4;
//	int M_direction;
//};




#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <bitset>
#include "FileManipulation.h"

#define G_GrayCodeMaxNum (1 << 10)
#define G_GrayCodeMaxBit 20    
#define pi CV_PI


enum DLPModel {
	DLP4710 = 0,
	DLP6500 = 1,
	DLP3010 = 2,
	DLP4500 = 3,
};

//--------------------------------------格雷码编码规则-------------------------------------------------//
/*******************************第一种方法**********************************/
//第零步，假设要创建n位数据的格雷码值，则第一值为 00000 共n个零
//第一步，改变最右边的位元值；
//第二步，改变右起第一个为1的位元的左边位元；
//第三步，第四步重复第一步和第二步，直到所有的格雷码产生完毕（换句话说，已经走了(2 ^ n) - 1 步）。


/********************************第二种方法******************************************/
//推荐使用这种方法
//根据二进制和格雷码之间的转换关系
//第一步，生成二进制码（这个其实很简单，计算机内存储的就是二进制的形式，所以有几位，就计算出 2^n 个数就可以了）
//第二步，根据二进制码，通过转换关系转换成格雷码
// （从左数第一位相同，随后用第一位的值与右侧相邻位做异或操作：相同为零，不同为一）
//第三步，依次生成格雷码码值


/**********************************第三种方法***************************************/
//镜像排列生成格雷码法
//当只有一位的时候，格雷码值要么是0，要么是1；
//当有两位时,格雷码首先是：0
//                         1
// 这个0和1就是只有一位时的格雷码分布
//随后镜像得到             1
//                         0
//然后再上半区最左侧加0，下半区最左侧加1
//               00
//               01
//               11
//               10
//当有三位时，右侧两位的格雷码就延续两位时的码值
//               00
//               01
//               11
//               10
//随后进行镜像-----------
//               10
//               11
//               01
//               00
// 然后在上半区最左侧加0，下半区最左侧加1  以此类推

/************************************第四种方法********************************************/
//正序逆序生成格雷码方法（对称法）





//--------------------------------------格雷码解码规则-------------------------------------------------//
//第一种方法
//第一步，针对一个格雷码值，左边第一位保持不变，作为解码后的二进制码的第一位
//第二步，对格雷码从左边第二位开始，与左侧第一位进行异或操作，得到第二位二进制码值
//第三步，重复上述操作，直至计算出最后一位二进制数据



//这里为了练习一下虚构函数以及子类和父类对象的操作，将上面的四种不同编码规则作为四种子类对象，重写父类虚函数
class GrayCodeOperation
{
public:
	GrayCodeOperation();

	virtual ~GrayCodeOperation();

	virtual std::vector<int> Generate_value() = 0;

	virtual std::vector<cv::Mat> Generate_pics() = 0;

	virtual std::vector<cv::Mat> Generate_pics_reverse() = 0;

	virtual std::vector<cv::Mat> GenerateAndSave_pics(std::string& file_path) = 0;

	virtual std::vector<cv::Mat> GenerateAndSave_pics_reverse(std::string& file_path) = 0;

	virtual std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) = 0;

	virtual std::vector<cv::Mat> GenerateAndSave_pics_block(std::string& file_path, int bit, int rowBlock, int colBlock) = 0;

	virtual cv::Mat Decode_Pics(std::string& file_path, int PicNum) = 0;

	virtual cv::Mat Decode_PicsPositiveAndNegative(std::string& file_path, int PicNum) = 0;

	virtual cv::Mat GeneratePatternWithValue(const std::string& path, int value) = 0;



protected:

	inline void ShowImages(const std::vector<cv::Mat>& images, const std::string& windowName = "GraycodePics") {
		cv::namedWindow(windowName, cv::WINDOW_NORMAL);
		for (size_t i = 0; i < images.size(); ++i) {
			if (images[i].empty()) {
				std::cerr << "Warning: Image at index " << i << " is empty, skipping." << std::endl;
				continue;
			}
			cv::imshow(windowName, images[i]);
			cv::waitKey(0);
		}
		cv::destroyWindow(windowName);
	}


};



class GrayMethod_1 : public GrayCodeOperation {

public:
	/*
	 * @ brief  子类方法一构造函数
	 * @ prarm  bit      生成格雷码的位数 
	 * @ prarm  width    生成图片的宽度
	 * @ param  height   生成图片的高度
	 * @ return 
	*/
	GrayMethod_1(int bit, int width, int height);


	/*
	 * @ brief   生成格雷码值
	 * @ return  通过int类型的容容器来存储生成的格雷码值
	*/
	std::vector<int> Generate_value() override;

	/*
	 * @ brief   生成格雷码图像
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> Generate_pics() override;

	/*
	 * @ brief   生成格雷码图像，按照划分块的方式生成，每一列的格雷码值不再是一样的。
	 * @ param   生成格雷码的位数（和构造函数中的格雷码位数是一个意思）
	 * @ param   生成格雷码区域的行块数
	 * @ param   生成格雷码区域的列块数
 	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) override;


private:
	int M_bit;//格雷码位数
	int M_width;//图像宽度
	int M_height;//图像高度

	char M_grayCodes[G_GrayCodeMaxNum][G_GrayCodeMaxBit];    //Gray's code sequence
};


//推荐使用的方法
class GrayMethod_2 : public GrayCodeOperation {

public:
	/*
	 * @ brief  子类方法二构造函数
	 * @ prarm  bit      生成格雷码的位数
	 * @ prarm  width    生成图片的宽度
	 * @ param  height   生成图片的高度
	 * @ return
	*/
	GrayMethod_2(int bit, enum DLPModel DLP, int dircetion);

	/*
	 * @ brief  子类方法二 析构函数
	 * @ return
	*/
	~GrayMethod_2();

	/*
	 * @ brief   生成格雷码值
	 * @ return  通过int类型的容容器来存储生成的格雷码值
	*/
	std::vector<int> Generate_value() override;


	/*
	 * @ brief   生成指定灰度值的整张图像
	 * @ param  file_path   保存图像的路径
	 * @ param  value       指定的灰度值
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	cv::Mat GeneratePatternWithValue(const std::string& path, int value);

	/*
	 * @ brief   生成一系列具有不同均匀灰度值的图像序列
	 * @ param   filepath    图片存放路径
	 * @ param   step        生成图像时灰度值的步长
	*/
	void generationDiffGrayValuePattern(const std::string& filepath, int step);


	/*
	 * @ brief   生成格雷码图像
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> Generate_pics() override;

	/*
	 * @ brief   生成格雷码图像——反码图像
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> Generate_pics_reverse() override;


	/*
	 * @ brief   生成并保存格雷码图像
	 * @ param  file_path   保存格雷码图像的路径
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> GenerateAndSave_pics(std::string& file_path) override;

	/*
	 * @ brief   生成并保存格雷码图像——反码图像
	 * @ param  file_path   保存格雷码图像的路径
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> GenerateAndSave_pics_reverse(std::string& file_path) override;

	/*
	 * @ brief   生成格雷码图像，按照划分块的方式生成，每一列的格雷码值不再是一样的。
	 * @ param  bit        生成格雷码的位数（和构造函数中的格雷码位数是一个意思）
	 * @ param  rowBlock   生成格雷码区域的行块数
	 * @ param  colBlock   生成格雷码区域的列块数
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) override;

	/*
	 * @ brief   生成格雷码图像，按照划分块的方式生成，每一列的格雷码值不再是一样的。
	 * @ param  file_path  存放图像的路径
	 * @ param  bit        生成格雷码的位数（和构造函数中的格雷码位数是一个意思）
	 * @ param  rowBlock   生成格雷码区域的行块数
	 * @ param  colBlock   生成格雷码区域的列块数
	 * @ return  通过Mat类型的容容器来存储生成的格雷码图像
	*/
	std::vector<cv::Mat> GenerateAndSave_pics_block(std::string& file_path, int bit, int rowBlock, int colBlock) override;

	/*
	 * @ brief   解码格雷码图像：需要是二值化后的图像
	 * @ param   格雷码图像存放的路径
	 * @ param   格雷码图像的数量
	 * @ return  通过Mat类型的数据来存储解码出的格雷码图像
	*/
	cv::Mat Decode_Pics(std::string& file_path, int PicNum) override;

	/*
	 * @ brief   解码格雷码图像：正反码图像,无需对图像进行二值化处理
	 * @ param   格雷码图像存放的路径
	 * @ param   格雷码图像的数量
	 * @ return  通过Mat类型的数据来存储解码出的格雷码图像
	*/
	cv::Mat Decode_PicsPositiveAndNegative(std::string& file_path, int PicNum) override;


private:
	//文件操作成员变量
	FileManipulation M_ImageOperation;

	int M_bit;//格雷码位数
	int M_Width;//图像宽度
	int M_Height;//图像高度
	int M_Direction;//图像的方向
};



class GrayMethod_3 : public GrayCodeOperation {

public:

	GrayMethod_3();

	std::vector<int> Generate_value() override;

	std::vector<cv::Mat> Generate_pics() override;

};



class PhasePatternOperation {

public:
	/*
	 * @ brief  条纹图案生成类的构造函数
	 * @ prarm  f         相移的图案的频率
	 * @ prarm  bit       格雷码图像的位数
	 * @ param  DLP       DLP投影仪的型号，主要用来设置生成图片的高度和宽度
	 * @ param  direction 生成相移图案的方向，水平或者竖直
	 * @ param  num       生成相移图案的步数，默认为四步相移
 	 * @ return
	*/
	PhasePatternOperation(float f, enum DLPModel DLP, int direction = 0, int num = 4);
	
	~PhasePatternOperation();

	//Generate phase-shifted stripe patterns
	std::vector<cv::Mat> CreatPhasePatterns();


	/*
	 * @ brief  生成并保存相移条纹图案
	 * @ param  filepath     保存相移条纹图案的路径
	 * @ return 返回存放生成相移条纹图案的图像序列
	*/
	std::vector<cv::Mat> BuildAndSave_PhaseShiftingPattern(const std::string& filepath);

private:

	inline void ShowPhaseImages(const std::vector<cv::Mat>& images, const std::string& windowName = "PhaseShiftingPics") {
		cv::namedWindow(windowName, cv::WINDOW_NORMAL);
		for (size_t i = 0; i < images.size(); ++i) {
			if (images[i].empty()) {
				std::cerr << "Warning: Image at index " << i << " is empty, skipping." << std::endl;
				continue;
			}
			cv::imshow(windowName, images[i]);
			cv::waitKey(0);
		}
		cv::destroyWindow(windowName);
	}


	//文件操作成员变量
	FileManipulation M_ImageOperation;

	//For storing phase-shifted stripe patterns
	std::vector<cv::Mat>M_PhasePicture;

	int M_PhaseShiftingNum;//相依图像的步数
	int M_GrayCodeBit;//格雷码图像的位数
	float M_F;//相依图像的频率
	int M_Direction;//相依图像的方向
	int M_Width;//生成图像的宽度
	int M_Height;//生成图像的高度


};
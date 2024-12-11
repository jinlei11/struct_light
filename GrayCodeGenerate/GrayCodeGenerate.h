/**************************本套算法库用于实现多种方法的格雷码图像生成*******************************/
#pragma once
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <bitset>

#define G_GrayCodeMaxNum (1 << 10)
#define G_GrayCodeMaxBit 20

//--------------------------------------格雷码编码规则-------------------------------------------------//
/*******************************第一种方法**********************************/
//第零步，假设要创建n位数据的格雷码值，则第一值为 00000 共n个零
//第一步，改变最右边的位元值；
//第二步，改变右起第一个为1的位元的左边位元；
//第三步，第四步重复第一步和第二步，直到所有的格雷码产生完毕（换句话说，已经走了(2 ^ n) - 1 步）。


/********************************第二种方法******************************************/
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

	virtual std::vector<cv::Mat> Generate_pics_block(int bit,int rowBlock, int colBlock) = 0;

	virtual cv::Mat Decode_Pics(std::string file_path, int PicNum) = 0;

private:

	

};



class GrayMethod_1 : public GrayCodeOperation {

public:
	
	GrayMethod_1(int bit, int width , int height);


	std::vector<int> Generate_value() override;

	std::vector<cv::Mat> Generate_pics() override;

	std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) override;





private:
	int M_bit;
	int M_width;
	int M_height;

	char M_grayCodes[G_GrayCodeMaxNum][G_GrayCodeMaxBit];    //Gray's code sequence
};



class GrayMethod_2 : public GrayCodeOperation {

public:

	GrayMethod_2(int bit, int width, int height);

	std::vector<int> Generate_value() override;

	std::vector<cv::Mat> Generate_pics() override;

	std::vector<cv::Mat> Generate_pics_block(int bit, int rowBlock, int colBlock) override;


	cv::Mat Decode_Pics(std::string file_path, int PicNum) override;

private:
	int M_bit;
	int M_width;
	int M_height;
};



class GrayMethod_3 : public GrayCodeOperation {

public:
	GrayMethod_3();

	std::vector<int> Generate_value() override;

	std::vector<cv::Mat> Generate_pics() override;

};
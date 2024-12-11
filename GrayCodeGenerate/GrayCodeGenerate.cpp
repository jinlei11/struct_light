#include "GrayCodeGenerate.h"

//父类
GrayCodeOperation::GrayCodeOperation() {

};

GrayCodeOperation::~GrayCodeOperation() {

};


//子类1
GrayMethod_1::GrayMethod_1(int bit, int width, int height)
    :M_bit(bit),M_height(height),M_width(width)
{

};

/***********************方法1：改变最右侧位的值 + 改变右数第一个1的左侧位值*****************************/
//这种方法好像只能操作字符类型的数据。
std::vector<int> GrayMethod_1::Generate_value() {

    char digit[G_GrayCodeMaxBit];               //one-group Gray code

    std::vector<int> GrayCode_hex;
    GrayCode_hex.push_back(1);
    bool Operation_flag = true;

    //generate 00000
    for (int i = 0; i < M_bit; i++) {
        digit[i] = '0';
    }

    for (int count = 0; count < (1 << M_bit); count++) {  // 2^bits 个格雷码
        // 保存上一位格雷码值
        for (int i = 0; i < M_bit; i++) {
            //grayCodes[count][i] = digit[i];          // 从左到右存储
            M_grayCodes[count][M_bit - 1 - i] = digit[i];  // 从右到左存储
        }
        M_grayCodes[count][M_bit] = '\0';  // 结束符，确保字符串完整

        if (Operation_flag) {
            digit[0] = (digit[0] == '0') ? '1' : '0';  // 如果是奇数轮次，翻转第一个位
        }
        else {
            // 计算第一个为1的位置
            int j = 0;
            for ( ; j < M_bit && digit[j] == '0'; j++);
            // 如果已经是最后一个格雷码，跳出    
            if (j == M_bit - 1) {
                break;
            }  
            digit[j + 1] = (digit[j + 1] == '0') ? '1' : '0';
        }
        Operation_flag = !Operation_flag;  // 切换奇偶标志
    }

    printf("\n格雷码序列：\n");
    for (int i = 0; i < (1 << M_bit); i++) {
        printf("%s\n", M_grayCodes[i]);
    }
    return GrayCode_hex;
};


std::vector<cv::Mat> GrayMethod_1::Generate_pics() {

    std::vector<cv::Mat> Images;
    int region_width = M_width / (1 << M_bit);  // 每个区域的宽度
    Generate_value();

    for (int BitIndex = 0; BitIndex < M_bit; BitIndex++) {
        cv::Mat img(M_height, M_width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像

        for (int region = 0; region < (1<<M_bit); region++) {
            int grayValue = (M_grayCodes[region][BitIndex] == '1') ? 0 : 255;  // 根据格雷码位值设置颜色
            // 将当前区域的所有像素值设置为该值（黑或白）
            for (int row = 0; row < M_height; row++) {
                for (int col = region * region_width; col < (region + 1) * region_width; col++) {
                    img.at<uchar>(row, col) = grayValue;
                }
            }
        }
        Images.push_back(img);
    }
    for (int i = 0; i < M_bit; i++) {
        cv::namedWindow("GraycodePics",cv::WindowFlags::WINDOW_AUTOSIZE);
        cv::imshow("GraycodePics", Images[i]);
        cv::waitKey(0);
    }
        

    return Images;
};

std::vector<cv::Mat> GrayMethod_1::Generate_pics_block(int bit, int rowBlock, int colBlock) {

    //根据需要生成的区域传入bit值,生成格雷码序列
    std::vector<int> GrayCode_hex;
    for (int i = 0; i < (1 << bit); i++) {
        GrayCode_hex.push_back(i ^ (i >> 1));
    }

    std::vector<cv::Mat> Images;
    int region_width = M_width / colBlock;
    int region_height = M_height / rowBlock;
    int grayValue = 0;
    int initial_value = 1;//初始值不是从零开始，而是从一开始
    // 生成图像
    for (int BitIndex = 0; BitIndex < bit; BitIndex++) {
        cv::Mat img(M_height, M_width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像
        int Jump = 0;
        // 遍历所有区域，根据格雷码的当前位决定区域颜色
        for (int region_row = 0; region_row < rowBlock; region_row++) {
            for (int region_col = 0; region_col < colBlock; region_col++) {
                // 跳过第85个格雷码，就是按照二进制生成的格雷码的第85个（指的是顺序）
                if (region_row * colBlock + region_col + 1 == 85) {
                    Jump = 1;
                    grayValue = (GrayCode_hex[region_row * colBlock + region_col + initial_value + Jump] >> (bit - BitIndex - 1)) & 1;
                }
                else {
                    grayValue = (GrayCode_hex[region_row * colBlock + region_col + initial_value + Jump] >> (bit - BitIndex - 1)) & 1; //获取对应的位
                }
                grayValue = (grayValue == 1) ? 255 : 0;
                // 将当前区域的所有像素值设置为该值（黑或白）
                for (int row = region_row * region_height; row < (region_row + 1) * region_height; row++) {
                    for (int col = region_col * region_width; col < (region_col + 1) * region_width; col++) {
                        img.at<uchar>(row, col) = grayValue;
                        int k = 5;
                    }
                }
            }
        }
        Images.push_back(img);  // 添加当前生成的图像到结果中
    }

    // 显示生成的图像
    for (int i = 0; i < bit; i++) {
        cv::namedWindow("GraycodePics", cv::WindowFlags::WINDOW_NORMAL);
        cv::imshow("GraycodePics", Images[i]);
        cv::waitKey(0);
    }

    std::string SavePath = R"(D:\postgraduate\code\GrayCodeGenerate\GraycodePics\)";
    for (int i = 0; i < bit; i++) {
        cv::imwrite(SavePath + std::to_string(i) + ".bmp", Images[i]);
    }

    return Images;
};




/************************方法二：根据二进制码和格雷码之间的转换关系进行编码******************************/
GrayMethod_2::GrayMethod_2(int bit, int width, int height)
    :M_bit(bit), M_height(height), M_width(width)
{

};

std::vector<int> GrayMethod_2::Generate_value() {
    std::vector<int> GrayCode_hex;

    for (int i = 0; i < (1<<M_bit); i++) {
        GrayCode_hex.push_back(i ^ (i >> 1));
    }

    for (int i = 0; i < (1 << M_bit); i++) {
        std::cout << "常规十进制：" << i << " 对应的格雷码（二进制）："
            << std::bitset<10>(GrayCode_hex[i]) << std::endl;
    }
    return GrayCode_hex;
};

std::vector<cv::Mat> GrayMethod_2::Generate_pics() {
    std::vector<int> value = Generate_value();
    std::vector<cv::Mat> Images;

    int region_width = M_width / (1 << M_bit);  // 每个区域的宽度

    // 生成图像
    for (int BitIndex = 0; BitIndex < M_bit; BitIndex++) {
        cv::Mat img(M_height, M_width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像

        // 遍历所有区域，根据格雷码的当前位决定区域颜色
        for (int region = 0; region < (1<<M_bit); region++) {
            // 获取当前格雷码的指定位值
            int grayValue = (value[region] >> (M_bit - BitIndex - 1)) & 1;  // 获取对应的位
            grayValue = (grayValue == 1) ? 255 : 0;  // 根据格雷码位值设置颜色（1 -> 0, 0 -> 255）

            // 将当前区域的所有像素值设置为该值（黑或白）
            for (int row = 0; row < M_height; row++) {
                for (int col = region * region_width; col < (region + 1) * region_width; col++) {
                    img.at<uchar>(row, col) = grayValue;
                }
            }
        }

        Images.push_back(img);  // 添加当前生成的图像到结果中
    }

    // 显示生成的图像
    for (int i = 0; i < M_bit; i++) {
        cv::namedWindow("GraycodePics", cv::WindowFlags::WINDOW_AUTOSIZE);
        cv::imshow("GraycodePics", Images[i]);
        cv::waitKey(0);
    }
    return Images;
};

std::vector<cv::Mat> GrayMethod_2::Generate_pics_block(int bit, int rowBlock, int colBlock) {

    //根据需要生成的区域传入bit值,生成格雷码序列
    std::vector<int> GrayCode_hex;
    for (int i = 0; i < (1 << bit); i++) {
        GrayCode_hex.push_back(i ^ (i >> 1));
    }

    std::vector<cv::Mat> Images;
    int region_width = M_width / colBlock;
    int region_height = M_height / rowBlock;
    int grayValue = 0;
    int initial_value = 1;//初始值不是从零开始，而是从一开始
    // 生成图像
    for (int BitIndex = 0; BitIndex < bit; BitIndex++) {
        cv::Mat img(M_height, M_width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像
        int Jump = 0;
        // 遍历所有区域，根据格雷码的当前位决定区域颜色
        for (int region_row = 0; region_row < rowBlock; region_row++) {
            for (int region_col = 0; region_col < colBlock; region_col++) {
                // 跳过第85个格雷码，就是按照二进制生成的格雷码的第85个（指的是顺序）
                if (region_row * colBlock + region_col +1 == 85) {
                    Jump = 1;
                    grayValue = (GrayCode_hex[region_row * colBlock + region_col + initial_value  + Jump] >> (bit - BitIndex - 1)) & 1;
                }
                else {
                    grayValue = (GrayCode_hex[region_row * colBlock + region_col + initial_value  + Jump] >> (bit - BitIndex - 1)) & 1; //获取对应的位
                }
                grayValue = (grayValue == 1) ? 255 : 0;             
                // 将当前区域的所有像素值设置为该值（黑或白）
                for (int row = region_row * region_height; row < (region_row+1) * region_height; row++) {
                    for (int col = region_col * region_width; col < (region_col + 1) * region_width; col++) {
                        img.at<uchar>(row, col) = grayValue;
                        int k = 5;
                    }
                }
            }
        }
        Images.push_back(img);  // 添加当前生成的图像到结果中
    }

    // 显示生成的图像
    for (int i = 0; i < bit; i++) {
        cv::namedWindow("GraycodePics", cv::WindowFlags::WINDOW_NORMAL);
        cv::imshow("GraycodePics", Images[i]);
        cv::waitKey(0);
    }

    std::string SavePath = R"(D:\postgraduate\code\GrayCodeGenerate\GraycodePics\)";
    for (int i = 0; i < bit; i++) {
        cv::imwrite(SavePath + std::to_string(i) + ".bmp", Images[i]);
    }

    return Images;
};

//初始方法
cv::Mat GrayMethod_2::Decode_Pics(std::string file_path, int PicNum) {
    std::vector<cv::Mat> GrayCodeSequence;
    //这种读取所有
    for (int i = 0; i < PicNum; i++) {
        cv::Mat image= cv::imread(file_path + std::to_string(i) + ".bmp");
        GrayCodeSequence.push_back(image);
    }
    cv::Mat decodedImage = cv::Mat::zeros(M_height, M_width, CV_8UC1);
    for (int row = 0; row < M_height; row++) {
        for (int col = 0; col < M_width; col++) {
            int grayValue = 0;
            for (int BitIndex = 0; BitIndex < PicNum; BitIndex++) {
                //Guaranteed read-only is safer, code is more readable and maintainable, avoid double-counting subscripts
                const cv::Mat& img = GrayCodeSequence[BitIndex];
                int bitValue = (img.at<uchar>(row, col) == 255) ? 1 : 0;
                grayValue |= (bitValue << (M_bit - BitIndex - 1));           
            }
            int binaryValue = 0;
            //最高位不变，用最高位异或下一位，得到下一位值，以此类推
            //B=G⊕(G>>1)⊕(G>>2)⊕…
            for (int temp = grayValue; temp; temp >>= 1) {
                binaryValue ^= temp;
            }
            // 将解码后的二进制值存储到图像中
            decodedImage.at<uchar>(row, col) = static_cast<uchar>(binaryValue);
        }
    }
    return decodedImage;
};
////优化后方法
//cv::Mat GrayMethod_2::Decode_Pics(std::string file_path, int PicNum) {
//    cv::Mat decodedImage = cv::Mat::zeros(M_height, M_width, CV_8UC1);
//
//    // 提前构建 Gray-to-Binary 转换表
//    std::vector<int> grayToBinary(1 << PicNum, 0);
//    for (int i = 0; i < (1 << PicNum); i++) {
//        int binaryValue = 0;
//        for (int temp = i; temp; temp >>= 1) {
//            binaryValue ^= temp;
//        }
//        grayToBinary[i] = binaryValue;
//    }
//
//    // 解码逻辑
//    for (int BitIndex = 0; BitIndex < PicNum; BitIndex++) {
//        cv::Mat img = cv::imread(file_path + std::to_string(BitIndex) + ".bmp", cv::IMREAD_GRAYSCALE);
//        if (img.empty() || img.rows != M_height || img.cols != M_width) {
//            throw std::runtime_error("Error reading image: " + file_path + std::to_string(BitIndex) + ".bmp");
//        }
//
//        for (int row = 0; row < M_height; row++) {
//            const uchar* imgData = img.ptr<uchar>(row);
//            uchar* decodedData = decodedImage.ptr<uchar>(row);
//
//            for (int col = 0; col < M_width; col++) {
//                int bitValue = (imgData[col] == 255) ? 1 : 0;
//                decodedData[col] |= (bitValue << (M_bit - BitIndex - 1));
//            }
//        }
//    }
//
//    // 转换为二进制值
//    for (int row = 0; row < M_height; row++) {
//        uchar* decodedData = decodedImage.ptr<uchar>(row);
//        for (int col = 0; col < M_width; col++) {
//            decodedData[col] = static_cast<uchar>(grayToBinary[decodedData[col]]);
//        }
//    }
//
//    return decodedImage;
//}




/********************************方法三：镜像排列生成格雷码方法***************************************/
//镜像排列生成格雷码方法
GrayMethod_3::GrayMethod_3() {

};


std::vector<int> GrayMethod_3::Generate_value() {
    std::vector<int> GrayCode_hex;
	std::cout << "子类3" << std::endl;


    return GrayCode_hex;
};

std::vector<cv::Mat> GrayMethod_3::Generate_pics() {

    std::vector<cv::Mat> Images;


    return Images;

};


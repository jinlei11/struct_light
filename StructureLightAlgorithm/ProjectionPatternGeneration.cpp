#include"ProjectionPatternGeneration.h"



//
//BuildGrayPattern::BuildGrayPattern(int bit, enum DLPModel DLP, int direction) {
//    this->m_Bit = bit;
//    this->M_direction = direction;
//    switch (DLP)
//    {
//    case DLP4710:
//        this->m_Width = 1920;
//        this->m_Height = 1280;
//        break;
//    case DLP6500:
//        this->m_Width = 1920;
//        this->m_Height = 1080;
//        break;
//    case DLP3010:
//
//        break;
//    case DLP4500:
//        this->m_Width = 912;
//        this->m_Height = 1140;
//        break;
//    default:
//        break;
//    }
//};
//
//
//deque<string> BuildGrayPattern::BuildGrayCode() {
//    char gray[10] = { 0 };
//    char temp[10] = { 0 };
//    for (int i = 0; i < m_Bit; i++) {
//        gray[i] = '0';
//    }
//    m_GrayArr.push_back(gray);
//    flag = true;
//    while (1) {
//        if (flag) {
//            change_bit(gray[m_Bit - 1]);
//        }
//        else {
//            //Find the location of the first code point counting from the right that is a 1
//            for (b = m_Bit - 1; b > 0 && gray[b] == '0'; b--);
//            if (b == 0) {
//                break;
//            }
//            change_bit(gray[b - 1]);
//        }
//        m_GrayArr.push_back(gray);
//        //Change flag status
//        flag = !flag;
//    }
//    cout << m_GrayArr.size() << endl;
//    for (int i = 0; i < m_GrayArr.size(); i++) {
//        cout << m_GrayArr[i] << endl;
//    }
//    return m_GrayArr;
//}
//
//void BuildGrayPattern::qufan(char* a, char* b, int n) {
//    for (int i = 0; i < n; i++) {
//        b[i] = a[n - 1 - i];
//    }
//}
//
//int BuildGrayPattern::jugevalue(char a) {
//    switch (a) {
//    case '0':
//        b = 0;
//        break;
//    case '1':
//        b = 255;
//        break;
//    default:
//        break;
//    }
//    return b;
//}
//
//
//
//
//vector<cv::Mat> BuildGrayPattern::BuildGrayP() {
//    cv::Mat image = cv::Mat::zeros(m_Height, m_Width, CV_8UC1);
//    int a = m_Height;
//    int b = m_Width;
//    float c, d;
//    int value;
//    //计算出一个位元所占的像素值
//    //一共需要生成m_bit张图像
//    for (int k = 0; k < m_Bit; k++) {
//        //一共需要判断值的个数，就是2的N次幂
//        for (int m = 0; m < pow(2, m_Bit); m++) {
//            //因为将格雷码竖着向下读，所以先改变第一位
//            value = jugevalue(m_GrayArr[m][k]);
//            d = b / pow(2, m_Bit);
//            //修改像素值
//            for (int j = m * d; j < (m + 1) * d/*列数*/; j++) {
//                for (int i = 0; i < image.rows/*行数*/; i++) {
//                    image.at<uchar>(i, j) = value;
//                }
//            }
//            //test
//            //char filename[50];
//            // sprintf_s(filename,"C:\\Users\\86187\\Desktop\\Code\\1013buildgraypicture\\bmp\\%d.bmp",k);
//            // cv::imwrite(filename,image);
//            // cv::namedWindow("ouput",cv::WINDOW_NORMAL);
//
//        }
//        //test
//        //cv::imshow("dsdd", image);
//        //cv::waitKey(0);
//        //cv::imwrite("C:\\Users\\86187\\source\\repos\\GrayBuildCMAKE\\Result\\" + to_string(k) + ".bmp", image);
//        //浅拷贝，需要使用深拷贝的方式，建立图片的副本
//        M_GrayPicture.push_back(image.clone());
//    }
//    return M_GrayPicture;
//}
//
//void BuildGrayPattern::buildBlack(cv::Size Image_size) {
//    cv::Mat Black = cv::Mat::zeros(Image_size.height, Image_size.width, CV_8UC1);
//    for (int i = 0; i < Black.rows; i++) {
//        for (int j = 0; j < Black.cols; j++) {
//            Black.at<uchar>(i, j) = 0;
//        }
//    }
//    cv::imwrite("..\Projection_Pattern\black.bmp", Black);
//}
//
//void BuildGrayPattern::buildWhite(cv::Size Image_size) {
//    cv::Mat White = cv::Mat::zeros(Image_size.height, Image_size.width, CV_8UC1);
//    for (int i = 0; i < White.rows; i++) {
//        for (int j = 0; j < White.cols; j++) {
//            White.at<uchar>(i, j) = 255;
//        }
//    }
//    cv::imwrite("..\Projection_Pattern\white.bmp", White);
//}
//
//void BuildGrayPattern::buildBlack(const std::string& path) {
//    cv::Mat Black = cv::Mat::zeros(m_Height, m_Width, CV_8UC1);
//    for (int i = 0; i < Black.rows; i++) {
//        for (int j = 0; j < Black.cols; j++) {
//            Black.at<uchar>(i, j) = 0;
//        }
//    }
//    cv::imwrite(path + "black.bmp", Black);
//};
//
//void BuildGrayPattern::buildWhite(const std::string& path) {
//    cv::Mat White = cv::Mat::zeros(m_Height, m_Width, CV_8UC1);
//    for (int i = 0; i < White.rows; i++) {
//        for (int j = 0; j < White.cols; j++) {
//            White.at<uchar>(i, j) = 255;
//        }
//    }
//    cv::imwrite(path + "white.bmp", White);
//};
//
//
//void BuildGrayPattern::buildBlack(const std::string& path, int value) {
//    cv::Mat Black = cv::Mat::zeros(m_Height, m_Width, CV_8UC1);
//    for (int i = 0; i < Black.rows; i++) {
//        for (int j = 0; j < Black.cols; j++) {
//            Black.at<uchar>(i, j) = value;
//        }
//    }
//    cv::imwrite(path, Black);
//};
//void BuildGrayPattern::buildWhite(const std::string& path, int value) {
//    cv::Mat White = cv::Mat::zeros(m_Height, m_Width, CV_8UC1);
//    for (int i = 0; i < White.rows; i++) {
//        for (int j = 0; j < White.cols; j++) {
//            White.at<uchar>(i, j) = value;
//        }
//    }
//    cv::imwrite(path, White);
//};
//
//
//
//void BuildGrayPattern::SaveGray(const std::string path) {
//    string s_path = path;
//
//    for (int i = 0; i < M_GrayPicture.size(); i++) {
//        cv::imwrite(s_path + to_string(i) + ".bmp", M_GrayPicture[i]);
//    }
//}
//
//void BuildGrayPattern::BuildAndSave_GrayPattren(const std::string path) {
//    cv::Mat image = cv::Mat::zeros(m_Height, m_Width, CV_8UC1);
//    int a = m_Height;
//    int b = m_Width;
//    float c, d;
//    int value;
//    //Calculate the value of pixels occupied by a single bit,
//    //A total of m_bit images need to be generated
//    if (M_direction == 0) {
//        for (int k = 0; k < m_Bit; k++) {
//            //The total number of values that need to be judged is the Nth power of 2
//            for (int m = 0; m < pow(2, m_Bit); m++) {
//                //Because the Gray code is read vertically downwards, the first bit is changed first
//                value = jugevalue(m_GrayArr[m][k]);
//                d = b / pow(2, m_Bit);
//                //Modify pixel value
//                for (int j = m * d; j < (m + 1) * d/*列数*/; j++) {
//                    for (int i = 0; i < image.rows/*行数*/; i++) {
//                        image.at<uchar>(i, j) = value;
//                    }
//                }
//            }
//            //Shallow copy, you need to use a deep copy to create a copy of the image
//            M_GrayPicture.push_back(image.clone());
//        }
//    }
//    else {
//        for (int k = 0; k < m_Bit; k++) {
//            // 需要判断的值的总数为 2 的 m_Bit 次幂
//            for (int m = 0; m < pow(2, m_Bit); m++) {
//                // 因为格雷码是垂直向下读取的，所以首先改变第一个比特位
//                value = jugevalue(m_GrayArr[m][k]);
//                c = a / pow(2, m_Bit);
//                // 修改像素值
//                for (int i = m * c; i < (m + 1) * c; i++) {
//                    for (int j = 0; j < image.cols; j++) {
//                        image.at<uchar>(i, j) = value;
//                    }
//                }
//            }
//            // 浅拷贝，需要使用深拷贝来创建图像的副本
//            M_GrayPicture.push_back(image.clone());
//        }
//    }
//
//
//
//
//
//    for (int i = 0; i < M_GrayPicture.size(); i++) {
//        //这种格式下s_path就不需要加“”
//        cv::imwrite(path + to_string(i) + ".bmp", M_GrayPicture[i]);
//    }
//
//};
//
//
//
//
////以下是生成条纹的类
//vector<cv::Mat> CreateTW::creatP() {
//    for (int k = 0; k < N; k++) {
//        //cv_8uc1不需要加cv：：貌似是定义在全局的
//        cv::Mat image = cv::Mat::zeros(m_height, m_width, CV_8UC1);
//        //判断生成竖直方向条纹还是水平方向条纹
//        int a = (M_direction == 0 ? m_height : m_width);
//        int b = (M_direction == 0 ? m_width : m_height);
//        for (int i = 0; i < a; i++) {
//            for (int j = 0; j < b; j++) {
//                if (M_direction == 0) {
//                    //
//                    image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * f * j / m_width + 2 * pi * k / N);
//                }
//                else {
//                    image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * f * j / m_width + 2 * pi * k / N);
//                }
//            }
//        }
//        M_PhasePicture.push_back(image);
//    }
//    //test
//    //for (int i = 0; i < N; i++) {
//        //cv::imwrite("C:\\Users\\86187\\source\\repos\\GrayBuildCMAKE\\TW\\" + to_string(i) + ".bmp", imgs[i]);
//        //使用这个函数，读取图片存放的路径,这个函数又不行了
//        //sprintf_s(filename, "C:\\Users\\86187\\source\\repos\\GrayBuildCMAKE\\TW\\%d.bmp", i);
//        //cv::imwrite(filename, imgs[i]);
//    //}
//    return M_PhasePicture;
//}
//
//CreateTW::CreateTW(float f, int bit, int width, int height, int Direction) {
//    this->f = f;
//    this->m_bit = bit;
//    this->m_width = width;
//    this->m_height = height;
//    this->M_direction = Direction;
//}
//
//CreateTW::CreateTW(float f, int bit, enum DLPModel DLP, int direction) {
//    this->f = f;
//    this->m_bit = bit;
//    this->M_direction = direction;
//
//    switch (DLP)
//    {
//    case DLP4710:
//        this->m_width = 1920;
//        this->m_height = 1280;
//        break;
//    case DLP6500:
//        this->m_width = 1920;
//        this->m_height = 1080;
//        break;
//    case DLP3010:
//
//        break;
//    case DLP4500:
//        this->m_width = 912;
//        this->m_height = 1140;
//        break;
//    default:
//        break;
//    }
//};
//
//void CreateTW::SavePhase(std::string& path) {
//    for (int i = 0; i < M_PhasePicture.size(); i++) {
//        cv::imwrite(path + to_string(i) + ".bmp", M_PhasePicture[i]);
//    }
//
//}
//
//
//
//
////Generate and save phase shift patterns
//void CreateTW::BuildAndSave_phase_shift(const std::string& path) {
//    for (int k = 0; k < N; k++) {
//        //cv_8uc1 doesn't need to add cv:: seems to be defined globally
//        cv::Mat image = cv::Mat::zeros(m_height, m_width, CV_8UC1);
//        //Determine whether a vertical or horizontal stripe is generated.
//        int a = (M_direction == 0 ? m_height : m_width);
//        int b = (M_direction == 0 ? m_width : m_height);
//        //生成竖条纹
//        if (M_direction == 0) {
//            for (int i = 0; i < a; i++) {
//                for (int j = 0; j < b; j++) {
//                    if (M_direction == 0) {
//                        image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * f * j / m_width + 2 * pi * k / N);
//                    }
//                }
//            }
//        }
//        //生成列条纹
//        else {
//            for (int j = 0; j < a; j++) {
//                for (int i = 0; i < b; i++) {
//                    image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * f * i / m_width + 2 * pi * k / N);
//                }
//            }
//        }
//        M_PhasePicture.push_back(image);
//    }
//    for (int i = 0; i < M_PhasePicture.size(); i++) {
//        cv::imwrite(path + to_string(i) + ".bmp", M_PhasePicture[i]);
//    }
//};



//父类
GrayCodeOperation::GrayCodeOperation() {

};

GrayCodeOperation::~GrayCodeOperation() {

};


//子类1
GrayMethod_1::GrayMethod_1(int bit, int width, int height)
    :M_bit(bit), M_height(height), M_width(width)
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
            for (; j < M_bit && digit[j] == '0'; j++);
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

        for (int region = 0; region < (1 << M_bit); region++) {
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
        cv::namedWindow("GraycodePics", cv::WindowFlags::WINDOW_AUTOSIZE);
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

    //std::string SavePath = R"(D:\postgraduate\code\GrayCodeGenerate\GraycodePics\)";
    //for (int i = 0; i < bit; i++) {
    //    cv::imwrite(SavePath + std::to_string(i) + ".bmp", Images[i]);
    //}

    return Images;
};
/***********************方法1：改变最右侧位的值 + 改变右数第一个1的左侧位值*****************************/




/************************方法二：根据二进制码和格雷码之间的转换关系进行编码******************************/
GrayMethod_2::GrayMethod_2(int bit, enum DLPModel DLP, int dircetion)
    :M_bit(bit), M_Direction(dircetion)
{
    switch (DLP)
    {
    case DLP4710:
        this->M_Width = 1920;
        this->M_Height = 1280;
        break;
    case DLP6500:
        this->M_Width = 1920;
        this->M_Height = 1080;
        break;
    case DLP3010:
        this->M_Width = 1280;
        this->M_Height = 720;
        break;
    case DLP4500:
        this->M_Width = 912;
        this->M_Height = 1140;
        break;
    default:
        break;
    }

};

GrayMethod_2::~GrayMethod_2() {
    // 如果有需要释放的资源可以写这里，否则空的也可以
}

std::vector<int> GrayMethod_2::Generate_value() {
    std::vector<int> GrayCode_hex;
    for (int i = 0; i < (1 << M_bit); i++) {
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

    int region_width = M_Width / (1 << M_bit);  // 每个区域的宽度

    // 生成图像
    for (int BitIndex = 0; BitIndex < M_bit; BitIndex++) {
        cv::Mat img(M_Height, M_Width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像

        // 遍历所有区域，根据格雷码的当前位决定区域颜色
        for (int region = 0; region < (1 << M_bit); region++) {
            // 获取当前格雷码的指定位值
            int grayValue = (value[region] >> (M_bit - BitIndex - 1)) & 1;  // 获取对应的位
            grayValue = (grayValue == 1) ? 255 : 0;  // 根据格雷码位值设置颜色（1 -> 255, 0 -> 0）

            // 将当前区域的所有像素值设置为该值（黑或白）
            for (int row = 0; row < M_Height; row++) {
                for (int col = region * region_width; col < (region + 1) * region_width; col++) {
                    img.at<uchar>(row, col) = grayValue;
                }
            }
        }
        Images.push_back(img);  // 添加当前生成的图像到结果中
    }
    // 显示生成的图像
    ShowImages(Images);
    return Images;
};



std::vector<cv::Mat>GrayMethod_2::Generate_pics_reverse() {
    std::vector<int> value = Generate_value();
    std::vector<cv::Mat> Images;

    int region_width = M_Width / (1 << M_bit);  // 每个区域的宽度

    // 生成图像
    for (int BitIndex = 0; BitIndex < M_bit; BitIndex++) {
        cv::Mat img(M_Height, M_Width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像

        // 遍历所有区域，根据格雷码的当前位决定区域颜色
        for (int region = 0; region < (1 << M_bit); region++) {
            // 获取当前格雷码的指定位值
            int grayValue = (value[region] >> (M_bit - BitIndex - 1)) & 1;  // 获取对应的位
            grayValue = (grayValue == 1) ? 0 : 255;  // 根据格雷码位值设置颜色（1 -> 0, 0 -> 255）

            // 将当前区域的所有像素值设置为该值（黑或白）
            for (int row = 0; row < M_Height; row++) {
                for (int col = region * region_width; col < (region + 1) * region_width; col++) {
                    img.at<uchar>(row, col) = grayValue;
                }
            }
        }
        Images.push_back(img);  // 添加当前生成的图像到结果中
    }
    // 显示生成的图像
    ShowImages(Images);
    return Images;
};


std::vector<cv::Mat> GrayMethod_2::GenerateAndSave_pics(std::string& file_path) {
    std::vector<int> value = Generate_value();
    std::vector<cv::Mat> Images;

    int region_width = M_Width / (1 << M_bit);  // 每个区域的宽度

    // 生成图像
    for (int BitIndex = 0; BitIndex < M_bit; BitIndex++) {
        cv::Mat img(M_Height, M_Width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像

        // 遍历所有区域，根据格雷码的当前位决定区域颜色
        for (int region = 0; region < (1 << M_bit); region++) {
            // 获取当前格雷码的指定位值
            int grayValue = (value[region] >> (M_bit - BitIndex - 1)) & 1;  // 获取对应的位
            grayValue = (grayValue == 1) ? 255 : 0;  // 根据格雷码位值设置颜色（1 -> 255, 0 -> 0）

            // 将当前区域的所有像素值设置为该值（黑或白）
            for (int row = 0; row < M_Height; row++) {
                for (int col = region * region_width; col < (region + 1) * region_width; col++) {
                    img.at<uchar>(row, col) = grayValue;
                }
            }
        }

        Images.push_back(img);  // 添加当前生成的图像到结果中
    }
    // 显示生成的图像
    ShowImages(Images);
    //保存生成的图像
    M_ImageOperation.SavePics(Images, file_path);
    return Images;
};


std::vector<cv::Mat> GrayMethod_2::GenerateAndSave_pics_reverse(std::string& file_path) {
    std::vector<int> value = Generate_value();
    std::vector<cv::Mat> Images;

    int region_width = M_Width / (1 << M_bit);  // 每个区域的宽度

    // 生成图像
    for (int BitIndex = 0; BitIndex < M_bit; BitIndex++) {
        cv::Mat img(M_Height, M_Width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像

        // 遍历所有区域，根据格雷码的当前位决定区域颜色
        for (int region = 0; region < (1 << M_bit); region++) {
            // 获取当前格雷码的指定位值
            int grayValue = (value[region] >> (M_bit - BitIndex - 1)) & 1;  // 获取对应的位
            grayValue = (grayValue == 1) ? 0 : 255;  // 根据格雷码位值设置颜色（1 -> 0, 0 -> 255）

            // 将当前区域的所有像素值设置为该值（黑或白）
            for (int row = 0; row < M_Height; row++) {
                for (int col = region * region_width; col < (region + 1) * region_width; col++) {
                    img.at<uchar>(row, col) = grayValue;
                }
            }
        }

        Images.push_back(img);  // 添加当前生成的图像到结果中
    }
    // 显示生成的图像
    ShowImages(Images);
    //保存生成的图像
    M_ImageOperation.SavePics(Images, file_path);
    return Images;
};





cv::Mat GrayMethod_2::GeneratePatternWithValue(const std::string& path, int value) {
    cv::Mat Image(M_Height, M_Width, CV_8UC1, cv::Scalar(value)); // 注意高是行数，宽是列数
    cv::imwrite(path + "Image.bmp", Image);
    return Image;
}


void GrayMethod_2::generationDiffGrayValuePattern(const std::string& filepath,
                                                  int step) {
    for (int i = 0; i < 255; i += step) {
        GeneratePatternWithValue(filepath, i);
    }
};


std::vector<cv::Mat> GrayMethod_2::Generate_pics_block(int bit, int rowBlock, int colBlock) {

    //根据需要生成的区域传入bit值,生成格雷码序列
    std::vector<int> GrayCode_hex;
    for (int i = 0; i < (1 << bit); i++) {
        GrayCode_hex.push_back(i ^ (i >> 1));
    }

    std::vector<cv::Mat> Images;
    int region_width = M_Width / colBlock;
    int region_height = M_Height / rowBlock;
    int grayValue = 0;
    int initial_value = 1;//初始值不是从零开始，而是从一开始
    // 生成图像
    for (int BitIndex = 0; BitIndex < bit; BitIndex++) {
        cv::Mat img(M_Height, M_Width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像
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
    ShowImages(Images);
    return Images;
};



std::vector<cv::Mat> GrayMethod_2::GenerateAndSave_pics_block(std::string& file_path, int bit, int rowBlock, int colBlock) {
    //根据需要生成的区域传入bit值,生成格雷码序列
    std::vector<int> GrayCode_hex;
    for (int i = 0; i < (1 << bit); i++) {
        GrayCode_hex.push_back(i ^ (i >> 1));
    }

    std::vector<cv::Mat> Images;
    int region_width = M_Width / colBlock;
    int region_height = M_Height / rowBlock;
    int grayValue = 0;
    int initial_value = 1;//初始值不是从零开始，而是从一开始
    // 生成图像
    for (int BitIndex = 0; BitIndex < bit; BitIndex++) {
        cv::Mat img(M_Height, M_Width, CV_8UC1, cv::Scalar(255));  // 初始化为白色图像
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
    ShowImages(Images);
    //保存图像
    M_ImageOperation.SavePics(Images, file_path);
    return Images;
};





//初始方法
cv::Mat GrayMethod_2::Decode_Pics(std::string& file_path, int PicNum) {
    //读取格雷码图像
    std::vector<cv::Mat> GrayCodeSequence  = M_ImageOperation.ReadPics(file_path);
    //解码图像
    cv::Mat decodedImage = cv::Mat::zeros(M_Height, M_Width, CV_8UC1);
    for (int row = 0; row < M_Height; row++) {
        for (int col = 0; col < M_Width; col++) {
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


cv::Mat GrayMethod_2::Decode_PicsPositiveAndNegative(std::string& file_path, int PicNum) {
    //读取格雷码图像
    std::vector<cv::Mat> GrayCodeSequence = M_ImageOperation.ReadPics(file_path);

    //解码图像
    cv::Mat decodedImage = cv::Mat::zeros(M_Height, M_Width, CV_8UC1);
    for (int row = 0; row < M_Height; row++) {
        for (int col = 0; col < M_Width; col++) {
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













/************************方法二：根据二进制码和格雷码之间的转换关系进行编码******************************/




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




/*****************************************  生成条纹类  *********************************************/

PhasePatternOperation::PhasePatternOperation(float f,
                                             enum DLPModel DLP, 
                                             int direction ,
                                             int num)
    :M_F(f),M_Direction(direction),M_PhaseShiftingNum(num){
    switch (DLP)
    {
    case DLP4710:
        this->M_Width = 1920;
        this->M_Height = 1280;
        break;
    case DLP6500:
        this->M_Width = 1920;
        this->M_Height = 1080;
        break;
    case DLP3010:
        this->M_Width = 1280;
        this->M_Height = 720;
        break;
    case DLP4500:
        this->M_Width = 912;
        this->M_Height = 1140;
        break;
    default:
        break;
    }
};


PhasePatternOperation::~PhasePatternOperation() {


};


std::vector<cv::Mat> PhasePatternOperation::CreatPhasePatterns() {
    for (int k = 0; k < M_PhaseShiftingNum; k++) {
        //cv_8uc1不需要加cv：：貌似是定义在全局的
        cv::Mat image = cv::Mat::zeros(M_Height, M_Width, CV_8UC1);
        //判断生成竖直方向条纹还是水平方向条纹
        int a = (M_Direction == 0 ? M_Height : M_Width);
        int b = (M_Direction == 0 ? M_Width : M_Height);
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < b; j++) {
                if (M_Direction == 0) {
                    image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * M_F * j / M_Width + 2 * pi * k / M_PhaseShiftingNum);
                }
                else {
                    image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * M_F * j / M_Width + 2 * pi * k / M_PhaseShiftingNum);
                }
            }
        }
        M_PhasePicture.push_back(image);
    }
    return M_PhasePicture;
};

//Generate and save phase shift patterns
std::vector<cv::Mat> PhasePatternOperation::BuildAndSave_PhaseShiftingPattern(const std::string& path) {
    for (int k = 0; k < M_PhaseShiftingNum; k++) {
        //cv_8uc1不需要加cv：：貌似是定义在全局的
        cv::Mat image = cv::Mat::zeros(M_Height, M_Width, CV_8UC1);
        //判断生成竖直方向条纹还是水平方向条纹
        int a = (M_Direction == 0 ? M_Height : M_Width);
        int b = (M_Direction == 0 ? M_Width : M_Height);
        //生成竖条纹
        if (M_Direction == 0) {
            for (int i = 0; i < a; i++) {
                for (int j = 0; j < b; j++) {
                    if (M_Direction == 0) {
                        image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * M_F * j / M_Width + 2 * pi * k / M_PhaseShiftingNum);
                    }
                }
            }
        }
        //生成列条纹
        else {
            for (int j = 0; j < a; j++) {
                for (int i = 0; i < b; i++) {
                    image.at<uchar>(i, j) = 127.5 + 127.5 * cos(2 * pi * M_F * i / M_Width + 2 * pi * k / M_PhaseShiftingNum);
                }
            }
        }
        M_PhasePicture.push_back(image);
    }
    ShowPhaseImages(M_PhasePicture);
    M_ImageOperation.SavePics(M_PhasePicture, path);
    return M_PhasePicture;
};




/*****************************************  生成条纹类  *********************************************/

#include "GrayCodeGenerate.h"


int main() {

	GrayCodeOperation* TEST = new GrayMethod_2(7,1024,768);
	//TEST->Generate_pics_block(7,10,12);
	std::string filepath = R"(D:\postgraduate\code\GrayCodeGenerate\GraycodePics\)";
	TEST->Decode_Pics(filepath, 7);
	return 0;

}
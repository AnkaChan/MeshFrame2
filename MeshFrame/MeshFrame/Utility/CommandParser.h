#ifndef __COMMAND_PARSER__
#define __COMMAND_PARSER__
#include "cxxopts.hpp"

#define REQUIRE_NUM_INPUTS(n)\
if (argc < n+1) {\
        std::cout << "[Error] Need at least " << n << " input parameters!" << std::endl; \
        return -1;\
}

struct BaseCmdArgsParser {
	//bool checkFileNamesMatch = true;
	//int maxNumFilesToProcess = -1;

	BaseCmdArgsParser() :
		options("CalibConfig", "Configs for calibration application.")
	{
		options.add_options()
			//("c,checkFileNamesMatch", "Sort each lists of files and only keep file with name exists in all path lists.", cxxopts::value<bool>(checkFileNamesMatch))
			//("m, maxNumFilesToProcess", "Maximum number of files to process by the calibrator.", cxxopts::value<int>(maxNumFilesToProcess))
			;
	}

	void parse(int argc, char** argv) {
		try
		{
			auto result = options.parse(argc, argv);
		}
		catch (const cxxopts::OptionException& e)
		{
			std::cout << "Error parsing options: " << e.what() << std::endl;
			std::cout << options.help();
			exit(1);
		}
	}

	cxxopts::Options options;


};

#endif // !__COMMAND_PARSER__

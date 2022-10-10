#ifndef __PARSER__
#define __PARSER__

#define PUT_TO_JSON(j, arg) j[#arg] = arg
#define EXTRACT_FROM_JSON(j, arg) parseJsonParameters(j, #arg, arg)

#include <iostream>
#include "../Json/json.hpp"
#define REQUIRE_NUM_INPUTS(n)\
if (argc < n+1) {\
        std::cout << "[Error] Need at least " << n << " input parameters!" << std::endl; \
        return -1;\
}

namespace MF {
	inline bool loadJson(std::string  filePath, nlohmann::json& j) {
		std::ifstream ifs(filePath);
		if (ifs.is_open())
		{
			try
			{
				ifs >> j;
				return true;
			}
			catch (nlohmann::json::exception& e)
			{
				std::cout << e.what() << '\n';
				return false;
			}
		}
		else
		{
			std::cout << "Fail to open: " << filePath << '\n';
			return false;
		}

	}

	inline nlohmann::json tryGetJson(nlohmann::json j, std::string  name) {
		nlohmann::json param;
		try
		{
			param = j.at(name);
		}
		catch (nlohmann::json::exception& e)
		{
			std::cout << e.what() << '\n';
			std::cout << name << " does not exist in json and will be set to default value." << std::endl;
			return false;
		}
		return param;
	}

	template<typename T>
	bool parseJsonParameters(nlohmann::json j, const std::string& name, T& param) {
		try
		{
			nlohmann::json paramJson = j.at(name);
			param = paramJson.get<T>();
		}
		catch (nlohmann::json::exception& e)
		{
			std::cout << e.what() << '\n';
			std::cout << name << " does not exist in json and will be set to default value." << std::endl;
			return false;
		}
		return true;
	}

	struct BaseJsonConfig {

		virtual bool fromJson(nlohmann::json& j) = 0;
		virtual bool toJson(nlohmann::json& j) = 0;

		bool loadFromJsonFile(const std::string& inJsonParamFile) {
			nlohmann::json j;
			if (loadJson(inJsonParamFile, j))
			{
				if (fromJson(j))
				{

					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				return false;
			}
		}

		bool writeToJsonFile(std::string filePath, int indent = -1) {
			std::ofstream ofs(filePath);
			if (ofs.is_open())
			{
				nlohmann::json j;
				toJson(j);

				ofs << j.dump(indent);
				return true;
			}
			else
			{
				return false;
			}
		}
	};

}




#endif // !__COMMAND_PARSER__

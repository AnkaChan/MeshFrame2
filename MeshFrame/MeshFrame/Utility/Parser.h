#ifndef __PARSER__
#define __PARSER__

#define PUT_TO_JSON(j, arg) j[#arg] = arg
#define EXTRACT_FROM_JSON(j, arg) MF::parseJsonParameters(j, #arg, arg)

#include <iostream>
#include <fstream>
#include "../Json/json.hpp"
#define REQUIRE_NUM_INPUTS(n)\
if (argc < n+1) {\
        std::cout << "[Error] Need at least " << n << " input parameters!" << std::endl; \
        return -1;\
}

#define PARSE_VEC3_INDEXED_BY_ROUND_BRACKET(jsonParam, vec3)\
std::array<double, 3> vec3##Arr;\
if (MF::parseJsonParameters(jsonParam, #vec3, vec3##Arr))\
{\
	vec3(0) = vec3##Arr[0];\
	vec3(1) = vec3##Arr[1];\
	vec3(2) = vec3##Arr[2];\
}

#define PARSE_VEC3_INDEXED_BY_SQUARE_BRACKET(jsonParam, vec3)\
std::array<double, 3> vec3##Arr;\
if (MF::parseJsonParameters(jsonParam, #vec3, vec3##Arr))\
{\
	vec3[0] = vec3##Arr[0];\
	vec3[1] = vec3##Arr[1];\
	vec3[2] = vec3##Arr[2];\
}\

#define PUT_TO_JSON_VEC3_INDEXED_BY_ROUND_BRACKET(jsonParam, vec3)\
std::array<double, 3> vec3##yArr;\
vec3##yArr[0] = vec3(0);\
vec3##yArr[1] = vec3(1);\
vec3##yArr[2] = vec3(2);\
jsonParam[#vec3] = vec3##yArr;

#define PUT_TO_JSON_VEC3_INDEXED_BY_SQUARE_BRACKET(jsonParam, vec3)\
std::array<double, 3> vec3##yArr;\
vec3##yArr[0] = vec3[0];\
vec3##yArr[1] = vec3[1];\
vec3##yArr[2] = vec3[2];\
jsonParam[#vec3] = vec3##yArr;

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

	inline bool saveJson(std::string  filePath, nlohmann::json& j, int indent = -1) {
		std::ofstream ofs(filePath);
		if (ofs.is_open())
		{
			try
			{
				ofs << j.dump(indent) << std::endl;
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
			std::cout << "When loading: " << name << "\n"
			    << e.what() << "\n"
				<< name << " does not exist in json and will be set to default value." << std::endl;
			return false;
		}
		return true;
	}

	template<typename T>
	bool convertJsonParameters(nlohmann::json j, T& param) {
		try
		{
			param = j.get<T>();
		}
		catch (nlohmann::json::exception& e)
		{
			std::cout<< e.what() << "\n";
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

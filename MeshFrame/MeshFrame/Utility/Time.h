#ifndef __MESHFRAME_TIME_H__
#define __MESHFRAME_TIME_H__

#include <ctime>
#include <chrono>

namespace MF {
	struct Time
	{
		static std::string getCurTimeString() {
			std::time_t curtime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
			std::string timeStr = std::ctime(&curtime);
			// remove the '\n' at end
			timeStr.pop_back();

			timeStr = STR::replace(timeStr, " ", "_");
			timeStr = STR::replace(timeStr, ":", "-");

			return timeStr;
		}
	};
}

#endif // !__MESHFRAME_TIME_H__

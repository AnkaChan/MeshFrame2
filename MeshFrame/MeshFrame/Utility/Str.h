#ifndef _MESHFRAME_STR_H_
#define _MESHFRAME_STR_H_
#include <vector>
#include <string>

namespace MF {
	struct STR
	{
		static void split(std::vector<std::string>& result, const char* str, char c = ' ')
		{
			do
			{
				const char* begin = str;

				while (*str != c && *str)
					str++;

				result.push_back(std::string(begin, str));
			} while (0 != *str++);

			return;
		}

		static void padTo(std::string& str, const size_t num, const char paddingChar = ' ')
		{
			if (num > str.size())
				str.insert(0, num - str.size(), paddingChar);
		}

		static std::string replace(std::string& inStr, const std::string subStrToReplace, const std::string replacement, bool replaceAll = true) {
			std::string outStr;
			size_t index = 0;
			outStr = inStr;
			do {
				/* Locate the substring to replace. */
				index = outStr.find(subStrToReplace, index);
				if (index == std::string::npos) break;

				/* Make the replacement. */
				outStr.replace(index, subStrToReplace.size(), replacement);

				/* Advance index forward so the next iteration doesn't pick it up as well. */
				index += 3;
			} while (replaceAll);

			return outStr;
		};
	};
}
#endif
#pragma once

#include <algorithm>
#include <vector>
#include <string>

//#ifdef _WIN32
//#include <io.h>
//#endif

namespace MF {
	struct IO
	{
		struct FileParts
		{
			FileParts() {}
			FileParts(std::string filename) {
				std::replace(filename.begin(), filename.end(), '\\', '/'); // replace all '\' to '/', 

				int idx0 = filename.rfind("/");

				int idx1 = filename.rfind(".");

				path = filename.substr(0, idx0 + 1);
				name = filename.substr(idx0 + 1, idx1 - idx0 - 1);
				if (idx1 != -1) {
					ext = filename.substr(idx1);
				}
				else
				{
					ext = "";
				}
			}
			std::string path;
			std::string name;
			std::string ext;
		};

		static FileParts fileparts(std::string filename)
		{
			std::replace(filename.begin(), filename.end(), '\\', '/'); // replace all '\' to '/', 

			int idx0 = filename.rfind("/");

			int idx1 = filename.rfind(".");

			FileParts fp;
			fp.path = filename.substr(0, idx0 + 1);
			fp.name = filename.substr(idx0 + 1, idx1 - idx0 - 1);
			if (idx1 != -1) {
				fp.ext = filename.substr(idx1);
			}
			else
			{
				fp.ext = "";
			}

			return fp;
		}

		// Returns false on success, true on error
		static bool createFolder(std::string folderName);
		static bool folderExists(const char* folderName);

		static bool copyFile(std::string src, std::string dst);

	};
}
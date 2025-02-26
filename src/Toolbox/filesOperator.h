#pragma once
#include <string>
#include <vector>
#include <io.h>
void getFiles(std::string path, std::vector<std::string>& files)
{
	//文件句柄
	intptr_t  hFile = 0;
	//文件信息
	struct _finddata_t fileinfo;
	std::string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//如果是目录,迭代之
			//如果不是,加入列表
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

#include <qstring.h>
void truncateFilePath(std::string &file)
{
	QString fileName = QString::fromStdString(file);
	int id = fileName.lastIndexOf("/");
	if (id == -1)
	{
		id = fileName.lastIndexOf("\\");
		if (id == -1)
			return;
	}
	file = fileName.right(fileName.length() - id - 1).toLatin1().data();
}

void truncateFileExtension(std::string &file)
{
	QString fileName = QString::fromStdString(file);
	int id = fileName.lastIndexOf(".");
	if (id == -1)
		return;
	fileName.truncate(id);
	file = fileName.toLatin1().data();
}
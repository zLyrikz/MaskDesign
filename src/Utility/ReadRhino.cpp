#include "ReadRhino.h"
#include <fstream>
#include <sstream>
#include <iostream>
using std::fstream;
using std::ifstream;
using std::stringstream;
using std::string;
using std::cout;
using std::endl;

ReadRhino::ReadRhino()
{
}

void ReadRhino::Read3DPoints(const std::string _file, std::vector<Eigen::VectorXd>& _points)
{
	ifstream file(_file);
	char line[1024] = { 0 };
	if (file.is_open())
	{
		_points.clear();

		while (file.getline(line, sizeof(line)))
		{
			Eigen::VectorXd point(3);
			stringstream word(line);
			word >> point(0);
			word >> point(1);
			word >> point(2);
			_points.push_back(point);
		}
		file.close();
	}
}

void ReadRhino::Read3DPoints_Brackets(const std::string _file, std::vector<Eigen::VectorXd>& _points, bool _rotate)
{
	ifstream file(_file);
	string line;
	if (file.is_open())
	{
		_points.clear();

		while (getline(file, line))
		{
			// line = {x, y, z}
			std::istringstream words(line);
			string x;// {x,
			string y;// y,
			string z;// z}
			words >> x;
			words >> y;
			words >> z;
			x = x.substr(1, x.length() - 2);
			y = y.substr(0, y.length() - 1);
			z = z.substr(0, z.length() - 1);

			Eigen::VectorXd point(3);
			std::istringstream word_x(x);
			std::istringstream word_y(y);
			std::istringstream word_z(z);
			word_x >> point(0);
			word_y >> point(1);
			word_z >> point(2);

			//cout << point << endl;
			if (_rotate)
			{
				double old_y = point(1);
				point(1) = point(2);
				point(2) = -old_y;
			}

			_points.push_back(point);
		}


		file.close();
	}
}

void ReadRhino::Read3DPoints_Brackets(const std::string _file, std::vector<Mesh::Point>& _points, bool _rotate)
{
	ifstream file(_file);
	string line;
	if (file.is_open())
	{
		_points.clear();

		while (getline(file, line))
		{
			// line = {x, y, z}
			std::istringstream words(line);
			string x;// {x,
			string y;// y,
			string z;// z}
			words >> x;
			words >> y;
			words >> z;
			x = x.substr(1, x.length() - 2);
			y = y.substr(0, y.length() - 1);
			z = z.substr(0, z.length() - 1);
			Mesh::Point point;
			std::istringstream word_x(x);
			std::istringstream word_y(y);
			std::istringstream word_z(z);
			word_x >> point[0];
			word_y >> point[1];
			word_z >> point[2];

			//cout << point << endl;
			if (_rotate)
			{
				double old_y = point[1];
				point[1] = point[2];
				point[2] = -old_y;
			}
			_points.push_back(point);
		}


		file.close();
	}
}

void ReadRhino::ReadVector(const std::string& _file, double* _x, int _dim) const
{
	ifstream file(_file);
	char line[10000] = { 0 };
	if (file.is_open() && _dim * 11 < 10000)
	{
		file.getline(line, sizeof(line));	

		stringstream word(line);

		for (int iDim = 0; iDim < _dim; ++iDim)
		{
			word >> _x[iDim];
		}
		file.close();
	}
	else
	{
		cout << "[WARNING ReadRhino::ReadVector]file cannot open or vector too large" << endl;
	}
}

void ReadRhino::ReversePoints(std::vector<Eigen::VectorXd>& _points) const
{
	std::vector<Eigen::VectorXd> points(_points.size());
	for (int iPoint = 0; iPoint < _points.size(); ++iPoint)
	{
		points[iPoint] = _points[_points.size() - 1 - iPoint];
	}
	_points = points;
}

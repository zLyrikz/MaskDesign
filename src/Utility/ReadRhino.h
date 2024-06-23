#pragma once
#include <string>
#include <vector>
#include <Eigen/Core>
#include "../MeshViewer/MeshDefinition.h"

// read data saved from Rhino exported files
// also for reading txt file data 

class ReadRhino
{
public:
	ReadRhino();

	// file.txt 
	// points separated by space
	// x1 y1 z1
	// x2 y2 z2
	// ...
	void Read3DPoints(const std::string _file, std::vector<Eigen::VectorXd>& _points);
	// file.txt 
	// points separated by comma and space and with a bracket
	// {x1, y1, z1}
	// {x2, y2, z2}
	// ...
	// set rotate parameter = true if rotate around x axis make z axis point to y axis
	void Read3DPoints_Brackets(const std::string _file, std::vector<Eigen::VectorXd>& _points, bool _rotate = false);
	void Read3DPoints_Brackets(const std::string _file, std::vector<Mesh::Point>& _points, bool _rotate = false);// TODO: learn how to combine Eigen::VectorXd, Mesh::Point with template type

	// file.txt
	// a vector each element separated by space
	// x0 x1 ... xn
	void ReadVector(const std::string& _file, double* _x, int _dim) const;

	void ReversePoints(std::vector<Eigen::VectorXd>& _points) const;
};


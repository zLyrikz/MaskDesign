#pragma once
#include "../MaskInterface/DiscreteCushionPiece.h"
#include <vector>

class FrameRelax
{
public:
	FrameRelax();
	~FrameRelax();

	// input the cushion surface wire frame, the intersection info, and the expand time to each side
	// can achieve global relaxation if choose large expand to make a interval cover the whole frames 
	void Relax(std::vector<DiscreteCushionPiece>& _cross_sections, std::vector<std::pair<int, int>>& _intervals, int _side_expand = 2, bool _extra_iteration = false) const;
	// global relaxation, i.e. intervals cover all cross sections
	// also, all cross sections are set as circles
	void GlobalCircleRelax(float _radius, std::vector<DiscreteCushionPiece>& _cross_sections, bool _extra_iteration = false) const;

public:
	//utility function

	// relax one time for one interval
	void OneTimeRelax(std::vector<DiscreteCushionPiece>& _cross_sections, const std::pair<int, int>& _interval) const;
	// global relax one time
	void OneTimeGlobalRelax_FixFrame0(std::vector<DiscreteCushionPiece>& _cross_sections) const;
	// global relax one time, periodically relax all frames, no frame is fixed
	void OneTimeGlobalRelax(std::vector<DiscreteCushionPiece>& _cross_sections) const;

	void RotateDiscreteCushionPiece(const Eigen::Vector3d& _new_tangent, DiscreteCushionPiece& _cross_section) const;

	void ExpandIntervals(int total_frame, std::vector<std::pair<int, int>>& _intervals, int _side_expand) const;
};


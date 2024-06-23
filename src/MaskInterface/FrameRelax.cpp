#include "FrameRelax.h"
#include "SelfIntersectionTest.h"
#include "../MeshProcessing/FindRotation.h"
#include "../Utility/DivisionWithReminder.h"
#include <iostream>
using std::cout;
using std::endl;
FrameRelax::FrameRelax()
{
}

FrameRelax::~FrameRelax()
{
}

void FrameRelax::Relax(std::vector<DiscreteCushionPiece>& _cross_sections, std::vector<std::pair<int, int>>& _intervals, int _side_expand, bool _extra_iteration) const
{
	SelfIntersectionTest check_self_intersection;
	//cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
	//cout << "relaxation starts..."  << endl;
	//cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
	
	ExpandIntervals(_cross_sections.size(), _intervals, _side_expand);

	int counter = 0;
	int max_iteration = 50000;
	std::vector<std::pair<int, int>> intervals = _intervals;
	do
	{
		for (int iPart = 0; iPart < _intervals.size(); ++iPart)
		{
			OneTimeRelax(_cross_sections, _intervals[iPart]);// always use the original input interval for relaxation
		}
		check_self_intersection.TestIntersection(_cross_sections, intervals);
		++counter;
	} while (intervals.size() > 0 && counter < max_iteration);

	// try to do some extra iteration to achieve better smoothness
	if (_extra_iteration && counter < max_iteration)
	{
		for (int iExtra = 0; iExtra < counter * 2.0; ++iExtra)
		{
			for (int iPart = 0; iPart < _intervals.size(); ++iPart)
			{
				OneTimeRelax(_cross_sections, _intervals[iPart]);
			}
		}
	}

	//cout << "iterate time =" << (_extra_iteration ? int(counter + counter * 2.0) : counter) << endl;

}

void FrameRelax::GlobalCircleRelax(float _radius, std::vector<DiscreteCushionPiece>& _cross_sections, bool _extra_iteration) const
{
	SelfIntersectionTest check_self_intersection;
	//cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
	//cout << "relaxation starts..." << endl;


	int counter = 0;
	int max_iteration = 2000;
	std::vector<std::pair<int, int>> intervals;
	check_self_intersection.TestCircleIntersection(_radius, _cross_sections, intervals);
	while (intervals.size() > 0 && counter < max_iteration)
	{
		OneTimeGlobalRelax(_cross_sections);
		
		// update intersect intervals
		check_self_intersection.TestCircleIntersection(_radius, _cross_sections, intervals);
		++counter;
	} 

	// try to do some extra iteration to achieve better smoothness
	if (_extra_iteration && counter < max_iteration)
	{
		for (int iExtra = 0; iExtra < counter * 2.0; ++iExtra)
		{
			OneTimeGlobalRelax(_cross_sections);

		}
	}

	//cout << "iterate time =" << (_extra_iteration ? int(counter + counter * 2.0) : counter) << endl;
	//cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;

}


void FrameRelax::OneTimeRelax(std::vector<DiscreteCushionPiece>& _cross_sections, const std::pair<int, int>& _interval) const
{
	int begin_idx = _interval.first;
	int end_idx = _interval.second;
	if (end_idx < begin_idx)// periodic case
	{
		end_idx += _cross_sections.size();
	}

	int num_frame = end_idx - begin_idx + 1;
	std::vector<Eigen::Vector3d> original_tangents(num_frame);
	for (int iFrame = 0; iFrame < num_frame; ++iFrame)
	{
		int current_idx = (begin_idx + iFrame) % _cross_sections.size();
		int next_idx = (current_idx + 1) % _cross_sections.size();

		// record the original tangent before change it
		original_tangents[iFrame] = _cross_sections[current_idx].rotated_tangent_;

		if (iFrame > 0)
		{
			// average of last and next
			Eigen::Vector3d new_tangent = ((original_tangents[iFrame - 1] + _cross_sections[next_idx].rotated_tangent_) / 2.0).normalized();
			RotateDiscreteCushionPiece(new_tangent, _cross_sections[current_idx]);
		}
		else
		{
			// average of last and next
			int last_idx = (current_idx - 1 >= 0) ? (current_idx - 1) : (_cross_sections.size() - 1);
			Eigen::Vector3d new_tangent = ((_cross_sections[last_idx].rotated_tangent_ + _cross_sections[next_idx].rotated_tangent_) / 2.0).normalized();
			RotateDiscreteCushionPiece(new_tangent, _cross_sections[current_idx]);
		}

	}
}

void FrameRelax::OneTimeGlobalRelax_FixFrame0(std::vector<DiscreteCushionPiece>& _cross_sections) const
{
	int num_frame = _cross_sections.size() - 1;// fix only one frame(frame 0), relax all others
	int start_frame = 1;// cannot be changed for the current implementation(due to the index thing)
	std::vector<Eigen::Vector3d> original_tangents(num_frame);
	for (int iFrame = 0; iFrame < num_frame; ++iFrame)
	{
		int current_idx = iFrame + start_frame;// ranging from 1 to _cross_sections.size() - 1
		int last_idx = current_idx - 1;
		int next_idx = (current_idx + 1) % _cross_sections.size();

		// record the original tangent before change it
		original_tangents[iFrame] = _cross_sections[current_idx].rotated_tangent_;
		Eigen::Vector3d new_tangent;// average of last and next
		if (iFrame == 0)
		{
			new_tangent = ((_cross_sections[last_idx].rotated_tangent_ + _cross_sections[next_idx].rotated_tangent_) / 2.0).normalized();
		}
		else
		{
			new_tangent = ((original_tangents[iFrame - 1] + _cross_sections[next_idx].rotated_tangent_) / 2.0).normalized();
		}


		RotateDiscreteCushionPiece(new_tangent, _cross_sections[current_idx]);
	}
}

void FrameRelax::OneTimeGlobalRelax(std::vector<DiscreteCushionPiece>& _cross_sections) const
{
	int num_frame = _cross_sections.size();
	std::vector<Eigen::Vector3d> original_tangents(num_frame);// record original tangents
	for (int iFrame = 0; iFrame < num_frame; ++iFrame)
	{
		original_tangents[iFrame] = _cross_sections[iFrame].rotated_tangent_;
	}
	for (int iFrame = 0; iFrame < num_frame; ++iFrame)
	{
		DivisionWithReminder find_reminder;
		int last_idx = find_reminder.getReminder(iFrame - 1, _cross_sections.size());
		int next_idx = (iFrame + 1) % _cross_sections.size();

		Eigen::Vector3d new_tangent;// average of last and next using the original tangents
		new_tangent = ((original_tangents[last_idx] + original_tangents[next_idx]) / 2.0).normalized();
		
		RotateDiscreteCushionPiece(new_tangent, _cross_sections[iFrame]);
	}
}


void FrameRelax::RotateDiscreteCushionPiece(const Eigen::Vector3d& _new_tangent, DiscreteCushionPiece& _cross_section) const
{
	Eigen::AngleAxisd rotation;
	FindRotation find_rotation;
	find_rotation.VectorA2DirectionB_3D(_cross_section.rotated_tangent_, _new_tangent, rotation);

	_cross_section.rotation_ = rotation * _cross_section.rotation_;
	_cross_section.rotated_tangent_ = _new_tangent;

	// rotate the discrete cross section
	for (auto& point : _cross_section.discrete_cross_section_)
	{
		point = rotation * (point - _cross_section.frame_value_) + _cross_section.frame_value_;
	}
}

void FrameRelax::ExpandIntervals(int total_frame, std::vector<std::pair<int, int>>& _intervals, int _side_expand) const
{
	int num_interval = _intervals.size();
	if (num_interval > 0 && _side_expand >= 1)
	{
		// this indicate one index that's not having intersection with others
		// will be used if some interval cover all indexs
		// i.e. return one interval: [free_index_begin, free_index_end]
		int free_index_begin = (_intervals[0].second + 1) % total_frame;
		int free_index_end = _intervals[0].second;

		// get expended begin and end index
		for (auto& interval : _intervals)
		{
			int num_frame = interval.second - interval.first + 1;
			if (num_frame < 0)
			{
				interval.first -= total_frame;// convert this to the "true index", not the index of the vetor
				num_frame += total_frame;
			}

			interval.first -= num_frame * _side_expand;
			interval.second += num_frame * _side_expand;
		}

		// merge overlap intervals
		for (std::vector<std::pair<int, int>>::iterator iter = _intervals.begin(); iter != (_intervals.end() - 1); )
		{
			if (iter->second >= (iter + 1)->first)// current one overlap with the next
			{
				// merge current one with the next
				if (iter->first > (iter + 1)->first)
				{
					iter->first = (iter + 1)->first;// choose the smaller one
				}
				if (iter->second < (iter + 1)->second)
				{
					iter->second = (iter + 1)->second;// choose the larger one
				}

				// if the interval covers all frames, choose a decent interval and return 
				if (iter->second - iter->first + 1 >= total_frame)
				{
					_intervals = std::vector<std::pair<int, int>>{ std::make_pair(free_index_begin, free_index_end) };
					return;
				}

				iter = _intervals.erase(iter + 1);// erase the next one
				--iter;// pull iter back to the current interval
			}
			else
			{
				++iter;
			}
		}
		// check if last interval and the first overlap
		if (_intervals.size() > 1)
		{
			// pull the last interval one period back, then compare it with the first
			int last_end = _intervals.back().second - total_frame;
			int last_begin = _intervals.back().first - total_frame;
			int first_begin = _intervals[0].first;

			if (last_end >= first_begin)// overlap
			{
				//merge last to the first
				if (_intervals[0].first > last_begin)
				{
					_intervals[0].first = last_begin;
				}
				if (_intervals[0].second < last_end)
				{
					_intervals[0].second = last_end;
				}
				_intervals.pop_back();
			}
			// if merged interval covers all frames, choose a decent interval and return 
			if (_intervals[0].second - _intervals[0].first + 1 >= total_frame)
			{
				_intervals = std::vector<std::pair<int, int>>{ std::make_pair(free_index_begin, free_index_end) };
				return;
			}
		}
		else if (_intervals.size() == 1)
		{
			if (_intervals[0].second - _intervals[0].first + 1 >= total_frame)
			{
				_intervals = std::vector<std::pair<int, int>>{ std::make_pair(free_index_begin, free_index_end) };
				return;
			}
		}

		// get interval index in [0,vector size-1)
		for (auto& interval : _intervals)
		{
			DivisionWithReminder get_index;
			interval.first = get_index.getReminder(interval.first, total_frame);
			interval.second = get_index.getReminder(interval.second, total_frame);
		}
	}
}
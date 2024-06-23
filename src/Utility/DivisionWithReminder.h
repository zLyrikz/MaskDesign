#pragma once

class DivisionWithReminder
{
public:
	DivisionWithReminder();

	// a = q*b + r; 0 <= r < |b|
	// here b > 0
	// usually, b taken as the vector size, a is the periodic index, r is the standard index of the vector
	// return r
	int getReminder(int _a, unsigned int _b);
};


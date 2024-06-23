#include "DivisionWithReminder.h"

DivisionWithReminder::DivisionWithReminder()
{
}

int DivisionWithReminder::getReminder(int _a, unsigned int _b)
{
    int r = _a;

    if (_a < 0)
    {
        while (r < 0)
        {
            r += _b;
        }
    }
    else
    {
        r = _a % _b;
    }

    return r;
}

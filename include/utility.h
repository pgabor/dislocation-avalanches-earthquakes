/*
Copyright (C) 2018-2021 Gábor Péterffy <peterffy95@gmail.com>

This file is part of sdddstEQ.

sdddstEQ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

sdddstEQ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sdddstEQ. If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef UTILITY_H
#define UTILITY_H

#include <time.h>
#include <sys/time.h>

inline void normalize(double &n)
{
    while (n < -0.5)
    {
        n += 1.0;
    }

    while (n >= 0.5)
    {
        n -= 1.0;
    }
}

inline double X(const double & x)
{
    return sin(2.0 * M_PI * x) * 0.5 / M_PI;
}

inline double X2(const double & x)
{
    return (1.0 - cos(2.0 * M_PI * x)) * 0.5 / M_PI / M_PI;
}

inline double E(const double & x, const double & y, const double & K)
{
    return exp(-K*(X2(x)+X2(y)));
}

inline double X_dx(const double & x)
{
    return cos(2.0*M_PI*x);
}

inline double X2_dx(const double & x)
{
    return sin(2.0*M_PI*x)/M_PI;
}

inline double E_dx(const double & x, const double & y, const double & K)
{
    return -E(x, y, K) * K * X2_dx(x);
}

inline double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

inline double get_cpu_time()
{
    return (double)clock() / CLOCKS_PER_SEC;
}


#endif // UTILITY_H

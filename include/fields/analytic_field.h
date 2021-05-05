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

#ifndef ANALYTIC_FIELD_H
#define ANALYTIC_FIELD_H

#include "fields/field.h"
#include "constants.h"

#include <cmath>


namespace sdddst {

class AnalyticField: public Field
{
public:
    AnalyticField();
    virtual ~AnalyticField();

    virtual double xy(double dx, double dy);
    virtual double xy_diff_x(double dx, double dy);
};

}

#endif

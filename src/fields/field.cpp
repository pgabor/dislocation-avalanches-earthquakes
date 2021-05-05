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

#include "fields/field.h"

using namespace sdddst;

Field::Field()
{
    // Nothing to do
}

Field::~Field()
{
    // Nothing to do
}

double Field::xy(double, double)
{
    return 0;
}

double Field::xy_diff_x(double, double)
{
    return 0;
}

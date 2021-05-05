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

#ifndef CONSTANTS
#define CONSTANTS

#define SCALE_FACTOR_AALTO 200.0 // 200 b sized system
#define EPS 1e-12
#define ANALYTIC_FIELD_N 4
#define DEFAULT_CUTOFF_MULTIPLIER 1.0
#define DEFAULT_CUTOFF 1.0
#define DEFAULT_PRECISION 1e-8
#define DEFAULT_ITERATION_COUNT 2
#define DEFAULT_TIME_LIMIT 0.0
#define DEFAULT_STEP_SIZE 1.0
#define DEFAULT_SIM_TIME 0.0
#define DEFAULT_KASQR 1.65*1.65*1e6 / 256.0
#define DEFAULT_A 1e-4 * 16.0
#define DEFAULT_EXTERNAL_FIELD 0.0

#endif // CONSTANTS

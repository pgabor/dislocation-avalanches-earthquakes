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

#include "fields/analytic_field.h"

using namespace sdddst;

double f(const double & dx, const double & cos2piy, const double & dy)
{
    double cosh2pix = 0;
    if (dx * dx + dy * dy  > 1e-6)
    {
        cosh2pix = cosh(dx * 2.0 * M_PI);
        return  dx * (cosh2pix * cos2piy - 1.0) / ((cosh2pix - cos2piy) * (cosh2pix - cos2piy)) * 2. * M_PI * M_PI;
    }

    double dx2pi = dx * 2.0 * M_PI;
    double dy2pi = dy * 2.0 * M_PI;
    double cosh2pixminus1 = pow(dx2pi, 6.0)/720.0 + pow(dx2pi, 4.0)/24.0 + pow(dx2pi, 2) * 0.5;
    double cos2piyminus1 = -pow(dy2pi, 6.0) / 720.0 + pow(dy2pi, 4) / 24.0 - pow(dy2pi, 2) * 0.5;
    double coshminuscos = (pow(dy2pi, 6.0) + pow(dx2pi, 6.0)) / 720. + (pow(dx2pi, 4.0) - pow(dy2pi, 4.0)) / 24.0 + (pow(dx2pi, 2.0) + pow(dy2pi, 2.0)) * 0.5;
    return ((cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1)*dx)/(pow(coshminuscos, 2.0)) * 2.0 * pow(M_PI, 2.0);
}

template<int images>
double g(const double & dx, const double & cos2piy, const double & dy)
{
    return g<images-1>(dx,cos2piy, dy)+f(dx-double(images), cos2piy, dy)+f(dx+double(images), cos2piy, dy);
}

template<>
double g<0>(const double & dx, const double & cos2piy, const double & dy)
{
    return f(dx, cos2piy, dy);
}

double f_dx(const double & dx, const double & cos2piy, const double & dy)
{
    double cosh2pix = cosh(M_PI * 2.0 * dx);
    double sinh2pix = sinh(M_PI * 2.0 * dx);

    if (dx * dx + dy * dy  > 1e-6)
    {
        return ((cosh2pix * cos2piy - 1.0) / pow(cosh2pix-cos2piy, 2.0) +
                dx * (sinh2pix * cos2piy * 2.0 * M_PI / pow(cosh2pix-cos2piy, 2.0) -
                      (cosh2pix*cos2piy-1.0)/pow(cosh2pix-cos2piy, 3.0)*4.0*M_PI*sinh2pix)) * M_PI * M_PI * 2.0;
    }

    double dx2pi = dx * 2.0 * M_PI;
    double dy2pi = dy * 2.0 * M_PI;
    double cosh2pixminus1 = pow(dx2pi, 6.0)/720.0 + pow(dx2pi, 4.0)/24.0 + pow(dx2pi, 2) * 0.5;
    double cos2piyminus1 = -pow(dy2pi, 6.0) / 720.0 + pow(dy2pi, 4) / 24.0 - pow(dy2pi, 2) * 0.5;
    double coshminuscos = (pow(dy2pi, 6.0) + pow(dx2pi, 6.0)) / 720. + (pow(dx2pi, 4.0) - pow(dy2pi, 4.0)) / 24.0 + (pow(dx2pi, 2.0) + pow(dy2pi, 2.0)) * 0.5;
    double cosysinhx = - (pow(dx2pi, 7.0) * pow(dy2pi, 6.0)) / 3628800.0
                        -(pow(dx2pi, 5.0) * pow(dy2pi, 6.0)) / 86400.0
                        -(pow(dx2pi, 3.0) * pow(dy2pi, 6.0)) / 4320.0
                        -(pow(dy2pi, 6.0) * dx2pi) / 720
                        +(pow(dx2pi, 7.0) * pow(dy2pi, 4.0)) / 120960.0
                        +(pow(dx2pi, 5.0) * pow(dy2pi, 4.0)) / 2880.0
                        +(pow(dx2pi, 3.0) * pow(dy2pi, 4.0)) / 144.0
                        +(pow(dy2pi, 4.0) * dx2pi) / 24.0
                        -(pow(dx2pi, 7.0) * pow(dy2pi, 2.0)) / 1080.0
                        -(pow(dx2pi, 5.0) * pow(dy2pi, 2.0)) / 240.0
                        -(pow(dx2pi, 3.0) * pow(dy2pi, 2.0)) / 12.0
                        -(pow(dy2pi, 2.0) * dx2pi) * 0.5
                        +pow(dx2pi, 7.0) / 5040.0
                        +pow(dx2pi, 5.0) / 120.0
                        +pow(dx2pi, 3.0) / 6.0
                        +dx2pi;
    return (
                (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1)/pow(coshminuscos, 2.0) +
                ((cosysinhx*dx)/(pow(coshminuscos, 2.0)) -
                (cosh2pixminus1 * cos2piyminus1 + cos2piyminus1 + cosh2pixminus1)/(pow(coshminuscos, 3.0))*dx*sinh(dx2pi) * 2.0)* M_PI * 2.0
           ) * 2.0 * pow(M_PI, 2.0);


}

template<int images>
double g_dx(const double & dx, const double & cos2piy, const double & dy)
{
    return g_dx<images-1>(dx, cos2piy, dy) + f_dx(dx-double(images), cos2piy, dy) + f_dx(dx+double(images), cos2piy, dy);
}

template<>
double g_dx<0>(const double & dx, const double & cos2piy, const double & dy)
{
    return f_dx(dx, cos2piy, dy);
}


AnalyticField::AnalyticField() :
    Field()
{
    // Nothing to do
}

AnalyticField::~AnalyticField()
{
    // Nothing to do
}

double AnalyticField::xy(double dx, double dy)
{
    double cos2piy = cos(M_PI * 2.0 * dy);
    if (dx < 0)
    {
        return - dx * f(dx + double(ANALYTIC_FIELD_N+1), cos2piy, dy) +
                (1.0 + dx) * f(dx - double(ANALYTIC_FIELD_N), cos2piy, dy) +
                f(dx + double(ANALYTIC_FIELD_N), cos2piy, dy) +
                g<ANALYTIC_FIELD_N-1>(dx, cos2piy, dy);
    }

    return dx * f(dx - (double(ANALYTIC_FIELD_N+1)), cos2piy, dy) +
            (1.0 - dx) * f(dx + double(ANALYTIC_FIELD_N), cos2piy, dy) +
            f(dx - double(ANALYTIC_FIELD_N), cos2piy, dy) +
            g<ANALYTIC_FIELD_N-1>(dx, cos2piy, dy);
}

double AnalyticField::xy_diff_x(double dx, double dy)
{
    double cos2piy = cos(M_PI * 2.0 * dy);
    if (dx < 0)
    {
        return -f(dx+double(ANALYTIC_FIELD_N)+1.0, cos2piy, dy)-
                dx*f_dx(dx+double(ANALYTIC_FIELD_N)+1.0, cos2piy, dy)+
                f_dx(dx-double(ANALYTIC_FIELD_N), cos2piy, dy)*(1.0+dx) +
                f(dx-double(ANALYTIC_FIELD_N), cos2piy, dy) +
                f_dx(dx + double(ANALYTIC_FIELD_N), cos2piy, dy) +
                g_dx<ANALYTIC_FIELD_N-1>(dx, cos2piy, dy);
    }

    return  f(dx-double(ANALYTIC_FIELD_N)-1.0, cos2piy, dy)+
            dx*f_dx(dx-double(ANALYTIC_FIELD_N)-1.0, cos2piy, dy)+
            f_dx(dx+double(ANALYTIC_FIELD_N), cos2piy, dy)*(1.0-dx) -
            f(dx+double(ANALYTIC_FIELD_N), cos2piy, dy) +
            f_dx(dx - double(ANALYTIC_FIELD_N), cos2piy, dy) +
            g_dx<ANALYTIC_FIELD_N-1>(dx, cos2piy, dy);
}

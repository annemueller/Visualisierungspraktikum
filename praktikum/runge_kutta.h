//
// Created by yves on 19.01.17.
//

#ifndef VISPRAKTIKUM_RUNGE_KUTTA_H
#define VISPRAKTIKUM_RUNGE_KUTTA_H

#include <fantom/fields.hpp>
using namespace fantom;

std::vector<Point3>
runge_kutta(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, size_t steps, double h,
            Point3 start, bool vor, Point3 normale, std::shared_ptr<const Grid<3 >> grid, Cell c);

#endif //VISPRAKTIKUM_RUNGE_KUTTA_H

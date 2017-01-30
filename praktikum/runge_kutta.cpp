//
// Created by yves on 19.01.17.
//

#include "runge_kutta.h"

Point3 schnittpunkt(Point3 q, Point3 normale, Point3 start){

    double n1 = normale[0];
    double n2 = normale[1];
    double n3 = normale[2];

    double q1 = q[0];
    double q2 = q[1];
    double q3 = q[2];

    double b = n1*start[0] + n2*start[1] + n3*start[2];

    double lambda = (b-(n1*q1)-(n2*q2)-(n3*q3)) / (pow(n1,2) + pow(n2, 2) +  pow(n3, 2));

    Point3 x;
    x[0] = q1 + (lambda * n1);
    x[1] = q2 + (lambda * n2);
    x[2] = q3 + (lambda * n3);

    return x;
}

std::vector<Point3>
runge_kutta(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, size_t steps, double h,
            Point3 start, bool vor, Point3 normale, std::shared_ptr<const Grid<3 >> grid, Cell c) {
    auto evaluator = field->makeEvaluator();

    std::vector<Point3> points;
    //return points;
    Point3 point = {start[0], start[1], start[2]};
    points.push_back(point);

    for (size_t i = 0; i < steps; i++) {
        if (evaluator->reset(point)) {
            // Berechne q1-q4
            auto v = evaluator->value();
            Point3 q1 = h * v;
            q1 = point + 0.5 * q1;
            if (vor == false){
                q1 = -q1;
            }
            q1 = schnittpunkt(q1, normale, start);
            if (evaluator->reset(q1)) {
                v = evaluator->value();
                Point3 q2 = h * v;
                q2 = point + 0.5 * q2;
                if (vor ==false){
                    q2 = -q2;
                }
                q2 = schnittpunkt(q2, normale, start);
                if (evaluator->reset(q2)) {
                    v = evaluator->value();
                    Point3 q3 = h * v;
                    q3 = point + 0.5 * q3;
                    if (vor ==false){
                        q3 = -q3;
                    }
                    q3 = schnittpunkt(q3, normale, start);
                    if (evaluator->reset(q3)) {
                        v = evaluator->value();
                        Point3 q4 = h * v;
                        q4 = point + 0.5 * q4;
                        if (vor ==false){
                            q4 = -q4;
                        }
                        q4 = schnittpunkt(q4, normale, start);
                        point = (1.0 / 6.0) * (q1 + 2 * q2 + 2 * q3 + q4);

                        if (grid->contains(c, point)) {
                            points.push_back(point);
                        } else{
                            break;
                        }
                    } else {
                        break; // nicht mehr in Domain
                    }

                } else {
                    break; // nicht mehr in Domain
                }
            } else {
                break; // nicht mehr in Domain
            }
        } else {
            break; // nicht mehr in Domain
        }
    }
    return points;
}
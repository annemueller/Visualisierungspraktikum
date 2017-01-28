//
// Created by yves on 19.01.17.
//

#include "runge_kutta.h"

std::vector<Point3>
runge_kutta(std::shared_ptr<const TensorFieldContinuous<3, Vector3> > &field, size_t steps, double h,
            Point3 start, bool vor) {
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
            if (vor == false){
                q1 = -q1;
            }
            if (evaluator->reset(point + 0.5 * q1)) {
                v = evaluator->value();
                Point3 q2 = h * v;
                if (vor ==false){
                    q2 = -q2;
                }
                if (evaluator->reset(point + 0.5 * q2)) {
                    v = evaluator->value();
                    Point3 q3 = h * v;
                    if (vor ==false){
                        q3 = -q3;
                    }
                    if (evaluator->reset(point + 0.5 * q3)) {
                        v = evaluator->value();
                        Point3 q4 = h * v;
                        if (vor ==false){
                            q4 = -q4;
                        }
                        // TODO: q4 vorher auf Dreieck/Oberflaeche druecken!
                        point = point + (1.0 / 6.0) * (q1 + 2 * q2 + 2 * q3 + q4);
                        // TODO: Punkt noch im Dreieck? sonst abbrechen und punkt nicht in Liste speichern!
                        points.push_back(point);
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
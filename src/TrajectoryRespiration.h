//
// Created by Hanna Maria Hanson on 07/08/2020.
//

#ifndef JEMRIS_TRAJECTORYRESPIRATION_H
#define JEMRIS_TRAJECTORYRESPIRATION_H

#include "TrajectoryInterface.h"
#include "NDData.h"

class TrajectoryRespiration : public TrajectoryInterface {
public:
    TrajectoryRespiration();
    virtual ~TrajectoryRespiration();

    virtual void GetValueDerived(double time, double *value);

    virtual void LoadFile(string filename);

    double dx;
    size_t xi;
    double dy;
    size_t yi;
    double dz;
    size_t zi;
    size_t xip;
    size_t yip;
    size_t zip;

    double int_part;
    double store_x;
    double store_y;
    double store_z;
    double m_ap_interpolated[3];
    double m_si_interpolated[3];
    double m_model_offset_interpolated[3];

protected:
    vector<double> m_res;
    vector<double> m_offset;
    NDData<double> m_field;
    NDData<double> m_ap;
    NDData<double> m_si;
    NDData<double> m_model_offset;
    NDData<double> m_breathing_trace;
};

#endif //JEMRIS_TRAJECTORYRESPIRATION_H

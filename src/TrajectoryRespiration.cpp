//
// Created by Hanna Maria Hanson on 07/08/2020.
//

#include "TrajectoryRespiration.h"

#include "BinaryContext.h"
#include <cmath>

TrajectoryRespiration:: TrajectoryRespiration() {
}

TrajectoryRespiration::~TrajectoryRespiration() {
}

void TrajectoryRespiration::LoadFile(string filename) {

    // file is in HDF5 format
    BinaryContext bc(filename, IO::IN);
    if (bc.Status() != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }

    if (bc.Read(m_field, "resolution", "/model") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }
    m_res = m_field.Data();

    if (bc.Read(m_field, "offset", "/model") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }
    m_offset = m_field.Data();

    if (bc.Read(m_ap, "ap", "/model") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }



    if (bc.Read(m_si, "si", "/model") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }


    if (bc.Read(m_model_offset, "model_offset", "/model") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }

    if (bc.Read(m_breathing_trace, "breathing_trace", "/model") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }


    // Samples automatically have their positions offset by half their size to
    // centre them around (0, 0, 0), so we do the same to the field
    for (size_t i = 0; i < 3; i++) {
        m_offset[i] -= 0.5 * (m_ap.Dim(i + 1) - 1) * m_res[i + 1];
    }

    // this is necessary to fool the TrajectoryInterface into calling
    // GetValueDerived() as it requires at least one time in the list
    m_time.push_back(100);

}



void TrajectoryRespiration::GetValueDerived(double time, double *value) {


    // time is cyclical so we calculate mod total duration
    // but space is clamped at the edges
    double int_part;
    double dt = modf(time / m_res[0], &int_part);
    size_t ti = int(int_part) % m_breathing_trace.Dim(1);
    double dx = modf((value[0] - m_offset[0]) / m_res[1], &int_part);
    size_t xi = min(max(size_t(int_part), size_t(0)), m_ap.Dim(1));
    double dy = modf((value[1] - m_offset[1]) / m_res[2], &int_part);
    size_t yi = min(max(size_t(int_part), size_t(0)), m_ap.Dim(2));
    double dz = modf((value[2] - m_offset[2]) / m_res[3], &int_part);
    size_t zi = min(max(size_t(int_part), size_t(0)), m_ap.Dim(3));



    // also calculate indices plus one for the interpolation
    size_t xip = min(xi + 1, m_ap.Dim(1) - 1);
    size_t yip = min(yi + 1, m_ap.Dim(2) - 1);
    size_t zip = min(zi + 1, m_ap.Dim(3) - 1);
    size_t tip = (ti + 1) % m_breathing_trace.Dim(1);


    NDData<double> DVF(3,2,2,2,2);

    // calculate DVF for indices defined above
    for (size_t i=0; i < 3; i++) {
        double ap = m_breathing_trace(0,ti);
        double si = m_breathing_trace(1,ti);


        DVF(i,0,0, 0, 0) = ap * m_ap(i,xi,yi,zi) + si * m_si(i,xi,yi,zi) + m_model_offset(i,xi,yi,zi);
        DVF(i,1,0, 0, 0) = ap * m_ap(i,xip,yi,zi) + si * m_si(i,xip,yi,zi) + m_model_offset(i,xip,yi,zi);
        DVF(i,0,1, 0, 0) = ap * m_ap(i,xi,yip,zi) + si * m_si(i,xi,yip,zi) + m_model_offset(i,xi,yip,zi);
        DVF(i,1,1, 0, 0) = ap * m_ap(i,xip,yip,zi) + si * m_si(i,xip,yip,zi) + m_model_offset(i,xip,yip,zi);



        DVF(i,0,0,1,0) = ap * m_ap(i,xi,yi,zip) + si * m_si(i,xi,yi,zip) + m_model_offset(i,xi,yi,zip);
        DVF(i,1,0,1,0) = ap * m_ap(i,xip,yi,zip) + si * m_si(i,xip,yi,zip) + m_model_offset(i,xip,yi,zip);
        DVF(i,0,1,1,0) = ap * m_ap(i,xi,yip,zip) + si * m_si(i,xi,yip,zip) + m_model_offset(i,xi,yip,zip);
        DVF(i,1,1,1,0) = ap * m_ap(i,xip,yip,zip) + si * m_si(i,xip,yip,zip) + m_model_offset(i,xip,yip,zip);

        ap = m_breathing_trace(0,tip);
        si = m_breathing_trace(1,tip);


        DVF(i,0,0, 0, 1) = ap * m_ap(i,xi,yi,zi) + si * m_si(i,xi,yi,zi) + m_model_offset(i,xi,yi,zi);
        DVF(i,1,0, 0, 1) = ap * m_ap(i,xip,yi,zi) + si * m_si(i,xip,yi,zi) + m_model_offset(i,xip,yi,zi);
        DVF(i,0,1, 0, 1) = ap * m_ap(i,xi,yip,zi) + si * m_si(i,xi,yip,zi) + m_model_offset(i,xi,yip,zi);
        DVF(i,1,1, 0, 1) = ap * m_ap(i,xip,yip,zi) + si * m_si(i,xip,yip,zi) + m_model_offset(i,xip,yip,zi);



        DVF(i,0,0,1,1) = ap * m_ap(i,xi,yi,zip) + si * m_si(i,xi,yi,zip) + m_model_offset(i,xi,yi,zip);
        DVF(i,1,0,1,1) = ap * m_ap(i,xip,yi,zip) + si * m_si(i,xip,yi,zip) + m_model_offset(i,xip,yi,zip);
        DVF(i,0,1,1,1) = ap * m_ap(i,xi,yip,zip) + si * m_si(i,xi,yip,zip) + m_model_offset(i,xi,yip,zip);
        DVF(i,1,1,1,1) = ap * m_ap(i,xip,yip,zip) + si * m_si(i,xip,yip,zip) + m_model_offset(i,xip,yip,zip);


    }

    for (size_t i=0; i < 3; i++) {
        value[i] += (1 - dt) * (1 - dx) * (1 - dy) * (1 - dz) * DVF(i, 0, 0, 0, 0) +
                    (dt) * (1 - dx) * (1 - dy) * (1 - dz) * DVF(i, 0, 0, 0, 1) +
                    (1 - dt) * (dx) * (1 - dy) * (1 - dz) * DVF(i, 1, 0, 0, 0) +
                    (1 - dt) * (1 - dx) * (dy) * (1 - dz) * DVF(i, 0, 1, 0, 0) +
                    (1 - dt) * (1 - dx) * (1 - dy) * (dz) * DVF(i, 0, 0, 1, 0) +
                    (dt) * (dx) * (1 - dy) * (1 - dz) * DVF(i, 1, 0, 0, 1) +
                    (dt) * (1 - dx) * (dy) * (1 - dz) * DVF(i, 0, 1, 0, 1) +
                    (dt) * (1 - dx) * (1 - dy) * (dz) * DVF(i, 0, 0, 1, 1) +
                    (1 - dt) * (dx) * (dy) * (1 - dz) * DVF(i, 1, 1, 0, 0) +
                    (1 - dt) * (dx) * (1 - dy) * (dz) * DVF(i, 1, 0, 1, 0) +
                    (1 - dt) * (1 - dx) * (dy) * (dz) * DVF(i, 0, 1, 1, 0) +
                    (dt) * (dx) * (dy) * (1 - dz) * DVF(i, 1, 1, 0, 1) +
                    (dt) * (dx) * (1 - dy) * (dz) * DVF(i, 1, 0, 1, 1) +
                    (dt) * (1 - dx) * (dy) * (dz) * DVF(i, 0, 1, 1, 1) +
                    (1 - dt) * (dx) * (dy) * (dz) * DVF(i, 1, 1, 1, 0) +
                    (dt) * (dx) * (dy) * (dz) * DVF(i, 1, 1, 1, 1);
    }


}
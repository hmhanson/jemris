//
// Created by bennyrowland on 3/20/18.
//

#include "TrajectoryDeformation.h"

#include "BinaryContext.h"
#include <cmath>

TrajectoryDeformation::TrajectoryDeformation() {

}

TrajectoryDeformation::~TrajectoryDeformation() {
}

void TrajectoryDeformation::LoadFile(string filename) {

    // file is in HDF5 format with 3 fields:
    // resolution and offset are 3-tuples of doubles and
    // field is a 5d shape (component, x, y, z, t) storing
    // the deformation field as a grid
    BinaryContext bc(filename, IO::IN);
    if (bc.Status() != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }

    if (bc.Read(m_field, "resolution", "/deformation") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }
    m_res = m_field.Data();

    if (bc.Read(m_field, "offset", "/deformation") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }
    m_offset = m_field.Data();

    if (bc.Read(m_field, "field", "/deformation") != IO::OK) {
        cout << bc.Status();
        exit(-1);
    }

    // Samples automatically have their positions offset by half their size to
    // centre them around (0, 0, 0), so we do the same to the field
    for (size_t i = 0; i < 3; i++) {
        m_offset[i] -= 0.5 * (m_field.Dim(i + 1) - 1) * m_res[i + 1];
    }

    // this is necessary to fool the TrajectoryInterface into calling
    // GetValueDerived() as it requires at least one time in the list
    m_time.push_back(100);
}

void TrajectoryDeformation::GetValueDerived(double time, double *value) {

    // time is cyclical so we calculate mod total duration
    // but space is clamped at the edges
    double int_part;
    double dt = modf(time / m_res[0], &int_part);
    size_t ti = int(int_part) % m_field.Dim(4);
    double dx = modf((value[0] - m_offset[0]) / m_res[1], &int_part);
    size_t xi = min(max(size_t(int_part), size_t(0)), m_field.Dim(1));
    double dy = modf((value[1] - m_offset[1]) / m_res[2], &int_part);
    size_t yi = min(max(size_t(int_part), size_t(0)), m_field.Dim(2));
    double dz = modf((value[2] - m_offset[2]) / m_res[3], &int_part);
    size_t zi = min(max(size_t(int_part), size_t(0)), m_field.Dim(3));

    // also calculate indices plus one for the interpolation
    size_t xip = min(xi + 1, m_field.Dim(1) - 1);
    size_t yip = min(yi + 1, m_field.Dim(2) - 1);
    size_t zip = min(zi + 1, m_field.Dim(3) - 1);
    size_t tip = (ti + 1) % m_field.Dim(4);

    for (size_t i=0; i < 3; i++) {
        value[i] += (1 - dt) * (1 - dx) * (1 - dy) * (1 - dz) * m_field(i, xi, yi, zi, ti) +
                    (dt) * (1 - dx) * (1 - dy) * (1 - dz) * m_field(i, xi, yi, zi, tip) +
                    (1 - dt) * (dx) * (1 - dy) * (1 - dz) * m_field(i, xip, yi, zi, ti) +
                    (1 - dt) * (1 - dx) * (dy) * (1 - dz) * m_field(i, xi, yip, zi, ti) +
                    (1 - dt) * (1 - dx) * (1 - dy) * (dz) * m_field(i, xi, yi, zip, ti) +
                    (dt) * (dx) * (1 - dy) * (1 - dz) * m_field(i, xip, yi, zi, tip) +
                    (dt) * (1 - dx) * (dy) * (1 - dz) * m_field(i, xi, yip, zi, tip) +
                    (dt) * (1 - dx) * (1 - dy) * (dz) * m_field(i, xi, yi, zip, tip) +
                    (1 - dt) * (dx) * (dy) * (1 - dz) * m_field(i, xip, yip, zi, ti) +
                    (1 - dt) * (dx) * (1 - dy) * (dz) * m_field(i, xip, yi, zip, ti) +
                    (1 - dt) * (1 - dx) * (dy) * (dz) * m_field(i, xi, yip, zip, ti) +
                    (dt) * (dx) * (dy) * (1 - dz) * m_field(i, xip, yip, zi, tip) +
                    (dt) * (dx) * (1 - dy) * (dz) * m_field(i, xip, yi, zip, tip) +
                    (dt) * (1 - dy) * (dy) * (dz) * m_field(i, xi, yip, zip, tip) +
                    (1 - dt) * (dx) * (dy) * (dz) * m_field(i, xip, yip, zip, ti) +
                    (dt) * (dx) * (dy) * (dz) * m_field(i, xip, yip, zip, tip);
    }
}

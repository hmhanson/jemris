//
// Created by bennyrowland on 3/20/18.
//

#ifndef TRAJECTORYDEFORMATION_H
#define TRAJECTORYDEFORMATION_H

#include "TrajectoryInterface.h"
#include "NDData.h"

/**
 * @brief Deformation Trajectory
 */
class TrajectoryDeformation : public TrajectoryInterface {
public:
    TrajectoryDeformation();
    virtual ~TrajectoryDeformation();

    virtual void GetValueDerived(double time, double *value);

    virtual void LoadFile(string filename);

protected:
    vector<double> m_res;
    vector<double> m_offset;
    NDData<double> m_field;
};


#endif //TRAJECTORYDEFORMATION_H

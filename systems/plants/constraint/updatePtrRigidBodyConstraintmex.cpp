#include "mex.h"
#include <iostream>
#include "drakeUtil.h"
#include <Eigen/Dense>
#include "RigidBodyConstraint.h"
#include "RigidBodyManipulator.h"
#include "constructPtrRigidBodyConstraint.h"
#include <cstdio>

using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  RigidBodyConstraint* constraint = (RigidBodyConstraint*) getDrakeMexPointer(prhs[0]);
  int constraint_type = constraint->getType();
  mwSize strlen = mxGetNumberOfElements(prhs[1])+1;
  char* field = new char[strlen];
  mxGetString(prhs[1],field,strlen);
  string field_str(field);
  switch(constraint_type)
  {
    case RigidBodyConstraint::QuasiStaticConstraintType:
      {
        QuasiStaticConstraint* cnst = (QuasiStaticConstraint*) constraint;
        if(field_str=="active")
        {
          // setActive(qsc_ptr,flag)
          if(mxGetNumberOfElements(prhs[2]) != 1)
          {
            mexErrMsgIdAndTxt("Drake:updatePtrRigidBodyConstraintmex:BadInputs","QuasiStaticConstraint:flag must be a single boolean");
          }
          bool* flag = mxGetLogicals(prhs[2]);
          QuasiStaticConstraint* cnst_new = new QuasiStaticConstraint(*cnst);
          cnst_new->setActive(*flag);
          plhs[0] = createDrakeConstraintMexPointer((void*) cnst_new,"deleteRigidBodyConstraintmex","QuasiStaticConstraint");
        }
        else if(field_str=="factor")
        {// setShrinkFactor(qsc_ptr,factor)
          if(mxGetNumberOfElements(prhs[2]) != 1)
          {
            mexErrMsgIdAndTxt("Drake:updatePtrRigidBodyConstraintmex:BadInputs","QuasiStaticConstraint:shrink factor must be a double scalar");
          }
          double factor = mxGetScalar(prhs[2]);
          if(factor<=0.0)
          {
            mexErrMsgIdAndTxt("Drake:updatePtrRigidBodyConstraintmex:BadInputs","QuasiStaticConstraint:shrink factor should be a positive scalar");
          }
          QuasiStaticConstraint* cnst_new = new QuasiStaticConstraint(*cnst);
          cnst_new->setShrinkFactor(factor);
          plhs[0] = createDrakeConstraintMexPointer((void*) cnst_new,"deleteRigidBodyConstraintmex","QuasiStaticConstraint");
        }
        else if(field_str=="contact")
        {
          // addContact(qsc_ptr,body1, body1_pts, body2, body2_pts,...)
          int num_new_bodies = (nrhs-2)/2;
          int* new_bodies = new int[num_new_bodies];
          MatrixXd* new_body_pts = new MatrixXd[num_new_bodies];
          for(int idx = 0;idx<num_new_bodies;idx++)
          {
            new_bodies[idx] = (int) mxGetScalar(prhs[2+idx*2])-1;
            int npts = mxGetN(prhs[3+idx*2]);
            MatrixXd new_body_pts_tmp(3,npts);
            memcpy(new_body_pts_tmp.data(),mxGetPr(prhs[3+idx*2]),sizeof(double)*3*npts);
            new_body_pts[idx].resize(4,npts);
            new_body_pts[idx].block(0,0,3,npts) = new_body_pts_tmp;
            new_body_pts[idx].row(3) = MatrixXd::Ones(1,npts);
          }
          QuasiStaticConstraint* cnst_new = new QuasiStaticConstraint(*cnst);
          cnst_new->addContact(num_new_bodies,new_bodies,new_body_pts);
          plhs[0] = createDrakeConstraintMexPointer((void*) cnst_new,"deleteRigidBodyConstraintmex","QuasiStaticConstraint");
          delete[] new_bodies;
          delete[] new_body_pts;
        }
        else if(field_str=="robot")
        {
          RigidBodyManipulator* robot = (RigidBodyManipulator*) getDrakeMexPointer(prhs[2]);
          QuasiStaticConstraint* cnst_new = new QuasiStaticConstraint(*cnst);
          cnst_new->updateRobot(robot);
          plhs[0] = createDrakeConstraintMexPointer((void*) cnst_new,"deleteRigidBodyConstraintmex","QuasiStaticConstraint");
        }
        else if(field_str=="robotnum")
        {
          int num_robot = mxGetNumberOfElements(prhs[2]);
          double* robotnum_tmp = new double[num_robot];
          int* robotnum = new int[num_robot];
          memcpy(robotnum_tmp,mxGetPr(prhs[2]),sizeof(double)*num_robot);
          for(int i = 0;i<num_robot;i++)
          {
            robotnum[i] = (int) robotnum_tmp[i]-1;
          }
          set<int> robotnumset(robotnum,robotnum+num_robot);
          QuasiStaticConstraint* cnst_new = new QuasiStaticConstraint(*cnst);
          cnst_new->updateRobotnum(robotnumset);
          plhs[0] = createDrakeConstraintMexPointer((void*) cnst_new,"deleteRigidBodyConstraintmex","QuasiStaticConstraint");
          delete[] robotnum_tmp;
          delete[] robotnum;
        }
        else
        {
          mexErrMsgIdAndTxt("Drake:updatePtrRigidBodyConstraintmex:BadInputs","QuasiStaticConstraint:argument 2 is not accepted");
        }
      }
      break;
    case RigidBodyConstraint::PostureConstraintType:
      {
        PostureConstraint* pc = (PostureConstraint*) constraint;
        if(field_str=="bounds")
        { // setJointLimits(pc,joint_idx,lb,ub)
          int num_idx = mxGetM(prhs[1]);
          if(!mxIsNumeric(prhs[1]) || mxGetN(prhs[1]) != 1 || !mxIsNumeric(prhs[2]) || mxGetM(prhs[2]) != num_idx || mxGetN(prhs[2]) != 1 || !mxIsNumeric(prhs[3]) || mxGetM(prhs[3]) != num_idx || mxGetN(prhs[3]) != 1)
          {
            mexErrMsgIdAndTxt("Drake:constructPtrRIgidBodyConstraint:BadInputs","PostureConstraint:joint_idx, lb and ub must be of the same length numerical vector");
          }
          double* joint_idx_tmp = new double[num_idx];
          int* joint_idx = new int[num_idx];
          memcpy(joint_idx_tmp,mxGetPr(prhs[1]),sizeof(double)*num_idx);
          for(int i = 0;i<num_idx;i++)
          {
            joint_idx[i] = (int) joint_idx_tmp[i]-1;
          }
          double* lb = new double[num_idx];
          double* ub = new double[num_idx];
          memcpy(lb,mxGetPr(prhs[2]),sizeof(double)*num_idx);
          memcpy(ub,mxGetPr(prhs[3]),sizeof(double)*num_idx);
          PostureConstraint* pc_new = new PostureConstraint(*pc);
          pc_new->setJointLimits(num_idx,joint_idx,lb,ub);
          delete[] joint_idx_tmp; delete[] joint_idx; delete[] lb; delete[] ub;
          plhs[0] = createDrakeConstraintMexPointer((void*)pc_new,"deleteRigidBodyConstraintmex","PostureConstraint");
        }
      }
      break;
    default:
      mexErrMsgIdAndTxt("Drake:updatePtrRigidBodyConstraintmex:BadInputs","Unsupported constraint type");
      break;
  }
}

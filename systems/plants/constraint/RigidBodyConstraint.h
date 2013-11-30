#ifndef __CONSTRAINT_H__
#define __CONSTRAINT_H__
#include <iostream>
#include "../RigidBodyManipulator.h"
#include <vector>
#include <set>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../../../util/drakeQuatUtil.h"
#include <sstream>


namespace DrakeRigidBodyConstraint{
  extern Eigen::Vector4d com_pts;
  extern const int WorldCoMDefaultRobotNum[1];
  extern Eigen::Vector2d default_tspan;
}

void drakePrintMatrix(const Eigen::MatrixXd &mat);


class RigidBodyConstraint
{
  protected:
    int category;
    int type;
  public:
    RigidBodyConstraint(int category):category(category),type(0){};
    RigidBodyConstraint(const RigidBodyConstraint &rhs);
    int getType() {return this->type;};
    int getCategory() {return this->category;};
    virtual ~RigidBodyConstraint(void) = 0;
    /* In each category, constraint class share the same function interface, this value needs to be in consistent with that in MATLAB*/
    static const int SingleTimeKinematicConstraintCategory       = -1;
    static const int MultipleTimeKinematicConstraintCategory     = -2;
    static const int QuasiStaticConstraintCategory               = -3;
    static const int PostureConstraintCategory                   = -4;
    static const int MultipleTimeLinearPostureConstraintCategory = -5;
    static const int SingleTimeLinearPostureConstraintCategory   = -6;
    /* Each non-abstrac RigidBodyConstraint class has a unique type. Please stay in consistent with the value in MATLAB*/
    static const int QuasiStaticConstraintType                  = 1;
    static const int PostureConstraintType                      = 2;
    static const int SingleTimeLinearPostureConstraintType      = 3;
    static const int AllBodiesClosestDistanceConstraintType     = 4;
    static const int WorldEulerConstraintType                   = 5;
    static const int WorldGazeDirConstraintType                 = 6;
    static const int WorldGazeOrientConstraintType              = 7;
    static const int WorldGazeTargetConstraintType              = 8;
    static const int RelativeGazeTargetConstraintType           = 9;
    static const int WorldCoMConstraintType                     = 10;
    static const int WorldPositionConstraintType                = 11;
    static const int WorldPositionInFrameConstraintType         = 12;
    static const int WorldQuatConstraintType                    = 13;
    static const int Point2PointDistanceConstraintType          = 14;
    static const int Point2LineSegDistConstraintType            = 15;
    static const int WorldFixedPositionConstraintType           = 16;
    static const int WorldFixedOrientConstraintType             = 17;
    static const int WorldFixedBodyPoseConstraintType           = 18;
    static const int PostureChangeConstraintType                = 19;
};



class QuasiStaticConstraint: public RigidBodyConstraint
{
  protected:
    RigidBodyManipulator* robot;
    std::set<int> m_robotnumset;
    double tspan[2];
    double shrinkFactor;
    bool active;
    int num_bodies;
    int num_pts;
    std::vector<int> bodies;
    std::vector<int> num_body_pts;
    std::vector<Eigen::MatrixXd> body_pts;
    static const std::set<int> defaultRobotNumSet;
  public:
    QuasiStaticConstraint(RigidBodyManipulator* robot, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan,const std::set<int> &robotnumset = QuasiStaticConstraint::defaultRobotNumSet); 
    QuasiStaticConstraint(const QuasiStaticConstraint &rhs);
    bool isTimeValid(const double* t);
    int getNumConstraint(const double* t);
    void eval(const double* t,double* weights,Eigen::VectorXd &c, Eigen::MatrixXd &dc); 
    void bounds(const double* t,Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    void name(const double* t,std::vector<std::string> &name_str);
    virtual ~QuasiStaticConstraint(void);
    bool isActive(){return this->active;};
    RigidBodyManipulator* getRobotPointer(){return robot;};
    int getNumWeights() {return this->num_pts;};
    void addContact(int num_new_bodies, const int* body, const MatrixXd* body_pts);
    void setShrinkFactor(double factor);
    void setActive(bool flag){this->active = flag;};
    void updateRobot(RigidBodyManipulator *robot);
    void updateRobotnum(std::set<int> &robotnumset);
};

class PostureConstraint: public RigidBodyConstraint
{
  protected:
    double tspan[2];
    double* lb;
    double* ub;
    RigidBodyManipulator *robot;
  public:
    PostureConstraint(RigidBodyManipulator *model, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    PostureConstraint(const PostureConstraint& rhs);
    bool isTimeValid(const double* t);
    void setJointLimits(int num_idx, const int* joint_idx, const double* lb, const double* ub);
    void bounds(const double* t,double* joint_min, double* joint_max);
    RigidBodyManipulator* getRobotPointer(){return robot;};
    virtual ~PostureConstraint(void);
};

class MultipleTimeLinearPostureConstraint: public RigidBodyConstraint
{
  protected:
    RigidBodyManipulator* robot;
    int numValidTime(const std::vector<bool> &valid_flag);
    void validTimeInd(const std::vector<bool> &valid_flag, Eigen::VectorXi &valid_t_ind);
  public:
    double tspan[2];
    MultipleTimeLinearPostureConstraint(RigidBodyManipulator *model, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    MultipleTimeLinearPostureConstraint(const MultipleTimeLinearPostureConstraint& rhs);
    std::vector<bool> isTimeValid(const double* t, int n_breaks);
    RigidBodyManipulator* getRobotPointer(){return robot;};
    void eval(const double* t, int n_breaks, const Eigen::MatrixXd &q, Eigen::VectorXd &c,Eigen::SparseMatrix<double> &dc);
    virtual int getNumConstraint(const double* t, int n_breaks) = 0;
    virtual void feval(const double* t, int n_breaks, const Eigen::MatrixXd &q, Eigen::VectorXd &c) = 0;
    virtual void geval(const double* t, int n_breaks, Eigen::VectorXi &iAfun, Eigen::VectorXi &jAvar, Eigen::VectorXd &A) = 0;
    virtual void name(const double* t, int n_breaks, std::vector<std::string> &name_str) = 0;
    virtual void bounds(const double* t, int n_breaks, Eigen::VectorXd &lb, Eigen::VectorXd &ub) = 0;
    virtual ~MultipleTimeLinearPostureConstraint(){};
};

class SingleTimeLinearPostureConstraint: public RigidBodyConstraint
{
  protected:
    RigidBodyManipulator* robot;
    Eigen::VectorXi iAfun;
    Eigen::VectorXi jAvar;
    Eigen::VectorXd A;
    Eigen::VectorXd lb;
    Eigen::VectorXd ub;
    int num_constraint;
    Eigen::SparseMatrix<double> A_mat;
    double tspan[2];
  public:
    SingleTimeLinearPostureConstraint(RigidBodyManipulator* robot, const VectorXi &iAfun, const VectorXi &jAvar, const VectorXd &A, const VectorXd &lb, const VectorXd &ub, const Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    SingleTimeLinearPostureConstraint(const SingleTimeLinearPostureConstraint& rhs);
    bool isTimeValid(const double* t);
    int getNumConstraint(const double* t);
    RigidBodyManipulator* getRobotPointer(){return robot;};
    void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    void feval(const double* t, const Eigen::VectorXd &q, Eigen::VectorXd &c);
    void geval(const double* t, Eigen::VectorXi &iAfun, Eigen::VectorXi &jAvar, Eigen::VectorXd &A);
    void eval(const double* t, const Eigen::VectorXd &q, Eigen::VectorXd &c, Eigen::SparseMatrix<double> &dc);
    void name(const double* t, std::vector<std::string> &name_str);
};

class SingleTimeKinematicConstraint: public RigidBodyConstraint
{
  protected:
    RigidBodyManipulator *robot;
    int num_constraint;
  public:
    double tspan[2];
    SingleTimeKinematicConstraint(RigidBodyManipulator *model, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    SingleTimeKinematicConstraint(const SingleTimeKinematicConstraint &rhs);
    bool isTimeValid(const double* t);
    int getNumConstraint(const double* t);
    RigidBodyManipulator* getRobotPointer(){return robot;};
    virtual void eval(const double* t,Eigen::VectorXd &c, Eigen::MatrixXd &dc) = 0;
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub) = 0;
    virtual void name(const double* t, std::vector<std::string> &name_str) = 0;
    virtual void updateRobot(RigidBodyManipulator *robot);
    virtual ~SingleTimeKinematicConstraint(){};
};

class MultipleTimeKinematicConstraint : public RigidBodyConstraint
{
  protected:
    RigidBodyManipulator *robot;
    int numValidTime(const double* t,int n_breaks);
  public:
    double tspan[2];
    MultipleTimeKinematicConstraint(RigidBodyManipulator *model, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    MultipleTimeKinematicConstraint(const MultipleTimeKinematicConstraint &rhs);
    std::vector<bool> isTimeValid(const double* t,int n_breaks);
    RigidBodyManipulator* getRobotPointer(){return robot;};
    virtual int getNumConstraint(const double* t,int n_breaks) = 0;
    void eval(const double* t, int n_breaks,const Eigen::MatrixXd &q,Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void eval_valid(const double* valid_t, int num_valid_t,const Eigen::MatrixXd &valid_q,Eigen::VectorXd &c, Eigen::MatrixXd &dc_valid) = 0;
    virtual void bounds(const double* t, int n_breaks, Eigen::VectorXd &lb, Eigen::VectorXd &ub) = 0;
    virtual void name(const double* t, int n_breaks, std::vector<std::string> &name_str) = 0;
    virtual void updateRobot(RigidBodyManipulator *robot);
    virtual ~MultipleTimeKinematicConstraint(){};
};

class PositionConstraint : public SingleTimeKinematicConstraint
{
  protected:
    double* lb;
    double* ub;
    bool* null_constraint_rows;
    Eigen::MatrixXd pts; 
    int n_pts;
    virtual void evalPositions(Eigen::MatrixXd &pos,Eigen::MatrixXd &J) = 0;
  public:
    PositionConstraint(RigidBodyManipulator *model, const Eigen::MatrixXd &pts,Eigen::MatrixXd lb, Eigen::MatrixXd ub, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    PositionConstraint(const PositionConstraint& rhs);
    virtual void eval(const double* t,Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~PositionConstraint();
};

class WorldPositionConstraint: public PositionConstraint
{
  protected:
    int body;
    std::string body_name;
    virtual void evalPositions(Eigen::MatrixXd &pos, Eigen::MatrixXd &J);
  public:
    WorldPositionConstraint(RigidBodyManipulator *model, int body, const Eigen::MatrixXd &pts, Eigen::MatrixXd lb, Eigen::MatrixXd ub, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual ~WorldPositionConstraint();
};

class WorldCoMConstraint: public PositionConstraint
{
  protected:
    std::set<int> m_robotnum;
    int body;
    std::string body_name;
    virtual void evalPositions(Eigen::MatrixXd &pos, Eigen::MatrixXd &J);
    static const std::set<int> defaultRobotNumSet;
  public:
    WorldCoMConstraint(RigidBodyManipulator *model, Eigen::Vector3d lb, Eigen::Vector3d ub, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan, const std::set<int> &robotnum = WorldCoMConstraint::defaultRobotNumSet);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    void updateRobotnum(const std::set<int> &robotnum);
    virtual ~WorldCoMConstraint(); 
};
/*
class RelativePositionConstraint: public PositionConstraint
{
  protected:
    int bodyA_idx;
    int bodyB_idx;
    std::string bodyA_name;
    std::string bodyB_name;
    Matrix<double,4,4> bodyB_bpTb;
    virtual void evalPositions(Eigen::MatrixXd &pos, Eigen::MatrixXd &J);
  public:
    RelativePositionConstraint(RigidBodyManipulator *model, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &lb, const Eigen::MatrixXd &ub, int bodyA_idx, int bodyB_idx, const Matrix<double,4,4> &bTbp, const Vector2d &tspan);
    virtual void name(const double* t, std::vector<std::string> &name-str);
    virtual ~RelativePositionConstraint();
};*/

class QuatConstraint: public SingleTimeKinematicConstraint
{
  protected:
    double tol;
    virtual void evalOrientationProduct(double &prod, Eigen::MatrixXd &dprod) = 0;
  public:
    QuatConstraint(RigidBodyManipulator *model, double tol, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~QuatConstraint();
};

class WorldQuatConstraint: public QuatConstraint
{
  protected:
    int body;
    std::string body_name;
    Eigen::Vector4d quat_des;
    virtual void evalOrientationProduct(double &prod, Eigen::MatrixXd &dprod);
  public:
    WorldQuatConstraint(RigidBodyManipulator *model, int body, Eigen::Vector4d quat_des, double tol, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual ~WorldQuatConstraint();
};

class EulerConstraint: public SingleTimeKinematicConstraint
{
  protected:
    double* ub;
    double* lb;
    bool null_constraint_rows[3];
    double* avg_rpy;
    virtual void evalrpy(Eigen::Vector3d &rpy, Eigen::MatrixXd &J) = 0;
  public:
    EulerConstraint(RigidBodyManipulator *model, Vector3d lb, Vector3d ub, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    EulerConstraint(const EulerConstraint &rhs);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~EulerConstraint();
};

class WorldEulerConstraint: public EulerConstraint
{
  protected:
    int body;
    std::string body_name;
    virtual void evalrpy(Eigen::Vector3d &rpy, Eigen::MatrixXd &J);
  public:
    WorldEulerConstraint(RigidBodyManipulator *model, int body, Vector3d lb, Vector3d ub, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual ~WorldEulerConstraint();
};

class GazeConstraint : public SingleTimeKinematicConstraint
{
  protected:
    Eigen::Vector3d axis;
    double conethreshold;
  public:
    GazeConstraint(RigidBodyManipulator *model, Eigen::Vector3d axis, double conethreshold = 0.0, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual ~GazeConstraint(void){};
};

class GazeOrientConstraint : public GazeConstraint
{
  protected:
    double threshold;
    Eigen::Vector4d quat_des;
    virtual void evalOrientation(Eigen::Vector4d &quat, Eigen::MatrixXd &dquat_dq) = 0;
  public:
    GazeOrientConstraint(RigidBodyManipulator* model, Eigen::Vector3d axis, Eigen::Vector4d quat_des, double conethreshold, double threshold, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void bounds(const double* t,Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~GazeOrientConstraint(void){};
};

class WorldGazeOrientConstraint: public GazeOrientConstraint
{
  protected:
    int body;
    std::string body_name;
    virtual void evalOrientation(Eigen::Vector4d &quat, Eigen::MatrixXd &dquat_dq);
  public:
    WorldGazeOrientConstraint(RigidBodyManipulator* model, int body, Eigen::Vector3d axis, Eigen::Vector4d quat_des,double conethreshold, double threshold, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void name(const double* t,std::vector<std::string> &name_str);
    virtual ~WorldGazeOrientConstraint(){};
};

class GazeDirConstraint: public GazeConstraint
{
  protected:
    Eigen::Vector3d dir;
  public:
    GazeDirConstraint(RigidBodyManipulator* model, Eigen::Vector3d axis, Eigen::Vector3d dir,double conethreshold, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~GazeDirConstraint(void){};
};

class WorldGazeDirConstraint: public GazeDirConstraint
{
  protected:
    int body;
    std::string body_name;
  public:
    WorldGazeDirConstraint(RigidBodyManipulator* model, int body,Eigen::Vector3d axis, Eigen::Vector3d dir, double conethreshold, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual ~WorldGazeDirConstraint(void){};
};

class GazeTargetConstraint: public GazeConstraint
{
  protected:
    Eigen::Vector3d target;
    Eigen::Vector4d gaze_origin;
  public:
    GazeTargetConstraint(RigidBodyManipulator* model, Eigen::Vector3d axis, Eigen::Vector3d target, Eigen::Vector4d gaze_origin, double conethreshold, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~GazeTargetConstraint(void){};
};

class WorldGazeTargetConstraint: public GazeTargetConstraint
{
  protected:
    int body;
    std::string body_name;
  public:
    WorldGazeTargetConstraint(RigidBodyManipulator* model, int body, Eigen::Vector3d axis, Eigen::Vector3d target, Eigen::Vector4d gaze_origin, double conethreshold, Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual ~WorldGazeTargetConstraint(void){};
};

class RelativeGazeTargetConstraint: public GazeTargetConstraint
{
  protected:
    int bodyA_idx;
    int bodyB_idx;
    std::string bodyA_name;
    std::string bodyB_name;
  public:
    RelativeGazeTargetConstraint(RigidBodyManipulator* model, int bodyA_idx, int bodyB_idx, const Eigen::Vector3d &axis, const Vector3d &target, const Vector4d &gaze_origin,  double conethreshold, Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual ~RelativeGazeTargetConstraint(void){};
};

class Point2PointDistanceConstraint: public SingleTimeKinematicConstraint
{
  protected:
    int bodyA;
    int bodyB;
    Eigen::MatrixXd ptA;
    Eigen::MatrixXd ptB;
    Eigen::VectorXd dist_lb;
    Eigen::VectorXd dist_ub;
  public:
    Point2PointDistanceConstraint(RigidBodyManipulator* model, int bodyA, int bodyB, const Eigen::MatrixXd &ptA, const Eigen::MatrixXd &ptB, const VectorXd &lb, const Eigen::VectorXd &ub, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~Point2PointDistanceConstraint(void){};
};

class Point2LineSegDistConstraint: public SingleTimeKinematicConstraint
{
  protected:
    int pt_body;
    int line_body;
    Eigen::Vector4d pt;
    Eigen::MatrixXd line_ends;
    double dist_lb;
    double dist_ub;
  public:
    Point2LineSegDistConstraint(RigidBodyManipulator* model, int pt_body, const Eigen::Vector4d &pt, int line_body, const Eigen::Matrix<double,4,2> &line_ends,double dist_lb, double dist_ub, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t, Eigen::VectorXd &c, Eigen::MatrixXd &dc);
    virtual void name(const double* t, std::vector<std::string> &name_str);
    virtual void bounds(const double* t, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual ~Point2LineSegDistConstraint(void){};
};

class WorldFixedPositionConstraint: public MultipleTimeKinematicConstraint
{
  protected:
    int body;
    std::string body_name;
    Eigen::MatrixXd pts;
  public:
    WorldFixedPositionConstraint(RigidBodyManipulator* model, int body, const Eigen::MatrixXd &pts,const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual int getNumConstraint(const double* t, int n_breaks);
    virtual void eval_valid(const double* valid_t, int num_valid_t,const Eigen::MatrixXd &valid_q,Eigen::VectorXd &c, Eigen::MatrixXd &dc_valid);
    virtual void bounds(const double* t,int n_breaks, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual void name(const double* t, int n_breaks,std::vector<std::string> &name_str);
    virtual ~WorldFixedPositionConstraint(void){};
};

class WorldFixedOrientConstraint: public MultipleTimeKinematicConstraint
{
  protected:
    int body;
    std::string body_name;
  public:
    WorldFixedOrientConstraint(RigidBodyManipulator* model, int body, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual int getNumConstraint(const double* t, int n_breaks);
    virtual void eval_valid(const double* valid_t, int num_valid_t,const Eigen::MatrixXd &valid_q,Eigen::VectorXd &c, Eigen::MatrixXd &dc_valid);
    virtual void bounds(const double* t,int n_breaks, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual void name(const double* t, int n_breaks,std::vector<std::string> &name_str);
    virtual ~WorldFixedOrientConstraint(void){};
};

class WorldFixedBodyPoseConstraint: public MultipleTimeKinematicConstraint
{
  protected:
    int body;
    std::string body_name;
  public:
    WorldFixedBodyPoseConstraint(RigidBodyManipulator* model, int body, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual int getNumConstraint(const double* t, int n_breaks);
    virtual void eval_valid(const double* valid_t, int num_valid_t,const Eigen::MatrixXd &valid_q,Eigen::VectorXd &c, Eigen::MatrixXd &dc_valid);
    virtual void bounds(const double* t,int n_breaks, Eigen::VectorXd &lb, Eigen::VectorXd &ub);
    virtual void name(const double* t, int n_breaks,std::vector<std::string> &name_str);
    virtual ~WorldFixedBodyPoseConstraint(void){};
};

class AllBodiesClosestDistanceConstraint : public SingleTimeKinematicConstraint
{
  protected:
    double ub;
    double lb;
  public:
    AllBodiesClosestDistanceConstraint(RigidBodyManipulator* model, 
                                       double lb, double ub,
                                       Eigen::Vector2d tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual void eval(const double* t,Eigen::VectorXd& c, Eigen::MatrixXd& dc);
    virtual void name(const double* t, std::vector<std::string> &name);
    virtual void bounds(const double* t, Eigen::VectorXd& lb, Eigen::VectorXd& ub);
    virtual ~AllBodiesClosestDistanceConstraint(){};
};

class WorldPositionInFrameConstraint: public WorldPositionConstraint
{
  protected:
    Eigen::Matrix4d T_world_to_frame;
    Eigen::Matrix4d T_frame_to_world;
    virtual void evalPositions(Eigen::MatrixXd &pos, Eigen::MatrixXd &J);
  public:
    WorldPositionInFrameConstraint(RigidBodyManipulator *model, int body, 
        const Eigen::MatrixXd &pts, const Eigen::Matrix4d& T_world_to_frame, 
        Eigen::MatrixXd lb, Eigen::MatrixXd ub, const Eigen::Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    //virtual void updateRobot(RigidBodyManipulator *robot);
    virtual ~WorldPositionInFrameConstraint();
};

class PostureChangeConstraint: public MultipleTimeLinearPostureConstraint
{
  protected:
    Eigen::VectorXi joint_ind;
    Eigen::VectorXd lb_change;
    Eigen::VectorXd ub_change;
    virtual void setJointChangeBounds(const Eigen::VectorXi &joint_ind, const Eigen::VectorXd &lb_change, const Eigen::VectorXd &ub_change);
  public:
    PostureChangeConstraint(RigidBodyManipulator* model, const Eigen::VectorXi &joint_ind, const Eigen::VectorXd &lb_change, const Eigen::VectorXd &ub_change, const Vector2d &tspan = DrakeRigidBodyConstraint::default_tspan);
    virtual int getNumConstraint(const double* t, int n_breaks);
    virtual void feval(const double* t, int n_breaks, const MatrixXd &q, VectorXd &c);
    virtual void geval(const double* t, int n_breaks, VectorXi &iAfun, VectorXi &jAvar, VectorXd &A);
    virtual void name(const double* t, int n_breaks, std::vector<std::string> &name_str);
    virtual void bounds(const double* t, int n_breaks, VectorXd &lb, VectorXd &ub);
    virtual ~PostureChangeConstraint(){};
};
#endif


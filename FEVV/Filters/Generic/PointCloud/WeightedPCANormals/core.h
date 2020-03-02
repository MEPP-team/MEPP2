// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

// -------------------------------------------------------------------------------------
//        Iterative weighted PCA for robust and edge-aware normal vector estimation
//--------------------------------------------------------------------------------------
// Julia Sanchez, Florence Denis, David Coeurjolly, Florent dupont, Laurent Trassoudaine, Paul Checchin
// Liris (Lyon), Institut Pascal (Clermont Ferrand)
// Région Auvergne Rhône Alpes ARC6
// Private use for reviewers only
// --------------------------------------------------------------------------------------


#include <math.h>
#include <algorithm>
#include <chrono>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <flann/flann.hpp>

using namespace std;

//Constants with effect on result
const float mu_max2 = 0.35;         // initial mu when starting the optimization for selection 2
                                    // increase if too many normals stuck in medium position/ decrease if too many face confusions when sampling anisotropy
const float div_fact = 1.01;

//Usual constants DO NOT MODIFY
const float epsilon = 1e-10;        // small value for initializations of comparisons
const float lim_error = 1e-5;       // to stop optimize when not moving
const float lim_diff = 1e-7;        // to stop optimize when not moving

//Detail constants
const int itr_min = 5;                      // minimum number of iterations to perform
const float mu_max = 1.0;                   // coefficient to compute the initial mu when starting optimize for the first time (emax = er_max*mu_max)
const int itr_per_mu = 1;                   // number of iterations to perform for each mu
const float noise_min = 0.000001;           // minimum noise to make it not infinite
const float likelihood_threshold = 0.95;    // value for evaluation/comparison between vectors
const float min_points_fact = 0.1;          // limit number of points of the neighborhood which must have an impact on the normal computation
const float thresh_weight = 0.25;           // threshold to select neighbors


class CApp{
public:
    CApp()
    {
        std::cout<<"add one pointcloud + its kdtree + the reference point you want to compute the normal on"<<std::endl<<std::endl;
        normalFirst0_ = NULL;
        normalFirst1_ = NULL;
        normalFirst2_ = NULL;
        normalSecond0_ = NULL;
        normalSecond1_ = NULL;
        normalSecond2_ = NULL;
        pointFirst_ = NULL;
        pointSecond_ = NULL;
    };

    CApp(Eigen::Matrix<float, Eigen::Dynamic, 3> *pc, flann::Index<flann::L2<float>>* tree, float noise)
    {
        setPc(pc);
        setTree(tree);
        noise_ = noise;
        normalFirst0_ = NULL;
        normalFirst1_ = NULL;
        normalFirst2_ = NULL;
        normalSecond0_ = NULL;
        normalSecond1_ = NULL;
        normalSecond2_ = NULL;
        pointFirst_ = NULL;
        pointSecond_ = NULL;
    };

    CApp(Eigen::Matrix<float, Eigen::Dynamic, 3> *pc, flann::Index<flann::L2<float>>* tree, int ref, float noise)
    {
        setPc(pc);
        setTree(tree);
        setRef(ref);
        noise_ = noise;
        normalFirst0_ = NULL;
        normalFirst1_ = NULL;
        normalFirst2_ = NULL;
        normalSecond0_ = NULL;
        normalSecond1_ = NULL;
        normalSecond2_ = NULL;
        pointFirst_ = NULL;
        pointSecond_ = NULL;
    };

    CApp(float divFact, float limMu, float limMuPos, Eigen::Matrix<float, Eigen::Dynamic, 3> *pc, flann::Index<flann::L2<float>>* tree, int ref, float noise)
    {
        divFact_ = divFact;
        limMu_ = limMu;
        limMuPos_ = limMuPos;
        noise_ = noise;
        setPc(pc);
        setTree(tree);
        setRef(ref);
        normalFirst0_ = NULL;
        normalFirst1_ = NULL;
        normalFirst2_ = NULL;
        normalSecond0_ = NULL;
        normalSecond1_ = NULL;
        normalSecond2_ = NULL;
        pointFirst_ = NULL;
        pointSecond_ = NULL;

    };

    ~CApp()
    {
        if(normalFirst0_!= NULL)
            delete normalFirst0_;
        if(normalFirst1_!= NULL)
            delete normalFirst1_;
        if(normalFirst2_!= NULL)
            delete normalFirst2_;
        if(normalSecond0_!= NULL)
            delete normalSecond0_;
        if(normalSecond1_!= NULL)
            delete normalSecond1_;
        if(normalSecond2_!= NULL)
            delete normalSecond2_;
        if(pointFirst_!= NULL)
            delete pointFirst_;
        if(pointSecond_!= NULL)
            delete pointSecond_;
    };

    void setTree(flann::Index<flann::L2<float>> *t){tree_ = t;}
    void setPc(Eigen::Matrix<float, Eigen::Dynamic, 3> *pc){pointcloud_ = pc;};
    void setRef (int ref);
    void setNormal(Eigen::Vector3f norm){normal = norm;};
    void setPoint(Eigen::Vector3f point);
    void setLimMu ( float limMu){limMu_ = limMu;}
    void setNbrNeighbors(int N){n_neigh_ = N;}
    void setLimMuPos ( float limMuPos){limMuPos_ = limMuPos;}
    void setDivFact ( float divFact){divFact_= divFact;}
    void setNoise ( float noise){noise_= noise;}
    void setParams(float divFact, float curvature)
    {
        divFact_= divFact;
        curvature_ = curvature;
    };
    Eigen::Vector3f getNormal(){return normal;}
    Eigen::Vector3f getPoint (){return pt;}
    Eigen::Matrix<float, Eigen::Dynamic, 3> getNeighborhood(){return neighborhood_;}
    int get_N_neigh(){return n_neigh_;}

    void reinitPoint();
    void ComputeDist();
    void SearchFLANNTreeKnn(flann::Index<flann::L2<float>>* index,
                                Eigen::Vector3f& input,
                                std::vector<int>& indices,
                                std::vector<float>& dists,
                                int nn);
    void selectNeighborsKnn(int N);
    void selectNeighborsRadius(float r);

    void init1();
    void reinitFirst0();
    void init2();
    void Optimize(bool first);
    void OptimizePos(bool first, float thresh_weigh);
    Eigen::Vector3f getEdgeDirection(int it);
    void evaluate(int *impact, float *moyError, float impacter_weigh);
    void select_normal();
    void addSymmetricDist();
    bool isOnEdge();
    bool isSecondOption();



    bool isNan();
    Eigen::Vector3f finalNormal_;
    Eigen::Vector3f finalPos_;

    Eigen::Vector3f normal;
    bool SuspectedOnEdge_;


private:
    int n_neigh_;
    float mu_;
    float noise_;
    float limMu_;
    float limMuPos_;
    float curvature_;
    float divFact_;
    float thresh_proj_;
    Eigen::Matrix<float, Eigen::Dynamic, 3> *pointcloud_;
    flann::Index<flann::L2<float>>* tree_;
    Eigen::Matrix<float, Eigen::Dynamic, 3> neighborhood_;
    Eigen::Matrix<float, Eigen::Dynamic, 3> dist_;
    int ref_;
    Eigen::Vector3f pt;
    Eigen::Vector3f ptRef_;
    Eigen::Vector3f centroid_;
    Eigen::Vector3f* normalFirst0_;
    Eigen::Vector3f* normalFirst1_;
    Eigen::Vector3f* normalFirst2_;
    Eigen::Vector3f* normalSecond0_;
    Eigen::Vector3f* normalSecond1_;
    Eigen::Vector3f* normalSecond2_;
    Eigen::Vector3f* pointFirst_;
    Eigen::Vector3f* pointSecond_;
    float theta;
    float phi;
    Eigen::VectorXf weights_;
    std::vector<float> error_;
    std::vector<float> diff_error_;
    void ComputeWeights();
    float GetMaxError();
    void ComputeNormProjError(std::vector<float>& er_proj);
    void orient();
    void actualizeMu();
    void init_weight();
    Eigen::Vector3f getThirdEigenVector(Eigen::Matrix3f& C);
    void save_itr(int itr);
    int impactFirst_;
    int impactSecond_;
    float moyErrorFirst_;
    float moyErrorSecond_;
    bool jamais_fait;
    float er_max;
    float mu_init2_;
};

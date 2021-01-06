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


#if defined _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4267 4244)
#endif

//Extract neighbors
inline void CApp::selectNeighborsKnn(int N)
{
    n_neigh_ = N;
    std::vector<float> dis;
    std::vector<int> neigh;

    int Required_neighbors_nbr = n_neigh_+1;
    SearchFLANNTreeKnn(tree_, pt, neigh, dis, Required_neighbors_nbr);
    neighborhood_.resize (n_neigh_,3);
    std::vector<int> to_erase;

    for (int i = 1; i < neigh.size(); ++i)
    {
        if( dis[i] != 0)
            neighborhood_.row(i-1) = pointcloud_->row(neigh[i]);
        else
        {
            ++Required_neighbors_nbr;
            SearchFLANNTreeKnn(tree_, pt, neigh, dis, Required_neighbors_nbr);
            to_erase.push_back(i);
            for(int k = to_erase.size()-1; k>=0; --k)
            {
                neigh.erase(neigh.begin() + to_erase[k]);
                dis.erase(dis.begin() + to_erase[k]);
            }
            --i;
        }
    }
    weights_.resize(neighborhood_.rows());
    dist_.resize (neighborhood_.rows(), 3);
    ComputeDist();
    jamais_fait = 1;
}

//Flann search (internally used in selectNeighborsKnn )
inline void CApp::SearchFLANNTreeKnn(flann::Index<flann::L2<float>>* index,
                            Eigen::Vector3f& input,
                            std::vector<int>& indices,
                            std::vector<float>& dists,
                            int nn)
{
    int dim = input.size();

    std::vector<float> query;
    query.resize(dim);
    for (int i = 0; i < dim; i++)
        query[i] = input(i);
    flann::Matrix<float> query_mat(&query[0], 1, dim);

    indices.resize(nn);
    dists.resize(nn);
    flann::Matrix<int> indices_mat(&indices[0], 1, nn);
    flann::Matrix<float> dists_mat(&dists[0], 1, nn);

  /**
   * Find the "num_closest" nearest neighbors to the \a query_point[0:dim-1].
   * Their indices are stored inside the result object. \sa radiusSearch,
   * findNeighbors \note nChecks_IGNORED is ignored but kept for compatibility
   * with the original FLANN interface. \return Number `N` of valid points in
   * the result set. Only the first `N` entries in `out_indices` and
   * `out_distances_sq` will be valid. Return may be less than `num_closest`
   * only if the number of elements in the tree is less than `num_closest`.
   */
  /*
    size_t knnSearch(
        const ElementType *query_point, 
        const size_t num_closest,
        IndexType *out_indices, 
        DistanceType *out_distances_sq,
        )
     */

    index->knnSearch(query_mat, indices_mat, dists_mat, nn, flann::SearchParams(128));
}


//Fills dist_ with distances of the neighbors from current PCA reference point
inline void CApp::ComputeDist()
{
    dist_ = neighborhood_ - Eigen::VectorXf::Ones(n_neigh_) * pt.transpose();
}

//Computes 3rd eigen vector of a covariance matrix
inline Eigen::Vector3f CApp::getThirdEigenVector(Eigen::Matrix3f& C)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(C);
    return es.eigenvectors().col(0);
}

//first initialization with PCA + initialize mu_lim and tau
inline void CApp::init1()
{
    //PCA
    Eigen::MatrixXf centered_points (n_neigh_, 3);
    Eigen::Vector3f mean_point = neighborhood_.colwise().mean();
    centered_points = neighborhood_ - Eigen::VectorXf::Ones(n_neigh_)*(mean_point.transpose());

    Eigen::Matrix3f C = centered_points.transpose()*centered_points/n_neigh_;
    normal = getThirdEigenVector(C);

    normalFirst0_ = new Eigen::Vector3f();
    *normalFirst0_ = normal;
    finalNormal_ = normal;
    finalPos_ = pt;

    orient();

    //mu_lim and tau
    float r = 0;
    int N_boundary_neihborhood = (int)(0.2*neighborhood_.rows());
    for(int i = 1; i<=N_boundary_neihborhood; ++i)
        r += (neighborhood_.row(neighborhood_.rows()-i)-pt.transpose()).norm();
    r /= N_boundary_neihborhood;
    float Xmax = r*r/(2.0*curvature_);
    float theta_m = 2*acos(1-Xmax/curvature_);
    float mean;
    float thresh;
    if(abs(theta_m)>epsilon)
    {
        mean = 2*curvature_*sin(theta_m/2)/theta_m; // mean of distances to plane
        thresh = sqrt((pow(curvature_,2)/2)*(1+(sin(theta_m)/theta_m)) - pow(mean,2));//sqrt( (1/theta_m)*( pow(curvature_,2)*( (theta_m/4) - (sin(theta_m)/4) ) - 2*mean*curvature_*sin(theta_m/2) + pow(mean,2)*theta_m/2 ) );
    }
    else
    {
        mean = 0;
        thresh = 0;
    }
    float emax = Xmax+noise_/2;

    thresh_proj_ = (thresh + noise_);
    limMuPos_ = emax*emax;
}

inline void CApp::reinitFirst0()
{
    normal = *normalFirst0_;
}

//second initialization -> perpendicular to edge and n01
inline void CApp::init2()
{
    normal = *normalFirst2_;
    if(abs(normalFirst0_->dot(*normalFirst2_))<0.996)
    {
        Eigen::Vector3f edge_direction = normalFirst0_->cross(*normalFirst2_);
        edge_direction /= edge_direction.norm();
        normal = edge_direction.cross(*normalFirst2_);
        normal /= normal.norm();
        normalSecond0_ = new Eigen::Vector3f();
        *normalSecond0_ = normal;
        SuspectedOnEdge_ = 1;
    }
    else
    {
        SuspectedOnEdge_ = 0;
        finalNormal_ = *normalFirst2_;
        finalPos_ = *pointFirst_;
    }
}

//First weighting -> all weights in 1
inline void CApp::init_weight()
{
    for (int i = 0; i<weights_.size(); ++i)
        weights_[i] = 1;
}

//Checks if point has neigborhood on edge (use tau)
inline bool CApp::isOnEdge()
{
    Eigen::Vector3f pt_moy = Eigen::Vector3f::Zero();
    for (int c = 0; c < n_neigh_; c++)
        pt_moy += neighborhood_.row(c);

    pt_moy /= n_neigh_;

    float proj_moy = 0;
    for (int c = 0; c < n_neigh_; c++)
        proj_moy += pow(normal.dot(neighborhood_.row(c).transpose()-pt_moy),2);

    proj_moy /= n_neigh_;
    proj_moy = sqrt(proj_moy);
    return ( proj_moy>thresh_proj_ );
}

//Actualizes weights
inline void CApp::ComputeWeights()
{
    for (int c = 0; c < n_neigh_; c++)
        weights_(c) = pow(mu_/(mu_+pow(normal.dot(dist_.row(c)),2)),2); //Scaled version of Geman McClure estimator
}

//Compute maximal residual
inline float CApp::GetMaxError()
{
    float er_max = 0;
    float er_proj;
    for (int c = 0; c < n_neigh_; c++)
    {
        er_proj = normal.dot(dist_.row(c));
        if(abs(er_proj)>er_max)
            er_max = abs(er_proj);
    }
    return er_max;
}

//Rough estimation
inline void CApp::Optimize(bool first)
{
    ComputeDist();

    if (first)
    {
        er_max = GetMaxError();
        mu_ = pow(mu_max*er_max,2);
        mu_ = std::max(mu_,limMuPos_);
    }
    else
    {
        std::vector<float> er_proj_squared(n_neigh_);
        for (int c = 0; c < n_neigh_; c++)
            er_proj_squared[c] = pow(normal.dot(dist_.row(c)),2);
        std::sort(er_proj_squared.begin(), er_proj_squared.end());
        mu_ = er_proj_squared[(int)(n_neigh_*mu_max2)];
        mu_ = std::max(mu_,limMuPos_);

        mu_init2_ = mu_;
    }

    int itr = 0;

    while(itr < itr_min || mu_ > limMuPos_)
    {
        if( mu_ > limMuPos_ && (itr % itr_per_mu) == 0)
            mu_ /= divFact_;

        ComputeWeights();
        Eigen::Matrix3f C = dist_.transpose()*weights_.asDiagonal()*dist_;// /weights_.sum();

        normal = getThirdEigenVector(C);
        ++itr;
    }

    if (first)
    {
        normalFirst1_ = new Eigen::Vector3f();
        *normalFirst1_ = normal;
    }
    else
    {
        normalSecond1_ = new Eigen::Vector3f();
        *normalSecond1_ = normal;
    }
}

//Refinement
inline void CApp::OptimizePos(bool first, float thresh_weight)
{
    if(!first)
    {
        if(!isSecondOption())
        {
            finalNormal_ = *normalFirst2_;
            finalPos_ = *pointFirst_;
            return;
        }
          mu_ = 0.5*mu_init2_;
          mu_ = std::max(mu_,limMuPos_);
    }
    else
    {
        mu_ = pow(0.5*er_max,2);
        mu_ = std::max(mu_,limMuPos_);
    }

    //-----------------------------------------------------------------------------------------
    float sum_poids = 0;
    float moy_proj = 0;

    int it = itr_per_mu*(int)(log(mu_/(limMuPos_) )/log(divFact_));

    if(it<itr_min)
        it=itr_min;

    int itr = 0;

    while( itr< it)
    {
        if( mu_ > limMuPos_ && (itr % itr_per_mu) == 0)
            mu_ /= divFact_;

        //actualize position of the PCA reference

        sum_poids = 0;
        moy_proj = 0;

        for (int c = 0; c < n_neigh_; c++)
        {
            if(weights_(c) > thresh_weight)
            {
                sum_poids += weights_(c);
                moy_proj += weights_(c)*normal.dot(dist_.row(c));
            }
        }
        moy_proj /= sum_poids;

        pt = pt + moy_proj*normal;

        ComputeDist();
        ComputeWeights();

        //actualize normal
        Eigen::Matrix3f C = dist_.transpose()*weights_.asDiagonal()*dist_;
        normal = getThirdEigenVector(C);

         ++itr;
    }

     orient();

     if(first)
     {
         normalFirst2_ = new Eigen::Vector3f();
         *normalFirst2_ = normal;
         pointFirst_ = new Eigen::Vector3f();
         *pointFirst_ = pt;

         evaluate(&impactFirst_, &moyErrorFirst_, thresh_weight);
     }
     else
     {
         normalSecond2_ = new Eigen::Vector3f();
         *normalSecond2_ = normal;
         pointSecond_ = new Eigen::Vector3f();
         *pointSecond_ = pt;
         evaluate(&impactSecond_, &moyErrorSecond_, thresh_weight);
         select_normal();
     }
}

//orient normal to the exterior of the edge
inline void CApp::orient()
{
    float moy_err = 0;
    for (int c = 0; c < n_neigh_; c++)
        moy_err += normal.dot(dist_.row(c));

    if(moy_err> epsilon)
        normal = -normal;
}

inline bool CApp::isSecondOption()
{
    return (abs(normalFirst1_->dot(*normalSecond1_))<likelihood_threshold);
}

inline void CApp::setRef( int ref)
{
    ref_ = ref;
    pt = pointcloud_->row(ref);
    ptRef_ = pt;
}

inline void CApp::reinitPoint()
{
    setPoint(ptRef_);
}

inline void CApp::evaluate(int *impact, float *mean, float weight_tresh)
{
    float imp = 0;
    for(int c = 0; c<n_neigh_; ++c)
    {
        if(weights_(c)>weight_tresh)
            ++imp;
    }

    *impact = imp;
    *mean = (pt-ptRef_).dot(normal);
}

inline void CApp::setPoint(Eigen::Vector3f point)
{
    pt = point;
    ComputeDist();
}

//Check if result contains nan
inline bool CApp::isNan()
{
    return (finalNormal_(0) != finalNormal_(0) || finalPos_(0) != finalPos_(0));
}

inline void CApp::select_normal()
{
    int min_points = (int)(min_points_fact*(float)neighborhood_.rows());
    if( moyErrorFirst_<moyErrorSecond_ && impactFirst_>min_points )
    {
        finalNormal_ = *normalFirst2_;
        finalPos_ = *pointFirst_;
    }
    else if(impactSecond_>min_points)
    {
        finalNormal_  = *normalSecond2_;
        finalPos_  = *pointSecond_;
    }
    else if(impactFirst_>min_points)
    {
        finalNormal_  = *normalFirst2_;
        finalPos_  = *pointFirst_;
    }
    else if( moyErrorFirst_<moyErrorSecond_ )
    {
        finalNormal_  = *normalFirst2_;
        finalPos_  = *pointFirst_;
    }
    else
    {
        finalNormal_  = *normalSecond2_;
        finalPos_  = *pointSecond_;
    }
}


#if defined _MSC_VER
#pragma warning(pop)
#endif

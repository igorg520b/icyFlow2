#ifndef MODELPRMS_H
#define MODELPRMS_H

#include <QObject>
#include <cmath>
#include "numbercrunching.h"

namespace icy {
class ModelPrms;
class NumberCrunching;
}

class icy::ModelPrms : public QObject
{
    Q_OBJECT
    Q_PROPERTY(double sim_InitialTimeStep MEMBER InitialTimeStep NOTIFY propertyChanged)
    Q_PROPERTY(int sim_MaxSteps MEMBER MaxSteps NOTIFY propertyChanged)
    Q_PROPERTY(double sim_IndentationVelocity MEMBER IndentationVelocity NOTIFY propertyChanged)

    // integration
    Q_PROPERTY(double intg_NewmarkBeta MEMBER NewmarkBeta NOTIFY propertyChanged)
    Q_PROPERTY(double intg_NewmarkGamma MEMBER NewmarkGamma NOTIFY propertyChanged)
    Q_PROPERTY(double intg_ConvergenceEpsilon MEMBER ConvergenceEpsilon NOTIFY propertyChanged)
    Q_PROPERTY(double intg_ConvergenceCutoff MEMBER ConvergenceCutoff NOTIFY propertyChanged)
    Q_PROPERTY(double intg_maxDamagePerStep MEMBER maxDamagePerStep NOTIFY propertyChanged)
    Q_PROPERTY(double intg_maxFailPerStep MEMBER maxFailPerStep NOTIFY propertyChanged)
    Q_PROPERTY(double intg_maxIterations MEMBER maxIterations NOTIFY propertyChanged)
    Q_PROPERTY(double intg_minIterations MEMBER minIterations NOTIFY propertyChanged)
//    Q_PROPERTY(double intg_gravity MEMBER gravity NOTIFY propertyChanged)

    // collisions
    Q_PROPERTY(double coll_penaltyK MEMBER penaltyK NOTIFY propertyChanged)
    Q_PROPERTY(double coll_DistanceEpsilon MEMBER DistanceEpsilon NOTIFY propertyChanged)
    Q_PROPERTY(double coll_ReconstructBVH MEMBER ReconstructBVH NOTIFY propertyChanged)

    // material
    Q_PROPERTY(double mat_Y MEMBER Y NOTIFY propertyChanged)
    Q_PROPERTY(double mat_rho WRITE set_rho READ get_rho)
    Q_PROPERTY(double mat_nu WRITE set_nu READ get_nu)
    Q_PROPERTY(double mat_dampingMass MEMBER dampingMass NOTIFY propertyChanged)
    Q_PROPERTY(double mat_dampingStiffness MEMBER dampingStiffness NOTIFY propertyChanged)

    // cz
    Q_PROPERTY(double cz_alpha WRITE setCZAlpha READ getCZAlpha)
    Q_PROPERTY(double cz_beta WRITE setCZBeta READ getCZBeta)
    Q_PROPERTY(double cz_lambda_n WRITE set_lambda_n READ get_lambda_n)
    Q_PROPERTY(double cz_lambda_t WRITE set_lambda_t READ get_lambda_t)
    Q_PROPERTY(double cz_phi_n WRITE set_phi_n READ get_phi_n)
    Q_PROPERTY(double cz_phi_t WRITE set_phi_t READ get_phi_t)
    Q_PROPERTY(double cz_sigma_max WRITE set_sigma_max READ get_sigma_max)
    Q_PROPERTY(double cz_tau_max WRITE set_tau_max READ get_tau_max)
    Q_PROPERTY(double cz_deln READ get_del_n)
    Q_PROPERTY(double cz_delt READ get_del_t)

public:
    double IndentationVelocity = 0.1;//0.0001; // 0.1 mm/s
    double InitialTimeStep = 0.05;
    int MaxSteps = 100;
    double nThreshold = 0, tThreshold = 0; // CZ peak traction values

    // material
    double Y = 6000000;//5e9;
    double rho = 916.2;
    double dampingMass = 0.0005;
    double dampingStiffness = 0.0005;
    double nu = 0.3;

    // integration
    double NewmarkBeta = 0.5;
    double NewmarkGamma = 1.0;
    double ConvergenceEpsilon = 0.005;
    double ConvergenceCutoff = 1E-8;
    double maxDamagePerStep = 0.01;
    double maxFailPerStep = 0.01;
    int maxIterations = 10;
    int minIterations = 3;
//    double gravity = -9.8;

    // collisions
    double penaltyK = 30000;
    double DistanceEpsilon = 1E-15;
    int ReconstructBVH = 10;

    // cz parameters
    double alpha = 4, beta = 4, lambda_n = 0.015, lambda_t = 0.015;
    double phi_n = 3, phi_t = 3; // fracture energy
    double sigma_max = 120000, tau_max = 140000;
    double del_n = 0, del_t = 0;

    // computed variables
    double totalVolume;
    double M[12][12] = {};   // mass distribution matrix (needs to be multiplied by mass of individual elements)
    double p_m, p_n;
    double pMtn, pMnt; // < phi_t - phi_n >, < phi_n - phi_t >
    double gam_n, gam_t;

    double E[6][6]={}; // recompute when nu or rho changes

    ModelPrms() {
        // initialize M

        double coeff = 1.0 / 20.0;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int m = 0; m < 3; m++)
                {
                    int col = i * 3 + m;
                    int row = j * 3 + m;
                    M[col][row] = (col == row) ? 2 * coeff : coeff;
                }




        Recompute();
        RecomputeE();
    }

private:
    void setCZAlpha(double value) { alpha=value; Recompute(); emit propertyChanged(); }
    double getCZAlpha() {return alpha;}
    void setCZBeta(double value) { beta=value; Recompute(); emit propertyChanged(); }
    double getCZBeta() {return beta;}
    void set_lambda_n(double value) {lambda_n = value; Recompute(); emit propertyChanged(); }
    double get_lambda_n() {return lambda_n;}
    void set_lambda_t(double value) {lambda_t = value; Recompute(); emit propertyChanged(); }
    double get_lambda_t() {return lambda_t;}
    void set_phi_n(double value) {phi_n = value; Recompute(); emit propertyChanged(); }
    double get_phi_n() {return phi_n;}
    void set_phi_t(double value) {phi_t = value; Recompute(); emit propertyChanged(); }
    double get_phi_t() {return phi_t;}
    void set_sigma_max(double value) {sigma_max = value; Recompute(); emit propertyChanged(); }
    double get_sigma_max() {return sigma_max;}
    void set_tau_max(double value) {tau_max = value; Recompute(); emit propertyChanged(); }
    double get_tau_max() {return tau_max;}
    double get_del_n() {return del_n;}
    double get_del_t() {return del_t;}

    double Macaulay(double a, double b) { if (a > b) return a - b; else return 0; }

    void Recompute()
    {
        pMnt = Macaulay(phi_n, phi_t);
        pMtn = Macaulay(phi_t, phi_n);

        double rn_sq = lambda_n * lambda_n;
        double rt_sq = lambda_t * lambda_t;
        p_m = (alpha * (alpha - 1.0) * rn_sq) / (1.0 - alpha * rn_sq);
        p_n = (beta * (beta - 1.0) * rt_sq) / (1.0 - beta * rt_sq);

        if (phi_n < phi_t)
        {
            gam_n = pow(alpha / p_m, p_m);
            gam_t = -phi_t * pow(beta / p_n, p_n);
        }
        else
        {
            gam_n = -phi_n * pow(alpha / p_m, p_m);
            gam_t = pow(beta / p_n, p_n);
        }

        del_n = (phi_n / sigma_max) * alpha * lambda_n * pow((1.0 - lambda_n), (alpha - 1.0)) * ((alpha / p_m) + 1.0) * pow(((alpha / p_m) * lambda_n + 1.0), (p_m - 1.0));
        del_t = (phi_t / tau_max) * beta * lambda_t * pow((1.0 - lambda_t), (beta - 1.0)) * ((beta / p_n) + 1.0) * pow(((beta / p_n) * lambda_t + 1.0), (p_n - 1.0));

        nThreshold = del_n * lambda_n;
        tThreshold = del_t * lambda_t;
    }

    void RecomputeE()
    {
        double coeff1 = Y / ((1.0 + nu) * (1.0 - 2.0 * nu));
        E[0][0] = E[1][1] = E[2][2] = (1.0 - nu) * coeff1;
        E[0][1] = E[0][2] = E[1][2] = E[1][0] = E[2][0] = E[2][1] = nu * coeff1;
        E[3][3] = E[4][4] = E[5][5] = (0.5 - nu) * coeff1;
    }

    void set_rho(double value) {rho = value; RecomputeE(); emit propertyChanged(); }
    double get_rho() {return rho;}
    void set_nu(double value) {nu = value; RecomputeE(); emit propertyChanged(); }
    double get_nu() {return nu;}

signals:
    void propertyChanged();
};

#endif // MODELPRMS_H

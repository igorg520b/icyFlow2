#ifndef MODELPRMS_H
#define MODELPRMS_H

#include <QObject>
#include <cmath>

namespace icy { class ModelPrms; }

class icy::ModelPrms : public QObject
{
    Q_OBJECT
    Q_PROPERTY(double sim_InitialTimeStep MEMBER InitialTimeStep NOTIFY propertyChanged)
    Q_PROPERTY(int sim_MaxSteps MEMBER MaxSteps NOTIFY propertyChanged)
    Q_PROPERTY(double sim_IndentationVelocity MEMBER IndentationVelocity NOTIFY propertyChanged)

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
    double IndentationVelocity = 0.0001; // 0.1 mm/s
    double InitialTimeStep = 0.01;
    int MaxSteps = 300;

    double nThreshold = 0, tThreshold = 0; // CZ peak traction values

private:
    // cz parameters
    double _alpha = 3, _beta = 3, _lambda_n = 0.02, _lambda_t = 0.02;
    double _phi_n = 30, _phi_t = 30;
    double _sigma_max = 4e5, _tau_max = 15e5;
    double _del_n = 0, _del_t = 0;

    // computed variables
    double totalVolume;
    double E[6][6];
    double M[12][12];
    double G_fn, G_ft; // fracture energy
    double f_tn, f_tt;
    double rn, rt;       // lambda_n, lambda_t
    double p_m, p_n;
    double pMtn, pMnt; // < phi_t - phi_n >, < phi_n - phi_t >
    double gam_n, gam_t;
    double sf[3][3];
    double B[3][18];

    void setCZAlpha(double value) { _alpha=value; Recompute(); emit propertyChanged(); }
    double getCZAlpha() {return _alpha;}
    void setCZBeta(double value) { _beta=value; Recompute(); emit propertyChanged(); }
    double getCZBeta() {return _beta;}
    void set_lambda_n(double value) {_lambda_n = value; Recompute(); emit propertyChanged(); }
    double get_lambda_n() {return _lambda_n;}
    void set_lambda_t(double value) {_lambda_t = value; Recompute(); emit propertyChanged(); }
    double get_lambda_t() {return _lambda_t;}
    void set_phi_n(double value) {_phi_n = value; Recompute(); emit propertyChanged(); }
    double get_phi_n() {return _phi_n;}
    void set_phi_t(double value) {_phi_t = value; Recompute(); emit propertyChanged(); }
    double get_phi_t() {return _phi_t;}
    void set_sigma_max(double value) {_sigma_max = value; Recompute(); emit propertyChanged(); }
    double get_sigma_max() {return _sigma_max;}
    void set_tau_max(double value) {_tau_max = value; Recompute(); emit propertyChanged(); }
    double get_tau_max() {return _tau_max;}
    double get_del_n() {return _del_n;}
    double get_del_t() {return _del_t;}

    double Macaulay(double a, double b) { if (a > b) return a - b; else return 0; }

    void Recompute()
    {
        pMnt = Macaulay(_phi_n, _phi_t);
        pMtn = Macaulay(_phi_t, _phi_n);

        double rn_sq = _lambda_n * _lambda_n;
        double rt_sq = _lambda_t * _lambda_t;
        p_m = (_alpha * (_alpha - 1.0) * rn_sq) / (1.0 - _alpha * rn_sq);
        p_n = (_beta * (_beta - 1.0) * rt_sq) / (1.0 - _beta * rt_sq);

        if (_phi_n < _phi_t)
        {
            gam_n = pow(_alpha / p_m, p_m);
            gam_t = -_phi_t * pow(_beta / p_n, p_n);
        }
        else
        {
            gam_n = -_phi_n * pow(_alpha / p_m, p_m);
            gam_t = pow(_beta / p_n, p_n);
        }

        _del_n = (_phi_n / _sigma_max) * _alpha * _lambda_n * pow((1.0 - _lambda_n), (_alpha - 1.0)) * ((_alpha / p_m) + 1.0) * pow(((_alpha / p_m) * _lambda_n + 1.0), (p_m - 1.0));
        _del_t = (_phi_t / _tau_max) * _beta * _lambda_t * pow((1.0 - _lambda_t), (_beta - 1.0)) * ((_beta / p_n) + 1.0) * pow(((_beta / p_n) * _lambda_t + 1.0), (p_n - 1.0));

        nThreshold = _del_n * _lambda_n;
        tThreshold = _del_t * _lambda_t;
    }

signals:
    void propertyChanged();
};

#endif // MODELPRMS_H

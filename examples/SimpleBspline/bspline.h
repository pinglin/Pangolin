#ifndef BSPLIN_H
#define BSPLIN_H

#include <vector>
#include <array>

#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

static Matrix4d CubicBsplineMatrix()
{
    Matrix4d C;
    C << 1.0, 4.0, 1.0, 0.0, -3.0, 0.0, 3.0, 0.0, 3.0, -6.0, 3.0, 0.0, -1.0, 3.0, -3.0, 1.0;
    return C / 6.0;
}

//! @brief B-Spline template
template<typename _Tp,int dim>
class Bspline
{
public:
    enum BsplinType {OPEN, CLOSED} type;

    Bspline() : type(OPEN)
    {
    }

    ~Bspline()
    {
        knot_pts.clear();
        ctrl_pts.clear();
    }

    void clear()
    {
        ctrl_pts.clear();
        knot_pts.clear();
    }

    bool is_ready()
    {
        return knot_pts.size() > 3;
    }

    void add_knot_pt(array<_Tp,dim> const& pt)
    {
        knot_pts.push_back(pt);
        cvt_knot_to_ctrl_cubic();
    }

    void add_knot_pts(vector<array<_Tp,dim> > const& pts)
    {
        knot_pts = pts;
        cvt_knot_to_ctrl_cubic();
    }

    void add_ctrl_pt(array<_Tp,dim> const& pt)
    {
        ctrl_pts.push_back(pt);
        cvt_ctrl_to_knot_cubic();
    }

    void add_ctrl_pts(vector<array<_Tp,dim> > const& pts)
    {
        ctrl_pts = pts;
        cvt_ctrl_to_knot_cubic();
    }

    array<_Tp,dim> get_ctrl_pt(size_t const p_idx) const
    {
        return ctrl_pts.at(p_idx);
    }

    array<_Tp,dim> get_ctrl_first_pt() const
    {
        return ctrl_pts.front();
    }

    array<_Tp,dim> get_ctrl_last_pt() const
    {
        return ctrl_pts.back();
    }

    array<_Tp,dim> get_knot_pt(size_t const p_idx) const
    {
        return knot_pts.at(p_idx);
    }

    array<_Tp,dim> get_knot_first_pt() const
    {
        return knot_pts.front();
    }

    array<_Tp,dim> get_knot_last_pt() const
    {
        return knot_pts.back();
    }

    size_t get_pt_idx(int const pt_idx) const
    {
        size_t num_ctrl_pts = ctrl_pts.size();

        if(type == CLOSED)
        {
            /* Cycling */
            if(pt_idx < 0)
                return num_ctrl_pts+pt_idx;

            if(pt_idx >= ctrl_pts.size())
                return pt_idx-num_ctrl_pts;
        }
        else
        {
            /* Truncating */
            if(pt_idx < 0)
                return 0;

            if(pt_idx >= ctrl_pts.size())
                return num_ctrl_pts-1;
        }

        return pt_idx;
    }

    size_t get_num_knot_pts() const { return knot_pts.size(); }
    size_t get_num_ctrl_pts() const { return ctrl_pts.size(); }

    void set_bspline_type(BsplinType const type)
    {
        this->type = type;
        cvt_knot_to_ctrl_cubic();
    }

    string get_bspline_type() const
    {
        switch(type)
        {
        case OPEN:
            return "Open B-spline";
        case CLOSED:
            return "Closed B-spline";
        }
    }

    array<_Tp,dim> cubic_intplt(int const pt_idx, _Tp const t, size_t const d_order = 0) const
    {

        if(ctrl_pts.size() < 4)
            return array<_Tp,dim>();

        if(t > 1.0 || t < 0)
        {
            cerr << "Knot is not normalised!" << endl;
            return array<_Tp,dim>();
        }

        /* Basis function coefficients */
        Vector4d B;
        if(d_order == 0)
        {
            Vector4d t_vec(1, t, t*t, t*t*t);
            B = t_vec.transpose()*CubicBsplineMatrix();
        }
        else if(d_order == 1)
        {
            Vector4d t_vec(0, 1, 2*t, 3*t*t);
            B = t_vec.transpose()*CubicBsplineMatrix();
        }
        else if(d_order == 2)
        {
            Vector4d t_vec(0, 0, 2, 6*t);
            B = t_vec.transpose()*CubicBsplineMatrix();
        }
        else if(d_order == 3)
        {
            Vector4d t_vec(0, 0, 0, 6);
            B = t_vec.transpose()*CubicBsplineMatrix();
        }

        array<_Tp,dim> pt;
        for(int i = 0; i < dim; ++i)
        {
            pt[i] = ctrl_pts.at(get_pt_idx(pt_idx-1))[i]*B[0] +
                    ctrl_pts.at(get_pt_idx(pt_idx))[i]*B[1] +
                    ctrl_pts.at(get_pt_idx(pt_idx+1))[i]*B[2] +
                    ctrl_pts.at(get_pt_idx(pt_idx+2))[i]*B[3];
        }

        return pt;
    }

private:
    void cvt_ctrl_to_knot_cubic()
    {
        size_t num_ctrl_pts = get_num_ctrl_pts();
        if(num_ctrl_pts > 3)
        {
            MatrixXd B = Eigen::MatrixXd::Zero(num_ctrl_pts, num_ctrl_pts);

            if(type == OPEN)
            {
                B(0, 0) = 1.0;
                B(num_ctrl_pts-1, num_ctrl_pts-1) = 1.0;

                for(int c = 1; c < num_ctrl_pts-1; ++c)
                {
                    B(c, (c-1)%num_ctrl_pts) = 1.0/6.0;
                    B(c, (c)%num_ctrl_pts) = 2.0/3.0;
                    B(c, (c+1)%num_ctrl_pts) = 1.0/6.0;
                }
            }

            if(type == CLOSED)
            {
                for(int c = 0; c < num_ctrl_pts; ++c)
                {
                    B(c, c%num_ctrl_pts) = 1.0/6.0;
                    B(c, (c+1)%num_ctrl_pts) = 2.0/3.0;
                    B(c, (c+2)%num_ctrl_pts) = 1.0/6.0;
                }
            }

            knot_pts.clear();

            array<_Tp,dim> pts[num_ctrl_pts];
            for(int d = 0; d < dim; ++d)
            {
                Eigen::VectorXd ctrl(num_ctrl_pts);
                for(int c = 0; c < num_ctrl_pts; ++c)
                    ctrl[c] = ctrl_pts[c][d];

                Eigen::VectorXd knot = B*ctrl;

                for(int k = 0; k < num_ctrl_pts; ++k)
                    pts[k][d] = knot(k);
            }

            for(int c = 0; c < num_ctrl_pts; ++c)
                knot_pts.push_back(pts[c]);
        }
    }


    void cvt_knot_to_ctrl_cubic()
    {
        size_t num_knot_pts = get_num_knot_pts();
        if(num_knot_pts > 3)
        {
            MatrixXd B = Eigen::MatrixXd::Zero(num_knot_pts, num_knot_pts);
            MatrixXd inv_B = Eigen::MatrixXd::Zero(num_knot_pts, num_knot_pts);

            if(type == OPEN)
            {
                B(0, 0) = 1.0;
                B(num_knot_pts-1, num_knot_pts-1) = 1.0;

                for(int i = 1; i < num_knot_pts-1; ++i)
                {
                    B(i, (i-1)%num_knot_pts) = 1.0/6.0;
                    B(i, (i)%num_knot_pts) = 2.0/3.0;
                    B(i, (i+1)%num_knot_pts) = 1.0/6.0;
                }

                inv_B = B.inverse();
            }

            if(type == CLOSED)
            {
                for(int i = 0; i < num_knot_pts; ++i)
                {
                    B(i, i%num_knot_pts) = 1.0/6.0;
                    B(i, (i+1)%num_knot_pts) = 2.0/3.0;
                    B(i, (i+2)%num_knot_pts) = 1.0/6.0;
                }

                inv_B = B.inverse();
            }

            ctrl_pts.clear();

            array<_Tp,dim> pts[num_knot_pts];
            for(int d = 0; d < dim; ++d)
            {
                Eigen::VectorXd knot(num_knot_pts);
                for(int k = 0; k < num_knot_pts; ++k)
                    knot[k] = knot_pts[k][d];

                Eigen::VectorXd ctrl = inv_B*knot;

                for(int k = 0; k < num_knot_pts; ++k)
                    pts[k][d] = ctrl(k);
            }

            for(int k = 0; k < num_knot_pts; ++k)
                ctrl_pts.push_back(pts[k]);

        }
    }

    /* Knot points */
    vector<array<_Tp,dim> > knot_pts;

    /* Control points */
    vector<array<_Tp,dim> > ctrl_pts;

};

#endif // BSPLIN_H

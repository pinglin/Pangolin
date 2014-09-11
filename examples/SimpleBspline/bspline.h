#ifndef BSPLIN_H
#define BSPLIN_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

//! @brief B-Spline template
template<typename _Tp,int dim>
class Bspline
{

public:
    enum BsplinType {OPEN, CLOSED} type;

    Bspline() : type(OPEN)
    {
        CubicBsplineMatrix << 1.0, 4.0, 1.0, 0.0, -3.0, 0.0, 3.0, 0.0, 3.0, -6.0, 3.0, 0.0, -1.0, 3.0, -3.0, 1.0;
        CubicBsplineMatrix /= 6.0;
    }

    ~Bspline()
    {
        knot_pts.conservativeResize(NoChange, 0);
        ctrl_pts.conservativeResize(NoChange, 0);
    }

    void clear()
    {
        knot_pts.conservativeResize(NoChange, 0);
        ctrl_pts.conservativeResize(NoChange, 0);
    }

    bool is_ready()
    {
        return knot_pts.size() > 3;
    }

    void add_knot_pt(Matrix<_Tp,dim,1> const& pt)
    {
        knot_pts.conservativeResize(NoChange, knot_pts.cols()+1);
        knot_pts.rightCols(1) = pt;
        cvt_knot_to_ctrl_cubic();
    }

    void add_knot_pts(Matrix<_Tp,dim,Dynamic> const& pts)
    {
        knot_pts = pts;
        cvt_knot_to_ctrl_cubic();
    }

    void add_ctrl_pt(Matrix<_Tp,dim,1> const& pt)
    {
        ctrl_pts.conservativeResize(NoChange, ctrl_pts.cols()+1);
        ctrl_pts.rightCols(1) = pt;
        cvt_ctrl_to_knot_cubic();
    }

    void add_ctrl_pts(Matrix<_Tp,dim,Dynamic> const& pts)
    {
        ctrl_pts = pts;
        cvt_ctrl_to_knot_cubic();
    }

    Array<_Tp,dim,1> get_knot_pt(size_t const p_idx) const
    {
        return knot_pts.col(p_idx);
    }

    Array<_Tp,dim,1> get_knot_first_pt() const
    {
        return knot_pts.leftCols(1);
    }

    Array<_Tp,dim,1> get_knot_last_pt() const
    {
        return knot_pts.rightCols(1);
    }

    Matrix<_Tp,dim,1> get_ctrl_pt(size_t const p_idx) const
    {
        return ctrl_pts.col(p_idx);
    }

    Array<_Tp,dim,1> get_ctrl_first_pt() const
    {
        return ctrl_pts.leftCols(1);
    }

    Array<_Tp,dim,1> get_ctrl_last_pt() const
    {
        return ctrl_pts.rightCols(1);
    }

    size_t get_pt_idx(int const pt_idx) const
    {
        size_t num_ctrl_pts = get_num_ctrl_pts();

        if(type == CLOSED)
        {
            /* Cycling */
            if(pt_idx < 0)
                return num_ctrl_pts+pt_idx;

            if(pt_idx >= num_ctrl_pts)
                return pt_idx-num_ctrl_pts;
        }
        else
        {
            /* Truncating */
            if(pt_idx < 0)
                return 0;

            if(pt_idx >= num_ctrl_pts)
                return num_ctrl_pts-1;
        }

        return pt_idx;
    }

    size_t get_num_knot_pts() const { return knot_pts.cols(); }
    size_t get_num_ctrl_pts() const { return ctrl_pts.cols(); }

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

    Matrix<_Tp,dim,1> cubic_intplt(int const pt_idx, _Tp const t, size_t const d_order = 0) const
    {

        if(get_num_ctrl_pts() < 4)
            return Matrix<_Tp,dim,1>();

        if(t > 1.0 || t < 0)
        {
            cerr << "Knot is not normalised!" << endl;
            return Matrix<_Tp,dim,1>();
        }

        /* Basis function coefficients */
        Matrix<_Tp,4,1> B;
        if(d_order == 0)
        {
            Matrix<_Tp,4,1> t_vec(1, t, t*t, t*t*t);
            B = t_vec.transpose()*CubicBsplineMatrix;
        }
        else if(d_order == 1)
        {
            Matrix<_Tp,4,1> t_vec(0, 1, 2*t, 3*t*t);
            B = t_vec.transpose()*CubicBsplineMatrix;
        }
        else if(d_order == 2)
        {
            Matrix<_Tp,4,1> t_vec(0, 0, 2, 6*t);
            B = t_vec.transpose()*CubicBsplineMatrix;
        }
        else if(d_order == 3)
        {
            Matrix<_Tp,4,1> t_vec(0, 0, 0, 6);
            B = t_vec.transpose()*CubicBsplineMatrix;
        }

        Matrix<_Tp,dim,1> pt;
        for(int i = 0; i < dim; ++i)
        {
            pt[i] = ctrl_pts(i, get_pt_idx(pt_idx-1))*B[0] +
                    ctrl_pts(i, get_pt_idx(pt_idx))*B[1] +
                    ctrl_pts(i, get_pt_idx(pt_idx+1))*B[2] +
                    ctrl_pts(i, get_pt_idx(pt_idx+2))*B[3];
        }

        return pt;
    }

private:
    void cvt_ctrl_to_knot_cubic()
    {
        size_t num_ctrl_pts = get_num_ctrl_pts();
        if(num_ctrl_pts > 3)
        {
            Matrix<_Tp,Dynamic,Dynamic> B = Matrix<_Tp,Dynamic,Dynamic>::Zero(num_ctrl_pts, num_ctrl_pts);

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

            knot_pts = (B*ctrl_pts.transpose()).transpose();

        }
    }


    void cvt_knot_to_ctrl_cubic()
    {
        size_t num_knot_pts = get_num_knot_pts();
        if(num_knot_pts > 3)
        {
            Matrix<_Tp,Dynamic,Dynamic> B = Matrix<_Tp,Dynamic,Dynamic>::Zero(num_knot_pts, num_knot_pts);
            Matrix<_Tp,Dynamic,Dynamic> inv_B;

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

            ctrl_pts = (inv_B*knot_pts.transpose()).transpose();

        }
    }

    /* Cubic Bspline basis matrix */
    Matrix<_Tp,4,4> CubicBsplineMatrix;

    /* Knot points */
    Matrix<_Tp,dim,Dynamic> knot_pts;

    /* Control points */
    Matrix<_Tp,dim,Dynamic> ctrl_pts;

};

#endif // BSPLIN_H

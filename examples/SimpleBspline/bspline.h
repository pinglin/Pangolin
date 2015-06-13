#include <iostream>
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

    Bspline() : lod(30), type(OPEN)
    {
        CubicBsplineMatrix << 1.0, 4.0, 1.0, 0.0, -3.0, 0.0, 3.0, 0.0, 3.0, -6.0, 3.0, 0.0, -1.0, 3.0, -3.0, 1.0;
        CubicBsplineMatrix /= 6.0;
    }

    ~Bspline()
    {
        knot_pts.conservativeResize(NoChange, 0);
        ctrl_pts.conservativeResize(NoChange, 0);
    }

    void Reset()
    {
        knot_pts.conservativeResize(NoChange, 0);
        ctrl_pts.conservativeResize(NoChange, 0);
    }

    bool IsReady() const
    {
        return knot_pts.cols() > 3 && ctrl_pts.cols() > 3;
    }

    void AddFrontKnotPt(Matrix<_Tp,dim,1> const& pt)
    {
        Matrix<_Tp,dim,Dynamic> pre_knot_pts = knot_pts;
        knot_pts.conservativeResize(NoChange, knot_pts.cols()+1);
        knot_pts.block(0, 1, knot_pts.rows(), knot_pts.cols()-1) = pre_knot_pts;
        knot_pts.leftCols(1) = pt;
        CvtKnotToCtrlCubic();
    }

    void AddFrontCtrlPt(Matrix<_Tp,dim,1> const& pt)
    {
        Matrix<_Tp,dim,Dynamic> pre_ctrl_pts = ctrl_pts;
        ctrl_pts.conservativeResize(NoChange, ctrl_pts.cols()+1);
        ctrl_pts.block(0, 1, ctrl_pts.rows(), ctrl_pts.cols()-1) = pre_ctrl_pts;
        ctrl_pts.leftCols(1) = pt;
        CvtCtrlToKnotCubic();
    }

    void AddBackKnotPt(Matrix<_Tp,dim,1> const& pt)
    {
        knot_pts.conservativeResize(NoChange, knot_pts.cols()+1);
        knot_pts.rightCols(1) = pt;
        CvtKnotToCtrlCubic();
    }

    void AddKnotPt(Matrix<_Tp,dim,1> const& pt, int const idx_insert)
    {
        Matrix<_Tp,dim,Dynamic> pre_knot_pts = knot_pts;
        knot_pts.conservativeResize(NoChange, knot_pts.cols()+1);

        int idx = idx_insert;
        if(idx < 1) { // insert in front
            knot_pts.block(0, 1, knot_pts.rows(), knot_pts.cols()-1) = pre_knot_pts;
            idx = 0; // make sure idx is valid
        } else if(idx > knot_pts.cols() - 1) { // insert in back
            idx = knot_pts.cols() - 1;
        } else // insert between
            knot_pts.block(0, idx+1, knot_pts.rows(), knot_pts.cols()-1-idx) = pre_knot_pts.block(0, idx, knot_pts.rows(), knot_pts.cols()-1-idx);

        knot_pts.col(idx) = pt;

        CvtKnotToCtrlCubic();
    }

    void AddBackKnotPts(Matrix<_Tp,dim,Dynamic> const& pts)
    {
        knot_pts = pts;
        CvtKnotToCtrlCubic();
    }

    void AddBackCtrlPt(Matrix<_Tp,dim,1> const& pt)
    {
        ctrl_pts.conservativeResize(NoChange, ctrl_pts.cols()+1);
        ctrl_pts.rightCols(1) = pt;
        CvtCtrlToKnotCubic();
    }

    void AddBackCtrlPts(Matrix<_Tp,dim,Dynamic> const& pts)
    {
        ctrl_pts = pts;
        CvtCtrlToKnotCubic();
    }

    void RemoveFrontKnotPt()
    {
        if(knot_pts.cols() > 0) {
            Matrix<_Tp,dim,Dynamic> pre_knot_pts = knot_pts.block(0, 1, knot_pts.rows(), knot_pts.cols()-1);
            knot_pts = pre_knot_pts;
            CvtKnotToCtrlCubic();
        }
    }

    void RemoveFrontCtrlPt()
    {
        if(ctrl_pts.cols() > 0) {
            Matrix<_Tp,dim,Dynamic> pre_ctrl_pts = ctrl_pts.block(0, 1, ctrl_pts.rows(), ctrl_pts.cols()-1);
            ctrl_pts = pre_ctrl_pts;
            CvtCtrlToKnotCubic();
        }
    }

    void RemoveBackKnotPt()
    {
        if(knot_pts.cols() > 0) {
            knot_pts.conservativeResize(NoChange, knot_pts.cols()-1);
            CvtKnotToCtrlCubic();
        }
    }

    void RemoveBackCtrlPt()
    {
        if(ctrl_pts.cols() > 0) {
            ctrl_pts.conservativeResize(NoChange, ctrl_pts.cols()-1);
            CvtCtrlToKnotCubic();
        }
    }

    void RemoveKnotPt(int const idx) {
        if(knot_pts.cols() > 0) {
            if(idx <= 0)
                RemoveFrontKnotPt();
            else if(idx >= knot_pts.cols() - 1)
                RemoveBackKnotPt();
            else {
                Matrix<_Tp,dim,Dynamic> pre_knot_pts = knot_pts;
                knot_pts.conservativeResize(NoChange, knot_pts.cols() - 1);
                knot_pts.block(0, idx, knot_pts.rows(), knot_pts.cols() - idx) = pre_knot_pts.block(0, idx+1, pre_knot_pts.rows(), pre_knot_pts.cols() - idx - 1);
            }

            CvtKnotToCtrlCubic();
        }
    }

    void SetKnotPt(size_t const p_idx, Matrix<_Tp,dim,Dynamic> const& pt)
    {
        knot_pts.col(GetPtIdx(p_idx)) = pt;
        CvtKnotToCtrlCubic();
    }

    void SetCtrlPt(size_t const p_idx, Matrix<_Tp,dim,Dynamic> const& pt)
    {
        ctrl_pts.col(GetPtIdx(p_idx)) = pt;
        CvtCtrlToKnotCubic();
    }

    Matrix<_Tp,dim,Dynamic> GetKnotPts() const
    {
        return knot_pts;
    }

    Matrix<_Tp,dim,1> GetKnotPt(size_t const p_idx) const
    {
        return knot_pts.col(GetPtIdx(p_idx));
    }

    Matrix<_Tp,dim,1> GetFrontKnotPt() const
    {
        return knot_pts.leftCols(1);
    }

    Matrix<_Tp,dim,1> GetBackKnotPt() const
    {
        return knot_pts.rightCols(1);
    }

    Matrix<_Tp,dim,Dynamic> GetCtrlPts() const
    {
        return ctrl_pts;
    }

    Matrix<_Tp,dim,1> GetCtrlPt(size_t const p_idx) const
    {
        return ctrl_pts.col(GetPtIdx(p_idx));
    }

    Matrix<_Tp,dim,1> GetFrontCtrlPt() const
    {
        return ctrl_pts.leftCols(1);
    }

    Matrix<_Tp,dim,1> GetBackCtrlPt() const
    {
        return ctrl_pts.rightCols(1);
    }

    size_t GetNumKnotPts() const { return knot_pts.cols(); }
    size_t GetNumCtrlPts() const { return ctrl_pts.cols(); }

    void SetBsplineType(BsplinType const type)
    {
        this->type = type;
        CvtKnotToCtrlCubic();
    }

    string GetBsplineType() const {
        switch(type)
        {
        case OPEN:
            return "Open B-spline";
        case CLOSED:
            return "Closed B-spline";
        }
    }


    Matrix<_Tp,dim,1> GetDerivative(int const knot_idx, int const order) {

        Matrix<_Tp,dim,1> deriv = Matrix<_Tp,dim,1>::Zero();

        if(IsReady()) {

            if(knot_idx < 0 || knot_idx >= GetNumKnotPts()) {
                cerr << "[Bspline] GetDerivative(knot_idx, order) Knot index out of range." << endl;
                return deriv;
            }

            if(order < 0 || order > 2) {
                cerr << "[Bspline] GetDerivative(knot_idx, order) Order out of range." << endl;
                return deriv;
            }

            int pt;
            if(knot_idx == 0)
                pt = knot_idx-1;
            else
                pt = knot_idx;

            deriv = CubicIntplt(pt, 0, order);
        }

        return deriv;
    }

    _Tp GetLength(int const k_i, int const k_j) {

        _Tp arc_length = 0.0;

        if(IsReady()) {

            if(k_i < 0 || k_j >= GetNumKnotPts()) {
                cerr << "[Bspline]: GetLength(k_i, k_j) out of range." << endl;
                return arc_length;
            }

            int p0, p1;
            if(k_i == 0 && k_j == 0 || k_i == GetNumKnotPts()-1 && k_j == GetNumKnotPts()-1) {
                p0 = k_i;
                p1 = k_j;
            } else if(k_i == 0 && k_j == GetNumKnotPts()-1) {
                p0 = k_i-1;
                p1 = k_j+1;
            } else if(k_i == 0) {
                p0 = k_i-1;
                p1 = k_j;
            } else if(k_j == GetNumKnotPts()-1) {
                p0 = k_i;
                p1 = k_j+1;
            } else {
                p0 = k_i;
                p1 = k_j;
            }

            for(int pt_idx = p0, seg_idx = p0+1; seg_idx <= p1; ++pt_idx, ++seg_idx) {
                for(int d = 0; d < lod; ++d) {
                    Matrix<_Tp,dim,1> pt1 = CubicIntplt(pt_idx, (d)/static_cast<_Tp>(lod));
                    Matrix<_Tp,dim,1> pt2 = CubicIntplt(pt_idx, (d+1)/static_cast<_Tp>(lod));

                    arc_length += (pt2 - pt1).norm();
                }
            }
        }

        return arc_length;
    }

    _Tp GetLength() {

        _Tp arc_length = 0.0;

        if(IsReady()) {
            for(int pt_idx = -1, seg_idx = 0; seg_idx <= GetNumCtrlPts(); ++pt_idx, ++seg_idx) {
                for(int d = 0; d < lod; ++d) {
                    Matrix<_Tp,dim,1> pt1 = CubicIntplt(pt_idx, (d)/static_cast<_Tp>(lod));
                    Matrix<_Tp,dim,1> pt2 = CubicIntplt(pt_idx, (d+1)/static_cast<_Tp>(lod));

                    arc_length += (pt2 - pt1).norm();
                }
            }
        }

        return arc_length;
    }

    _Tp GetAvgArcLength() {

        return GetLength() / _Tp(GetNumKnotPts()-1);;

    }

    // Applies a transformation matrix to knot points
    void TransformKnots(Matrix<_Tp,dim+1,dim+1> const& M)
    {
        for(size_t col = 0; col < knot_pts.cols(); ++col)
            knot_pts.col(col) = M.topLeftCorner(dim,dim) * knot_pts.col(col) + M.topRightCorner(dim,1);
        CvtKnotToCtrlCubic();
    }


    void KnotEquidist()
    {

        if(IsReady()) {

            bool has_converged = false;

            while(!has_converged) {

                _Tp const avg_arc_len = GetAvgArcLength();

                Matrix<_Tp,dim,Dynamic> new_knot_pts;
                new_knot_pts.conservativeResize(NoChange, GetNumKnotPts());
                new_knot_pts.leftCols(1) = knot_pts.leftCols(1);

                _Tp begin_arc_len = 0.0;
                _Tp end_arc_len = 0.0;
                _Tp knot_arc_len = avg_arc_len;

                size_t new_knot_idx = 1;
                for(int pt_idx = -1, seg_idx = 0; seg_idx <= GetNumKnotPts(); ++pt_idx, ++seg_idx) {

                    begin_arc_len = end_arc_len;

                    for(int d = 0; d < lod; ++d) {

                        Matrix<_Tp,dim,1> pt1 = CubicIntplt(pt_idx, _Tp(d)/_Tp(lod));
                        Matrix<_Tp,dim,1> pt2 = CubicIntplt(pt_idx, _Tp(d+1)/_Tp(lod));

                        end_arc_len += (pt2 - pt1).norm();
                    }

                    while(end_arc_len - knot_arc_len > 1e-1) {

                        _Tp t = _Tp(knot_arc_len - begin_arc_len) / _Tp(end_arc_len - begin_arc_len);

                        if(t > 1e-2)
                            new_knot_pts.col(new_knot_idx) = CubicIntplt(pt_idx, t);
                        else
                            new_knot_pts.col(new_knot_idx) = knot_pts.col(new_knot_idx);

                        knot_arc_len += avg_arc_len;
                        new_knot_idx++;
                    }
                }

                new_knot_pts.rightCols(1) = knot_pts.rightCols(1);

                if((new_knot_pts - knot_pts).norm() < 1e-1)
                    has_converged = true;

                knot_pts = new_knot_pts;
                CvtKnotToCtrlCubic();
            }
        }

    }


    /* 2D only */
    vector<Vector2i> GetIntegerBsplinePts() const {

        vector<Vector2i> continuous_pts;
        continuous_pts.reserve(1e4);

        if(IsReady()) {
            for(int pt_idx = -1, seg_idx = 0; seg_idx <= GetNumCtrlPts(); ++pt_idx, ++seg_idx) {

                for(int d = 0; d <= GetLOD(); ++d) {

                     Matrix<_Tp,dim,1> pt = CubicIntplt(pt_idx, d/float(GetLOD()));
                     Vector2i int_pt(pt[0], pt[1]);

                    if(continuous_pts.size() == 0)
                        continuous_pts.push_back(int_pt);
                    else if(continuous_pts.back() != int_pt) {
                        int x0 = continuous_pts.back()[0];
                        int y0 = continuous_pts.back()[1];
                        int x1 = int_pt[0];
                        int y1 = int_pt[1];

                        int d_x = x1 - x0;
                        int d_y = y1 - y0;

                        if(d_x != 0) {
                            for(int i = 1; i <= abs(d_x); ++i) {
                                int inter_x = x0 + (d_x < 0 ? -i : i);
                                int inter_y = y0 + round(d_y*(float(inter_x-x0)/float(d_x)));
                                Vector2i inter_pt(inter_x, inter_y);
                                if(continuous_pts.back() != inter_pt)
                                    continuous_pts.push_back(inter_pt);
                            }
                        }

                        if(d_y != 0) {
                            for(int i = 1; i <= abs(d_y); ++i) {
                                int inter_y = y0 + (d_y < 0 ? -i : i);
                                int inter_x = x0 + round(d_x*(float(inter_y-y0)/float(d_y)));
                                Vector2i inter_pt(inter_x, inter_y);
                                if(continuous_pts.back() != inter_pt)
                                    continuous_pts.push_back(inter_pt);
                            }
                        }

                    }
                }
            }
        }

        return continuous_pts;
    }

    Matrix<_Tp,dim,1> CubicIntplt(int const pt_idx, _Tp const t, size_t const d_order = 0) const {

        if(GetNumCtrlPts() < 4)
            return Matrix<_Tp,dim,1>();

        if(t > 1.0 || t < 0) {
            cerr << "Curve parameter is not normalised!" << endl;
            return Matrix<_Tp,dim,1>();
        }

        /* Basis function coefficients */
        Matrix<_Tp,4,1> B;
        if(d_order == 0) {
            Matrix<_Tp,4,1> t_vec(1, t, t*t, t*t*t);
            B = t_vec.transpose()*CubicBsplineMatrix;
        } else if(d_order == 1) {
            Matrix<_Tp,4,1> t_vec(0, 1, 2*t, 3*t*t);
            B = t_vec.transpose()*CubicBsplineMatrix;
        } else if(d_order == 2) {
            Matrix<_Tp,4,1> t_vec(0, 0, 2, 6*t);
            B = t_vec.transpose()*CubicBsplineMatrix;
        } else if(d_order == 3) {
            Matrix<_Tp,4,1> t_vec(0, 0, 0, 6);
            B = t_vec.transpose()*CubicBsplineMatrix;
        }

        Matrix<_Tp,dim,1> pt;
        for(int i = 0; i < dim; ++i) {
            pt[i] = ctrl_pts(i, GetPtIdx(pt_idx-1))*B[0] +
                    ctrl_pts(i, GetPtIdx(pt_idx))*B[1] +
                    ctrl_pts(i, GetPtIdx(pt_idx+1))*B[2] +
                    ctrl_pts(i, GetPtIdx(pt_idx+2))*B[3];
        }

        return pt;
    }

    void SetLOD(size_t const lod) { this->lod = lod; }
    size_t GetLOD() const { return lod; }

    size_t GetPtIdx(int const pt_idx) const
    {
        size_t num_ctrl_pts = GetNumCtrlPts();
        size_t num_knot_pts = GetNumKnotPts();

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

            if(pt_idx >= max(num_ctrl_pts, num_knot_pts))
                return max(num_ctrl_pts, num_knot_pts)-1;
        }

        return pt_idx;
    }

private:
    void CvtCtrlToKnotCubic()
    {
        size_t num_ctrl_pts = GetNumCtrlPts();

        if(num_ctrl_pts > 3)
        {

            knot_pts.conservativeResize(NoChange, num_ctrl_pts);
            knot_pts = Matrix<_Tp,dim,Dynamic>::Zero(dim, num_ctrl_pts);

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

    void CvtKnotToCtrlCubic()
    {
        size_t num_knot_pts = GetNumKnotPts();
        if(num_knot_pts > 3)
        {
            ctrl_pts.conservativeResize(NoChange, num_knot_pts);
            ctrl_pts = Matrix<_Tp,dim,Dynamic>::Zero(dim, num_knot_pts);

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

    _Tp Abs(_Tp const& x) {
        return (x >= 0) ? x : -x;
    }

    /* Cubic Bspline basis matrix */
    Matrix<_Tp,4,4> CubicBsplineMatrix;

    /* Knot points */
    Matrix<_Tp,dim,Dynamic> knot_pts;

    /* Control points */
    Matrix<_Tp,dim,Dynamic> ctrl_pts;

    /* Level of details */
    size_t lod;

};

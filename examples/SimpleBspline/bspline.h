#ifndef CASCADE_BSPLINE
#define CASCADE_BSPLINE

#include <boost/array.hpp>
#include <vector>

namespace cascade {

using namespace boost;
using namespace std;

enum BsplineType {OPEN, CLOSED};

struct Quadratic {};
struct Cubic {};

template<typename order,typename T,int dim>
class Bspline
{
private:
    Bspline();
};


//! @brief Cubic B-Spline
template<typename T,int dim>
class Bspline<Cubic,T,dim>
{

public:
    Bspline() : type(OPEN)
    {        
    }

    ~Bspline() { ctrl_pts.clear(); }

    void ClearCtrlPoints() { ctrl_pts.clear(); }

    void AddCtrlPoint(T const pt[dim])
    {
        ctrl_pts.push_back(boost::array<T,dim>());
        copy(pt, pt+dim, ctrl_pts.back().elems);
    }

    void AddCtrlPoint(boost::array<T,dim> const pt)
    {
        ctrl_pts.push_back(pt);
    }

    boost::array<T,dim> GetFirstCtrlPoint() const
    {
        return ctrl_pts.front();
    }

    boost::array<T,dim> GetLastCtrlPoint() const
    {
        return ctrl_pts.back();
    }

    boost::array<T,dim> GetCtrlPoint(int const p_idx) const
    {

        if(type == CLOSED)
        {
            size_t num_ctrl_pts = ctrl_pts.size();

            if(p_idx < 0)
                return ctrl_pts.at(num_ctrl_pts+p_idx);

            if(p_idx >= ctrl_pts.size())
                return ctrl_pts.at(p_idx-num_ctrl_pts);

            return ctrl_pts.at(p_idx);
        }

        if(p_idx < 0)
            return ctrl_pts.at(0);

        if(p_idx >= ctrl_pts.size())
            return ctrl_pts.back();

        return ctrl_pts.at(p_idx);

    }

    void SetCtrlPoint(int const p_idx, boost::array<T,dim> const pt)
    {
        if(p_idx < 0)
            ctrl_pts.at(0) = pt;

        if(p_idx >= ctrl_pts.size())
            ctrl_pts.back() = pt;

        ctrl_pts.at(p_idx) = pt;
    }

    size_t GetNumCtrlPoints() const { return ctrl_pts.size(); }

    boost::array<T,dim> Interpolate(int const p_idx, T const t) const
    {

        if(ctrl_pts.size() < 4)
        {
            cerr << "The number of the rest of control points is less than 4 for the point: " << p_idx << "!" << endl;
            return boost::array<T,dim>();
        }

        if(t > 1.0)
        {
            cerr << "Knot is not normalised!" << endl;
            return boost::array<T,dim>();
        }

        /* Basis function coefficients */
        T B[4];
        T it = 1.0 - t;
        B[0] = it*it*it/6.0;
        B[1] = (3*t*t*t - 6*t*t +4)/6.0;
        B[2] = (-3*t*t*t +3*t*t + 3*t + 1)/6.0;
        B[3] =  t*t*t/6.0;

        boost::array<T,dim> pt;
        for(int i = 0; i < dim; ++i)
            pt[i] = GetCtrlPoint(p_idx-1)[i]*B[0] +
                    GetCtrlPoint(p_idx)[i]*B[1] +
                    GetCtrlPoint(p_idx+1)[i]*B[2] +
                    GetCtrlPoint(p_idx+2)[i]*B[3];

        return pt;

    }

    int GetNumCtrlPts() const { return ctrl_pts.size(); }

    void SetBsplineType(BsplineType const type) { this->type = type; }

    string GetBsplineType() const
    {
        switch(type)
        {
        case OPEN:
            return "Open B-spline";
        case CLOSED:
            return "Closed B-spline";
        }
    }

private:

    BsplineType type;

    /* control points */
    vector<boost::array<T,dim> > ctrl_pts;

};

} // namespace cascade {

#endif // CASCADE_BSPLINE

#include "..\src\MeshViewer\MeshDefinition.h"
#include "scalar_def.h"
#include "constant_data.h"
#include <Eigen/Core>
namespace CABT
{
    class tet2
    {
    public:
        tet2() { }
        ~tet2() { }

    public:
        static tet2_constant_data constant_data;
        static subdivide_tree sub_tree;
        mat3_10 control_val;
        mat3_10 control_dir;
        void init(mat3_10& val, mat3_10& dir)
        {
            control_val = val; control_dir = dir;
        }
        void calc_sampling_jacobi_det(scalar t, mat20_4 &cof);
        bool search_cubic_root(scalar a, scalar b, scalar c, scalar d, scalar& time);
        bool tet_status(mat20_4& cof, scalar& time);
        void get_min_time(mat20_4& cof, scalar& time);
        void run(scalar& time)
        {
            mat20_4 cof;
            calc_sampling_jacobi_det(time, cof);
            get_min_time(cof, time);
        }
    };
    
}
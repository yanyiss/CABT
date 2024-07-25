#include "btet.h"
namespace CABT
{
	void tet2::calc_sampling_jacobi_det(scalar t, mat20_4& det)
	{
		mat3_20 knot;
		scalar inverse3 = 1.0 / 3.0;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4 - i; ++j)
			{
				for (int k = 0; k < 4 - i - j; ++k)
				{
					knot(0, o3tensor2vec(i, j, k)) = i * inverse3;
					knot(1, o3tensor2vec(i, j, k)) = j * inverse3;
					knot(2, o3tensor2vec(i, j, k)) = k * inverse3;
				}
			}
		}
		mat3_3 jacobi;
		for (int n = 0; n < 4; ++n)
		{
			mat3_10 val = control_val + n / 3.0 * t * control_dir;
			for (int i = 0; i < 4; ++i)
			{
				for (int j = 0; j < 4 - i; ++j)
				{
					for (int k = 0; k < 4 - i - j; ++k)
					{
						int id4ijk = o3tensor2vec(i, j, k);
						jacobi.col(0) = 2.0 * ((val.col(o2tensor2vec(1, 0, 0)) - val.col(o2tensor2vec(0, 0, 0))) * (1.0 - knot(0, id4ijk) - knot(1, id4ijk) - knot(2, id4ijk))
							+ (val.col(o2tensor2vec(2, 0, 0)) - val.col(o2tensor2vec(1, 0, 0))) * knot(0, id4ijk)
							+ (val.col(o2tensor2vec(1, 1, 0)) - val.col(o2tensor2vec(0, 1, 0))) * knot(1, id4ijk)
							+ (val.col(o2tensor2vec(1, 0, 1)) - val.col(o2tensor2vec(0, 0, 1))) * knot(2, id4ijk));
						jacobi.col(1) = 2.0 * ((val.col(o2tensor2vec(0, 1, 0)) - val.col(o2tensor2vec(0, 0, 0))) * (1.0 - knot(0, id4ijk) - knot(1, id4ijk) - knot(2, id4ijk))
							+ (val.col(o2tensor2vec(1, 1, 0)) - val.col(o2tensor2vec(1, 0, 0))) * knot(0, id4ijk)
							+ (val.col(o2tensor2vec(0, 2, 0)) - val.col(o2tensor2vec(0, 1, 0))) * knot(1, id4ijk)
							+ (val.col(o2tensor2vec(0, 1, 1)) - val.col(o2tensor2vec(0, 0, 1))) * knot(2, id4ijk));
						jacobi.col(2) = 2.0 * ((val.col(o2tensor2vec(0, 0, 1)) - val.col(o2tensor2vec(0, 0, 0))) * (1.0 - knot(0, id4ijk) - knot(1, id4ijk) - knot(2, id4ijk))
							+ (val.col(o2tensor2vec(1, 0, 1)) - val.col(o2tensor2vec(1, 0, 0))) * knot(0, id4ijk)
							+ (val.col(o2tensor2vec(0, 1, 1)) - val.col(o2tensor2vec(0, 1, 0))) * knot(1, id4ijk)
							+ (val.col(o2tensor2vec(0, 0, 2)) - val.col(o2tensor2vec(0, 0, 1))) * knot(2, id4ijk));
						det(id4ijk, n) = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1)
							- jacobi(2, 0) * jacobi(1, 1) * jacobi(0, 2) - jacobi(1, 0) * jacobi(0, 1) * jacobi(2, 2) - jacobi(0, 0) * jacobi(2, 1) * jacobi(1, 2);
					}
				}
			}
		}
		for (int i = 0; i < 20; ++i)
		{
			det(i, 1) = -5.5 * det(i, 0) + 9 * det(i, 1) - 4.5 * det(i, 2) + det(i, 3);
			det(i, 2) = 9 * det(i, 0) - 22.5 * det(i, 1) + 18 * det(i, 2) - 4.5 * det(i, 3);
			det(i, 3) = -4.5 * det(i, 0) + 13.5 * det(i, 1) - 13.5 * det(i, 2) + 4.5 * det(i, 3);
		}
	}

	void tet2::calc_bezier_cof(mat20_4& cof)
	{
		cof = transform * cof;
	}
}
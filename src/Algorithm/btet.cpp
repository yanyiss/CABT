#include "btet.h"
namespace CABT
{
	void tet2::calc_sampling_jacobi_det(scalar t, mat20_4& cof)
	{
		auto& o3tensor2vec = constant_data.o3tensor2vec;
		auto& o2tensor2vec = constant_data.o2tensor2vec;
		mat3_20 knot;
		scalar inverse3 = scalar(1) / scalar(3);
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4 - i; ++j)
			{
				for (int k = 0; k < 4 - i - j; ++k)
				{
					knot(0, o3tensor2vec(i, j, k)) = scalar(i) * inverse3;
					knot(1, o3tensor2vec(i, j, k)) = scalar(j) * inverse3;
					knot(2, o3tensor2vec(i, j, k)) = scalar(k) * inverse3;
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
						cof(id4ijk, n) = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1)
							- jacobi(2, 0) * jacobi(1, 1) * jacobi(0, 2) - jacobi(1, 0) * jacobi(0, 1) * jacobi(2, 2) - jacobi(0, 0) * jacobi(2, 1) * jacobi(1, 2);
					}
				}
			}
		}
		for (int i = 0; i < 20; ++i)
		{
			cof(i, 1) = -5.5 * cof(i, 0) + 9 * cof(i, 1) - 4.5 * cof(i, 2) + cof(i, 3);
			cof(i, 2) = 9 * cof(i, 0) - 22.5 * cof(i, 1) + 18 * cof(i, 2) - 4.5 * cof(i, 3);
			cof(i, 3) = -4.5 * cof(i, 0) + 13.5 * cof(i, 1) - 13.5 * cof(i, 2) + 4.5 * cof(i, 3);
			cof(i, 1) = scalar(-9) / scalar(2) * cof(i, 0) + scalar(9) * cof(i, 1) - scalar(9) / scalar(2) * cof(i, 2) + cof(i, 3);
			cof(i, 2) = scalar(9) * cof(i, 0) - scalar(45) / scalar(2) * cof(i, 1) + scalar(18) * cof(i, 2) - scalar(9) / scalar(2) * cof(i, 3);
			cof(i, 3) = scalar(-9) / scalar(2) * cof(i, 0) + scalar(27) / scalar(2) * cof(i, 1) - scalar(27) / scalar(2) * cof(i, 2) + scalar(9) / scalar(2) * cof(i, 3);
		}
	}


	bool tet2::tet_status(mat20_4& cof, scalar &time)
	{
		//if no root in (0, time) return false;
		//else set time = the min root t, and return true;

	}

	void tet2::get_min_time(mat20_4 &cof, scalar &time)
	{
		sub_tree.tree[0][0] = cof;
		if (!tet_status(cof, time))
			return;
		int index[4] = { 0,8,64,256 };
		for (int i = 1; i < 4; ++i)
		{
			bool flag = false;
			for (int j = 0; j < index[i]; ++j)
			{
				sub_tree.tree[i][j] = constant_data.son_transform[j % 8] * sub_tree.tree[i - 1][j / 8];
				flag = flag || tet_status(sub_tree.tree[i][j], time);
			}
			if (!flag)
			{
				break;
			}
		}
	}
}
//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_PBC_Boundary_normal.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)

TFDTD_PBC_Boundary_normal::TFDTD_PBC_Boundary_normal()
{
	isFDTD = 0;
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1 / sqrt(eps_0 * mu_0);

	_pi = acos(-1);
}

TFDTD_PBC_Boundary_normal::~TFDTD_PBC_Boundary_normal() {}

void TFDTD_PBC_Boundary_normal::FDTD_PBC_Ex_cal()
{
	for (i = 0; i < txsize; i++)
		for (j = 1; j < tysize; j++)
			for (k = 1; k < tzsize; k++)
			{
				mx = Media_Ex[i][j][k];
				Ex[i][j][k] = C1[mx] * Ex[i][j][k]
					+ Cb[mx] / E_dy[j] * (Hz[i][j][k] - Hz[i][j - 1][k])
					- Cb[mx] / E_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);
			}

	// y-axis PBC
	//j==0 PBC
	j = 0;
	for (i = 0; i < txsize; i++)
		for (k = 1; k < tzsize; k++)
		{
			mx = Media_Ex[i][j][k];
			Ex[i][j][k] = Ex[i][j][k]
				+ Cb[mx] / E_dy[j] * (Hz[i][j][k] - Hz[i][tysize - 1][k])
				- Cb[mx] / E_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);

		}

	//j=tysize;
	for (i = 0; i < txsize; i++)
		for (k = 1; k < tzsize; k++)
		{
			Ex[i][tysize][k] = Ex[i][0][k];
		}
}

void TFDTD_PBC_Boundary_normal::FDTD_PBC_Ey_cal()
{
	for (i = 1; i < txsize; i++)
		for (j = 0; j < tysize; j++)
			for (k = 1; k < tzsize; k++)
			{
				my = Media_Ey[i][j][k];
				Ey[i][j][k] = Ey[i][j][k]
					+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
					- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[i - 1][j][k]);
			}

	//i==0 PBC
	i = 0;
	for (j = 0; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			my = Media_Ey[i][j][k];
			Ey[i][j][k] = Ey[i][j][k]
				+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
				- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[txsize - 1][j][k]);
		}


	//	i=txsize;
	for (j = 0; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			Ey[txsize][j][k] = Ey[0][j][k]; //*(*(*(*Dy + i) + j) + k)//
		}

}

void TFDTD_PBC_Boundary_normal::FDTD_PBC_Ez_cal()
{
	for (i = 1; i < txsize; i++)
		for (j = 1; j < tysize; j++)
			for (k = 1; k < tzsize; k++)
			{
				mz = Media_Ez[i][j][k];
				Ez[i][j][k] = Ez[i][j][k]
					+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
					- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j - 1][k]);
			}

	//i==0 PBC
	i = 0;
	for (j = 1; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			mz = Media_Ez[i][j][k];
			Ez[i][j][k] = Ez[i][j][k]
				+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 1][j][k])
				- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j - 1][k]);
		}

	// j==0 PBC
	j = 0;
	for (i = 1; i < txsize; i++)
		for (k = 1; k < tzsize; k++)
		{
			mz = Media_Ez[i][j][k];
			Ez[i][j][k] = Ez[i][j][k]
				+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
				- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][tysize - 1][k]);
		}

	i = 0;
	j = 0;
	for (k = 1; k < tzsize; k++)
	{
		mz = Media_Ez[i][j][k];
		Ez[i][j][k] = Ez[i][j][k]
			+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 1][j][k])
			- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][tysize - 1][k]);
	}

	// i=txsize;
	for (j = 0; j < tysize + 1; j++)
		for (k = 0; k < tzsize + 1; k++)
			Ez[txsize][j][k] = Ez[0][j][k];

	// j=tysize;
	for (i = 0; i < txsize + 1; i++)
		for (k = 0; k < tzsize + 1; k++)
			Ez[i][tysize][k] = Ez[i][0][k];

}

void TFDTD_PBC_Boundary_normal::FDTD_PBC_H_cal()
{
	for (i = 0; i < txsize; i++)
		for (j = 0; j < tysize; j++)
			for (k = 0; k < tzsize; k++)
			{
				Hx[i][j][k] = Hx[i][j][k]
					- Db / H_dy[j] * (Ez[i][j + 1][k] - Ez[i][j][k])
					+ Db / H_dz[k] * (Ey[i][j][k + 1] - Ey[i][j][k]);

				Hy[i][j][k] = Hy[i][j][k]
					- Db / H_dz[k] * (Ex[i][j][k + 1] - Ex[i][j][k])
					+ Db / H_dx[i] * (Ez[i + 1][j][k] - Ez[i][j][k]);

				Hz[i][j][k] = Hz[i][j][k]
					- Db / H_dx[i] * (Ey[i + 1][j][k] - Ey[i][j][k])
					+ Db / H_dy[j] * (Ex[i][j + 1][k] - Ex[i][j][k]);
			}
}

//void TFDTD_PBC_Boundary_normal::FDTD_PBC_D_Chiral_Cal()
//{
//    for(i=0;i<txsize;i++)
//	for(j=1;j<tysize;j++)
//	for(k=1;k<tzsize;k++)
//    {
//        *(*(*(*Dx+i)+j)+k) = *(*(*(*Dx + i) + j) + k)
//                                + Cb_dy * ( *(*(*(*Hz + i) + j) + k) - *(*(*(*Hz + i) + j-1) + k))
//                                - Cb_dz * ( *(*(*(*Hy + i) + j) + k) - *(*(*(*Hy + i) + j) + k-1));
//    }
//
//    for(i=1;i<txsize;i++)
//    for(j=0;j<tysize;j++)
//    for(k=1;k<tzsize;k++)
//    {
//        *(*(*(*Dy + i) + j) + k) = *(*(*(*Dy + i) + j) + k)
//                                 + Cb_dz * (*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx + i) + j) + k-1))
//                                 - Cb_dx * (*(*(*(*Hz + i) + j) + k) - *(*(*(*Hz + i-1) + j) + k));
//    }
//
//	for(i=1;i<txsize;i++)
//	for(j=1;j<tysize;j++)
//	for(k=1;k<tzsize;k++)
//    {
//        *(*(*(*Dz + i) + j) + k) = *(*(*(*Dz + i) + j) + k)
//                                + Cb_dx*(*(*(*(*Hy + i) + j) + k) - *(*(*(*Hy + i - 1) + j) + k))
//                                - Cb_dy*(*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx + i) + j - 1) + k));
//    }
//
//
//}

//
//// Old
//// x = 1:56, y = 0, X = 1:56, Y = 134
//void TFDTD_PBC_Boundary_normal::FDTD_PBC_Kiral_BC1()
//{
//    j = 0;
//    for(i=0;i<txsize;i++) {
//        for(k=0;k<tzsize;k++) {
//        *(*(*(*Dx+i)+j)+k) = *(*(*(*Dx + i) + j) + k)
//                            + Cb_dy * ( *(*(*(*Hz + i) + j) + k) - *(*(*(*Hz + i) + tysize-1) + k))
//                            - Cb_dz * ( *(*(*(*Hy + i) + j) + k) - *(*(*(*Hy + i) + j       ) + k-1));
//    }
//	}
//
//	for(i=1;i<txsize;i++) {
//		for(k=0;k<tzsize;k++) {
//	    	*(*(*(*Dz + i) + j) + k) = *(*(*(*Dz + i) + j) + k)
//                            + Cb_dx*(*(*(*(*Hy + i) + j) + k) - *(*(*(*Hy + i - 1) + j) + k))
//                            - Cb_dy*(*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx + i) + tysize - 1) + k));
//		}
//	}
//
//    j=tysize;
//    for (i = 0; i < txsize; i++) {
//    for (k = 0; k < tzsize; k++) {
//        *(*(*(*Dx+i)+j)+k) = *(*(*(*Dx + i) + 0) + k);
//		}
//	}
//
//    j=tysize;
//	for (i=0;i<txsize+1;i++) {
//		for (k=0;k<tzsize;k++) {
//
//			*(*(*(*Dz + i) + j) + k) =*(*(*(*Dz + i) + 0) + k); //*(*(*(*Dx + i) + j) + k)//
//        }
//	}
//}
////
//
//void TFDTD_PBC_Boundary_normal::FDTD_PBC_Kiral_BC2()
//{
//    int j_offset = 0;
//     //i==0 PBC
//    i = 0;
//	for (j = 0; j<aa1m1; j++) {
//        j_offset = j+aa1m2;
//		for (k = 0; k<tzsize; k++) {
//			*(*(*(*Dy + i) + j) + k) = *(*(*(*Dy + i) + j) + k)
//                                        + Cb_dz * (*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx + i) + j) + k-1))
//                                        - Cb_dx * (*(*(*(*Hz + i) + j) + k) - *(*(*(*Hz + txsize-1) + j_offset) + k));
//		}
//	}
//    for (j = aa1m1; j<tysize; j++){
//        j_offset = j - aa1m1;
//		for (k = 0; k<tzsize; k++) {
//			*(*(*(*Dy + i) + j) + k) = *(*(*(*Dy + i) + j) + k)
//                                        + Cb_dz * (*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx + i) + j) + k-1))
//                                        - Cb_dx * (*(*(*(*Hz + i) + j) + k) - *(*(*(*Hz + txsize-1) + j_offset) + k));
//		}
//	}
//    for (j = 1; j < aa1m1; j++) {
//        j_offset = j + aa1m2;
//        for (k = 0; k<tzsize; k++) {
//            *(*(*(*Dz+ i) + j) + k) = *(*(*(*Dz+ i) + j) + k)
//                                        + Cb_dx*(*(*(*(*Hy + i) + j) + k) - *(*(*(*Hy + txsize-1) + j_offset) + k))
//                                        - Cb_dy*(*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx + i       ) + j - 1    ) + k));
//        }
//        j_offset++;
//    }
//    for (j = aa1m1; j<tysize; j++){
//        j_offset = j - aa1m1;
//        for (k = 0; k<tzsize; k++) {
//            *(*(*(*Dz+ i) + j) + k) = *(*(*(*Dz+ i) + j) + k)
//                                        + Cb_dx*(*(*(*(*Hy + i) + j) + k) - *(*(*(*Hy + txsize-1) + j_offset) + k))
//                                        - Cb_dy*(*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx + i       ) + j - 1    ) + k));
//        }
//    }
//
//    i = 0;
//	j = 0;
//	for (k = 0; k<tzsize; k++)
//	{
//		*(*(*(*Dz + i) + j) + k) = *(*(*(*Dz + i) + j) + k)
//                    + Cb_dx*(*(*(*(*Hy + i) + j) + k) - *(*(*(*Hy+txsize - 1)+aa1m2)+k))
//					- Cb_dy*(*(*(*(*Hx + i) + j) + k) - *(*(*(*Hx+txsize )+aa1m2-1)+k));
//	}
//
//	//i==txsize PBC
//    i = txsize;
//	for (j = 0; j<aa1m2; j++) {
//        j_offset = j + aa1m1;
//		for (k = 0; k<tzsize; k++) {
//			*(*(*(*Dy + i) + j) + k) = *(*(*(*Dy + 0) + j_offset) + k);
//		}
//	}
//    for (j = aa1m2; j<tysize; j++){
//        j_offset = j - aa1m2;
//		for (k = 0; k<tzsize; k++) {
//			*(*(*(*Dy + i) + j) + k) = *(*(*(*Dy + 0) + j_offset) + k);
//		}
//	}
//    for (j = 0; j<aa1m2-1; j++) {
//        j_offset = j + aa1m1;
//        for (k = 1; k<tzsize; k++) {
//            *(*(*(*Dz + i) + j) + k) = *(*(*(*Dz + 0) + j_offset) + k);
//        }
//    }
//    for (j = aa1m2; j<tysize-1; j++){
//        j_offset = j - aa1m2;
//        for (k = 1; k<tzsize; k++) {
//            *(*(*(*Dz + i) + j) + k) = *(*(*(*Dz + 0) + j_offset) + k);
//        }
//    }
//}

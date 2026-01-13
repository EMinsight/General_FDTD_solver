//---------------------------------------------------------------------------
#pragma hdrstop
//#pragma boost
#include "TFDTDparameters.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

__fastcall TFDTDparameters::TFDTDparameters()
{
	m_pi = acos(-1.0);
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1 / sqrt(eps_0*mu_0);
    strcpy(JsonFileName, "./IEM_solver_input_format.json");
}

TFDTDparameters::~TFDTDparameters()
{
    if(isFDTD==1)
    {
        memfree(&Media_Ex, txsize, tysize + 1);
        memfree(&Media_Ey, txsize + 1, tysize);
        memfree(&Media_Ez, txsize + 1, tysize + 1);

        memfree(&Media_Hx, txsize + 1, tysize);
        memfree(&Media_Hy, txsize, tysize + 1);
        memfree(&Media_Hz, txsize, tysize);
    }
}

void TFDTDparameters::Set_Kappa_coefficient(double alphaT, double kappaT)
{
    this->kappa_max = (double)kappaT;
    m_kappa = 4;
    printf("kappa set\n");
}

void TFDTDparameters::Initial_Setting()
{
    Read_Media_Data();
    dt = 0.99 / c_0 / sqrt((1 / min1D(dx)) * (1 / min1D(dx)) + (1 / min1D(dy)) * (1 / min1D(dy)) + (1 / min1D(dz)) * (1 / min1D(dz)));
    Read_Material_Data();
    Set_Kappa_coefficient(0, 1);
    Spatial_Size_With_Kappa();
    
    Initiate_MediaSpace();
    Modeling_YeeCell();
}
void TFDTDparameters::Initial_Setting(struct_Field_carrier *pFc)
{
    PMLX = pFc->PML[0];
    PMLY = pFc->PML[1];
    PMLZ = pFc->PML[2];


    Read_Media_Data();
    dt = 0.99 / c_0 / sqrt((1 / min1D(dx)) * (1 / min1D(dx)) + (1 / min1D(dy)) * (1 / min1D(dy)) + (1 / min1D(dz)) * (1 / min1D(dz)));
    Read_Material_Data();
    Set_Kappa_coefficient(0, 1);
    Spatial_Size_With_Kappa();

    Initiate_MediaSpace();
    Modeling_YeeCell();

    pFc->txsize = this->txsize;
    pFc->tysize = this->tysize;
    pFc->tzsize = this->tzsize;

    pFc->media = this->Media;
    pFc->Media_Ex = this->Media_Ex;
    pFc->Media_Ey = this->Media_Ey;
    pFc->Media_Ez = this->Media_Ez;

    pFc->dt = this->dt;
    pFc->dx = this->dx;
    pFc->dy = this->dy;
    pFc->dz = this->dz;

    pFc->E_dx = this->E_dx;
    pFc->E_dy = this->E_dy;
    pFc->E_dz = this->E_dz;

    pFc->H_dx = this->H_dx;
    pFc->H_dy = this->H_dy;
    pFc->H_dz = this->H_dz;

    pFc->total_pol = this->total_pol;
}
void TFDTDparameters::set_PMLXYZ(int PMLX, int PMLY, int PMLZ)
{
    this->PMLX = PMLX;
    this->PMLY = PMLY;
    this->PMLZ = PMLZ;
}

void TFDTDparameters::Read_Media_Data()
{
    ifstream ifs(JsonFileName);
    IStreamWrapper isw(ifs);
    Jdocument.ParseStream(isw);;
    double Xlen, Ylen, Zlen;
    double dxT, dyT, dzT;
    double temp1, temp2;
    txsize = Jdocument["MESH"]["Voxelsize"]["txsize"].GetDouble();
    tysize = Jdocument["MESH"]["Voxelsize"]["tysize"].GetDouble();
    tzsize = Jdocument["MESH"]["Voxelsize"]["tzsize"].GetDouble();
	printf("Mesh size: %d x %d x %d\n", txsize, tysize, tzsize);

    dx = Tptr1<double>(txsize + 1); // �ּҰ� ���ϴ� �Լ� ������ ���� �迪�� �� ǥ�⸦ ���� +1�� ����.
    dy = Tptr1<double>(tysize + 1); // �ּҰ� ���ϴ� �Լ� ������ ���� �迪�� �� ǥ�⸦ ���� +1�� ����.
    dz = Tptr1<double>(tzsize + 1); // �ּҰ� ���ϴ� �Լ� ������ ���� �迪�� �� ǥ�⸦ ���� +1�� ����.
    dx[txsize] = -1;
    dy[tysize] = -1;
    dz[tzsize] = -1;
    E_dx = Tptr1<double>(txsize);
    E_dy = Tptr1<double>(tysize);
    E_dz = Tptr1<double>(tzsize);
    H_dx = Tptr1<double>(txsize);
    H_dy = Tptr1<double>(tysize);
    H_dz = Tptr1<double>(tzsize);

    for (int ii = 0; ii < txsize; ii++)
    {
        temp1 = Jdocument["MESH"]["CellLocation"]["xcoordinate"][ii].GetDouble();
        temp2 = Jdocument["MESH"]["CellLocation"]["xcoordinate"][ii + 1].GetDouble();
        dx[ii] = temp2 - temp1;
    }
    for (int jj = 0; jj < tysize; jj++)
    {
        temp1 = Jdocument["MESH"]["CellLocation"]["ycoordinate"][jj].GetDouble();
        temp2 = Jdocument["MESH"]["CellLocation"]["ycoordinate"][jj + 1].GetDouble();
        dy[jj] = temp2 - temp1;
    }
    for (int kk = 0; kk < tzsize; kk++)
    {
        temp1 = Jdocument["MESH"]["CellLocation"]["zcoordinate"][kk].GetDouble();
        temp2 = Jdocument["MESH"]["CellLocation"]["zcoordinate"][kk + 1].GetDouble();
        dz[kk] = temp2 - temp1;
    }
    Media = Tptr3<int>(txsize, tysize, tzsize);
    strcpy(MediaFileLocation, Jdocument["MESH"]["Voxel"].GetString());
    FILE* pVmesh; 
    pVmesh = fopen(MediaFileLocation, "rb");
    for (int i = 0; i < txsize; i++)
        for (int j = 0; j < tysize; j++)
            fread(Media[i][j], sizeof(int), tzsize, pVmesh);
    ifs.close();
    fclose(pVmesh);
}

void TFDTDparameters::Read_Material_Data()
{
    ifstream ifs(JsonFileName);
    IStreamWrapper isw(ifs);
    Jdocument.ParseStream(isw);;
    num_material = Jdocument["Material"].Size();
    pMp = Tptr1<Material_properties>(num_material);
    // find maximum pol
    for (int mm = 0; mm < num_material; mm++)
    {
        if (total_pol < Jdocument["Material"][mm]["num_pol"].GetDouble())
            total_pol = Jdocument["Material"][mm]["num_pol"].GetDouble();
    }
    C1 = Tptr1<double>(num_material);
    C2 = Tptr1<double>(num_material);
    eps_r_inf = Tptr1<double>(num_material);
    Cb = Tptr1<double>(num_material);
    if (total_pol != 0)
    {
        Ca_temp = Tptr2<double>(num_material, total_pol);
        Cb_temp = Tptr2<double>(num_material, total_pol);
    }
    for (int mm = 0; mm < num_material; mm++)
    {
        strcpy(pMp[mm].Name, Jdocument["Material"][mm]["Name"].GetString());
        if (strcmp(Jdocument["Material"][mm]["Type_of_material"].GetString(), "Dispersion") == 0)
        {
            double temp1 = 0, temp2 = 0, temp3 = 0;
            pMp[mm].TypeOfMaterials = 'D';

            pMp[mm].Ca = Tptr1<double>(total_pol);
            pMp[mm].Cb = Tptr1<double>(total_pol);
            pMp[mm].Cc = Tptr1<double>(total_pol);
            pMp[mm].Cd = Tptr1<double>(total_pol);
            pMp[mm].Ce = Tptr1<double>(total_pol);

            for (int pp = 0; pp < total_pol; pp++)
            {
                pMp[mm].Ca[pp] = Jdocument["Material"][mm]["Ca"][pp].GetDouble();
                pMp[mm].Cb[pp] = Jdocument["Material"][mm]["Cb"][pp].GetDouble();
                pMp[mm].Cc[pp] = Jdocument["Material"][mm]["Cc"][pp].GetDouble();
                pMp[mm].Cd[pp] = Jdocument["Material"][mm]["Cd"][pp].GetDouble();
                pMp[mm].Ce[pp] = Jdocument["Material"][mm]["Ce"][pp].GetDouble();
                pMp[mm].Epsr = Jdocument["Material"][mm]["permittivity_inf"].GetDouble();
                eps_r_inf[mm] = eps_0 * pMp[mm].Epsr;
                temp1 += pMp[mm].Cc[pp];
                temp2 += pMp[mm].Cd[pp];
                temp3 += pMp[mm].Ce[pp];
            }
            for (int pp = 0; pp < total_pol; pp++)
            {
                Ca_temp[mm][pp] = (pMp[mm].Ca[pp] + 1.) * 0.5 / (eps_r_inf[mm] / dt + temp1 * 0.5);
                Cb_temp[mm][pp] = (pMp[mm].Cb[pp]) * 0.5 / (eps_r_inf[mm] / dt + temp1 * 0.5);
            }
            C1[mm] = (eps_r_inf[mm] / dt - temp2 * 0.5) / (eps_r_inf[mm] / dt + temp1 * 0.5); // Ep
            C2[mm] = temp3 * 0.5 / (eps_r_inf[mm] / dt + temp1 * 0.5); // Epp
            Cb[mm] = 1. / (eps_r_inf[mm] / dt + temp1 * 0.5);// Curl H
            Db = dt / mu_0;
        }
        else if (strcmp(Jdocument["Material"][mm]["Type_of_material"].GetString(), "Nondispersion") == 0)
        {
            pMp[mm].TypeOfMaterials = 'N';
            pMp[mm].Epsr = Jdocument["Material"][mm]["permittivity"].GetDouble();
            pMp[mm].Sigma = Jdocument["Material"][mm]["conductivity"].GetDouble();
            eps_r_inf[mm] = eps_0 * pMp[mm].Epsr;
            //pMp[mm].loss_tan = (pMp[mm].Sigma / (pMp[mm].Epsr * eps_0 * 2.0 * m_pi * f0));
            C1[mm] = (((2.0 * eps_r_inf[mm]) - (pMp[mm].Sigma * dt)) /
                ((2.0 * eps_r_inf[mm]) + (pMp[mm].Sigma * dt))); //Ep
            C2[mm] = 0; // Epp
            Cb[mm] = ((2.0 * dt) / ((2.0 * eps_r_inf[mm]) + (pMp[mm].Sigma * dt))); // Curl H
            Db = dt / mu_0;
        }
    }
    ifs.close();
	printf("Material properties read completed.\n");
}

void TFDTDparameters::Initiate_MediaSpace()
{
	Media_Ex = Tptr3<int>(txsize, tysize + 1, tzsize + 1);
	Media_Ey = Tptr3<int>(txsize + 1, tysize, tzsize + 1);
	Media_Ez = Tptr3<int>(txsize + 1, tysize + 1, tzsize);

	Media_Hx = Tptr3<int>(txsize + 1, tysize, tzsize);
	Media_Hy = Tptr3<int>(txsize, tysize + 1, tzsize);
	Media_Hz = Tptr3<int>(txsize, tysize, tzsize + 1);
}

void TFDTDparameters::Initiate_Fields()
{
    
}

void TFDTDparameters::H_field_update()
{

#pragma omp parallel sections
    {
#pragma omp section
        Hx_field_update();

#pragma omp section
        Hy_field_update();

#pragma omp section
        Hz_field_update();
    }
//    int i = 0, j = 0, k = 0;
//#pragma omp parallel for collapse(3) schedule(dynamic) private(i,j,k)
//	for (i = 0; i<txsize; i++)
//        for (j = 0; j<tysize; j++)
//            for (k = 0; k<tzsize; k++)
//            {
//
//				Hx[i][j][k] = Hx[i][j][k] - Db/H_dy[j] * (Ez[i][j + 1][k] - Ez[i][j][k])
//										  + Db/H_dz[k] * (Ey[i][j][k + 1] - Ey[i][j][k]);
//				Hy[i][j][k] = Hy[i][j][k] - Db/H_dz[k] * (Ex[i][j][k + 1] - Ex[i][j][k])
//										  + Db/H_dx[i] * (Ez[i + 1][j][k] - Ez[i][j][k]);
//				Hz[i][j][k] = Hz[i][j][k] - Db/H_dx[i] * (Ey[i + 1][j][k] - Ey[i][j][k])
//										  + Db/H_dy[j] * (Ex[i][j + 1][k] - Ex[i][j][k]);
//
//            }
}
void TFDTDparameters::Hx_field_update()
{
    int i=0, j=0, k=0;
#pragma omp parallel for collapse(3) schedule(dynamic) private(i,j,k)
    for (i = 0; i < txsize; i++)
    for (j = 0; j < tysize; j++)
    for (k = 0; k < tzsize; k++)
    {
        Hx[i][j][k] = Hx[i][j][k] - Db / H_dy[j] * (Ez[i][j + 1][k] - Ez[i][j][k])
            + Db / H_dz[k] * (Ey[i][j][k + 1] - Ey[i][j][k]);
    }
}
void TFDTDparameters::Hy_field_update()
{
    int i=0, j=0, k=0;
#pragma omp parallel for collapse(3) schedule(dynamic) private(i,j,k)
    for (i = 0; i < txsize; i++)
    for (j = 0; j < tysize; j++)
    for (k = 0; k < tzsize; k++)
    {
        Hy[i][j][k] = Hy[i][j][k] - Db / H_dz[k] * (Ex[i][j][k + 1] - Ex[i][j][k])
            + Db / H_dx[i] * (Ez[i + 1][j][k] - Ez[i][j][k]);
    }
}
void TFDTDparameters::Hz_field_update()
{
    int i=0, j=0, k=0;
#pragma omp parallel for collapse(3) schedule(dynamic) private(i,j,k)
    for (i = 0; i < txsize; i++)
    for (j = 0; j < tysize; j++)
    for (k = 0; k < tzsize; k++)
    {
        Hz[i][j][k] = Hz[i][j][k] - Db / H_dx[i] * (Ey[i + 1][j][k] - Ey[i][j][k])
            + Db / H_dy[j] * (Ex[i][j + 1][k] - Ex[i][j][k]);
    }
}
void TFDTDparameters::Spatial_Size_With_Kappa()
{
    for (int ii = 1; ii < txsize; ii++)
    {
        E_dx[ii] = ((dx[ii] + dx[ii - 1]) / 2);
        H_dx[ii] = dx[ii];
    }
    for (int jj = 1; jj < tysize; jj++)
    {
        E_dy[jj] = ((dy[jj] + dy[jj - 1]) / 2);
        H_dy[jj] = dy[jj];
    }
    for (int kk = 1; kk < tzsize; kk++)
    {
        E_dz[kk] = ((dz[kk] + dz[kk - 1]) / 2);
        H_dz[kk] = dz[kk];
    }

    E_dx[0] = E_dx[1];
    H_dx[0] = H_dx[1];
    E_dy[0] = E_dy[1];
    H_dy[0] = H_dy[1];
    E_dz[0] = E_dz[1];
    H_dz[0] = H_dz[1];

    kappa_e_xprof = Tptr1<double>(txsize); // profile in x-direction
    kappa_h_xprof = Tptr1<double>(txsize);
    kappa_e_yprof = Tptr1<double>(tysize); // profile in x-direction
    kappa_h_yprof = Tptr1<double>(tysize);
    kappa_e_zprof = Tptr1<double>(tzsize); // profile in x-direction
    kappa_h_zprof = Tptr1<double>(tzsize);

    for (int i = 0; i < txsize; i++)
    {
        kappa_e_xprof[i] = 1;
        kappa_h_xprof[i] = 1;
    }
    for (int i = 0; i < tysize; i++)
    {
        kappa_e_yprof[i] = 1;
        kappa_h_yprof[i] = 1;
    }
    for (int i = 0; i < tzsize; i++)
    {
        kappa_e_zprof[i] = 1;
        kappa_h_zprof[i] = 1;
    }

    Kappa_parameter_initiateX();
    Kappa_parameter_initiateY();
    Kappa_parameter_initiateZ();

    for (int ii = 0; ii < txsize; ii++)
    {
        E_dx[ii] = E_dx[ii] * kappa_e_xprof[ii];
        H_dx[ii] = H_dx[ii] * kappa_h_xprof[ii];
    }
    for (int jj = 0; jj < tysize; jj++)
    {
        E_dy[jj] = E_dy[jj] * kappa_e_yprof[jj];
        H_dy[jj] = H_dy[jj] * kappa_h_yprof[jj];
    }
    for (int kk = 0; kk < tzsize; kk++)
    {
        E_dz[kk] = E_dz[kk] * kappa_e_zprof[kk];
        H_dz[kk] = H_dz[kk] * kappa_h_zprof[kk];
    }
}

void TFDTDparameters::Kappa_parameter_initiateX()
{
    double* kappa_e_xprof_cal;
    double* kappa_h_xprof_cal;
    kappa_e_xprof_cal = Tptr1<double>(PMLX + 1);
    kappa_h_xprof_cal = Tptr1<double>(PMLX + 1);

    for (i = 0; i < PMLZ + 1; i++)
    {
        kappa_e_xprof_cal[i] = 1.;
        kappa_h_xprof_cal[i] = 1.;
    }

    for (i = 0; i < PMLX + 1; i++)
    {
        kappa_e_xprof_cal[i] = 1. + (kappa_max - 1.) * pow((i + 0.0) / PMLX, m_kappa);
        kappa_h_xprof_cal[i] = 1. + (kappa_max - 1.) * pow((i + 1.0 / 2.0) / PMLX, m_kappa);
    }

    int ai, bi;
    for (i = 0; i < PMLX; i++)
    {
        bi = PMLX - 1 - i; //recursive index
        ai = PMLX - 1 - i; //recursive index
        kappa_e_xprof[i + 1] = kappa_e_xprof_cal[ai];
        kappa_e_xprof[txsize - PMLX + i] = kappa_e_xprof_cal[i];
        kappa_h_xprof[i] = kappa_h_xprof_cal[bi];
        kappa_h_xprof[txsize - PMLX + i] = kappa_h_xprof_cal[i];
    }

    memfree(&kappa_e_xprof_cal);
    memfree(&kappa_h_xprof_cal);
}

void TFDTDparameters::Kappa_parameter_initiateY()
{
    double* kappa_e_yprof_cal;
    double* kappa_h_yprof_cal;
    kappa_e_yprof_cal = Tptr1<double>(PMLY + 1);
    kappa_h_yprof_cal = Tptr1<double>(PMLY + 1);

    for (i = 0; i < PMLY + 1; i++)
    {
        kappa_e_yprof_cal[i] = 1.;
        kappa_h_yprof_cal[i] = 1.;
    }

    for (i = 0; i < PMLY + 1; i++)
    {
        kappa_e_yprof_cal[i] = 1. + (kappa_max - 1.) * pow((i + 0.0) / PMLY, m_kappa);
        kappa_h_yprof_cal[i] = 1. + (kappa_max - 1.) * pow((i + 1.0 / 2.0) / PMLY, m_kappa);
    }

    int aj, bj;
    for (j = 0; j < PMLY; j++)
    {
        bj = PMLY - 1 - j; //recursive index
        aj = PMLY - 1 - j; //recursive index
        kappa_e_yprof[j + 1] = kappa_e_yprof_cal[aj];
        kappa_e_yprof[tysize - PMLY + j] = kappa_e_yprof_cal[j];
        kappa_h_yprof[j] = kappa_h_yprof_cal[bj];
        kappa_h_yprof[tysize - PMLY + j] = kappa_h_yprof_cal[j];
    }

    memfree(&kappa_e_yprof_cal);
    memfree(&kappa_h_yprof_cal);
}

void TFDTDparameters::Kappa_parameter_initiateZ()
{
    double* kappa_e_zprof_cal;
    double* kappa_h_zprof_cal;
    kappa_e_zprof_cal = Tptr1<double>(PMLZ + 1);
    kappa_h_zprof_cal = Tptr1<double>(PMLZ + 1);

    for (k = 0; k < PMLZ + 1; k++)
    {
        kappa_e_zprof_cal[k] = 1.;
        kappa_h_zprof_cal[k] = 1.;
    }

    for (k = 0; k < PMLZ + 1; k++)
    {
        kappa_e_zprof_cal[k] = 1. + (kappa_max - 1.) * pow((k + 0.0) / PMLZ, m_kappa);
        kappa_h_zprof_cal[k] = 1. + (kappa_max - 1.) * pow((k + 0.5) / PMLZ, m_kappa);
    }

    int ak, bk;
    for (k = 0; k < PMLZ; k++)
    {
        bk = PMLZ - 1 - k; //recursive index
        ak = PMLZ - 1 - k; //recursive index
        kappa_e_zprof[k + 1] = kappa_e_zprof_cal[ak];
        kappa_e_zprof[tzsize - PMLZ + k] = kappa_e_zprof_cal[k];
        kappa_h_zprof[k] = kappa_h_zprof_cal[bk];
        kappa_h_zprof[tzsize - PMLZ + k] = kappa_h_zprof_cal[k];
    }

    FILE* pf_ke;
    FILE* pf_kh;
    pf_ke = fopen("kappa_e.txt", "w");
    pf_kh = fopen("kappa_h.txt", "w");
    for (int i = 0; i < tzsize; i++)
    {
        fprintf(pf_ke, "%lf\n", kappa_e_zprof[i]);
        fprintf(pf_kh, "%lf\n", kappa_h_zprof[i]);
    }

    fclose(pf_ke);
    fclose(pf_kh);
    memfree(&kappa_e_zprof_cal);
    memfree(&kappa_h_zprof_cal);
}
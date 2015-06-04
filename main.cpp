#include <iostream>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>
#include <N_PDS_Manager.h>
#include <N_PDS_ParMap.h>
#include <N_IO_CmdParse.h>
#include <N_DEV_ADMSbsimcmg.h>

#include <N_UTL_Misc.h>

#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_DEV_RegisterDevices.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_System.h>
#include <N_LAS_LAFactory.h>

#include "Epetra_SerialDenseMatrix.h"
using namespace std;

/*
A = [
0  0  1  0   0  -1
1  0  0  0  -1   0
0 -1  0  1   0   0
1  0 -1  0   0   0
-1 1  0  0   0   0
0  1 -1  0   0   0
-1 0  0  0   1   0
0  0 -1  0   0   1
0  1  0 -1   0   0
0  0  0  1  -1   0
0  0  0  1   0  -1
0  1  0  0  -1   0
0  0  0  0   1  -1
0  1  0  0   0  -1
]

*/

static const int admsNodeID_d = 0;
static const int admsNodeID_g = 1;
static const int admsNodeID_s = 2;
static const int admsNodeID_e = 3;
static const int admsNodeID_di = 4;
static const int admsNodeID_si = 5;

static const int admsProbeID_V_s_si = 0;
static const int admsProbeID_V_d_di = 1;
static const int admsProbeID_V_e_g = 2;
static const int admsProbeID_V_d_s = 3;
static const int admsProbeID_V_g_d = 4;
static const int admsProbeID_V_g_s = 5;
static const int admsProbeID_V_di_d = 6;
static const int admsProbeID_V_si_s = 7;
static const int admsProbeID_V_g_e = 8;
static const int admsProbeID_V_e_di = 9;
static const int admsProbeID_V_e_si = 10;
static const int admsProbeID_V_g_di = 11;
static const int admsProbeID_V_di_si = 12;
static const int admsProbeID_V_g_si = 13;

void SetupModelParamsForNMOS(N_DEV_DeviceMgr* pDevMgr)
{
  N_DEV_ModelBlock MB("nmos1", "nmos", 107);
  Xyce::Device::Param TP("L", 16e-9);
  MB.params.push_back(TP);
  pDevMgr->addDeviceModel(MB);
}

void SetupModelParamsForPMOS(N_DEV_DeviceMgr* pDevMgr)
{
  N_DEV_ModelBlock MB("pmos1", "pmos", 107);
  Xyce::Device::Param TP("L", 16e-9);
  MB.params.push_back(TP);
  pDevMgr->addDeviceModel(MB);
}

void loadDAEdQdx(Xyce::Device::ADMSbsimcmg::Instance* pInst, Epetra_SerialDenseMatrix& Q)
{
  Q(admsNodeID_d,admsNodeID_d ) +=  -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_d) +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_d_s);
  Q(admsNodeID_d,admsNodeID_g ) +=  +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_d) +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_di) +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_e) +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_si);
  Q(admsNodeID_d,admsNodeID_s ) +=  -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_d_s);
  Q(admsNodeID_d,admsNodeID_e ) +=  +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_e_si) +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_e);
  Q(admsNodeID_d,admsNodeID_di) +=  -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_di) -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_e_di) +pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_di_si);
  Q(admsNodeID_d,admsNodeID_si) +=  -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_e_si) -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_g_si) -pInst->dynamicContributions[admsNodeID_d].dx(admsProbeID_V_di_si);
                              
  Q(admsNodeID_g,admsNodeID_d ) +=  -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_d);
  Q(admsNodeID_g,admsNodeID_g ) +=  +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_d) +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_s) -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_e_g) +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_e) +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_si) +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_di);
  Q(admsNodeID_g,admsNodeID_s ) +=  -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_s);
  Q(admsNodeID_g,admsNodeID_e ) +=  +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_e_g) +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_e_si) +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_e);
  Q(admsNodeID_g,admsNodeID_di) +=  +pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_di_si) -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_di);
  Q(admsNodeID_g,admsNodeID_si) +=  -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_di_si) -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_e_si) -pInst->dynamicContributions[admsNodeID_g].dx(admsProbeID_V_g_si);
                              
  Q(admsNodeID_s,admsNodeID_d ) +=  +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_d_s);
  Q(admsNodeID_s,admsNodeID_g ) +=  +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_s) +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_di) +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_e) +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_si);
  Q(admsNodeID_s,admsNodeID_s ) +=  -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_s) -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_d_s);
  Q(admsNodeID_s,admsNodeID_e ) +=  +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_e_si) +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_e);
  Q(admsNodeID_s,admsNodeID_di) +=  -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_di) -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_e_di) +pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_di_si);
  Q(admsNodeID_s,admsNodeID_si) +=  -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_e_si) -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_g_si) -pInst->dynamicContributions[admsNodeID_s].dx(admsProbeID_V_di_si);
                              
  Q(admsNodeID_e,admsNodeID_d ) +=  -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_di_d);
  Q(admsNodeID_e,admsNodeID_g ) +=  -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_e_g) +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_g_di) +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_g_e) +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_g_si);
  Q(admsNodeID_e,admsNodeID_s ) +=  -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_si_s);
  Q(admsNodeID_e,admsNodeID_e ) +=  +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_e_g) +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_e_si) +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_g_e);
  Q(admsNodeID_e,admsNodeID_di) +=  +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_di_d) -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_g_di) -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_e_di) +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_di_si);
  Q(admsNodeID_e,admsNodeID_si) +=  +pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_si_s) -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_e_si) -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_g_si) -pInst->dynamicContributions[admsNodeID_e].dx(admsProbeID_V_di_si);
                              
  Q(admsNodeID_di,admsNodeID_d ) +=  -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_d) -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_di_d);
  Q(admsNodeID_di,admsNodeID_g ) +=  +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_d) +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_e) +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_si) +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_di);
  Q(admsNodeID_di,admsNodeID_s ) +=  -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_si_s);
  Q(admsNodeID_di,admsNodeID_e ) +=  +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_e_si) +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_e);
  Q(admsNodeID_di,admsNodeID_di) +=  +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_di_d) +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_di_si) -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_di);
  Q(admsNodeID_di,admsNodeID_si) +=  +pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_si_s) -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_di_si) -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_e_si) -pInst->dynamicContributions[admsNodeID_di].dx(admsProbeID_V_g_si);
  
  Q(admsNodeID_si,admsNodeID_d ) +=  -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_di_d);
  Q(admsNodeID_si,admsNodeID_g ) +=  +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_s) +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_e) +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_si) +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_di);
  Q(admsNodeID_si,admsNodeID_s ) +=  -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_s) -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_si_s);
  Q(admsNodeID_si,admsNodeID_e ) +=  +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_e_si) +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_e);
  Q(admsNodeID_si,admsNodeID_di) +=  +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_di_d) +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_di_si) -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_e_di) -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_di);
  Q(admsNodeID_si,admsNodeID_si) +=  +pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_si_s) -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_di_si) -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_e_si) -pInst->dynamicContributions[admsNodeID_si].dx(admsProbeID_V_g_si);

}

void loadDAEdFdx(Xyce::Device::ADMSbsimcmg::Instance* pInst, Epetra_SerialDenseMatrix& M)
{

  M(admsNodeID_d,admsNodeID_d ) +=  +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_d_di) -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_di_d);
  M(admsNodeID_d,admsNodeID_g ) +=  +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_g_e) +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_g_si) +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_g_di);
  M(admsNodeID_d,admsNodeID_s ) +=  -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_si_s);
  M(admsNodeID_d,admsNodeID_e ) +=  +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_e_si) +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_e_di) -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_g_e);
  M(admsNodeID_d,admsNodeID_di) +=  -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_d_di) +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_di_d) -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_e_di) +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_di_si) -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_g_di);
  M(admsNodeID_d,admsNodeID_si) +=  +pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_si_s) -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_e_si) -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_g_si) -pInst->staticContributions[admsNodeID_d].dx(admsProbeID_V_di_si);

  M(admsNodeID_g,admsNodeID_d ) +=  -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_di_d);
  M(admsNodeID_g,admsNodeID_g ) +=  +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_g_di) +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_g_e) +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_g_si);
  M(admsNodeID_g,admsNodeID_s ) +=  -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_si_s);
  M(admsNodeID_g,admsNodeID_e ) +=  +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_e_si) +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_e_di) -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_g_e);
  M(admsNodeID_g,admsNodeID_di) +=  -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_g_di) +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_di_si) -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_e_di) +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_di_d);
  M(admsNodeID_g,admsNodeID_si) +=  -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_di_si) -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_e_si) -pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_g_si) +pInst->staticContributions[admsNodeID_g].dx(admsProbeID_V_si_s);

  M(admsNodeID_s,admsNodeID_d ) +=  -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_di_d);
  M(admsNodeID_s,admsNodeID_g ) +=  +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_g_di) +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_g_e) +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_g_si);
  M(admsNodeID_s,admsNodeID_s ) +=  +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_s_si) -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_si_s);
  M(admsNodeID_s,admsNodeID_e ) +=  +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_e_si) +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_e_di) -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_g_e);
  M(admsNodeID_s,admsNodeID_di) +=  +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_di_d) -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_g_di) -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_e_di) +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_di_si);
  M(admsNodeID_s,admsNodeID_si) +=  -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_s_si) +pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_si_s) -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_e_si) -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_g_si) -pInst->staticContributions[admsNodeID_s].dx(admsProbeID_V_di_si);

  M(admsNodeID_e,admsNodeID_d ) +=  -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_di_d);
  M(admsNodeID_e,admsNodeID_g ) +=  +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_g_di) +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_g_si) +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_g_e);
  M(admsNodeID_e,admsNodeID_s ) +=  -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_si_s);
  M(admsNodeID_e,admsNodeID_e ) +=  +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_e_si) +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_e_di) -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_g_e);
  M(admsNodeID_e,admsNodeID_di) +=  -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_g_di) +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_di_d) +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_di_si) -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_e_di);
  M(admsNodeID_e,admsNodeID_si) +=  +pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_si_s) -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_g_si) -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_di_si) -pInst->staticContributions[admsNodeID_e].dx(admsProbeID_V_e_si);

  M(admsNodeID_di,admsNodeID_d ) +=  +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_d_di) -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_di_d);
  M(admsNodeID_di,admsNodeID_g ) +=  +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_g_e) +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_g_si) +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_g_di);
  M(admsNodeID_di,admsNodeID_s ) +=  -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_si_s);
  M(admsNodeID_di,admsNodeID_e ) +=  +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_e_si) +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_e_di) -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_g_e);
  M(admsNodeID_di,admsNodeID_di) +=  -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_d_di) +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_di_d) -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_e_di) +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_di_si) -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_g_di);
  M(admsNodeID_di,admsNodeID_si) +=  +pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_si_s) -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_e_si) -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_g_si) -pInst->staticContributions[admsNodeID_di].dx(admsProbeID_V_di_si);
  
  M(admsNodeID_si,admsNodeID_d ) +=  -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_di_d);
  M(admsNodeID_si,admsNodeID_g ) +=  +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_g_e) +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_g_si) +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_g_di);
  M(admsNodeID_si,admsNodeID_s ) +=  +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_s_si) -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_si_s);
  M(admsNodeID_si,admsNodeID_e ) +=  +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_e_si) +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_e_di) -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_g_e);
  M(admsNodeID_si,admsNodeID_di) +=  +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_di_d) -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_e_di) +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_di_si) -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_g_di);
  M(admsNodeID_si,admsNodeID_si) +=  -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_s_si) +pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_si_s) -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_e_si) -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_g_si) -pInst->staticContributions[admsNodeID_si].dx(admsProbeID_V_di_si);
}

void GetJacobianMatrix(Xyce::Device::ADMSbsimcmg::Instance* pInst, double Vd, double Vg, double Vs, double Ve)
{
  double Vdi = 0.29999963187823330824;
  double Vsi = 3.6812176753204769947e-07;
  int nNumNode = pInst->numIntVars + pInst->numExtVars;
  Epetra_SerialDenseMatrix Q(nNumNode, nNumNode);
  Epetra_SerialDenseMatrix M(nNumNode, nNumNode);
  pInst->updateIntermediateVarsMy(Vd, Vg, Vs, Ve, Vdi, Vsi);
  loadDAEdFdx(pInst, M);
  loadDAEdQdx(pInst, Q);

  cout << "M =\n" << M;
  cout << "Q =\n" << Q;

}

int main(int nArgc, char** pArgv) {  

  //Xyce::Device::Config<Xyce::Device::ADMSbsimcmg::Traits>& rConfig = Xyce::Device::Config<Xyce::Device::ADMSbsimcmg::Traits>::addConfiguration()
  //  .registerDevice("m", 107)
  //  .registerModelType("nmos", 107)
  //  .registerModelType("pmos", 107);

  Xyce::Device::registerDevices();
 
  N_DEV_InstanceBlock IB2;
  Xyce::IO::CmdParse cp; 
  N_DEV_DeviceMgr* pDevMgr = N_DEV_DeviceMgr::factory(cp);  // Allocate device manager:

  SetupModelParamsForPMOS(pDevMgr);
  SetupModelParamsForNMOS(pDevMgr);

  Xyce::Device::InstanceName instance("M1");
  IB2.setInstanceName(instance);
  IB2.setModelName("nmos1");
  cout << __LINE__ << "\n";
  Xyce::Device::ADMSbsimcmg::Instance* pInst = (Xyce::Device::ADMSbsimcmg::Instance*) pDevMgr->addDeviceInstance(IB2);

  pInst->updateTemperature(25);
  
  cout << __LINE__ << "\n";
  
  double Vd = 0.3;
  double Vg = 1.0;
  double Vs = 0.0;
  double Ve = 0.0;
 
  GetJacobianMatrix(pInst, Vd, Vg, Vs, Ve);

  return 0;
}

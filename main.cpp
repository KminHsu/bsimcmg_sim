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

using namespace std;

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

int main(int nArgc, char** pArgv) {  
  cout << "Hello !!\n";

  //Xyce::Device::Config<Xyce::Device::ADMSbsimcmg::Traits>& rConfig = Xyce::Device::Config<Xyce::Device::ADMSbsimcmg::Traits>::addConfiguration()
  //  .registerDevice("m", 107)
  //  .registerModelType("nmos", 107)
  //  .registerModelType("pmos", 107);

  Xyce::Device::registerDevices();
 
  N_DEV_InstanceBlock IB2;
  Xyce::IO::CmdParse cp; 
  N_DEV_DeviceMgr* pDevMgr = N_DEV_DeviceMgr::factory(cp);  // Allocate device manager:

  //pDevMgr->initializeAll();

  SetupModelParamsForPMOS(pDevMgr);
  SetupModelParamsForNMOS(pDevMgr);

  // cout << __LINE__ << "\n";

  // N_DEV_ModelBlock MB("nmos1", "nmos", 107);
  // MB.params.push_back(Xyce::Device::Param("L", 16e-9));
  // MB.params.push_back(Xyce::Device::Param("BULKMOD", 1));
  // pDevMgr->addDeviceModel(MB);

  // cout << __LINE__ << "\n";

  // N_DEV_ModelBlock MB1("pmos1", "pmos", 107);
  // MB1.params.push_back(Xyce::Device::Param("L", 16e-9));
  // MB1.params.push_back(Xyce::Device::Param("BULKMOD", 1));
  // pDevMgr->addDeviceModel(MB1);
  // 
  // cout << __LINE__ << "\n";

  Xyce::Device::InstanceName instance("M1");
  IB2.setInstanceName(instance);
  IB2.setModelName("nmos1");
  cout << __LINE__ << "\n";
  Xyce::Device::ADMSbsimcmg::Instance* pInst = (Xyce::Device::ADMSbsimcmg::Instance*) pDevMgr->addDeviceInstance(IB2);
  
  cout << __LINE__ << "\n";
  
  double Vd = 0.3;
  double Vg = 1.0;
  double Vs = 0.0;
  double Ve = 0.0;
  double Vdi = 0.0;
  double Vsi = 0.0;
 
  pInst->updateIntermediateVarsMy(Vd, Vg, Vs, Ve, Vdi, Vsi);


  return 0;
}

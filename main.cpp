#include <iostream>
#include <algorithm>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>
#include <N_PDS_Manager.h>
#include <N_PDS_ParMap.h>
#include <N_IO_CmdParse.h>
#include <N_DEV_ADMSbsimcmg.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_CompositeParam.h>
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

#include <QFile>
#include <QString>
#include <QRegExp>

#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Xyce;
using namespace Device;

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

static const Xyce::Device::Configuration* pGlobalConfig = NULL;

static QRegExp rxComment("[ ]*[*]+.*");
static QRegExp rxModel("[ ]*[.]model[ ]*([a-zA-Z0-9]+)[ ]*([a-zA-Z0-9]+)[ ](level|LEVEL)=([0-9]+).*");
static QRegExp rxKeyValue("[+][ ]*([A-Za-z0-9]+)[ ]*=[ ]*([0-9eE+-.]+).*");
 
void IncludeModelCard(N_DEV_DeviceMgr* pDevMgr, char * pFileName)
{
  QString qsModelFileName(pFileName);
  QFile modelFile( qsModelFileName );

  if( ! modelFile.open( QIODevice::ReadOnly | QIODevice::Text ) )
  {
    cout << "Can not open file for read : " << qsModelFileName.toLocal8Bit().data() ;
  }

 QString qsLine;
 QString modelName;
 QString modelType;
 int modelLevel = 107;
 QString key;
 double value;
 N_DEV_ModelBlock* pMB = NULL;
 while( ! modelFile.atEnd() )
  {
    qsLine = modelFile.readLine();
    if( qsLine == "\n" )
    {
      continue;
    }
    if( rxComment.exactMatch( qsLine) )
      continue;
    else if( rxModel.exactMatch( qsLine ) )
    {
      cout << qsLine.toLocal8Bit().data();
      cout << rxModel.cap(1).toLocal8Bit().data() << endl;
      cout << rxModel.cap(2).toLocal8Bit().data() << endl;
      cout << rxModel.cap(3).toLocal8Bit().data() << endl;
      cout << rxModel.cap(4).toLocal8Bit().data() << endl;
      
      modelName = rxModel.cap(1).trimmed();
      modelType = rxModel.cap(2).trimmed();
      modelLevel = rxModel.cap(4).trimmed().toInt();
    }
    else if( rxKeyValue.exactMatch( qsLine ) )
    {
      //cout << qsLine.toLocal8Bit().data();
      if( pMB == NULL )
      {
        if( !modelName.isEmpty() && !modelType.isEmpty() )
        {
          pMB = new N_DEV_ModelBlock(modelName.toLocal8Bit().data(), modelType.toLocal8Bit().data(), modelLevel);
        }
      }

      if( pMB != NULL )
      {
        cout << rxKeyValue.cap(1).toLocal8Bit().data() << " = "
             << rxKeyValue.cap(2).toLocal8Bit().data() << endl;
        key = rxKeyValue.cap(1).trimmed();
        value = rxKeyValue.cap(2).trimmed().toDouble();
        Xyce::Device::Param TP(key.toLocal8Bit().data(), value);
        //Xyce::Device::Param TP(key.toLocal8Bit().data(), (const char*)rxKeyValue.cap(2).trimmed().toLocal8Bit().data(), true);
 
        pMB->params.push_back(TP);
      }
    }
    else 
    {
      cout << "Undefined : " << qsLine.toLocal8Bit().data();
    }
  }
  if( NULL != pMB )
  {
    const Xyce::Device::Configuration* pConfig = Xyce::Device::Configuration::findConfiguration("m", modelLevel);
    if( NULL == pGlobalConfig )
      pGlobalConfig = pConfig;
    //cout << "pConfig : " << pConfig << "\n";
    if( NULL != pConfig )
    {
      cout << pConfig << endl;
      Xyce::Device::ParametricData<void> &model_parameters = pConfig->getModelParameters();
      
      const Xyce::Device::ParameterMap& parameter_map = model_parameters.getMap();
      for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
        cout << "model " << (*it).first << " = "; 
        string key = (*it).first;

	// for( pMB->params::const_iterator it2 = pMB->params.begin(); it2 != pMB->params.end(); ++it2) {
        //   if( 
        // }

        const Descriptor &descriptor = *(*it).second;
        if( descriptor.isType<int>() )
        {
           //cout << "Integer type : " << endl;
           int value = (int)Xyce::Device::getDefaultValue<int>(descriptor);
           cout << value << endl;
           if( pMB->params.end() != find( pMB->params.begin(), pMB->params.end(), Xyce::Device::Param(key, value) ) )
           {
             cout << "skip int !!" << endl;
             continue;
           }
           pMB->params.push_back(Xyce::Device::Param(key, value));
        }
        else if( descriptor.isType<double>() )
        {
           //cout << "Double type : " << endl;
           double value = (double)Xyce::Device::getDefaultValue<double>(descriptor);
           cout << value << endl;
           if( pMB->params.end() != find( pMB->params.begin(), pMB->params.end(), Xyce::Device::Param(key, value) ) )
           {
             cout << "skip double !!" << endl;
             continue;
           }
           pMB->params.push_back(Xyce::Device::Param(key, value));
        }
        else
        {
          cout << "Unknow Type !!" << endl;
        }
      }
    }

    cout << "Add Model : " << modelName.toLocal8Bit().data() << " " << modelType.toLocal8Bit().data() << " " << modelLevel << endl;
    pDevMgr->addDeviceModel(*pMB);
  }
}

//void SetupModelParamsForNMOS(N_DEV_DeviceMgr* pDevMgr)
//{
//  N_DEV_ModelBlock MB("nmos1", "nmos", 107);
//  Xyce::Device::Param TP("L", 16e-9);
//  MB.params.push_back(TP);
//  pDevMgr->addDeviceModel(MB);
//}
//
//void SetupModelParamsForPMOS(N_DEV_DeviceMgr* pDevMgr)
//{
//  N_DEV_ModelBlock MB("pmos1", "pmos", 107);
//  Xyce::Device::Param TP("L", 16e-9);
//  MB.params.push_back(TP);
//  pDevMgr->addDeviceModel(MB);
//}

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

  IncludeModelCard(pDevMgr, "modelcard.nmos_xyce");
  //IncludeModelCard(pDevMgr, "modelcard.pmos_xyce");
  //SetupModelParamsForPMOS(pDevMgr);
  //SetupModelParamsForNMOS(pDevMgr);

  Xyce::Device::InstanceName instance("M1");
  IB2.setInstanceName(instance);
  IB2.setModelName("nmos1");
  cout << __LINE__ << "\n";

  IB2.params.push_back( Xyce::Device::Param("L"   , 3e-8, true) );
  IB2.params.push_back( Xyce::Device::Param("TFIN", 1.5e-8, true) );
  IB2.params.push_back( Xyce::Device::Param("NFIN", 10.0, true) );
  IB2.params.push_back( Xyce::Device::Param("NRS" , 1.0, true) );
  IB2.params.push_back( Xyce::Device::Param("NRD" , 1.0, true) );
  
  Xyce::Device::ParametricData<void>& inst_parameters = pGlobalConfig->getInstanceParameters();
  const Xyce::Device::ParameterMap& parameter_map = inst_parameters.getMap();
  for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
    cout << "inst " << (*it).first << " = "; 
    string key = (*it).first;

    const Descriptor &descriptor = *(*it).second;
    if( descriptor.isType<int>() )
    {
       //cout << "Integer type : " << endl;
       int value = (int)Xyce::Device::getDefaultValue<int>(descriptor);
       cout << value << endl;
       if( IB2.params.end() != find( IB2.params.begin(), IB2.params.end(), Xyce::Device::Param(key, value) ) )
       {
         cout << "skip int !!" << endl;
         continue;
       }
       IB2.params.push_back(Xyce::Device::Param(key, value));
    }
    else if( descriptor.isType<double>() )
    {
       //cout << "Double type : " << endl;
       double value = (double)Xyce::Device::getDefaultValue<double>(descriptor);
       cout << value << endl;
       if( IB2.params.end() != find( IB2.params.begin(), IB2.params.end(), Xyce::Device::Param(key, value) ) )
       {
         cout << "skip double !!" << endl;
         continue;
       }
       IB2.params.push_back(Xyce::Device::Param(key, value));
    }
    else
    {
      cout << "Unknow Type !!" << endl;
    }
  }

  Xyce::Device::ADMSbsimcmg::Instance* pInst = (Xyce::Device::ADMSbsimcmg::Instance*) pDevMgr->addDeviceInstance(IB2);
  pInst->processParams();
  pDevMgr->updateTemperature(25.0);
  
  double Vd = 0.2999999999999999889;
  double Vg = 1.0;
  double Vs = 0.0;
  double Ve = 0.0;
 
  GetJacobianMatrix(pInst, Vd, Vg, Vs, Ve);

  return 0;
}
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
#include <QStringList>
#include <QRegExp>
#include <QMap>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
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

GG = [
   1.0000000000000001137e+03   0.0000000000000000000e+00   0.0000000000000000000e+00   0.0000000000000000000e+00  -1.0000000000000001137e+03   0.0000000000000000000e+00;
   0.0000000000000000000e+00   4.9922949716673672339e-12   0.0000000000000000000e+00   0.0000000000000000000e+00  -1.9650936694031349837e-12  -3.0272013022642322503e-12;
   0.0000000000000000000e+00   0.0000000000000000000e+00   1.0000000000000001137e+03   0.0000000000000000000e+00   0.0000000000000000000e+00  -1.0000000000000001137e+03;
   0.0000000000000000000e+00   0.0000000000000000000e+00   0.0000000000000000000e+00   9.5004850870021704238e-89  -9.5004850869942244848e-89  -7.9458140349868174397e-101;
  -1.0000000000000001137e+03   3.1071226013801485095e-04   0.0000000000000000000e+00   4.2344138429478520918e-10   1.0000008713331211538e+03  -1.1820458046115692695e-03;
   0.0000000000000000000e+00  -3.1071226513030983981e-04  -1.0000000000000001137e+03  -4.2344138429478520918e-10  -8.7133311906707636239e-04   1.0000011820458076954e+03
]

VV = [
   0.3;
   1.0;
   0.0;
   0.0;
   2.9999963187823330824e-01;
   3.6812176753204769947e-07
]

I = [
   3.6812176665534934727e-04;
   4.4027664798614864832e-12;
  -3.6812176753204774648e-04;
  -2.8501420287629128450e-89;
  0;
  0
]


G = [
   -GG(1,1:6) 1 0;
   -GG(2,1:6) 0 1;
   -GG(3,1:6) 0 0;
   -GG(4,1:6) 0 0;
   -GG(5,1:6) 0 0;
   -GG(6,1:6) 0 0;
   1 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0
]

Gn = [
   -GG(1,1:6) 1 0;
   -GG(2,1:6) 0 1;
   1 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0
]


f3 = [
0.00036812176668066243224
8.2468941300689627635e-13
-0.00036812176753204769227
-1.0013301663122570181e-89
2.6695876412241581477e-14
0]

Vn = [
   0.3;
   1.0;
   0;
   0;
   2.9999963187823330824e-01;
   3.6812176753204769947e-07;
   3.6812176665534934727e-04;
   4.4027664798614864832e-12
]

V = [
   0.3;
   1.0;
   0.0;
   0.0;
   2.9999963187823330824e-01;
   3.6812176753204769947e-07;
   3.6812176665534934727e-04;
   4.4027664798614864832e-12
]

II = [
   3.6812176665534934727e-04;
   4.4027664798614864832e-12;
  -3.6812176753204774648e-04;
  -2.8501420287629128450e-89;
  0;
  0;
  0;
  0
]

F3 = [
f3;
0.3;
1.0
]

disi = [
   2.9999963187823330824e-01;
   3.6812176753204769947e-07
]

Mdisi = [
-1.0000000000000001137e+03   0.0000000000000000000e+00 ;
-1.9650936694031349837e-12  -3.0272013022642322503e-12 ;
 0.0000000000000000000e+00  -1.0000000000000001137e+03 ;
-9.5004850869942244848e-89  -7.9458140349868174397e-101 ;
 0 0;
 0 0;
]

f = [
0.00036812176668066243224;
8.2468941300689627635e-13;
-0.00036812176753204769227;
-1.0013301663122570181e-89;
0;
0
]


MM = [
   -GG(1,1:2) 1 0;
   -GG(2,1:2) 0 1;
   1 0 0 0 ;
   0 1 0 0 
]

XX = [
   0.3;
   1.0;
   0.00036812176668066243224;
   8.2468941300689627635e-13
]

f3 = [
0.00036812176668066243224
8.2468941300689627635e-13
-0.00036812176753204769227
-1.0013301663122570181e-89
2.6695876412241581477e-14
0]

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


static QMap<QString, Xyce::Device::ADMSbsimcmg::Instance*> g_theInstMap;
static N_DEV_DeviceMgr* g_pDevMgr = NULL;
 
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
        //value = rxKeyValue.cap(2).trimmed().toDouble();
        //Xyce::Device::Param TP(key.toLocal8Bit().data(), value);
        Xyce::Device::Param TP(key.toLocal8Bit().data(), (const char*)rxKeyValue.cap(2).trimmed().toLocal8Bit().data(), true);
 
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

void updateIntermediateVarsMy(Xyce::Device::ADMSbsimcmg::Instance* pInst, 
                              double Vd, 
                              double Vg, 
                              double Vs, 
                              double Ve,
                              double Vdi,
                              double Vsi,
                              Epetra_SerialDenseMatrix& M,
                              Epetra_SerialDenseMatrix& Q,
                              Epetra_SerialDenseMatrix& F,
                              Epetra_SerialDenseMatrix& I
)
{
  pInst->updateIntermediateVarsMy(Vd, Vg, Vs, Ve, Vdi, Vsi);
  loadDAEdFdx(pInst, M);
  loadDAEdQdx(pInst, Q);
  I(admsNodeID_d,0) = pInst->staticContributions[admsNodeID_d].val();
  I(admsNodeID_g,0) = pInst->staticContributions[admsNodeID_g].val();
  I(admsNodeID_s,0) = pInst->staticContributions[admsNodeID_s].val();
  I(admsNodeID_e,0) = pInst->staticContributions[admsNodeID_e].val();
  I(admsNodeID_di,0) = pInst->staticContributions[admsNodeID_di].val();
  I(admsNodeID_si,0) = pInst->staticContributions[admsNodeID_si].val();

  F(admsNodeID_d,0) = pInst->staticContributions[admsNodeID_d].val();
  F(admsNodeID_g,0) = pInst->staticContributions[admsNodeID_g].val();
  F(admsNodeID_s,0) = pInst->staticContributions[admsNodeID_s].val();
  F(admsNodeID_e,0) = pInst->staticContributions[admsNodeID_e].val();
  F(admsNodeID_di,0) = pInst->staticContributions[admsNodeID_di].val();
  F(admsNodeID_si,0) = pInst->staticContributions[admsNodeID_si].val();
  cout << "M =\n" << M;
  cout << "Q =\n" << Q;
}

void GetJacobianMatrix(Xyce::Device::ADMSbsimcmg::Instance* pInst, 
double Vd,
double Vg,
double Vs,
double Ve,
double* pM,
double* pQ,
double* pF,
double* pI)
{
  unsigned int nIter = 100;
  double Vdi = 0.0;
  double Vsi = 0.0;
  int nNumExtNode = pInst->numExtVars;
  int nNumIntNode = pInst->numExtVars;
  int nNumNode = pInst->numIntVars + pInst->numExtVars;
  Epetra_SerialDenseMatrix Q_(nNumNode, nNumNode);
  Epetra_SerialDenseMatrix M_(nNumNode, nNumNode);
  Epetra_SerialDenseMatrix I_(nNumNode, 1);
  Epetra_SerialDenseMatrix F_(nNumNode, 1);
  Epetra_SerialDenseMatrix Vdisi(pInst->numIntVars, 1);
  Epetra_SerialDenseMatrix Fdisi(pInst->numIntVars, 1);
  Epetra_SerialDenseMatrix B_(pInst->numIntVars, pInst->numIntVars);
  Epetra_SerialDenseMatrix invB_(pInst->numIntVars, pInst->numIntVars);
  Epetra_SerialDenseSolver solver;

  cout << setprecision(20); 
  for(int i = 0; i < nIter; i++ )
  {
    M_.Scale(0.0);
    Q_.Scale(0.0);
    I_.Scale(0.0);
    F_.Scale(0.0);
 
    if( i == 0 )
    {
      updateIntermediateVarsMy(pInst, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, M_, Q_, F_, I_);
      Fdisi(0,0) = -(Vd*M_(4,0)+Vg*M_(4,1)+Vs*M_(4,2)+Ve*M_(4,3));
      Fdisi(1,0) = -(Vd*M_(5,0)+Vg*M_(5,1)+Vs*M_(5,2)+Ve*M_(5,3));
    }
    else
    {
      updateIntermediateVarsMy(pInst, Vd, Vg, Vs, Ve, Vdi, Vsi, M_, Q_, F_, I_);
      Fdisi(0,0) = -F_(admsNodeID_di,0);
      Fdisi(1,0) = -F_(admsNodeID_si,0);
    }
    cout << "M_ = " << endl << M_; 
    cout << "Q_ = " << endl << Q_; 
    cout << "F_ = " << endl << F_; 
    cout << "I_ = " << endl << I_; 

    B_(0,0) = M_(admsNodeID_di,admsNodeID_di); 
    B_(0,1) = M_(admsNodeID_di,admsNodeID_si);
    B_(1,0) = M_(admsNodeID_si,admsNodeID_di);
    B_(1,1) = M_(admsNodeID_si,admsNodeID_si);
    invB_ = B_;
    cout << "B_ = " << endl << B_; 

    cout << "line = " << __LINE__ << endl;
    
    //Crash in Matlab
    //solver.SetMatrix(invB_);
    //solver.Invert();

    invB_(0,0) =  B_(1,1); 
    invB_(0,1) = -B_(0,1); 
    invB_(1,0) = -B_(1,0);
    invB_(1,1) =  B_(0,0);
    invB_.Scale( 1.0/(B_(0,0)*B_(1,1) - B_(1,0)*B_(0,1)) ); 

    cout << "line = " << __LINE__ << endl;
    cout << "invB_ = " << endl << invB_; 
    cout << "Fdisi = " << endl << Fdisi; 

    cout << "line = " << __LINE__ << endl;
    //Not work in Matlab
    //Vdisi.Multiply('N', 'N', 1.0, invB_, Fdisi, 0.0); 
    Vdisi(0,0) = invB_(0,0)*Fdisi(0,0) + invB_(0,1)*Fdisi(1,0);
    Vdisi(1,0) = invB_(1,0)*Fdisi(0,0) + invB_(1,1)*Fdisi(1,0);
    cout << "Vdisi = " << endl << Vdisi; 

    if( i != 0 && abs(Vdisi(0,0)) < 1e-12 && abs(Vdisi(1,0)) < 1e-12 )
    {
      cout << "Device convergence!!" << endl;
      cout << "Iter = " << i+1 << endl;
      for(int m = 0; m < nNumExtNode; m++)
      {
        //pF[m] = F_(m,0);
        //pF[m] = -(M_(m,4)*Vdi+M_(m,5)*Vsi);
        //pI[m] = I_(m,0);
        for(int n = 0; n < nNumExtNode; n++)
        {
          pM[m*nNumExtNode+n] = M_(m,n);
          pQ[m*nNumExtNode+n] = Q_(m,n);
        }
      }

      for(int m = 0; m < nNumExtNode; m++)
      {
        pF[m] = -(M_(m,5)*Vdi+M_(m,6)*Vsi);
        pI[m] = I_(m,0);
      }
      Vdi = Vdi + Vdisi(0,0);
      Vsi = Vsi + Vdisi(1,0);
      break;
    }
    Vdi = Vdi + Vdisi(0,0);
    Vsi = Vsi + Vdisi(1,0);
  } 
  cout << "Vd = " << Vd << endl
       << "Vg = " << Vg << endl
       << "Vs = " << Vs << endl
       << "Ve = " << Ve << endl
       << "Vdi = " << Vdi << endl
       << "Vsi = " << Vsi << endl;
}

Xyce::Device::ADMSbsimcmg::Instance* CreateInst(char* psInstName, char* psModelName, char* psInstParams)
{
  QString qsInstParams(psInstParams);
  QStringList paramList = qsInstParams.split(" ");

  cout << psInstName << " " << psModelName << endl;
  QStringList keyValue;
  N_DEV_InstanceBlock IB2;
  Xyce::Device::InstanceName instance(psInstName);
  IB2.setInstanceName(instance);
  IB2.setModelName(psModelName);

  for(int i = 0; i < paramList.size(); i++)
  {
    if( paramList.at(i).trimmed() == "" )
     continue;
    keyValue = paramList.at(i).trimmed().split("=");
    cout << keyValue.at(0).toLocal8Bit().data() << " = " << keyValue.at(1).toLocal8Bit().data() << endl;
    IB2.params.push_back( Xyce::Device::Param((char*)keyValue.at(0).toLocal8Bit().data() , keyValue.at(1).toDouble(), true) );
  }

  //IB2.params.push_back( Xyce::Device::Param("L"   , 3e-8, true) );
  //IB2.params.push_back( Xyce::Device::Param("TFIN", 1.5e-8, true) );
  //IB2.params.push_back( Xyce::Device::Param("NFIN", 10.0, true) );
  //IB2.params.push_back( Xyce::Device::Param("NRS" , 1.0, true) );
  //IB2.params.push_back( Xyce::Device::Param("NRD" , 1.0, true) );
  Xyce::Device::ADMSbsimcmg::Instance* pInst = (Xyce::Device::ADMSbsimcmg::Instance*) g_pDevMgr->addDeviceInstance(IB2);

  cout << "Line:" << __LINE__ << "\n";
  pInst->processParams();
  cout << "Line:" << __LINE__ << "\n";
  pInst->getModel().processParams();
  cout << "Line:" << __LINE__ << "\n";

  g_theInstMap.insert( QString(psInstName), pInst);
  return pInst;
}

void Initialize()
{
  //  .registerDevice("m", 107)
  //  .registerModelType("nmos", 107)
  //  .registerModelType("pmos", 107);

  Xyce::Device::registerDevices();
  cout << "Line:" << __LINE__ << "\n";
 
  Xyce::IO::CmdParse cp; 
  static N_DEV_DeviceMgr* pDevMgr = N_DEV_DeviceMgr::factory(cp);  // Allocate device manager:
  g_pDevMgr = pDevMgr;
  pDevMgr->updateTemperature(25.0);
  cout << "Line:" << __LINE__ << "\n";

  IncludeModelCard(pDevMgr, "modelcard.nmos_xyce");
  IncludeModelCard(pDevMgr, "modelcard.pmos_xyce");
  cout << "Line:" << __LINE__ << "\n";
}

void BSIMCMG(char* psInstName, double Vd, double Vg, double Vs, double Ve, double* M, double* Q, double* F, double* I)
{
  //Xyce::Device::Config<Xyce::Device::ADMSbsimcmg::Traits>& rConfig = Xyce::Device::Config<Xyce::Device::ADMSbsimcmg::Traits>::addConfiguration()
  //Xyce::Device::ParametricData<void>& inst_parameters = pGlobalConfig->getInstanceParameters();
  //const Xyce::Device::ParameterMap& parameter_map = inst_parameters.getMap();
  //for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
  //  cout << "inst " << (*it).first << " = "; 
  //  string key = (*it).first;
  //  const Descriptor &descriptor = *(*it).second;
  //  if( descriptor.isType<int>() )
  //  {
  //     //cout << "Integer type : " << endl;
  //     int value = (int)Xyce::Device::getDefaultValue<int>(descriptor);
  //     cout << value << endl;
  //     if( IB2.params.end() != find( IB2.params.begin(), IB2.params.end(), Xyce::Device::Param(key, value) ) )
  //     {
  //       cout << "skip int !!" << endl;
  //       continue;
  //     }
  //     IB2.params.push_back(Xyce::Device::Param(key, value));
  //  }
  //  else if( descriptor.isType<double>() )
  //  {
  //     //cout << "Double type : " << endl;
  //     double value = (double)Xyce::Device::getDefaultValue<double>(descriptor);
  //     cout << value << endl;
  //     if( IB2.params.end() != find( IB2.params.begin(), IB2.params.end(), Xyce::Device::Param(key, value) ) )
  //     {
  //       cout << "skip double !!" << endl;
  //       continue;
  //     }
  //     IB2.params.push_back(Xyce::Device::Param(key, value));
  //  }
  //  else
  //  {
  //    cout << "Unknow Type !!" << endl;
  //  }
  //}

  //Xyce::Device::ADMSbsimcmg::Instance* pInst = CreateInst("M1", "nmos1", "L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0");
  Xyce::Device::ADMSbsimcmg::Instance* pInst = g_theInstMap[QString(psInstName)];
  if( NULL != pInst )
  {
    GetJacobianMatrix(pInst, Vd, Vg, Vs, Ve, M, Q, F, I);
  }
  else
    cout << "Can not find instance :" << psInstName << "\n";
}

int main(int nArgc, char** pArgv) {  

  double Vd = 0.2999999999999999889;
  double Vg = 1.0;
  double Vs = 0.0;
  double Ve = 0.0;
 
  Initialize();
  CreateInst("M1", "nmos1", "L=3e-8 TFIN=1.5e-8 NFIN=10.0 NRS=1.0 NRD=1.0");

  double* M;
  double* Q;
  double* F;
  double* I;

  M = new double[4*4];
  Q = new double[4*4];
  F = new double[4];
  I = new double[4];
  
  BSIMCMG("M1", Vd, Vg, Vs, Ve, M, Q, F, I);
 
  return 0;
}

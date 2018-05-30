// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TUnfoldV17Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TUnfold.h"
#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfoldSys.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TUnfoldV17(void *p = 0);
   static void *newArray_TUnfoldV17(Long_t size, void *p);
   static void delete_TUnfoldV17(void *p);
   static void deleteArray_TUnfoldV17(void *p);
   static void destruct_TUnfoldV17(void *p);
   static void streamer_TUnfoldV17(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldV17*)
   {
      ::TUnfoldV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldV17", ::TUnfoldV17::Class_Version(), "UHH2/MTopJet/Unfolding/TUnfold.h", 105,
                  typeid(::TUnfoldV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldV17::Dictionary, isa_proxy, 16,
                  sizeof(::TUnfoldV17) );
      instance.SetNew(&new_TUnfoldV17);
      instance.SetNewArray(&newArray_TUnfoldV17);
      instance.SetDelete(&delete_TUnfoldV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldV17);
      instance.SetDestructor(&destruct_TUnfoldV17);
      instance.SetStreamerFunc(&streamer_TUnfoldV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldBinningV17(void *p = 0);
   static void *newArray_TUnfoldBinningV17(Long_t size, void *p);
   static void delete_TUnfoldBinningV17(void *p);
   static void deleteArray_TUnfoldBinningV17(void *p);
   static void destruct_TUnfoldBinningV17(void *p);
   static void streamer_TUnfoldBinningV17(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldBinningV17*)
   {
      ::TUnfoldBinningV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldBinningV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldBinningV17", ::TUnfoldBinningV17::Class_Version(), "UHH2/MTopJet/Unfolding/TUnfoldBinning.h", 54,
                  typeid(::TUnfoldBinningV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldBinningV17::Dictionary, isa_proxy, 16,
                  sizeof(::TUnfoldBinningV17) );
      instance.SetNew(&new_TUnfoldBinningV17);
      instance.SetNewArray(&newArray_TUnfoldBinningV17);
      instance.SetDelete(&delete_TUnfoldBinningV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldBinningV17);
      instance.SetDestructor(&destruct_TUnfoldBinningV17);
      instance.SetStreamerFunc(&streamer_TUnfoldBinningV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldBinningV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldBinningV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldBinningXMLV17(void *p = 0);
   static void *newArray_TUnfoldBinningXMLV17(Long_t size, void *p);
   static void delete_TUnfoldBinningXMLV17(void *p);
   static void deleteArray_TUnfoldBinningXMLV17(void *p);
   static void destruct_TUnfoldBinningXMLV17(void *p);
   static void streamer_TUnfoldBinningXMLV17(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldBinningXMLV17*)
   {
      ::TUnfoldBinningXMLV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldBinningXMLV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldBinningXMLV17", ::TUnfoldBinningXMLV17::Class_Version(), "TUnfoldBinningXML.h", 56,
                  typeid(::TUnfoldBinningXMLV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldBinningXMLV17::Dictionary, isa_proxy, 16,
                  sizeof(::TUnfoldBinningXMLV17) );
      instance.SetNew(&new_TUnfoldBinningXMLV17);
      instance.SetNewArray(&newArray_TUnfoldBinningXMLV17);
      instance.SetDelete(&delete_TUnfoldBinningXMLV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldBinningXMLV17);
      instance.SetDestructor(&destruct_TUnfoldBinningXMLV17);
      instance.SetStreamerFunc(&streamer_TUnfoldBinningXMLV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldBinningXMLV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldBinningXMLV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldDensityV17(void *p = 0);
   static void *newArray_TUnfoldDensityV17(Long_t size, void *p);
   static void delete_TUnfoldDensityV17(void *p);
   static void deleteArray_TUnfoldDensityV17(void *p);
   static void destruct_TUnfoldDensityV17(void *p);
   static void streamer_TUnfoldDensityV17(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldDensityV17*)
   {
      ::TUnfoldDensityV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldDensityV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldDensityV17", ::TUnfoldDensityV17::Class_Version(), "UHH2/MTopJet/Unfolding/TUnfoldDensity.h", 49,
                  typeid(::TUnfoldDensityV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldDensityV17::Dictionary, isa_proxy, 16,
                  sizeof(::TUnfoldDensityV17) );
      instance.SetNew(&new_TUnfoldDensityV17);
      instance.SetNewArray(&newArray_TUnfoldDensityV17);
      instance.SetDelete(&delete_TUnfoldDensityV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldDensityV17);
      instance.SetDestructor(&destruct_TUnfoldDensityV17);
      instance.SetStreamerFunc(&streamer_TUnfoldDensityV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldDensityV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldDensityV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldSysV17(void *p = 0);
   static void *newArray_TUnfoldSysV17(Long_t size, void *p);
   static void delete_TUnfoldSysV17(void *p);
   static void deleteArray_TUnfoldSysV17(void *p);
   static void destruct_TUnfoldSysV17(void *p);
   static void streamer_TUnfoldSysV17(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldSysV17*)
   {
      ::TUnfoldSysV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldSysV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldSysV17", ::TUnfoldSysV17::Class_Version(), "UHH2/MTopJet/Unfolding/TUnfoldSys.h", 56,
                  typeid(::TUnfoldSysV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldSysV17::Dictionary, isa_proxy, 16,
                  sizeof(::TUnfoldSysV17) );
      instance.SetNew(&new_TUnfoldSysV17);
      instance.SetNewArray(&newArray_TUnfoldSysV17);
      instance.SetDelete(&delete_TUnfoldSysV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldSysV17);
      instance.SetDestructor(&destruct_TUnfoldSysV17);
      instance.SetStreamerFunc(&streamer_TUnfoldSysV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldSysV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldSysV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldV17::Class_Name()
{
   return "TUnfoldV17";
}

//______________________________________________________________________________
const char *TUnfoldV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldBinningV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldBinningV17::Class_Name()
{
   return "TUnfoldBinningV17";
}

//______________________________________________________________________________
const char *TUnfoldBinningV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldBinningV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldBinningV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldBinningV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldBinningXMLV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldBinningXMLV17::Class_Name()
{
   return "TUnfoldBinningXMLV17";
}

//______________________________________________________________________________
const char *TUnfoldBinningXMLV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldBinningXMLV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldBinningXMLV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldBinningXMLV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldDensityV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldDensityV17::Class_Name()
{
   return "TUnfoldDensityV17";
}

//______________________________________________________________________________
const char *TUnfoldDensityV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldDensityV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldDensityV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldDensityV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldSysV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldSysV17::Class_Name()
{
   return "TUnfoldSysV17";
}

//______________________________________________________________________________
const char *TUnfoldSysV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldSysV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldSysV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldSysV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TUnfoldV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldV17.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fA;
      R__b >> fL;
      R__b >> fY;
      R__b >> fVyy;
      R__b >> fBiasScale;
      R__b >> fX0;
      R__b >> fTauSquared;
      fXToHist.Streamer(R__b);
      fHistToX.Streamer(R__b);
      fSumOverY.Streamer(R__b);
      void *ptr_fConstraint = (void*)&fConstraint;
      R__b >> *reinterpret_cast<Int_t*>(ptr_fConstraint);
      void *ptr_fRegMode = (void*)&fRegMode;
      R__b >> *reinterpret_cast<Int_t*>(ptr_fRegMode);
      R__b >> fIgnoredBins;
      R__b >> fEpsMatrix;
      R__b >> fX;
      R__b >> fVxx;
      R__b >> fVxxInv;
      R__b >> fVyyInv;
      R__b >> fAx;
      R__b >> fChi2A;
      R__b >> fLXsquared;
      R__b >> fRhoMax;
      R__b >> fRhoAvg;
      R__b >> fNdf;
      int R__i;
      for (R__i = 0; R__i < 2; R__i++)
         R__b >> fDXDAM[R__i];
      for (R__i = 0; R__i < 2; R__i++)
         R__b >> fDXDAZ[R__i];
      R__b >> fDXDtauSquared;
      R__b >> fDXDY;
      R__b >> fEinv;
      R__b >> fE;
      R__b.CheckByteCount(R__s, R__c, TUnfoldV17::IsA());
   } else {
      R__c = R__b.WriteVersion(TUnfoldV17::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << fA;
      R__b << fL;
      R__b << fY;
      R__b << fVyy;
      R__b << fBiasScale;
      R__b << fX0;
      R__b << fTauSquared;
      fXToHist.Streamer(R__b);
      fHistToX.Streamer(R__b);
      fSumOverY.Streamer(R__b);
      R__b << (Int_t)fConstraint;
      R__b << (Int_t)fRegMode;
      R__b << fIgnoredBins;
      R__b << fEpsMatrix;
      R__b << fX;
      R__b << fVxx;
      R__b << fVxxInv;
      R__b << fVyyInv;
      R__b << fAx;
      R__b << fChi2A;
      R__b << fLXsquared;
      R__b << fRhoMax;
      R__b << fRhoAvg;
      R__b << fNdf;
      int R__i;
      for (R__i = 0; R__i < 2; R__i++)
         R__b << fDXDAM[R__i];
      for (R__i = 0; R__i < 2; R__i++)
         R__b << fDXDAZ[R__i];
      R__b << fDXDtauSquared;
      R__b << fDXDY;
      R__b << fEinv;
      R__b << fE;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldV17(void *p) {
      return  p ? new(p) ::TUnfoldV17 : new ::TUnfoldV17;
   }
   static void *newArray_TUnfoldV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldV17[nElements] : new ::TUnfoldV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldV17(void *p) {
      delete ((::TUnfoldV17*)p);
   }
   static void deleteArray_TUnfoldV17(void *p) {
      delete [] ((::TUnfoldV17*)p);
   }
   static void destruct_TUnfoldV17(void *p) {
      typedef ::TUnfoldV17 current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TUnfoldV17(TBuffer &buf, void *obj) {
      ((::TUnfoldV17*)obj)->::TUnfoldV17::Streamer(buf);
   }
} // end of namespace ROOT for class ::TUnfoldV17

//______________________________________________________________________________
void TUnfoldBinningV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldBinningV17.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> parentNode;
      R__b >> childNode;
      R__b >> nextNode;
      R__b >> prevNode;
      R__b >> fAxisList;
      R__b >> fAxisLabelList;
      R__b >> fHasUnderflow;
      R__b >> fHasOverflow;
      R__b >> fDistributionSize;
      R__b >> fFirstBin;
      R__b >> fLastBin;
      R__b >> fBinFactorFunction;
      R__b >> fBinFactorConstant;
      R__b.CheckByteCount(R__s, R__c, TUnfoldBinningV17::IsA());
   } else {
      R__c = R__b.WriteVersion(TUnfoldBinningV17::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << parentNode;
      R__b << childNode;
      R__b << nextNode;
      R__b << prevNode;
      R__b << fAxisList;
      R__b << fAxisLabelList;
      R__b << fHasUnderflow;
      R__b << fHasOverflow;
      R__b << fDistributionSize;
      R__b << fFirstBin;
      R__b << fLastBin;
      R__b << fBinFactorFunction;
      R__b << fBinFactorConstant;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldBinningV17(void *p) {
      return  p ? new(p) ::TUnfoldBinningV17 : new ::TUnfoldBinningV17;
   }
   static void *newArray_TUnfoldBinningV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldBinningV17[nElements] : new ::TUnfoldBinningV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldBinningV17(void *p) {
      delete ((::TUnfoldBinningV17*)p);
   }
   static void deleteArray_TUnfoldBinningV17(void *p) {
      delete [] ((::TUnfoldBinningV17*)p);
   }
   static void destruct_TUnfoldBinningV17(void *p) {
      typedef ::TUnfoldBinningV17 current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TUnfoldBinningV17(TBuffer &buf, void *obj) {
      ((::TUnfoldBinningV17*)obj)->::TUnfoldBinningV17::Streamer(buf);
   }
} // end of namespace ROOT for class ::TUnfoldBinningV17

//______________________________________________________________________________
void TUnfoldBinningXMLV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldBinningXMLV17.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TUnfoldBinningV17::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TUnfoldBinningXMLV17::IsA());
   } else {
      R__c = R__b.WriteVersion(TUnfoldBinningXMLV17::IsA(), kTRUE);
      TUnfoldBinningV17::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldBinningXMLV17(void *p) {
      return  p ? new(p) ::TUnfoldBinningXMLV17 : new ::TUnfoldBinningXMLV17;
   }
   static void *newArray_TUnfoldBinningXMLV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldBinningXMLV17[nElements] : new ::TUnfoldBinningXMLV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldBinningXMLV17(void *p) {
      delete ((::TUnfoldBinningXMLV17*)p);
   }
   static void deleteArray_TUnfoldBinningXMLV17(void *p) {
      delete [] ((::TUnfoldBinningXMLV17*)p);
   }
   static void destruct_TUnfoldBinningXMLV17(void *p) {
      typedef ::TUnfoldBinningXMLV17 current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TUnfoldBinningXMLV17(TBuffer &buf, void *obj) {
      ((::TUnfoldBinningXMLV17*)obj)->::TUnfoldBinningXMLV17::Streamer(buf);
   }
} // end of namespace ROOT for class ::TUnfoldBinningXMLV17

//______________________________________________________________________________
void TUnfoldDensityV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldDensityV17.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TUnfoldSysV17::Streamer(R__b);
      R__b >> fConstOutputBins;
      R__b >> fConstInputBins;
      R__b >> fOwnedOutputBins;
      R__b >> fOwnedInputBins;
      R__b >> fRegularisationConditions;
      R__b.CheckByteCount(R__s, R__c, TUnfoldDensityV17::IsA());
   } else {
      R__c = R__b.WriteVersion(TUnfoldDensityV17::IsA(), kTRUE);
      TUnfoldSysV17::Streamer(R__b);
      R__b << fConstOutputBins;
      R__b << fConstInputBins;
      R__b << fOwnedOutputBins;
      R__b << fOwnedInputBins;
      R__b << fRegularisationConditions;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldDensityV17(void *p) {
      return  p ? new(p) ::TUnfoldDensityV17 : new ::TUnfoldDensityV17;
   }
   static void *newArray_TUnfoldDensityV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldDensityV17[nElements] : new ::TUnfoldDensityV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldDensityV17(void *p) {
      delete ((::TUnfoldDensityV17*)p);
   }
   static void deleteArray_TUnfoldDensityV17(void *p) {
      delete [] ((::TUnfoldDensityV17*)p);
   }
   static void destruct_TUnfoldDensityV17(void *p) {
      typedef ::TUnfoldDensityV17 current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TUnfoldDensityV17(TBuffer &buf, void *obj) {
      ((::TUnfoldDensityV17*)obj)->::TUnfoldDensityV17::Streamer(buf);
   }
} // end of namespace ROOT for class ::TUnfoldDensityV17

//______________________________________________________________________________
void TUnfoldSysV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldSysV17.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TUnfoldV17::Streamer(R__b);
      R__b >> fDAinRelSq;
      R__b >> fDAinColRelSq;
      R__b >> fAoutside;
      R__b >> fSysIn;
      R__b >> fBgrIn;
      R__b >> fBgrErrUncorrInSq;
      R__b >> fBgrErrScaleIn;
      R__b >> fDtau;
      R__b >> fYData;
      R__b >> fVyyData;
      R__b >> fEmatUncorrX;
      R__b >> fEmatUncorrAx;
      R__b >> fDeltaCorrX;
      R__b >> fDeltaCorrAx;
      R__b >> fDeltaSysTau;
      R__b.CheckByteCount(R__s, R__c, TUnfoldSysV17::IsA());
   } else {
      R__c = R__b.WriteVersion(TUnfoldSysV17::IsA(), kTRUE);
      TUnfoldV17::Streamer(R__b);
      R__b << fDAinRelSq;
      R__b << fDAinColRelSq;
      R__b << fAoutside;
      R__b << fSysIn;
      R__b << fBgrIn;
      R__b << fBgrErrUncorrInSq;
      R__b << fBgrErrScaleIn;
      R__b << fDtau;
      R__b << fYData;
      R__b << fVyyData;
      R__b << fEmatUncorrX;
      R__b << fEmatUncorrAx;
      R__b << fDeltaCorrX;
      R__b << fDeltaCorrAx;
      R__b << fDeltaSysTau;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldSysV17(void *p) {
      return  p ? new(p) ::TUnfoldSysV17 : new ::TUnfoldSysV17;
   }
   static void *newArray_TUnfoldSysV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldSysV17[nElements] : new ::TUnfoldSysV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldSysV17(void *p) {
      delete ((::TUnfoldSysV17*)p);
   }
   static void deleteArray_TUnfoldSysV17(void *p) {
      delete [] ((::TUnfoldSysV17*)p);
   }
   static void destruct_TUnfoldSysV17(void *p) {
      typedef ::TUnfoldSysV17 current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TUnfoldSysV17(TBuffer &buf, void *obj) {
      ((::TUnfoldSysV17*)obj)->::TUnfoldSysV17::Streamer(buf);
   }
} // end of namespace ROOT for class ::TUnfoldSysV17

namespace {
  void TriggerDictionaryInitialization_TUnfoldV17Dict_Impl() {
    static const char* headers[] = {
"TUnfold.h",
"TUnfoldBinning.h",
"TUnfoldBinningXML.h",
"TUnfoldDensity.h",
"TUnfoldSys.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed5/include",
"/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/Unfolding/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TUnfoldV17Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Unfolding with support for L-curve analysis)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfold.h")))  TUnfoldV17;
class __attribute__((annotate(R"ATTRDUMP(Complex binning schemes for TUnfoldDensity)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfoldBinning.h")))  TUnfoldBinningV17;
class __attribute__((annotate(R"ATTRDUMP(Complex binning schemes for TUnfoldDensity)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfoldBinningXML.h")))  TUnfoldBinningXMLV17;
class __attribute__((annotate(R"ATTRDUMP(Unfolding with density regularisation)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfoldDensity.h")))  TUnfoldDensityV17;
class __attribute__((annotate(R"ATTRDUMP(Unfolding with support for systematic error propagation)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfoldDensity.h")))  TUnfoldSysV17;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TUnfoldV17Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TUnfold.h"
#include "TUnfoldBinning.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include "TUnfoldSys.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TUnfoldBinningV17", payloadCode, "@",
"TUnfoldBinningXMLV17", payloadCode, "@",
"TUnfoldDensityV17", payloadCode, "@",
"TUnfoldSysV17", payloadCode, "@",
"TUnfoldV17", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TUnfoldV17Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TUnfoldV17Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TUnfoldV17Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TUnfoldV17Dict() {
  TriggerDictionaryInitialization_TUnfoldV17Dict_Impl();
}

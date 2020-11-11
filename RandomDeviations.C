// program to measure random deviations

#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TH3.h"
#include "TH2.h"
#include "Math/MixMaxEngine.h"
#ifdef HAVE_RANLUX_ENGINE
#include "Math/RanLuxEngine.h"
#endif
#include <vector>
#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TClass.h"

#include <cstdint>

// needed if testing other engines than mixmax
// #if ROOT_VERSION_CODE >= ROOT_VERSION(6,7,7)
// #include "TestEngines.h"
// #endif

//const int N = 256;

int iRBegin = 1;   // must be 1 for MIXMAX , 0 for the others!!!
/// define a wrapper class for testing the random engines
namespace Test {

   
   template<class Engine>
   class TestEngine : public Engine {

   public:

      typedef typename Engine::StateInt_t StateInt_t;
      int Counter() const { return 0; }//Engine::Counter(); }

      void GetState(std::vector<StateInt_t> & state) { //const {
         Engine::GetState(state);
      }
      void SetState(const std::vector<StateInt_t> & state) {
         Engine::SetState(state);
      }
      //void SetCounter(int val) { fCount624 = val; }
      
      //static void SetSkipNumber(int /*nskip */) {  }

      // void MinCounter() {
      //    return 0; 
      // }
      // void MaxCounter() {
      //    return Engine::Size();
      // }
   };

   template<class Engine>
   class TestRanLux : public Engine {

   public:

      typedef typename Engine::StateInt_t StateInt_t;
      int Counter() const { return 0; }//Engine::Counter(); }

      void GetState(std::vector<StateInt_t> & state) { //const {
         // needs 96
         if (state.size() < 24) state.resize(24);
         std::vector<StateInt_t> tmp( Engine::Size());
         Engine::GetState(tmp);
         std::copy(tmp.begin(), tmp.begin()+24, state.begin());
      }
      void SetState(const std::vector<StateInt_t> & state) {
         std::vector<StateInt_t> tmp( state.begin(), state.end());
         for (int i = 0; i < 3; ++i) 
            tmp.insert(tmp.end(), state.begin(), state.end());
         
         Engine::SetState(tmp);
      }

      static size_t Size()  { return 24; } 
      //void SetCounter(int val) { fCount624 = val; }
      
      //static void SetSkipNumber(int /*nskip */) {  }

      // void MinCounter() {
      //    return 0; 
      // }
      // void MaxCounter() {
      //    return Engine::Size();
      // }
   };
//#endif

}

unsigned long M61 = 2305843009213693951ULL;  // 2^61-1


int debug = 0;

int nstate = 0;

int gIndexChange = 0;

int firstGen = 1;

using namespace ROOT::Math;

int normtype = 1; // 0 : Eucledian, 1 Chebyshev (max norm)

int engineSubType = 0; 

enum EDiffType {
   kOld,
   kSingle1,    // change one element in the matrix randomly by a distance of +/- 1
   kSingleN,   // change one element in the matrix randomly by a distance of +/- N= ndiffSize
   kFixed1,    // change one element in the matrix at position j = changeElement  by a distance of +/- 1
   kFixed1All,    // change one element in the matrix at position j = changeElement (varyin all time) by a distance of +/- 1
   kFixedNAll,    // change one element in the matrix at position j = changeElement (varyin all time) by a distance of +/- N
   kFixedN,    // change one element in the matrix at position j = changeElement  by a distance of +/- N= ndiffSize
   kDiff1,     // change randomly in the matrix one or more by keeping distance <= +/- 1
   kDiffN,     // change randomly in the matrix one or more by keeping distance <= +/- N = ndiffSize
   kUnique,    // set a given different state
};


EDiffType diffType =  kDiff1; // kFixed1All;
//unsigned long ndiffSize = 1;   // size of the difference (used in kDiffN)
unsigned long ndiffSize = 1UL << 54;   // size of the difference (used in kDiffN)
int changeElement = 0;


bool excludeChangeInElement0 = true;  // do not makke trajectory which have changes only in element 0. 

//#define USEDIFF

//#define USEDIFF_FIXED
//#define  USE_UNIQUE_STATE
bool generateInitialState = false;  // generate a random initial state
bool zeroState = false; // use a zero initial state

#ifdef USE_UNIQUE_STATE
//#include "special_mixmax_state.h"
const int N = 256;
#include "just_two_vectors.c"
std::vector<uint64_t> state_special1(Y1,Y1+256);
std::vector<uint64_t> state_special2(Y2,Y2+256);
#endif

std::vector< std::vector<int> > allFoundCombinations;
std::vector<int> currentCombination;
std::vector<int> currentCombination2;
int nDuplicateFound = 0;
bool allCombinationsFound = false; 

// Functions to make the initial state
template<class IntType> 
int MakeDiffStateOld(std::vector<IntType> & s2, TRandom& rndm)
{ 
   
   int nchange = 0; 
      
   int nstate = s2.size(); 
   std::vector<int> changedValues(nstate); 
   //int nchange = rndm.Integer(nstate)+1;
   //nchange = 1;
   // change a random number of elements 
   nchange = rndm.Integer(nstate)+1;   
   for (int i = 0; i < nchange; ++i) {
      int ichange = 0;
      int jj = 0;
      while (jj<1E9) { 
         ichange = rndm.Integer(nstate);
         if (changedValues[ichange] == 0) {
            changedValues[ichange] = 1; 
            break;
         }
         jj++;
      }

      int diff = 2*rndm.Integer(2)-1;  // generate between 1 and -1 
      if (s2[ichange] <= 1 ) diff = 1; 
      if (s2[ichange] >= std::numeric_limits<uint64_t>::max() ) diff = -1;
      IntType sold = s2[ichange];
      s2[ichange] += diff;
         
      if (debug>1) std::cout << "changing element " << ichange  << " from " << sold << " to " << s2[ichange] << " by " << diff << std::endl;
   }
   return nchange; 
}
   

// make the second state
// take a random change betwwwn a shere of +/- r

template<class IntType> 
int MakeDiffState(std::vector<IntType> & s2, TRandom& rndm, int64_t r = 1 )
{
#ifdef OLD // old code before Nov 2018
   int nstate = s2.size(); 
   int nchange = 0;
   int rmax = 2*r +1;
   while (nchange == 0) {
      ntry++; 
      int istart = rndm.Integer(nstate);
      if (excludeChangeInElement0 && istart == 0) continue; 
      for (int i = 0; i < nstate; ++i) {
         int ichange = istart + i;
         // to make sure we start randomly in the state
         if (ichange >= nstate) ichange -= nstate; 
         //int diff = rndm.Integer(3)-1;  // generate -1,0 and 1
         int jrndm =  rndm.Integer(rmax);
         int diff = 2*jrndm - r;  // generate -r,..0 and 1,2,3,4
         // if (diff == -1 && s2[ichange] == 0) diff = 1; 
         // if (s2[ichange] <= 1 ) diff = 1; 
         // if (s2[ichange] >= std::numeric_limits<uint64_t>::max() ) diff = -1;
         int64_t sold = s2[ichange] & M61; 
         // should take the modulus of p
         int64_t snew = sold + diff;
         if (diff < 0 && (-diff) > sold ) snew = M61 - (-diff -sold);
         if (diff > 0 && (sold + diff) >= M61 ) snew = (sold+diff)-M61; 
         s2[ichange] = snew;
         
         if (diff !=0) nchange++;

         if (debug>1) std::cout << "changing element " << ichange  << " from " << sold << " to " << s2[ichange] << " by " << diff << std::endl;
      }
   }
#endif
   // new code changing random combinations
   // assume r = 1
   if (r != 1) return 0; 
   int nstate = s2.size(); 
   int nchange = 0;
   int rmax = 2*r +1;

   std::vector<int> currentCombination(nstate );
   std::vector<int> currentCombination2(nstate );
   int ntry = 0; 
   // for each coordinate generate a diff = {-1,0,1} 
   while (nchange == 0 )  {
      ntry++;
      int istart =  0; //excludeChangeInElement0 ? 1 : 0; 
      for (int ichange = istart; ichange < nstate; ++ichange) {
         int jrndm =  rndm.Integer(3);
         int diff = jrndm - 1;  // generate -r,..0 and 1,2,3,4
         currentCombination[ichange] = diff;
      }
      
      bool isCombExisting = false; 
      // check if a similar combinations exist
      for (size_t kcomb = 0; kcomb < allFoundCombinations.size(); ++kcomb) {
         bool isNewComb = false;
         auto & comb =  allFoundCombinations[kcomb];
         //std::cout << "comb " << kcomb << std::endl;
         for (size_t k = 0 ; k < comb.size(); ++k) {
            if (comb[k] != currentCombination[k]) {
               isNewComb = true;
               break; 
            }
         }
         if (!isNewComb) {
            isCombExisting = true;
            break;
         }
      }
      if (isCombExisting) {
         nDuplicateFound++; 
         if (debug) std::cout << "existing comb found "   << nDuplicateFound << " and found " << allFoundCombinations.size() << std::endl;
         if (ntry > 10000) {
            // no sense continue searching for newe combination, stop
            allCombinationsFound = true; 
            return 0; 
         }
         continue;
      }
      

      // here if I have found new combinations
      // store it and it symmetrical
      if (debug) std::cout << "new combination found " << std::endl;
      nchange  = 1;

//    if (diff == 0) continue;
      for (int ichange = istart; ichange < nstate; ++ichange) {
         int diff = currentCombination[ichange];
         currentCombination2[ichange] = -diff;  // symmetric combination
         int64_t sold = s2[ichange] & M61; 
         // should take the modulus of p
         int64_t snew = sold + diff;
         if (diff < 0 && (-diff) > sold ) snew = M61 - (-diff -sold);
         if (diff > 0 && (sold + diff) >= M61 ) snew = (sold+diff)-M61; 
         s2[ichange] = snew;

         //if (diff !=0) nchange++;

         if (debug>1) std::cout << "changing element " << ichange  << " from " << sold << " to " << s2[ichange] << " by " << diff << std::endl;
      }
      allFoundCombinations.push_back(currentCombination);
      allFoundCombinations.push_back(currentCombination2);
   }
   if (debug>1) std::cout << "changing " << nchange << " elements " << std::endl;
   return nchange;
}


/// change only one element by +/- 1
template<class IntType> 
int MakeDiffStateFixed(std::vector<IntType> & s2, TRandom& rndm, int ichange = -1, int64_t r = 1)
{

   if (ichange < 0) 
      ichange = rndm.Integer(nstate);

   
   int64_t diff = 0L;
   if (r == 1) {
      // get diff randomly as +/- 1
      if ( rndm.Rndm() < 0.5) 
         diff = -1;
      else
         diff = 1; 
      //std::cout << "GenDifference is " << diff << std::endl;
   }
   else if (r < 0) {
      while (diff != 0) {
         r = -r;
         diff = static_cast<int64_t>( rndm.Integer(2*r+1) ) -r;
      }
   }
   else
      diff = r;

   // mask the top 3 bits
   int64_t sold = s2[ichange] & M61;
   int64_t snew = sold + diff;
   if (diff < 0 && (-diff) > sold ) snew = M61 - (-diff -sold);
   if (diff > 0 && (sold + diff) >= M61 ) snew = (sold+diff)-M61; 
   s2[ichange] = snew;

   if (debug>0) std::cout << "changing element " << ichange  << " from " << sold << " to " << s2[ichange] << " by " << diff << std::endl;
   return ichange;
}

template<class IntType> 
int MakeDiffSpecialState(std::vector<IntType> & s1, std::vector<IntType> & s2)
{
#ifdef USE_UNIQUE_STATE
   s1 =   state_special1; 
   s2 =   state_special2;
#endif
   return -1;
}

// FUNCTION COMPUTE DEVIATIONS
template<class Type> 
Type  ComputeDistance(const std::vector<Type> & v1, const std::vector<Type> & v2, uint64_t maxType, int & idistloc, std::vector<Type> & distVec) { 
      double distSum = 0;
      double maxProjDist = 0; 
      Type result =  -1;
      idistloc = -1;
      int nstate = v1.size();
      // to make sure we don't bias distances location we start randomly 
      int istart = gRandom->Integer(nstate);
      std::vector<Type> d(nstate);

      for (int j  = 0; j < nstate; ++j) {
         int i = istart + j;
         if (i>= nstate) i -= nstate;
         // skip computing distance from first and last element
         //if (i == 0 || i == nstate-1 ) continue;
         Type dist = 0;

         if (std::is_floating_point<Type>::value) {
            double d12 = std::abs(double(v1[i]-v2[i]) );
            double ddd = std::min(d12, 1.-d12);
            dist = ddd;
            d[i] = ddd;
         }
         else {  // here Type is an integer
            Type d12 = (v1[i] > v2[i] ) ? v1[i]-v2[i] : v2[i]-v1[i];
            std::cout << d12 ;
            uint64_t dist = std::min( (uint64_t) d12, (uint64_t) ( maxType - d12) );
            std::cout << "  - > " << dist << std::endl;
            d[i] = dist;
         }
         if (d[i] > maxProjDist) {
            idistloc = i;
            maxProjDist = d[i];
            d[i] = dist;
         }

         // compute using the metric
         if (normtype == 0) {
            distSum += dist*dist;
            if (debug > 1) std::cout << " i " << i << "  " << v1[i] << "  " << v2[i] << " diff " << dist << " cumdiff " << distSum << std::endl;
         }

         std::cout << "computed distance " << dist << "  " << d[i] << std::endl;
         
         if (normtype == 1) {
            if (debug > 1)
               std::cout << " i " << i << "  " << v1[i] << "  " << v2[i] << " distance is " << d[i] << " log2 " << TMath::Log2(d[i]) << " log10 " << TMath::Log10(d[i]) << std::endl;
         }
      }
      if (normtype == 0)  result = sqrt(distSum);
      else if (normtype == 1) result = maxProjDist; 

      if (debug>0) std::cout << "ComputeDistance: -------> result is " << result << " log2 " << TMath::Log2(result) << " log10 " << TMath::Log10(result) << std::endl;

      // normalize projected distances
      std::vector<double> dnorm(nstate); 
      if (result > 0) { 
         for (int i = 0; i < nstate; ++i)
            dnorm[i] = double(d[i])/double(result); 
      }
      if (debug > 0) {
         std::cout << "--------------------------------------> proj distances = ";
         for (int i = 0; i < nstate; ++i)
            std::cout << dnorm[i] << " , ";
         std::cout << std::endl;
      }
      distVec = d; 
      return result;
}


template<class Engine, class Type, class IntType> 
std::vector<Type> ComputeDeviations(int niter, int & iChangeIndex, std::vector<int> & locMaxPos, std::vector<std::vector<Type> > &distResult) { 

   // generate a state of random integer values

   // define the skipping number
   //Engine::SetSkipNumber(0);
   //Engine::SetFirstReturnElement(1);

   // to store the values
   distResult.clear();
   distResult.resize(niter+1);
 
   nstate = Engine::Size();
   const uint32_t imax32 = std::numeric_limits<uint32_t>::max();
   uint64_t maxType = Engine::MaxInt(); 
   // in case of runlux
   //maxType = std::pow(2,23); 

   if (debug) std::cout << "Running for size " << nstate << std::endl;

   std::vector<IntType> s1(nstate); 
   auto s2 = s1;

   Engine r1;
   Engine r2;

   //Engine::SetSpecialNumber(1.E12);
   
   TRandom2 rndm(0);

   if (generateInitialState) { 
   
      //MixMaxEngine r0;
   //r0.SetSeed( rndm.Integer(2000000000) );
      for (int i = 0; i < nstate; ++i) {
         uint64_t val1,val2=0;
         if (maxType > imax32) { 
            val1 = (uint64_t) rndm.Integer(imax32) << 32 | rndm.Integer(imax32);
            val2 = (uint64_t) rndm.Integer(imax32) << 32 | rndm.Integer(imax32);
         }
         else {
            val1 = rndm.Integer(maxType);
            val2 = rndm.Integer(maxType); 
         }
#ifdef SPECIAL_INITIAL_STATE
         //  do initial state with large bits (e.g 7)
         int imax = TMath::Power(2,61-54-1);
         uint64_t val1 = (uint64_t) rndm.Integer(imax) <<  54;
         uint64_t val2 = (uint64_t) rndm.Integer(imax) <<  54;
#endif
         s1[i] = val1;
         s2[i] = val2; 
      }
   }
   
   else if (!zeroState)  {
      // leave generator to get its random state

      // if generate random seed
      int seed = 1; 
      seed = rndm.Integer(imax32); 
      
      r1.SetSeed(seed);
      r2.SetSeed(seed);
      r1.GetState(s1);
      r2.GetState(s2);
   }

   if (debug > 1) std::cout << " GENERATED STATE " << s1[0] << " " << s1[1] << std::endl;;
   s2 = s1;
   int nchange = -1;
   switch (diffType) {
   case kOld:
      nchange = MakeDiffStateOld(s2,rndm);
      break;
   case kSingle1:
      nchange = MakeDiffStateFixed(s2,rndm,-1,1);
      break;
   case kSingleN:
      nchange = MakeDiffStateFixed(s2,rndm,-1,ndiffSize);
      break;
   case kFixed1:
      nchange = MakeDiffStateFixed(s2,rndm,changeElement,1);
      break;
   case kFixed1All:
      changeElement = iChangeIndex;
      nchange = MakeDiffStateFixed(s2,rndm,changeElement,1);
      break;
   case kFixedNAll:
      changeElement = iChangeIndex;
      nchange = MakeDiffStateFixed(s2,rndm,changeElement,ndiffSize);
      break;
   case kFixedN:
      nchange = MakeDiffStateFixed(s2,rndm,changeElement,ndiffSize);
      break;
   case kUnique:
      MakeDiffSpecialState(s1,s2);
      break;
   case kDiff1:
      MakeDiffState(s2,rndm,1);
      break;
   case kDiffN:
      MakeDiffState(s2,rndm,ndiffSize);
      break;
   default:
      MakeDiffState(s2,rndm);
   }
   
   iChangeIndex = nchange;
   
   // make the second state with a certain distance from the first one

    //if (debug) cout << "get counter " << r1.Counter() << endl;

   if (debug >=3)  { 
      std::cout << "S1" << std::endl;
      for (auto & e: s1) 
         std::cout << e << " , ";
      std::cout << std::endl;
      std::cout << "S2" << std::endl;
      for (auto & e: s2) 
         std::cout << e << " , ";
      std::cout << std::endl;
   }


   
   r1.SetState(s1);
   r2.SetState(s2);

   //if (debug) cout << "get counter " << r1.Counter() << endl;
   
   // now compute the distance

   // Type is the type used for the distance
   // IntType is the type used for the state
   std::vector<Type> v1(nstate); 
   std::vector<Type> v2(nstate); 

   // let's copy s1 in v1
   std::cout << "Get distance from initial state (no generation) - maxtype " << maxType << std::endl;
   for (int i = 0; i < nstate; ++i) {
      if (std::is_floating_point<Type>::value) {
         v1[i] = s1[i]/Type(maxType);
         v2[i] = s2[i]/Type(maxType);
      }
      else {
         v1[i] = s1[i] ;
         v2[i] = s2[i];
      }
      //std::cout << " i " << v1[i] << "  " << s1[i] << " -- " << v2[i] << "  " << s2[i] << std::endl;
   }

   std::vector<Type> result(niter+1);
   locMaxPos.clear();
   locMaxPos.reserve(niter+1);

   int idistloc = -1;
   Type dist0 = ComputeDistance(v1,v2,maxType,idistloc,distResult[0]);
   locMaxPos.push_back(idistloc); 
   
   if (debug)
      std::cout << "Initial distance is " << dist0 << " log2 " << TMath::Log2(dist0) << " log10 " << TMath::Log10(dist0) << std::endl;

   result[0] = dist0; 
   // if (std::is_floating_point<Type>::value) 
   //    initialLogDistance = TMath::Log10((double) dist0);
   // else
   //    initialLogDistance = TMath::Log2((double) dist0);

   // make iterations advancing the generators 
   
   for (int iter = 1; iter <= niter; ++iter) { 
      //if (iter >= 29 && iter <= 32) debug = 3;
      //else debug = 0; 

      if (debug > 1) std::cout << "\nPerforming iteration..... " << iter << std::endl; 

      // print state before
      if (debug >= 3) { 
         std::vector<IntType> s;
         std::cout << "Generator state before generating" << std::endl;
         std::cout << "counter is " << r1.Counter() << std::endl;
         r1.GetState(s);
         for (auto & e: s) 
            std::cout << e << " , ";
         std::cout << std::endl;
         std::cout << "State for second generator " << std::endl;
         r2.GetState(s);
         for (auto & e: s) 
            std::cout << e << " , ";
         std::cout << std::endl << std::endl;

      }

      //r1.Iterate();
      //r2.Iterate();
      // r1.SetCounter(nstate+1);
      // r2.SetCounter(nstate+1);

      int i0 = iRBegin;


      if (std::is_floating_point<Type>::value) { 
         for (int i = i0; i < nstate; ++i) {
            //for (int i = 1; i < nstate; ++i) {  // mixmax starts from 1 
            v1[i] = r1(); 
            v2[i] = r2();
            if (debug >=3) std::cout <<  i << " rndm " << v1[i] << "  " << v2[i] << std::endl;
         }
      }
      else {
         // for mix max iterations first element is not returned

         for (int i = i0; i < nstate; ++i) { 
            v1[i] = r1.IntRndm() ;
            v2[i] = r2.IntRndm(); 
            if (debug >=3) std::cout <<  i << " rndm " << v1[i] << "  " << v2[i] << " counter   " << r1.Counter() << "  " << r2.Counter() << std::endl;
            // if (i < 3) {
            //    // exclude first 2 elements
            //    v1[i] = 0;
            //    v2[i] = 0;
            // }
         }
         // get state 
         r1.GetState(s1);
         r2.GetState(s2);
         v1[0] =  s1[0];
         v2[0] =  s2[0];
         if (debug >=3) std::cout << " rndm  for i = 0 " << v1[0] << "  " << v2[0] << " counter   " << r1.Counter() << "  " << r2.Counter() << std::endl;
      }
      
      if (debug >= 3) { 
         std::cout << "Generator state after generating " << std::endl;
         std::cout << "counter is " << r1.Counter() << std::endl;
         std::vector<IntType> s;
         r1.GetState(s);
         for (auto & e: s) 
            std::cout << e << "  ";
         std::cout << std::endl;
         std::cout << "State for second generator " << std::endl;
         r2.GetState(s);
         for (auto & e: s) 
            std::cout << e << "  ";
         std::cout << std::endl << std::endl;
      }

      // we have at 0 the initial distance (typically 1)
      int idistloc = 0; 
      result[iter] = ComputeDistance(v1,v2,maxType,idistloc,distResult[iter]);
      locMaxPos.push_back(idistloc);

      if (debug) std::cout << "iteration " << iter << " result is " << result[iter] << "\t log2 " << TMath::Log2(result[iter])
                           << "\t log10 " << TMath::Log10(result[iter]) << " location " << locMaxPos.back() << std::endl;
   }
   return result; 

}

template<class Engine, class Type, class IntType>
void GetDeviations(int n = 10, int ntrial = 1000) {



   Engine etest;
   auto eng_class = TClass::GetClass(typeid(etest));
   TString engineName = "Unknown_Engine";
   if (eng_class) 
      engineName = eng_class->GetName(); 
   std::cout << "Running for " << engineName << " size " << Engine::Size() << std::endl;

   //if (engineName.Contains("MixMax")) generateInitialState = true; 
   
   TH2 * hresult = 0;
   TProfile * hprof = 0;
   /// 2D HISTOGRAM AND PROFILE OF ITERATION VS MAX DISTANCE 
   if (std::is_floating_point<Type>::value) {
      double min = TMath::Log10(std::numeric_limits<Type>::epsilon() ) -2;
      double max = 2; 
      hresult = new TH2D("result","result",n+1,-0.5,n+0.5,100,min,max);
      hprof = new TProfile("profResult","result",n+1,-0.5,n+0.5,min,max);
   }
   else {
      double min = -.5;
      double max = 60.5;//TMath::Log2(std::numeric_limits<Type>::max() ) + 1;
      int nbins = 2*int(max-min); 
      hresult = new TH2D("result","result",n+1,-0.5,n+0.5,nbins,min,max);
      hprof = new TProfile("profResult","result",n+1,-0.5,n+0.5,min,max,"s");
   }

   
   TH3 * hresult3 = 0;
   TH3 * hresult4 = 0;
   int n3 =  Engine::Size();
// #ifdef USEDIFF_FIXED
//    n3 = Engine::Size();
// #endif
// #ifdef USEDIFF   
//    n3 = Engine::Size();
// #endif
   if (n3 > 0)  {
      /// 3D HISTOGRAM OF NUMBER OF ITERATIONS VS NAX DISTANCE VS POSITION IN MATRIX OF MAXDIFF 
      if (std::is_floating_point<Type>::value) {
         double min = TMath::Log10(std::numeric_limits<Type>::epsilon() ) -2;
         double max = 2; 
         hresult3 = new TH3D("result3","result3",n+1,-0.5,n+0.5,100,min,max,n3,0,n3);
      }
      else {
         double min = -.5;
         double max = 60.5; //TMath::Log2(std::numeric_limits<Type>::max() ) + 1;
         int nbins = int(max-min);
         hresult3 = new TH3D("result_maxdiff","maxdiff results",n+1,-0.5,n+0.5,nbins,min,max,n3,0,n3);
         hresult4 = new TH3D("result_alldist","alldist results",n+1,-0.5,n+0.5,nbins,min,max,n3,0,n3);
      }
   }


 

   
   int nefftrials  = 0;

   int jstart = (excludeChangeInElement0 && (diffType == kFixed1All || diffType == kFixedNAll) ) ? 1 : 0;

   allFoundCombinations.reserve(2*ntrial);
   std::vector<int> currentCombination = std::vector<int>(Engine::Size() );
   std::vector<int> currentCombination2 = std::vector<int>(Engine::Size() );
   if (excludeChangeInElement0) {
      currentCombination[0] = 1;
      currentCombination2[0] = -1;
      allFoundCombinations.push_back(currentCombination); 
      allFoundCombinations.push_back(currentCombination2); 
   }
   
   for (int j = jstart; j < ntrial; ++j) {
      //for (int j = 5; j < 6; ++j) {
      nefftrials++;


      if (debug) std::cout << "Computing deviations - trial test " << j << std::endl;
      else if ( (j % (ntrial/10)) == 0) std::cout << "Computing deviations - trial test " << j << std::endl;
      int ichange = gIndexChange;
      if (diffType == kFixed1All || diffType == kFixedNAll) {
         ichange = j;
         if (j >=  Engine::Size()) break;
      }
      std::vector<int> locMaxPos(1);
      std::vector<std::vector<Type> > distResult; 
      auto r = ComputeDeviations<Engine,Type,IntType>(n, ichange,locMaxPos,distResult);
      if (allCombinationsFound) break; 
      assert((int) r.size() >= n+1);

      // fill histograms 
      for (int i = 0; i < n+1; ++i) {
         double y = 0;
         if (std::is_floating_point<Type>::value) 
            y = (r[i] > 0) ? TMath::Log10(r[i]) : -1; 
         else 
            y = (r[i] > 0) ? TMath::Log2(r[i]) : -1; 

         hresult->Fill( i, y );
         hprof->Fill( i, y );
         if (hresult3) {
            // USE INITIAL POSITION INSTEAD OF FINAL ONE FOR RECORDING MAX DEVIATIONS
            // THIS MIGHT NOT WORK FOR ALL TYPES of diff selected (BE CAREFUL)
            if (ichange == j)
               hresult3->Fill( i, y, ichange );
            else
               hresult3->Fill( i, y, locMaxPos[i] );
         }
         //std::cout << "dist result " << distResult.size() << std::endl;
         if (hresult4 && distResult.size() > i && distResult[i].size() > 0) {
            //std::cout << "iter " << i << " ... " << distResult[i].size() << std::endl;
            // this histogram plots all distances (we make an average?)
            for (size_t k = 0; k < distResult[i].size();  ++k) {
               //std::cout << "{ " << k << " , " << TMath::Log2(distResult[i][k]) << " } ";
               hresult4->Fill( i, TMath::Log2(distResult[i][k]), k );
            }
            //std::cout << std::endl;
         }
      }
   }



   
   TString title = TString::Format("Deviations for %s - N=%d;Number of Interations; log2(Deviation)",engineName.Data(), nstate);
   hresult->SetTitle(title);
   hprof->SetTitle(title);

   //if (std::is_floating_point<Type>::value) 
      

// #ifdef OLD   
//    // look at results
//    std::vector<double> x(r.size());   
//    std::vector<double> y(r.size()); 
//    for (int i = 0; i < n; ++i) {
//       x[i] = i+1;
//       if (std::is_floating_point<Type>::value) 
//          y[i] = TMath::Log10(r[i]);
//       else 
//          y[i] = TMath::Log2(r[i]);
//    }

//    TGraph * g = new TGraph(n, x.data(), y.data() );
//    g->SetPoint(n+1,0,0);
//    g->SetMarkerStyle(20);
//    g->Draw("AP");
// #endif

   gStyle->SetOptStat(0);
   gStyle->SetErrorX(0);
   //hresult->Draw("candle");
   hresult->Draw();
   TString fprefix;
   fprefix = TString::Format("Deviations_trials%d_%s_N%d",nefftrials,engineName.Data(),nstate);
   if (engineSubType != 0) fprefix = TString::Format("Deviations_trials%d_%s_N%d-%d",nefftrials,engineName.Data(),nstate,engineSubType);
   gPad->SaveAs(TString::Format("%s_1.pdf",fprefix.Data()));
   auto c2 = new TCanvas("c2","c2",700,700);
   //auto prof = hresult->ProfileX("profile",1, hresult->GetNbinsY(),"s");
  
   // std::cout << "Init distance  (log) = " << initialLogDistance << std::endl;

   
   auto g = new TGraphErrors(hprof);
//    for (int i = 0; i < g->GetN(); ++i) { 
//       g->SetPoint(i,i, hprof->GetBinContent(i+1));
// //      g->SetPoint(i,i, hprof->GetBinContent(i+1));
//    }
   // g->SetPoint(g->GetN(),0, initialLogDistance);
   

   g->SetMarkerStyle(20);
   std::cout << "style : " << g->GetMarkerStyle() << std::endl;
   g->GetXaxis()->SetRangeUser(0,g->GetN());
   g->GetYaxis()->SetRangeUser(-0.1,64);
   g->SetMarkerStyle(20);
   std::cout << "style : " << g->GetMarkerStyle() << std::endl;
   //gStyle->SetOptFit(11111);
   auto f1 = new TF1("f1","pol1",0,5);
   f1->SetNpx(1000);
   f1->SetRange(0,3.1);
   g->Draw("AP");
   std::cout << "style : " << g->GetMarkerStyle() << std::endl;
   //g->Fit(f1,"R");
   c2->Update();  // force the paint
   c2->SaveAs(TString::Format("%s_2.pdf",fprefix.Data()));
   //c2->SaveAs("c2.pdf");

   auto c3 = new TCanvas("c3","c3",800,800);

   std::vector<TH1*> ph(n+1);
   // start from second bin (iteration 1) 
   for (int i = 0 ; i < n+1; ++i) { 
      int ibin = i+1;
      ph[i] = hresult->ProjectionY(TString::Format("h%d",i),ibin,ibin);
   }

   int ncx = 2;
   int ncy = 2;
   if (n>9)     { ncx=4; ncy = 4; }
   else if (n > 6) { ncx=3; ncy = 3; }
   else if (n > 4) { ncx=2; ncy = 3; }
 

   c3->Divide(ncx,ncy);
//    if (n > 9) c3->Divide(4,4);
//    else if (n > 4) c3->Divide(3,3);
//    else if (n>2) c3->Divide(2,2);
//    else if (n>1) c3->Divide(1,2);
   
   for (int i = 1; i < n+1; ++i) {
      c3->cd(i);
      gPad->SetLogy(true);
      ph[i]->SetTitle(TString::Format("Distance for Iteration %d",i));
      //ph[i]->SetFillColor(kBlue-3);
      ph[i]->Draw();
   }
   //c3->SaveAs("c3.pdf");
   c3->SaveAs(TString::Format("%s_3.pdf",fprefix.Data()));


   if (hresult3) {
      
      auto c4 = new TCanvas("c4","c4",800,800);
      if (Engine::Size() > 62 )
         c4->Divide(ncx,ncy);
      else 
         c4->Divide(ncy,ncx);
      
      for (int i = 0; i < n+1; ++i) {
         c4->cd(i+1);
         hresult3->GetXaxis()->SetRange(i+1,i+1); 
         auto ph2 = hresult3->Project3D(TString::Format("h%d_ZY",i)); 
         ph2->SetTitle(TString::Format("Distance vs state maxdiff pos  for Iteration %d",i));
         ph2->Draw("COLZ");
      }
      auto c5 = new TCanvas("c5","c5",800,800);
      c5->Divide(ncx,ncy); 
      for (int i = 0; i < n+1; ++i) {
         auto p = c5->cd(i+1);
         p->SetLogy(1);
         hresult3->GetXaxis()->SetRange(i+1,i+1); 
         auto phz = hresult3->Project3D(TString::Format("h%d_Z",i)); 
         phz->SetTitle(TString::Format("State maxdiff pos for Iteration %d",i));
         phz->Draw();
      }
   }

   // plot histogram of all distances vs state size 
   if (hresult4) {
      std::cout << "create histo for all distances" << std::endl;
      
      auto c6 = new TCanvas("c6","c6",900,900);
      if (Engine::Size() > 62 )
         c6->Divide(ncx,ncy);
      else 
         c6->Divide(ncy,ncx);
 
      for (int i = 0; i < n+1; ++i) {
         std::cout << "plotting in  " << i+1 << std::endl;
         c6->cd(i+1);
         hresult4->GetXaxis()->SetRange(i+1,i+1); 
         auto ph2 = hresult4->Project3D(TString::Format("ph%d_ZY",i)); 
         ph2->SetTitle(TString::Format("Distance vs state pos  for Iteration %d",i));
         ph2->Draw("COLZ");
      }
   }

   


   TFile * file = TFile::Open(TString::Format("%s_test.root",fprefix.Data()),"RECREATE");
   //TFile * file = TFile::Open("RandomDeviations_test.root","UPDATE");
   hresult->Write(TString::Format("resultN%d",nstate) );
   hprof->Write(TString::Format("profResultN%d",nstate) );
   if (hresult3) hresult3->Write(TString::Format("result3N%d",nstate) );
   g->Write(TString::Format("deviations_N%d",nstate) );  
   file->Close();

   //c2->cd(0);
#ifdef G0   
   auto c6 = new TCanvas("C6"); 
   auto gg = (TGraph*) c2->GetListOfPrimitives()->FindObject("Graph_from_profResult")->Clone();
   gg->SetMarkerStyle(20);
   cout << gg->GetMarkerStyle() << endl;
   gg->Draw("AP");
//   c2->Update();
#endif


   std::cout << "number of found duplicates " << nDuplicateFound << std::endl;
}

// template<class Engine> 
// struct EngineWrapper {

//    Engine fEngine; 
//    EngineWrapper() : fEngine() {}

//    double operator() () { return fEngine(); }

//    double IntRndm() 
   
// }


void RandomDeviations(int niter = 10, int ntrial = 2000, int ichangeEl = -1) {

   printf("running random deviations\n");

   if (ichangeEl  >=  0 ) changeElement = ichangeEl; 

   //GetDeviations<Test::MT,uint32_t,uint32_t>(n,ntrial);

   //GetDeviations<Test::LCG,uint32_t,uint32_t>(n,ntrial);

   debug = 0;
   //ntrial = 2;
   engineSubType = 0;   // use different values if we have modfied some parameter of the engine
   //iRBegin = 1;  for mixmax

   GetDeviations<Test::TestEngine<MixMaxEngine<256,0>>,uint64_t,uint64_t>(niter,ntrial);
   
   //GetDeviations<Test::TestEngine<TRandomMixMax256>,uint64_t,uint64_t>(n,ntrial);

   //GetDeviations<Test::TestRanLux<RanLuxSEngine>,float,uint32_t>(niter,ntrial);
#ifdef HAVE_RANLUX_ENGINE
   GetDeviations<RanLuxMatrixEngine,uint32_t,uint32_t>(niter,ntrial);
#endif
   
}
int main() {
   RandomDeviations();
}

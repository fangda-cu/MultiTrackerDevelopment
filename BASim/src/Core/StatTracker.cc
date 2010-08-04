/**
 * \file StatTracker.cc
 *
 * \author smith@cs.columbia.edu
 * \date 06/23/2010
 */

#include "StatTracker.hh"

namespace BASim 
{

IntStatTracker::IntStatTrackerMap IntStatTracker::intTrackerMap;

IntStatTracker& IntStatTracker::getIntTracker(std::string name, int initval)
{
  IntStatTrackerMap::iterator iter = intTrackerMap.find(name);
  if (iter != intTrackerMap.end()) return *(iter->second);
  
  IntStatTracker* t = new IntStatTracker();
  t->m_value = initval;
  intTrackerMap[name] = t;
  return *t;
}

void IntStatTracker::deleteAllIntTrackers()
{
  for( IntStatTrackerMap::iterator iter = intTrackerMap.begin(); iter != intTrackerMap.end(); ++iter )
  {
    delete iter->second;
    iter->second = NULL;
  }

  intTrackerMap.clear();
}

void IntStatTracker::reportIntTrackers()
{
  IntStatTrackerMap::iterator iter;
  
  for( iter = intTrackerMap.begin(); iter != intTrackerMap.end(); ++iter )
  {
    std::cout << iter->first << " = " << *(iter->second) << std::endl;
  }
}

void IntStatTracker::clearIntTrackers()
{
  m_value = 0;
}

IntStatTracker& IntStatTracker::operator+=( const int& val )
{
  m_value += val;
  return *this;
}

IntStatTracker& IntStatTracker::operator=( const int& val )
{
  m_value = val;
  return *this;
}

IntStatTracker::IntStatTracker()
: m_value(0)
{}



DoubleStatTracker::DoubleStatTrackerMap DoubleStatTracker::doubleTrackerMap;

DoubleStatTracker& DoubleStatTracker::getDoubleTracker(std::string name, double initval)
{
  DoubleStatTrackerMap::iterator iter = doubleTrackerMap.find(name);
  if (iter != doubleTrackerMap.end()) return *(iter->second);
  
  DoubleStatTracker* t = new DoubleStatTracker();
  t->m_value = initval;
  doubleTrackerMap[name] = t;
  return *t;
}

void DoubleStatTracker::deleteAllDoubleTrackers()
{
  for( DoubleStatTrackerMap::iterator iter = doubleTrackerMap.begin(); iter != doubleTrackerMap.end(); ++iter )
  {
    delete iter->second;
    iter->second = NULL;
  }
  
  doubleTrackerMap.clear();
}

void DoubleStatTracker::reportDoubleTrackers()
{
  DoubleStatTrackerMap::iterator iter;
  
  for( iter = doubleTrackerMap.begin(); iter != doubleTrackerMap.end(); ++iter )
  {
    std::cout << iter->first << " = " << *(iter->second) << std::endl;
  }
}

void DoubleStatTracker::clearDoubleTrackers()
{
  m_value = 0;
}

DoubleStatTracker& DoubleStatTracker::operator+=( const double& val )
{
  m_value += val;
  return *this;
}

DoubleStatTracker& DoubleStatTracker::operator=( const double& val )
{
  m_value = val;
  return *this;
}

DoubleStatTracker::DoubleStatTracker()
: m_value(0)
{}
  
  
  

std::map<std::string,PairVectorBase*> PairVectorBase::pairVectorMap;

PairVectorBase::PairVectorBase( const std::pair<std::string,std::string>& label )
: m_label(label)
{}  

PairVectorBase::~PairVectorBase() 
{}
  
const std::pair<std::string,std::string>& PairVectorBase::getLabel() const
{
  return m_label;
}  

void PairVectorBase::saveToFiles( const std::string& file_prefix )
{
  for( std::map<std::string,PairVectorBase*>::const_iterator it = pairVectorMap.begin(); it != pairVectorMap.end(); ++it )
  {
    StringPair labels = it->second->getLabel();
    std::string file_name = file_prefix + "_" + labels.first + "_" + labels.second + ".txt";
    it->second->saveToFile(file_name);
  }
}
  
//template <class T, class U>
//void PairVectorBase::insertPair( const std::string& name, const std::pair<T,U>& inputpair, const StringPair& label )
//{
//  std::map<std::string,PairVectorBase*>::const_iterator it = pairVectorMap.find(name);
//  if( it == pairVectorMap.end() )
//  {
//    it = pairVectorMap.insert( std::pair<std::string,PairVector<T,U>*>(label) );
//  }
//  PairVector<T,U>* vecptr = dynamic_cast<PairVector<T,U>*>((*it).second);
//  if( vecptr == NULL )
//  {
//    std::cerr << "\033[31;1mERROR IN PAIRVECTORBASE:\033[m Requested container and pair type do not match." << std::endl;
//    exit(1);
//  }
//  vecptr->insert(inputpair);
//}

void PairVectorBase::clear()
{
  for( std::map<std::string,PairVectorBase*>::iterator it = pairVectorMap.begin(); it != pairVectorMap.end(); ++it )
  {
    assert( it->second != NULL );
    delete it->second;
    it->second = NULL;
  }
  pairVectorMap.clear();  
}
  
  
  
} // namespace BASim


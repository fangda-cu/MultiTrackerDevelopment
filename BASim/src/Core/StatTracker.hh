/**
 * \file StatTracker.hh
 *
 * \author smith@cs.columbia.edu
 * \date 06/23/2010
 */

#ifndef STATTRACKER_HH
#define STATTRACKER_HH

#include "Definitions.hh"

namespace BASim {
  
  class IntStatTracker
  {
  public:
    
    typedef std::map<std::string, IntStatTracker*> IntStatTrackerMap;
    
    /// returns the tracker for a given name; reuses the tracker if it exists
    static IntStatTracker& getIntTracker(std::string name, int initval = 0);
    
    /// deletes all trackers created via getTracker() function
    static void deleteAllIntTrackers();
    
    /// dumps information to stdout
    static void reportIntTrackers();
    
    /// reset any stats
    void clearIntTrackers();
    
    IntStatTracker& operator+=( const int& val );
    
    IntStatTracker& operator=( const int& val );
    
    friend std::ostream& operator<<( std::ostream& os, IntStatTracker& tracker )
    {
      os << tracker.m_value;
      return os;
    }
    
  private:
    static IntStatTrackerMap intTrackerMap; ///< static map of name->Tracker

    /// private tracker constructor; create trackers via static getTracker() function
    IntStatTracker();

    int m_value; ///< value of the tracker
  };

  class DoubleStatTracker
  {
  public:
    
    typedef std::map<std::string, DoubleStatTracker*> DoubleStatTrackerMap;
    
    /// returns the tracker for a given name; reuses the tracker if it exists
    static DoubleStatTracker& getDoubleTracker(std::string name, double initval = 0);
    
    /// deletes all trackers created via getTracker() function
    static void deleteAllDoubleTrackers();
    
    /// dumps information to stdout
    static void reportDoubleTrackers();
    
    /// reset any stats
    void clearDoubleTrackers();
    
    DoubleStatTracker& operator+=( const double& val );
    
    DoubleStatTracker& operator=( const double& val );
    
    friend std::ostream& operator<<( std::ostream& os, DoubleStatTracker& tracker )
    {
      os << tracker.m_value;
      return os;
    }
    
    double getVal() const { return m_value; }
    
  private:
    static DoubleStatTrackerMap doubleTrackerMap; ///< static map of name->Tracker
    
    /// private tracker constructor; create trackers via static getTracker() function
    DoubleStatTracker();
    
    double m_value; ///< value of the tracker
  };

  template <class T, class U>
  class PairVector;
  
  class PairVectorBase
  {
  public:
    
    //typedef std::map<std::string,PairVectorBase*> PairVectorBaseMap;
    typedef std::pair<std::string,std::string> StringPair;

    PairVectorBase( const std::pair<std::string,std::string>& label );
    virtual ~PairVectorBase();
    
    const std::pair<std::string,std::string>& getLabel() const;

    virtual void saveToFile( const std::string& file_name ) const = 0;

    static void saveToFiles( const std::string& file_prefix );
    
    template <class T, class U>
    static void insertPair( const std::string& name, const std::pair<T,U>& inputpair, const StringPair& label = StringPair("DEFAULT_LABEL","DEFAULT_LABEL") )
    {
      std::map<std::string,PairVectorBase*>::iterator it = pairVectorMap.find(name);
      if( it == pairVectorMap.end() )
      {
        PairVectorBase* pvb = new PairVector<T,U>(label);
        std::pair<std::string,PairVectorBase*> newpair(name,pvb); 
        pairVectorMap.insert(newpair);
        it = pairVectorMap.find(name);;
      }
      PairVector<T,U>* vecptr = dynamic_cast<PairVector<T,U>*>((*it).second);
      if( vecptr == NULL )
      {
        std::cerr << "\033[31;1mERROR IN PAIRVECTORBASE:\033[m Requested container and pair type do not match." << std::endl;
        exit(1);
      }
      vecptr->insert(inputpair);
    }
    
    static void clear();
    
  private:
    static std::map<std::string,PairVectorBase*> pairVectorMap;
    std::pair<std::string,std::string> m_label;
  };

  template <class T, class U>
  class PairVector : public PairVectorBase
  {
  public:

    PairVector( const std::pair<std::string,std::string>& label )
    : PairVectorBase( label )
    {}
    virtual ~PairVector() {}

    virtual void saveToFile( const std::string& file_name ) const
    {
      std::ofstream output(file_name.c_str());
      if( !output.is_open() )
      {
        std::cerr << "\033[31;1mWARNING IN PAIRVECTOR:\033[m Failed to open file for output: " << std::endl;
        exit(1);
      }
      
      output << "#" << getLabel().first << '\t' << getLabel().second << std::endl;
      
      for( int i = 0; i < (int) m_pairvec.size(); ++i )
        output << m_pairvec[i].first << '\t' << m_pairvec[i].second << std::endl;

      output.close();
      
      std::cout << "\033[35;1mPAIRVECTOR MESSAGE:\033[m Saving statistics to: " << file_name << std::endl;
    }
    
    void insert( std::pair<T,U> inputpair )
    {
      m_pairvec.push_back(inputpair);
    }

  private:    
    std::vector<std::pair<T,U> > m_pairvec;
  };  

} // namespace BASim

#endif // STATTRACKER_HH

#ifndef acceptance_db_included
#define acceptance_db_included

#include <vector>
#include <map>

class AcceptanceDB {
 public:
  AcceptanceDB(const char* which=0) { loadDB(which); }
  void add(int mw, int mn, double a2010, double a2011, double a2012=-1);
  double getBestEstimate(int mw, int mn, int year) const;

 private:
  struct AcceptPt {
    int mw, mn;
    double forYear(int year) const {
      if (year==2010) return a2010;
      if (year==2011) return a2011;
      if (year==2012) return a2012;

      return 0;
    }
    double a2010, a2011, a2012;
  };


  std::vector<AcceptPt> m_DB;
  std::map<int, std::map<int, AcceptPt*> > m_mmm;
  double interpol2d(int mw, int mn, 
		    int year,
		    const AcceptPt& a,const AcceptPt& b,
		    const AcceptPt& c,const AcceptPt& d) const;
  double dist(const AcceptPt& a, const AcceptPt& b) const;
  void loadDB(const char* which);
};

#endif // acceptance_db_included

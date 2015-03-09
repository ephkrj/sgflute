//note the new division (same for daytime/work):
//tracts<---------------->zones
//nHomeComm<------------->subzone_home
//nHomeNeighborhood<----->postcode_home
//family<---------------->houseid
#include<vector>
#include<list>
using namespace std;
enum {
  AG0,  // 0-4 years old
  AG1,  // 5-18 years old
  AG2,  // 19-29 years old
  AG3,  // 30-64 years old
  AG4,	// 65 and older
  TAG	// total age groups
};	// age groups


struct Person {
  int gender ,race ,citizen ,marital ,number_children ,nSchoolType ,qualification ,work ,income ,headhouse ,occupation ,religion ,partnerid ,houseid ,dwelling ,hasmaid ,mobility , nhours ,industry ,NS ,trans_mode ,trans_time ,oversea ,housesize ,f_nuclei ,generation ,youngest ,numwork , hincome ,momid ,dadid ,unit_id,exactAge;
  unsigned int id;	// agent id
  unsigned char age;	// age group **(change from char)
  unsigned char status;	// infection status
  unsigned char nWhichVload;    // which viral load trajectory (vload) to use?
  char iday;            // infected days
  unsigned char ibits;  // state info for infection
  unsigned char vbits;  // state info for vaccine
  unsigned char vday;   // vaccination timer (number of days elapsed)
  char nAVTimer;        // antiviral timer (number of tablets left)
  unsigned char nQuarantineTimer; // quarantine timer (number of days left)
  char nTravelTimer;	// travel time left
  char sourcetype;      // in which setting infection took place
  unsigned int sourceid;// keeping track of the source of infection
  
  //temporary variables for reading
  int postcode_home; //for neighborhood 
  int zone_work;
  int subzone_work;
  int postcode_work;
  //Home
  int zone_home;	      // (zone_home)
  unsigned int nHomeComm;     // ID of home community (subzone_home)
  unsigned int nHomeNeighborhood;    // (postcode_home) **(change from char)
  unsigned int family;	             // (family ID)
  unsigned int nFamilySize; // family size (housesize)
  unsigned char householdcluster;// household cluster **(change from char)
  double x_coor_home,y_coor_home; //home coordinates
  //Work
  unsigned int nDayTract;        // ID of daytime (work, school) tract (zone_work)
  unsigned int nDayComm;         // ID of daytime (work, school) community (subzone_work)
  int nDayNeighborhood; // ID of work neighborhood (postcode_work) **(change from char)
  int nWorkplace;	// work or school group (work_ID)
  double x_coor_work,y_coor_work; //work coordinates
  //School (need to create something equivalent to work for school)
  int nSchTract; //(zone_school)
  int nSchComm;  // (subzone_school)
  int nSchNeighborhood; //(postcode_school)
  int nSchplace; //(school_ID)
  double x_coor_school,y_coor_school; // school coordinates
  
  double fBaselineVES;    // baseline VES before vaccination
  unsigned char nVaccinePriority; // vaccine priority group (0 for no vaccine, 1 for highest priority, 2 for next-highest priority, etc). Also set to 0 if the person does not want vaccine.
  int nInfectedTime;	// keeping track of the time of infection
  unsigned int nVaccineRestrictionBits; // for keeping track of vaccine restriction categories
 // bool bVaccineEligible[NUMVACCINES];

  // temporary data
  double pri; // infectiousness multiplier
  double prs; // susceptibility multiplier
  bool bWantVac; // Should be vaccinated now (if available)
  bool bWantAV;  // Should get antiviral now (if available)
  /*friend ostream& operator<<(ostream& os, Person& p) {
    os << p.id << " " << (int)p.age << " " << p.nHomeComm << " " << p.nDayComm <<  " " << p.family << " " << (int)p.nWorkplace << " " << (int)p.nHomeNeighborhood << " " << (int)p.nFamilySize << " " << p.sourceid << " " << p.nInfectedTime << " " << (int)p.sourcetype << " ";
    if (isSusceptible(p)) {
      os << "s";
    }
    if (isInfected(p)) {
      os << "i";
    }
    if (isSymptomatic(p)) {
      os << "S";
    }
    return os;
  }

  // person age accessor functions
  friend inline bool isChild(const Person &p) { return (p.age<16); } // is a child?
  friend inline bool isWorkingAge(const Person &p) { return (p.age>=16 && p.age<=65); } // old enough to work?
  
  // person.status accessor functions (vaccination and illness state)
  friend inline bool isSusceptible(const Person &p) { return p.status&SUSCEPTIBLE; }
  friend inline void setSusceptible(Person &p) {p.status|=SUSCEPTIBLE;} 
  friend inline void clearSusceptible(Person &p) {p.status&=~SUSCEPTIBLE;}
  friend inline bool isInfected(const Person &p) { return p.status&INFECTED; }
  friend inline bool isInfectious(const Person &p) { return (p.status&INFECTED && p.iday>=0); }
  friend inline bool isVaccinated(const Person &p) { return p.status&VACCINATED; }
  friend inline bool isBoosted(const Person &p) { return p.status&BOOSTED; }
  friend inline bool isAntiviral(const Person &p) { return p.status&ANTIVIRAL; }
  friend inline bool isSymptomatic(const Person &p) { return p.status&SYMPTOMATIC; }
  friend inline void setSymptomatic(Person &p) {p.status|=SYMPTOMATIC;}
  friend inline void setInfected(Person &p) {p.status|=INFECTED;}
  friend inline void setVaccinated(Person &p) {p.status|=VACCINATED;}
  friend inline void setBoosted(Person &p) {p.status|=BOOSTED;}
  friend inline void setAntiviral(Person &p) {p.status|=ANTIVIRAL;}
  friend inline void clearAntiviral(Person &p) {p.status&=(~ANTIVIRAL);}
  friend inline bool needsBoost(const Person &p) { return p.status&NEEDSBOOST; }
  friend inline void setNeedsBoost(Person &p) {p.status|=NEEDSBOOST;}
  friend inline void clearNeedsBoost(Person &p) {p.status|=NEEDSBOOST;}

  // person.ibits accessor functions (response to infection)
  friend inline bool getWillBeSymptomatic(const Person &p) { return p.ibits&WILLBESYMPTOMATIC; }
  friend inline bool getWillBeAscertained(const Person &p) { return p.ibits&WILLBEASCERTAINED; }
  friend inline int getIncubationDays(const Person &p) { return p.ibits&0x7u; }
  friend inline unsigned int getWithdrawDays(const Person &p) { return (p.ibits>>3)&0x7u; } // if 0, then will not withdraw.  otherwise, withdraw this may days after infection
  friend inline void setWillBeSymptomatic(Person &p) {p.ibits|=WILLBESYMPTOMATIC; }
  friend inline void setWillBeAscertained(Person &p) {p.ibits|=WILLBEASCERTAINED; }
  friend inline void setIncubationDays(Person &p, unsigned int x) {p.ibits=(p.ibits&0xF8u)|x; }
  friend inline void setWithdrawDays(Person &p, unsigned int x) {p.ibits=(p.ibits&0xC7u)|(x<<3); }

  // person.vbits accessor functions (miscellaneous state)
  friend inline void setWhichVaccine(Person &p, unsigned char nWhich) { assert(nWhich<16); p.vbits|=nWhich; } // set which vaccine did this person got (0-15, where 0 is the default vaccine)
  friend inline unsigned char whichVaccine(const Person &p) { return p.vbits&0x0F; }	// which vaccine did this person get?
  friend inline bool isWithdrawn(const Person &p) { return p.vbits&WITHDRAWN; }
  friend inline void setWithdrawn(Person &p) {p.vbits|=WITHDRAWN;}
  friend inline bool isQuarantined(const Person &p) { return p.vbits&QUARANTINED; }
  friend inline void setQuarantined(Person &p) {p.vbits|=QUARANTINED;}
  friend inline void clearQuarantined(Person &p) {p.vbits&=(~QUARANTINED);}
  friend inline bool isAVProphylaxis(const Person &p) { return p.vbits&AVPROPHYLAXIS; }
  friend inline void setAVProphylaxis(Person &p) {p.vbits|=AVPROPHYLAXIS;}
  friend inline void clearAVProphylaxis(Person &p) {p.vbits&=~AVPROPHYLAXIS;}

  // vaccine restriction status accessor functions
  friend inline void setEssential(Person &p) { p.nVaccineRestrictionBits|=VESSENTIAL; }
  friend inline bool isEssential(const Person &p) { return p.nVaccineRestrictionBits&VESSENTIAL; }
  friend inline void setPregnant(Person &p) { p.nVaccineRestrictionBits|=VPREGNANT; }
  friend inline bool isPregnant(const Person &p) { return p.nVaccineRestrictionBits&VPREGNANT; }
  friend inline void setHighRisk(Person &p) { p.nVaccineRestrictionBits|=VHIGHRISK; }
  friend inline bool isHighRisk(const Person &p) { return p.nVaccineRestrictionBits&VHIGHRISK; }
  friend inline void setInfant(Person &p) { p.nVaccineRestrictionBits|=VINFANT; }
  friend inline bool isInfant(const Person &p) { return p.nVaccineRestrictionBits&VINFANT; }

  friend inline unsigned char getVaccinePriority(Person &p) { return p.nVaccinePriority; }
*/
};

struct Community {
  const static int TARGETCOMMUNITYSIZE=2000;
  const static int FAMILIESPERCLUSTER = 4;
  const static int WORKGROUPSIZE = 20;

  vector< unsigned int > workers;   // vector of non-resident worker ids
  // visitors is a copy of the original individuals, so we can assign it
  // new communities, families, workgroups, etc.
  // Just keep the ID and the home rank intact.
  list <Person> visitors;           // list of short-term travelers to this community
#ifdef PARALLEL
  vector <Person> immigrantworkers; // vector of workers from other nodes
#endif
  unsigned int id;	 // unique community ID (build from subzone_home from Person p)
  int nTractID;          // ID of its tract
  int nNumResidents;	 // number of residents (population size)
  unsigned int nFirstPerson, nLastPerson; // ID of first and last person(+1) in community
  int nNumWorkers;	 // number of workers (includes immigrants)
  int nNumNonWorkers;	 // number of resident non-workers
  int nNumWorkersLeaving;// number of residents who work outside the community
  int nNumWorkGroups;	 // number of work groups
  int nNumAge[TAG];	 // number of residents in different age groups
  int ninf[TAG];	 // number of residents currently infected
  int nsym[TAG];	 // number of residents currently symptomatic
  int nEverInfected[TAG];// number of residents ever infected
  int nEverSymptomatic[TAG]; // number of residents ever symptomatic
  int nEverAscertained[TAG]; // number of residents ever ascertained
  double cpcm[5];        // community-specific community contact rates
  double cpnh[5];        // community-specific neighborhood contact rates
  double daycpcm[5];     // community-specific daytime community contact rates
  double daycpnh[5];     // community-specific daytime neighborhood contact rates
  double cps[10];         // community-specific school contact rates
  double y;			//FOR COORD
  double x;			//FOR COORD
};

struct Zone {
	int id;         //zone id
	int pop;        //zone population
	int employable; //as the name suggest, ppl aged between 16-65 of a zone
	int employed;   //subset of employable that has a job (p.occupation!=-1) of a zone
	double x;
	double y;
	std::vector<int> workzones;
	std::vector<int> workpop;
};


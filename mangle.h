#ifndef	_MANGLE_H_
#define	_MANGLE_H_

#include	<cmath>
#include	<iostream>
#include	<fstream>
#include	<sstream>
#include	<iomanip>
#include	<vector>
#include	<list>
#include	<string>
#include	<exception>


// A class to read and look up quantities in Mangle masks.
// This is provided as a "header only library" for convenience and
// because the code is relatively light-weight.
// If compilation time becomes an issue, the main methods can be
// defined in a .cpp file.
//
// Author:	Martin White	(UCB)

// The Mangle package is available at:
// http://casa.colorado.edu/~ajsh/mangle/   or
// http://space.mit.edu/~molly/mangle/
// and the algorithms are described in a series of papers
// including Swanson et al. (2008; MNRAS, 387, 1391).





// All of the Mangle-specific routines are in the Mangle namespace.
namespace Mangle {

// These file-scope quantities are defined here so they can be easily
// changed, e.g. if an MPI version of the code is wanted.
inline void	myexit(int flag)
// A version of exit which can be easily replaced.
{
  exit(flag);
}

inline void	myexception(std::exception& e)
// A simple exception handler.
{
  std::cout<<"Exception "<<e.what()<<std::endl;
  std::cout.flush();
  exit(1);
}

class CapClass {
private:
  double	x,y,z,cm;
public:
  CapClass(const double x1, const double y1, const double z1, const double cm1){
    x = x1; y = y1; z = z1; cm = cm1;
  }
  bool incap(const double x0, const double y0, const double z0) const {
    double	cdot;
    bool	tmp;
    cdot = 1.0 - x*x0 - y*y0 - z*z0;
    if (cm < 0.0)
      tmp = cdot > (-cm);
    else
      tmp = cdot < cm;
    return(tmp);
  }
  bool incap(const double theta, const double phi) const {
    double	cdot,x0,y0,z0;
    bool	tmp;
    x0 = sin(theta)*cos(phi);
    y0 = sin(theta)*sin(phi);
    z0 = cos(theta); 
    cdot= 1.0 - x*x0 - y*y0 - z*z0;
    if (cm < 0.0)
      tmp = cdot > (-cm);
    else
      tmp = cdot < cm;
    return(tmp);
  }
};	// CapClass



class PolygonClass {
private:
  std::list<CapClass>	caps;
  long			pid;
  double		weight;
public:
  PolygonClass() {
    caps.clear();
    weight = 0;
    pid    = 0;
  }
  long getid() const {
    return(pid);
  }
  void setid(const long pid1) {
    pid = pid1;
  }
  double getwt() const {
    return(weight);
  }
  void setwt(const double wt) {
    weight = wt;
  }
  void addcap(CapClass cap1) {
    try {caps.push_back(cap1);} catch(std::exception& e) {myexception(e);}
  }
  bool inpoly(const double x, const double y, const double z) const {
  // Not used currently.
    bool tmp=true;
    for (std::list<CapClass>::const_iterator i=caps.begin();
         i!=caps.end() && tmp; ++i) {
      tmp = tmp && ( i->incap(x,y,z) );
    }
    return(tmp);
  }
  bool inpoly(const double theta, const double phi) const {
    bool tmp=true;
    for (std::list<CapClass>::const_iterator i=caps.begin();
         i!=caps.end() && tmp; ++i) {
      tmp = tmp && ( i->incap(theta,phi) );
    }
    return(tmp);
  }
};	// PolygonClass





class MaskClass {
// A class to handle masks.  When invoked it reads an ascii mask file.
// The main method returns a weight (usually completeness).
private:
  std::vector<PolygonClass>	polygons;
  std::vector< std::list<int> >	pixels;
  int				pixelres;
  char				pixeltype;
  double			totalmaskarea;
  long parsepoly(std::string sbuf, long& ncap, double& weight, long& pixel,
                 double& area) {
  // Reads the "polygon" line in a Mangle file, returning
  // #caps, weight, pixel and area.
  // This routine isn't pretty, but a future evolution should probably
  // move away from the ascii file format and for now this works.
    long	i,j,ipoly=-1;
    std::string	ss;
    if ( (i=sbuf.find("polygon")   )!=std::string::npos &&
         (j=sbuf.find("("))!=std::string::npos) {
      ss.assign(sbuf,i+std::string("polygon").size(),j);
      std::istringstream(ss) >> ipoly;
    }
    else {
      std::cerr << "Cannot parse " << sbuf << std::endl;
      myexit(1);
    }
    if ( (i=sbuf.find("(")   )!=std::string::npos &&
         (j=sbuf.find("caps"))!=std::string::npos) {
      ss.assign(sbuf,i+1,j);
      std::istringstream(ss) >> ncap;
    }
    else {
      std::cerr << "Cannot parse " << sbuf << std::endl;
      myexit(1);
    }
    if ( (i=sbuf.find("caps,") )!=std::string::npos &&
         (j=sbuf.find("weight"))!=std::string::npos) {
      ss.assign(sbuf,i+std::string("caps,").size(),j);
      std::istringstream(ss) >> weight;
    }
    else {
      std::cerr << "Cannot parse " << sbuf << std::endl;
      myexit(1);
    }
    if ( (i=sbuf.find("weight,"))!=std::string::npos &&
         (j=sbuf.find("pixel")  )!=std::string::npos) {
      ss.assign(sbuf,i+std::string("weight,").size(),j);
      std::istringstream(ss) >> pixel;
    }
    else {
      if (pixeltype=='u' || pixelres<0) {
        pixel = 0;
      }
      else {
        std::cerr << "Cannot parse " << sbuf << std::endl;
        myexit(1);
      }
    }
    // The area is between either "pixel" or "weight" and "str".
    if ( (i=sbuf.find("pixel,") )!=std::string::npos &&
         (j=sbuf.find("str"))!=std::string::npos) {
      ss.assign(sbuf,i+std::string("pixel,").size(),j);
      std::istringstream(ss) >> area;
    }
    if ( (i=sbuf.find("weight,") )!=std::string::npos &&
         (j=sbuf.find("str"))!=std::string::npos) {
      ss.assign(sbuf,i+std::string("weight,").size(),j);
      std::istringstream(ss) >> area;
    }
    return(ipoly);
  }
  long pixelnum(const double theta, const double phi) const {
  // For the "simple" pixelization we're just Cartesian in cos(theta) & phi.
    long ipix=0;
    if (pixelres>0) {
      long   ps=0,p2=1;
      for (int i=0; i<pixelres; ++i) {	// Work out # pixels/dim and start pix.
        p2  = p2<<1;
        ps += (p2/2)*(p2/2);
      }
      double cth = cos(theta);
      long   n   = (cth==1.0)?0:long( ceil( (1.0-cth)/2 * p2 )-1 );
      long   m   = long( floor( (phi/2./M_PI)*p2 ) );
      ipix= p2*n+m + ps;
    }
    return(ipix);
  }
public:
  void load(const char fname[]) {
  // Load the mask from fname -- this is provided this way so we can stagger
  // loading of large files in e.g. MPI calls.
  // This is hard-coded for the ascii .ply format.  Future evolution should
  // probably move away from this ascii format.
    std::ifstream fs(fname);
    if (!fs) {
      std::cerr << "Unable to open " << fname << std::endl;
      myexit(1);
    }
    totalmaskarea=0;
    // Read in the number of polygons, which should start the file.
    std::string sbuf;
    getline(fs,sbuf);
    if (fs.eof()) {std::cerr<<"Unexpected end-of-file."<<std::endl;myexit(1);}
    long npoly;
    std::istringstream(sbuf) >> npoly;
    try {
      polygons.resize(npoly);
    } catch(std::exception& e) {myexception(e);}
    // See if we're pixelized -- we only handle simplepix here, everything
    // else is treated as unpixelized.
    long i,j;
    getline(fs,sbuf);
    if ( (i=sbuf.find("pixelization"))!=std::string::npos &&
         (j=sbuf.find("s",i)         )!=std::string::npos) {
      std::string ss;
      ss.assign(sbuf,i+std::string("pixelization").size(),j);
      std::istringstream(ss) >> pixelres;
      pixeltype = 's';
    }
    else  {
      pixelres  = -1;
      pixeltype = 'u';
    }
    // For each polygon, create the appropriate "polygons" entry.
    long maxpix=-1;
    for (int ipoly=0; ipoly<npoly; ++ipoly) {
      long	ncap,pixel;
      double	weight,area;
      while(sbuf.find("polygon")==std::string::npos) {
        getline(fs,sbuf);
        if (fs.eof()) {
          std::cerr<<"Unexpected end-of-file."<<std::endl;
          myexit(1);
        }
      }
      j = parsepoly(sbuf,ncap,weight,pixel,area);
      polygons[ipoly].setid(pixel);
      polygons[ipoly].setwt(weight);
      if (pixel>maxpix) maxpix=pixel;
      if (weight>0.0) totalmaskarea += area;
      // and read in the caps.
      for (int icap=0; icap<ncap; ++icap) {
        getline(fs,sbuf);
        if (fs.eof()) {
          std::cerr<<"Unexpected end-of-file."<<std::endl;
          myexit(1);
        }
        double x1,y1,z1,cm1;
        std::istringstream(sbuf) >> x1 >> y1 >> z1 >> cm1;
        polygons[ipoly].addcap( Mangle::CapClass(x1,y1,z1,cm1) );
      }
    }
    fs.close();
    // If we're pixelized make a list of which polygons lie in each pixel.
    if (pixelres>=0) {
      try{
        pixels.resize(maxpix+1);
        for (int ipoly=0; ipoly<npoly; ++ipoly) {
          pixels[polygons[ipoly].getid()].push_back(ipoly);
        }
      } catch(std::exception& e) {myexception(e);}
    }
  }
  MaskClass() {
  // A light-weight constructor.
    totalmaskarea=  0;
    pixelres     = -1;
  }
  MaskClass(const char fname[]) {	// Load the mask from fname.
    totalmaskarea =  0;
    pixelres      = -1;
    load(fname);
  }
  long npolygons() const { // For general interest, how many polygons in mask
    return(polygons.size());
  }
  double totalarea() const { // The area in the mask.
    return(totalmaskarea);
  }
  double getweight(const double theta, const double phi) const {
  // This is the main method, returning the weight at (theta,phi).
  // Note, this is NOT ra and dec but "mathematical" theta and phi.
    bool	notfnd=true;
    double	wt=0;
    if (pixelres==-1) {	// Mask isn't pixelized.
      for (long ii=0; ii<polygons.size() && notfnd; ++ii)
        if (polygons[ii].inpoly(theta,phi)) {
          wt    = polygons[ii].getwt();
          notfnd= false;
        }
    }
    else {
      double phir=phi;
      while (phir<  0   ) phir += 2*M_PI;
      while (phir>2*M_PI) phir -= 2*M_PI;
      long ipix;
      if (pixeltype=='s')
        ipix= this->pixelnum(theta,phir);
      else {
        std::cerr << "Unknown pixelization scheme " << pixeltype << std::endl;
        myexit(1);
      }
      if (ipix<pixels.size()) {
        for (std::list<int>::const_iterator ii=pixels[ipix].begin();
             ii!=pixels[ipix].end() && notfnd; ++ii) {
          if (polygons[*ii].inpoly(theta,phi)) {
            wt    = polygons[*ii].getwt();
            notfnd= false;
          }
        }
      }
    }
    if (notfnd)
      return(0.0);
    else
      return(wt);
  }
};	// MaskClass

}	// namespace Mangle
#endif

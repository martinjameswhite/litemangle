#!/usr/bin/env python
#

import numpy  as N
import os
import sys
import re
import string


class LiteMangle:
    """
    LiteMangle:
    A Python class to implement some very basic mangle routines which
    allow us to manipulate a mangle mask.
    This is a subset of the more complicated Python-based Mangle package,
    though it contains more methods than strictly necessary.
    The class is initialized with 1 argument, the file name of an
    ascii string (in Mangle polygon format) containing the mask.
    """

    __author__ = "Martin White"
    __version__ = "1.1"
    __email__  = "mwhite@berkeley.edu"

    def incap_spam(self,cap,x0,y0,z0):
        """
        incap_spam(self,cap,x0,y0,z0):
        This is an internal routine which you shouldn't need to use.
        Returns True for each (theta,phi) that lies in the cap specified by the
        4-vector "cap" containing (x,y,z,cm), and False for the rest.
        """
        cd = 1.0-cap[0]*x0-cap[1]*y0-cap[2]*z0
        return(((cap[3]<0.)&(cd>N.fabs(cap[3])))|((cap[3]>0.)&(cd<cap[3])))

    def inpoly_spam(self,polygon,x0,y0,z0):
        """
        inpoly_spam(self,polygon,theta,phi):
        This is an internal routine which you shouldn't need to use.
        Returns True if (theta,phi) is in the polygon, i.e. if it is
        within all of the caps in the polygon.
        A polygon is a list giving the polygon number, the weight and
        then all of the caps (each cap is a 4-vector/4-list).
        """
        test = N.ones(len(x0),dtype=bool)
        for x in polygon[3:]:
            test &= self.incap_spam(x,x0,y0,z0)
        return(test)

    def get_polyids(self,ra,dec):
        """
        polyid(self,ra,dec):
        Return the ID number of polygons given arrays of RA and DEC
        in decimal degrees (assumed to be numpy arrays).
        """
        theta = N.pi/180. * (90.0-dec)
        phi   = N.pi/180. * ra
        sintheta = N.sin(theta)
        x0 = sintheta*N.cos(phi)
        y0 = sintheta*N.sin(phi)
        z0 = N.cos(theta)
        goodpolys = -N.ones(len(ra),dtype='i8')
        for poly in self.polylist:
            test = self.inpoly_spam(poly,x0,y0,z0)
            goodpolys[test] = poly[0]
        return(goodpolys)

    def get_areas(self,ra,dec):
        """
        get_areas(self,ra,dec):
        Return the areas of the polygons containing each RA/Dec pair.
        Result is in steradians.
        """
        theta = N.pi/180. * (90.0-dec)
        phi   = N.pi/180. * ra
        sintheta = N.sin(theta)
        x0 = sintheta*N.cos(phi)
        y0 = sintheta*N.sin(phi)
        z0 = N.cos(theta)
        goodpolys = -N.ones(len(ra))
        for poly in self.polylist:
            test = self.inpoly_spam(poly,x0,y0,z0)
            goodpolys[test] = poly[2]
        return(goodpolys)

    def get_all_areas(self):
        """
        get_all_areas(self):
        Return an array containing the areas of all polygons, indexed
        by polyid.
        Result is in steradians.
        """
        area = N.zeros(self.npoly)
        for poly in self.polylist:
            area[poly[0]] = poly[2]
        return(area)

    def get_weights(self,ra,dec):
        """
        get_weights(self,ra,dec):
        Return the weights of the polygons containing each RA/Dec pair,
        for arrays of RA and DEC in decimal degrees.
        """
        theta = N.pi/180. * (90.0-dec)
        phi   = N.pi/180. * ra
        sintheta = N.sin(theta)
        x0 = sintheta*N.cos(phi)
        y0 = sintheta*N.sin(phi)
        z0 = N.cos(theta)
        goodpolys = -N.ones(len(ra))
        for poly in self.polylist:
            test = self.inpoly_spam(poly,x0,y0,z0)
            goodpolys[test] = poly[1]
        return(goodpolys)

    def get_all_weights(self):
        """
        get_all_weights(self):
        Return an array containing the weights of all polygons, indexed
        by polyid.
        """
        weight = N.zeros(self.npoly)
        for poly in self.polylist:
            weight[poly[0]] = poly[1]
        return(weight)

    def total_area(self):
        """
        total_area(self):
        Returns the total area in the mask (i.e. the sum of the areas of
        each polygon) and the total "effective" area (i.e. the area weighted
        by the completeness).
        Returns (tot_area,eff_area).
        """
        tot_area,eff_area = 0.0,0.0
        for poly in self.polylist:
            tot_area += poly[2]
            eff_area += poly[2]*poly[1]
        return((tot_area,eff_area))

    def set_weights(self,weight):
        """
        set_weights(self,weight):
        Sets the weight entries of the polygons to the weight array.
        The weight array should be of length Npoly.
        """
        for i,poly in enumerate(self.polylist):
            self.polylist[i][1]=weight[poly[0]]

    def set_one_weight(self,polyid,weight):
        """
        set_one_weight(self,polyid,weight):
        Sets the weight entry of the single polygon "polyid" to weight.
        """
        for i,poly in enumerate(self.polylist):
            if poly[0]==polyid:
                self.polylist[i][1]=weight
                break

    def set_all_weights(self,weight=0.0):
        """
        set_all_weights(self,weight=0.0):
        Sets the weight entry of all polygons to "weight", a scalar,
        typically used to set everything to 0.
        """
        for i in range(len(self.polylist)):
            self.polylist[i][1]=weight

    def write_ply(self,fn):
        """
        write_ply(self,fn):
        Writes a Mangle-formatted polygon file containing the information
        in the class.
        This does not include pixelization information.  Any mask written
        in this manner should be pixelized after the fact using the Mangle
        command line tools, if necessary.
        """
        ff = open(fn,"w")
        ff.write("%d polygons\n"%len(self.polylist))
        for poly in self.polylist:
            str = "polygon %10d ( %d caps,"%(poly[0],len(poly[3:]))
            str+= " %.8f weight, %.15f str):\n"%(poly[1],poly[2])
            ff.write(str)
            for cap in poly[3:]:
                ff.write("%25.20f %25.20f %25.20f %25.20f\n"%\
                  (cap[0],cap[1],cap[2],cap[3]))
        ff.close()

    def __init__(self,fn):
        """
        __init__(self,fn):
        The class is initialized with the name of an ascii file containing
        the Mangle mask.
        """
        if not os.path.exists(fn):
            raise RuntimeError,"Can not find %s"%fn
        #
        # It's useful to pre-compile a regular expression for a mangle line
        # defining a polygon.
        ex1 = re.compile(r"polygon\s+(\d+)\s+\(\s*(\d+)\s+caps")
        ex2 = re.compile(r"(\d*\.?\d+)\s+weight")
        ex3 = re.compile(r"(\d*\.?\d+)\s+str")
        #
        ff = open(fn,"r")
        self.npoly = 0
        line = ff.readline()
        ss = re.match(r"(\d+)\s+polygons",line)
        if ss==None:
            raise RuntimeError,"Can not parse 1st line of %s"%fn
        else:
            self.npoly = int( ss.group(1) )
        #
        self.polylist = []
        #
        ss = ex1.match(line)
        while len(line)>0:
            while (ss==None)&(len(line)>0):
                line = ff.readline()
                ss   = ex1.match(line)
            if len(line)>0:
                ipoly= int(ss.group(1))
                ncap = int(ss.group(2))
                # Check to see if we have a weight.
                ss = ex2.search(line)
                if ss==None:
                    weight=0.0
                else:
                    weight=float(ss.group(1))
                # Check to see if we have an area.
                ss = ex3.search(line)
                if ss==None:
                    area= -1.0
                else:
                    area=float(ss.group(1))
                polyg= [ipoly,weight,area]
                for i in range(ncap):
                    line = ff.readline()
                    cap  = [float(x) for x in string.split(line)]
                    polyg.append(cap)
                self.polylist.append(polyg)
                ss=None
        ff.close()
        # Check whether the polyids are sequential and range from 0 to npoly
        # If they don't, then there may be a problem with the file.
        # NOTE: this should always be correct for current_boss_geometry.
        badcounter = 0
        for i,poly in enumerate(self.polylist):
            if i != poly[0]:
                badcounter += 1
        if badcounter > 0:
            print "WARNING!!!!"
            print "Found",badcounter,"polygons out of order."

        if len(self.polylist) != self.npoly:
            print "Got %d polygons, expecting %d."%\
              (len(self.polylist),self.npoly)

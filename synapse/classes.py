#
#    Module      : classes.py
#    Description : Various class implementations.
#
#    Copyright 2014 Max Larsson <max.larsson@liu.se>
#
#    This file is part of Synapse.
#
#    Synapse is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Synapse is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Synapse.  If not, see <http://www.gnu.org/licenses/>.

import sys
import exceptions
import os.path
import random
import itertools
import unicode_csv
import geometry
from fileIO import *
from stringconv import *


# Convenience functions
def dotProgress(x, linelength=80, char='.'):
    """ Simple progress indicator on sys.stdout
    """
    sys.stdout.write(char)
    if (x + 1) % linelength == 0:
        sys.stdout.write('\n')

#
# Classes
#
class Point(geometry.Point):
    def __init__(self, x=None, y=None, ptype="", profile=None):
        if isinstance(x, Point):
            geometry.Point.__init__(self, x.x, x.y)
        else:
            geometry.Point.__init__(self, x, y)
        self.profile = profile
        if self.profile != None:
            self.opt = self.profile.opt
        else:
            self.opt = None
        self.skipped = False
        self.type = ptype
        self.cluster = None
           

    def getDistToPosel(self, profile):
        """ Determines the particle's perpendicular distance to the 
            outline of the postprofileaptic element
        """
        distToPosel = self.perpendDist(profile.posel, posloc=profile.postsynloc)
        if distToPosel == None and self.type == "particle":
            self.skipped = True            
            if self.type == "particle":
                ProfileWarning(profile, "Unable to project on postsynaptic "
                                    "element\n   => skipping particle at"
                                    " %s" % self)
        return distToPosel

    def getDistToPrsel(self, profile):
        """ Determines the particle's perpendicular distance to the 
            outline of the presynaptic element
        """        
        if len(profile.prsel) > 0:
            return self.perpendDist(profile.prsel, posloc=profile.postsynloc)
        else:
            return None

    def getLateralLocation(self, profile):
        """ Determines the particle's lateral location, defined as the 
            location along the mediolateral synaptic axis of the projection of 
            the particle on the postsynaptic membrane.
        """
        # A particle is classified as synaptic if its projection on the 
        # postsynaptic membrane is within the spatial resolution limit 
        # of the postsynaptic membrane
        for psd in profile.psdli:
            if (self.lateralDistSyn(profile.posel, psd.posm)
                - profile.spatial_resolutionInPixels <= psd.posm.length() / 2.0):
                return "synaptic"
        # If not synaptic but still within the extent of the synaptic membrane
        if (self.lateralDistSyn(profile.posel, profile.totalPosm)
            < (profile.totalPosm.length() / 2.0)):
                return "within perforation"
        # Otherwise it may be perisynaptic, defined here as being outside the
        # synapse but within half a diameter of the nearest PSD... or perhaps
        # it should be a fixed distance, e g 200 nm?
        # I don't want to determine which PSD is closest and which edge of that
        # PSD faces the extrasynapse, so I just determine the distance to the 
        # nearest edge; if that is an internal edge, the particle will have been 
        # identified as within perforation above.
        for psd in profile.psdli:
            if self.lateralDistSyn(profile.posel, psd.posm) <= psd.posm.length():
                return "perisynaptic"
        # If none of the above
        return "extrasynaptic"

    def getStrictLateralLocation(self, profile):
        """ Determines the particle's lateral location, defined as the 
            location along the mediolateral synaptic axis of the projection of 
            the particle on the postsynaptic membrane. Does not take spatial
            resolution into account, such that a particle is considered 
            synaptic only when strictly projectable on the postsynaptic membrane.
        """
        # A particle is classified as synaptic if its projection on the 
        # postsynaptic membrane is within the postsynaptic membrane
        for psd in profile.psdli:
            if (self.lateralDistSyn(profile.posel, psd.posm)
                <= psd.posm.length() / 2.0):
                return "synaptic"
        # If not synaptic but still within the extent of the synaptic membrane
        if (self.lateralDistSyn(profile.posel, profile.totalPosm)
            < (profile.totalPosm.length() / 2.0)):
                return "within perforation"
        # Otherwise it may be perisynaptic, defined here as being outside the
        # synapse but within half a diameter of the nearest PSD... or perhaps
        # it should be a fixed distance, e g 200 nm?
        # I don't want to determine which PSD is closest and which edge of that
        # PSD faces the extrasynapse, so I just determine the distance to the 
        # nearest edge; if that is an internal edge, the particle will have been 
        # identified as within perforation above.
        for psd in profile.psdli:
            if self.lateralDistSyn(profile.posel, psd.posm) <= psd.posm.length():
                return "perisynaptic"
        # If none of the above
        return "extrasynaptic"

        
    def __axodendriticLocation(self):
        """ Determines the particle's location along the axodendritic axis.
        """        
        if self.distToPosel >= 0:
            return "postsynaptic"
        elif self.distToPrsel == None: # if there is no presynaptic membrane,
            return "not postsynaptic"  # that's about all we can say
        elif self.distToPrsel <= 0:
            return "presynaptic"
        else:
            return "neither pre- or postsynaptic" 


    def __lateralDistPSD(self, profile):
        mindist = sys.maxint
        closestPSD = None
        for psd in profile.psdli:
            d = self.lateralDistSyn(profile.posel, psd.posm) - (psd.posm.length() / 2)
            if d < mindist:
                mindist = d
                closestPSD = psd
        if not closestPSD:
            raise ProfileError(profile, "could not determine lateral distance "
                                    "to a PSD of particle at %s" % self)
        mindist = self.lateralDistSyn(profile.posel, closestPSD.posm)
        normdist = mindist / (closestPSD.posm.length() / 2)            
        return mindist, normdist, closestPSD

    def __psdAssociation(self, profile):
        isWithinPSD = False
        isAssociatedWithPSD = False
        associatedPSD = None
        mindist = sys.maxint
        for psd in profile.psdli:
            if self.isWithinPolygon(psd.psdposm):
                isWithinPSD = True
                isAssociatedWithPSD = True
                associatedPSD = psd
                break
            dist = self.perpendDistClosedPath(psd.psdposm,
                                              doNotCareIfOnOrOffSeg=True)
            if dist <= profile.spatial_resolutionInPixels:
                isAssociatedWithPSD = True
                if dist < mindist:
                    associatedPSD = psd
                    mindist = dist
        return isWithinPSD, isAssociatedWithPSD, associatedPSD

    def __isWithinPostsynapticMembraneShell(self, profile, opt):
        return (self.distToPosel != None 
                and abs(self.distToPosel) <= 
                        geometry.toPixelUnits(opt.shell_width, profile.pixelwidth))

    def __membraneAssociation(self, profile):
        posmAssociated = False
        presmAssociated = False       
        if (self.distToPosel != None and 
            abs(self.distToPosel) <= profile.spatial_resolutionInPixels):
            posmAssociated = True
        if (self.distToPrsel != None and  
            abs(self.distToPrsel) <= profile.spatial_resolutionInPixels):
            presmAssociated = True
        return posmAssociated, presmAssociated
         
    def isWithinHole(self, profile):
        """  Determine whether self is inside a profile hole
        """                
        for h in profile.holeli:
            if self.isWithinPolygon(h):
                return True
        return False

 
    def determineStuff(self, parentli, profile, opt):
        if self.isWithinHole(profile):
            if self.type == "particle": # no warning for random and grid points
                ProfileWarning(profile, "Discarding point at %s: Located "
                                    "within a profile hole" % self)
            self.skipped = True
            return
        self.distToPosel = self.getDistToPosel(profile)
        if self.distToPosel > geometry.toPixelUnits(opt.shell_width,
                                           profile.pixelwidth):
            if self.type == "particle":
                ProfileMessage(profile, "Discarding particle at %s: "
                                    "Located outside shell" % self)
            self.skipped = True
            return
        self.distToPrsel = self.getDistToPrsel(profile)
        self.lateralLocation = self.getLateralLocation(profile)
        self.strictLateralLocation = self.getStrictLateralLocation(profile)
        self.lateralDistPSD, self.normLateralDistPSD,\
        self.closestPSD = self.__lateralDistPSD(profile)
        self.axodendriticLocation = self.__axodendriticLocation()
        self.isWithinPSD, self.isAssociatedWithPSD,\
        self.associatedPSD  = self.__psdAssociation(profile)
        self.isWithinPostsynapticMembraneShell = self.__isWithinPostsynapticMembraneShell(profile, opt)
        self.isPostsynapticMembraneAssociated,\
        self.isPresynapticMembraneAssociated = self.__membraneAssociation(profile)



    def determineNearestNeighbour(self, pointli, opt):
        # Assumes that only valid (projectable, within shell etc) points are
        # in pointli
        nearestNeighbourDist = None
        #nearestNeighbourParticle = Particle()
        mindist = float(sys.maxint)
        for p in pointli:
            # I will instead exclude non-desired points from the supplied point
            # list *before* calling this function
            #if p is not self and p.isWithinProfile:
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
                    #minp = p
        if not mindist < float(sys.maxint):
            return None
        else:
            nearestNeighbourDist = mindist
            #nearestNeighbourParticle = minp
            return nearestNeighbourDist

    def determineNearestLateralNeighbour(self, pointli, profile, opt):
        # Assumes that only valid (projectable, within shell etc) points are
        # in pointli
        nearestLateralNeighbourDist = None
        #nearestLateralNeighbourParticle = Particle()
        mindist = float(sys.maxint)
        for p in pointli:
            if p is not self:
                d = self.lateralDistToPoint(p, profile.posel)
                if d < mindist:
                    mindist = d
                    #minp = p
        if not mindist < float(sys.maxint):
            return None
        else:
            nearestLateralNeighbourDist = mindist
            #nearestLateralNeighbourPoint = minp
            return nearestLateralNeighbourDist

    def lateralDistToPoint(self, p2, sel):
        """ Determine lateral distance to a point p2 along sel.
        """
        path = geometry.SegmentedPath()
        p2_project, p2_seg_project = p2.projectOnPathOrEndnode(sel)
        project, seg_project = self.projectOnPathOrEndnode(sel)
        path.extend([project, p2_project])
        if p2_seg_project < seg_project:
            path.reverse()
        for n in range(min(p2_seg_project, seg_project) + 1,
                       max(p2_seg_project, seg_project)):
            path.insert(len(path)-1, sel[n])
        L = path.length()
        return L


    def lateralDistSyn(self, sel, sm):
        """ Determine lateral distance to center of (post- or pre-)synaptic 
            membrane. 
        """
        return self.lateralDistToPoint(sm.centerPoint(), sel)
        
class PointList(list):
    def __init__(self, pointli, ptype):
        try:
            self.extend([Point(p.x, p.y, ptype) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError, 'not a list of Point elements'

class PSD(geometry.SegmentedPath):
    def __init__(self, pointli):
        geometry.SegmentedPath.__init__(self, pointli)
        
    def adjustPSD(self, profile):
        """ Adjust PSD coordinates so that the PSD ends exactly on 
            the postsynaptic element. 
        """
        # Partition psd into paths defined by the intersections with posel,
        # beginning and ending with the projections of the end nodes of psd.
        # Each path will then be completely on one side of posel.
        pathli = [] 
        pathli.append(geometry.SegmentedPath())
        pathli[0].append(self[0].projectOnPathOrEndnode(profile.posel)[0])
        for n in range(0, len(self)-1): 
            pathli[-1].append(self[n])
            for d in range(len(profile.posel)-1):
                p = geometry.segmentIntersection(self[n], self[n+1],
                                        profile.posel[d], profile.posel[d+1])
                if p:
                    break
            if p: 
                pathli[-1].append(p)   
                pathli.append(geometry.SegmentedPath())          
                pathli[-1].append(p)     
        pathli[-1].append(self[-1])
        pathli[-1].append(self[-1].projectOnPathOrEndnode(profile.posel)[0])
        # Now, look for the longest path. This is assumed to be the intended
        # psd. (Perhaps area is more relevant. However, this is messier because
        # we need to determine the part of dm enclosed by path for each path.)
        maxLength = 0
        for path in pathli:
            L = path.length()
            if L > maxLength:
                longestPath = path
                maxLength = L
        del self[:]
        self.extend(longestPath)
        return self
        
    def getPosm(self, profile):
        """ Return postsynaptic membrane (i e the part of the postsynaptic
            element which is bounded by the PSD). Assume that the PSD is 
            adjusted so as to end exactly on postsynaptic element, and
            that the PSD is oriented in the same direction as posm
            (ie node1 <= node2).
        """        
        posm = geometry.SegmentedPath()
        p1, node1 = self[0].projectOnPathOrEndnode(profile.posel)
        p2, node2 = self[-1].projectOnPathOrEndnode(profile.posel)
        if not None in (node1, node2):
            posm.extend(profile.posel[node1+1:node2+1])
            posm.insert(0, p1)
            posm.append(p2)
        if len(posm) == 0 or posm.length() == 0:
            raise ProfileError(profile, "Could not determine postsynaptic "
                                     "membrane")      
        # It appears that sometimes the first or last nodes get duplicated, 
        # resulting in a zero vector and division by zero when projecting on 
        # posm. Simply checking for the duplicated node and removing it if 
        # present obviously works, but perhaps I should delve deeper into the 
        # problem. The node should not be duplicated in the first place. 
        # - 12/15/2003
        try:
            if posm[0] == posm[1]:
                del posm[0]
        except IndexError:
            pass        
        try:
            if posm[-1] == posm[-2]:
                del posm[-1]
        except IndexError:
            pass
        return posm


    def getPrsm(self, profile):
        prsm = geometry.SegmentedPath()
        p1, node1 = self[0].projectOnPathOrEndnode(profile.prsel)
        p2, node2 = self[-1].projectOnPathOrEndnode(profile.prsel)
        if not None in (node1, node2):
            # Now we have the projections of the psd endnodes on prsel. But if
            # the prsel is far from a psd endnode, and the psd is very curved, 
            # the projection might not coincide with the prsel endnode, even 
            # if the prsel endnode is apposed to the psd. Therefore, see if the
            # prsel endnode projections are on the posm, and if these are  
            # closer than the projections of psd endnodes on prsel. If so, use 
            # them instead. Phew! (04/03/2003)
            p0, foo = profile.prsel[0].projectOnPath(self.posm) # should be self.posm? 10.11.07 / ML
            pN, foo = profile.prsel[-1].projectOnPath(self.posm)
            if p0 and p0.dist(p1) < self[0].dist(p1):
                p1 = profile.prsel[0]
                node1 = 0
            if pN and pN.dist(p2) < self[-1].dist(p2):
                p2 = profile.prsel[-1]
                node2 = len(profile.prsel)
            prsm.extend(profile.prsel[node1+1:node2+1])
            prsm.insert(0, p1)
            prsm.append(p2)
        if len(prsm) == 0 or prsm.length() == 0:
            raise ProfileError(profile, "Could not determine presynaptic "
                                     "membrane")            
        return prsm    

    def getPsdPosm(self, profile):
        """ Make a closed-path PSD by concatenating the PSD and the PoSM.
            Assume that PSD and posm share end nodes and that they are 
            oriented in the same direction. 
        """
        pol = geometry.SegmentedPath()
        pol.extend(self.posm)
        pol.reverse()
        pol.extend(self[1:-1])
        return pol 

    def getCleft(self, profile):
        """ Return the polygon outline of the synaptic cleft, defined as the 
            polygon outlining the area bounded by the prsm, the segments formed 
            by projecting the end nodes of prsm onto posm, and the part of posm
            between these projections.  
        """
        pol = geometry.SegmentedPath()
        if len(profile.prsel) == 0:
            return pol          
        proj0, seg0 = self.prsm[0].projectOnPathOrEndnode(self.posm)
        projN, segN = self.prsm[-1].projectOnPathOrEndnode(self.posm)
        pol.extend(self.prsm)
        if seg0 < segN:
            pol.reverse()
            pol.append(proj0)
            pol.extend(self.posm[seg0+1:segN+1])
            pol.append(projN)
        else:
            pol.append(projN)
            pol.extend(self.posm[segN+1:seg0+1])
            pol.append(proj0)
        return pol            
                          
    def cleftWidthAverage(self):
        """ Return synaptic cleft width, calculated by dividing the area of the 
            cleft outline by the length of posm. 
        """
        try:
            return self.cleft.area() / self.posm.length()
        except (TypeError, IndexError):
            return None
# End of class PSD

class ClusterData(list):
    def __init__(self, pli=[]):
        try:
            self.extend([Point(p.x, p.y) for p in pli])
        except (AttributeError, IndexError):
            raise TypeError, 'not a point list'
        self.convexHull = SegmentedPath()

    def lateralDistSyn(self, c2, mem):
        """ Determine lateral distance to a cluster c2 along postsynaptic
            element membrane.
        """
        centroid = Point(self.convexHull.centroid())
        centroid2 = Point(c2.convexHull.centroid())
        return centroid.lateralDistToPoint(centroid2, mem)


class ProfileData:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.outputfn = ""
        self.warnflag = False
        self.errflag = False             

    def process(self, opt):
        """ Parse profile data from a file and determine distances
        """
        try:
            self.__parse(opt)
            self.__checkPaths()
            sys.stdout.write("Processing profile...\n")
            self.spatial_resolutionInPixels = geometry.toPixelUnits(
                                                        opt.spatial_resolution,
                                                        self.pixelwidth)                 
            self.prsel.orientToPath(self.posel)
            for p in self.psdli:
                p.adjustPSD(self)
                p.orientToPath(self.posel)
                p.psdloc = p.centerPoint()
                p.posm = p.getPosm(self)
                p.prsm = p.getPrsm(self)
                p.psdposm = p.getPsdPosm(self)                
                p.cleft = p.getCleft(self)
                p.cleftWidth = p.cleftWidthAverage()
            self.totalPosm = self.__getTotalPosm()
            self.postsynloc = self.psdli[0].psdloc 
            for p in self.pli:
                p.determineStuff(self.pli, self, opt)
            self.pli = [p for p in self.pli if not p.skipped]
            if opt.useGrid:
                for g in self.gridli:
                    g.determineStuff(self.gridli, self, opt)
                self.gridli = [g for g in self.gridli if not g.skipped]
            if opt.useRandom:
                for r in self.randomli:
                    r.determineStuff(self.randomli, self, opt)
                self.randomli = [r for r in self.randomli if not r.skipped]
            self.getInterDistlis(opt)
            self.getClusters(opt)
            self.getMonteCarlo(opt)
            if opt.stop_requested:
                return            
            sys.stdout.write("Done.\n")            
            if opt.outputs['individual profiles']:
                self.__saveResults(opt)
        except ProfileError, (self, msg):
            sys.stdout.write("Error: %s\n" % msg)
            self.errflag = True

    def getMonteCarlo(self, opt):
        if opt.run_monte_carlo == False:
            self.mcli = []
        else:
            sys.stdout.write("Running Monte Carlo simulations...\n")
            self.mcli = self.runMonteCarlo(opt)

    def getInterDistlis(self, opt):
        if not opt.determine_interpoint_dists:
            return
        if not True in [val for key, val in opt.interpoint_relations.items()
                                if "simulated" not in key]:
            return
        sys.stdout.write("Determining interpoint distances...\n")
        if opt.interpoint_relations["particle - particle"]:
            self.pp_distli, self.pp_latdistli = self.getSameInterpointDistances(opt, self.pli)
        if opt.useRandom and opt.interpoint_relations["random - particle"]:
            self.rp_distli, self.rp_latdistli = \
                    self.getInterpointDistances2(opt, self.randomli, self.pli)


    def getClusters(self, opt):
        if not (opt.determine_clusters and opt.within_cluster_dist > 0):
            return
        sys.stdout.write("Determining clusters...\n")
        self.clusterli = self.determineClusters(self.pli, opt)
        self.processClusters(self.clusterli, opt)

    def getSameInterpointDistances(self, opt, pointli):
        dli = []
        latdli = []
        for i in range(0, len(pointli)):
            if opt.stop_requested:
                return [], []
            if opt.interpoint_dist_mode == 'all':
                for j in range(i + 1, len(pointli)):
                    if opt.interpoint_shortest_dist:
                        dli.append(pointli[i].dist(pointli[j]))
                    if opt.interpoint_lateral_dist:
                        latdli.append(pointli[i].lateralDistToPoint(
                                                                   pointli[j],
                                                                   self.posel))
            elif opt.interpoint_dist_mode == 'nearest neighbour':
                if opt.interpoint_shortest_dist:
                    dli.append(pointli[i].determineNearestNeighbour(pointli,
                                                                    opt))
                if opt.interpoint_lateral_dist:
                    latdli.append(pointli[i].determineNearestLateralNeighbour(
                                                            pointli, self, opt))
        dli = [d for d in dli if d != None]
        latdli = [d for d in latdli if d != None]
        return dli, latdli


    def getInterpointDistances2(self, opt, pointli, pointli2=[]):
        dli = []
        latdli = []
        for p in pointli:
            if opt.stop_requested:
                return [], []
            if opt.interpoint_dist_mode == 'all':
                for p2 in pointli2:
                    if opt.interpoint_shortest_dist:
                        dli.append(p.dist(p2))
                    if opt.interpoint_lateral_dist:
                        latdli.append(p.lateralDistToPoint(p2, self.posel))
            elif opt.interpoint_dist_mode == 'nearest neighbour':
                if opt.interpoint_shortest_dist:
                    dli.append(p.determineNearestNeighbour(pointli2, opt))
                if opt.interpoint_lateral_dist:
                    latdli.append(p.determineNearestLateralNeighbour(pointli2,
                                                                    self, opt))
        dli = [d for d in dli if d != None]
        latdli = [d for d in latdli if d != None]
        return dli, latdli


    def runMonteCarlo(self, opt):
        def isValid(p):
            d = p.perpendDist(self.posel)
            if (d == None or abs(d) > shell_width or p.isWithinHole(self) or
                p in mcli[n]["pli"]):
                return False
            if opt.monte_carlo_simulation_window == "whole profile":
                return True
            if opt.monte_carlo_strict_location == True:
                latloc = p.getStrictLateralLocation(self)
            else:
                latloc = p.getLateralLocation(self)
            if (opt.monte_carlo_simulation_window == "synapse" and
                latloc in ("synaptic", "within perforation")):
                    return True
            if (opt.monte_carlo_simulation_window == "synapse - perforations" and
                latloc == "synaptic"):
                    return True
            if (opt.monte_carlo_simulation_window == "synapse + perisynapse" and
                latloc in ("synaptic", "within perforation", "perisynaptic")):
                    return True
            # TODO : include PSD as a window option
            return False


        pli = self.pli
        if opt.monte_carlo_strict_location:
            locvar = "strictLateralLocation"
        else:
            locvar = "lateralLocation"
        if opt.monte_carlo_simulation_window == "whole profile":
            numpoints = len(pli)   # particles outside shell have been discarded
        elif opt.monte_carlo_simulation_window == "synapse":
            numpoints = len([p for p in pli if p.__dict__[locvar] in 
                                ("synaptic", "within perforation")])
        elif opt.monte_carlo_simulation_window == "synapse - perforations":
            numpoints = len([p for p in pli
                                if p.__dict__[locvar] == "synaptic"])
        elif opt.monte_carlo_simulation_window == "synapse + perisynapse":
            numpoints = len([p for p in pli if p.__dict__[locvar] in
                            ("synaptic", "within perforation", "perisynaptic")])
        else:
            return []
        allpaths = geometry.SegmentedPath([p for path in (self.posel, self.prsel,
                                        [q for psd in self.psdli for q in psd],
                                        [r for h in self.holeli for r in h])
                                  for p in path])
        box = allpaths.boundingBox()
        shell_width = geometry.toPixelUnits(opt.shell_width, self.pixelwidth)
        mcli = []
        for n in range(0, opt.monte_carlo_runs):
            if opt.stop_requested:
                return []
            dotProgress(n)
            mcli.append({"pli": [],
                         "simulated - simulated": {"dist": [], "latdist": []},
                         "simulated - particle": {"dist": [], "latdist": []},
                         "particle - simulated": {"dist": [], "latdist": []},
                         "clusterli": []})
            for _ in range(0, numpoints):
                while True:
                    x = random.randint(int(box[0].x - shell_width),
                                       int(box[1].x + shell_width) + 1)
                    y = random.randint(int(box[0].y - shell_width),
                                       int(box[2].y + shell_width) + 1)
                    p = Point(x, y)
                    if isValid(p):
                        break
                # escape the while loop when a valid simulated point is found
                mcli[n]["pli"].append(p)
            for p in mcli[n]["pli"]:
                p.determineStuff(mcli[n]["pli"], self, opt)
            if opt.interpoint_relations["simulated - simulated"]:
                    distlis = self.getSameInterpointDistances(opt,
                                                            mcli[n]["pli"])
                    mcli[n]["simulated - simulated"]["dist"].append(distlis[0])
                    mcli[n]["simulated - simulated"]["latdist"].append(distlis[1])
            if opt.interpoint_relations["simulated - particle"]:
                distlis = self.getInterpointDistances2(opt, mcli[n]["pli"],
                                                       pli)
                mcli[n]["simulated - particle"]["dist"].append(distlis[0])
                mcli[n]["simulated - particle"]["latdist"].append(distlis[1])
            if opt.interpoint_relations["particle - simulated"]:
                distlis = self.getInterpointDistances2(opt, pli,
                                                       mcli[n]["pli"])
                mcli[n]["particle - simulated"]["dist"].append(distlis[0])
                mcli[n]["particle - simulated"]["latdist"].append(distlis[1])
        if opt.determine_clusters:
            for n, li in enumerate(mcli):
                dotProgress(n)
                mcli[n]["clusterli"] = self.determineClusters(li["pli"], opt)
                self.processClusters(mcli[n]["clusterli"], opt)
        sys.stdout.write("\n")
        return mcli



    def processClusters(self, clusterli, opt):
        for c in clusterli:
            if opt.stop_requested:
                return            
            c.convexHull = geometry.convexHullGraham(c)
            c.distToPath = c.convexHull.centroid().perpendDist(self.posel)
        for c in clusterli:
            if opt.stop_requested:
                return            
            c.nearestCluster = ClusterData()
            if len(clusterli) == 1:
                c.distToNearestCluster = -1
                return
            c.distToNearestCluster = sys.maxint
            for c2 in clusterli:
                if c2 != c:
                    d = c.lateralDistSyn(c2, self.posel)
                    if  d < c.distToNearestCluster:
                        c.distToNearestCluster = d
                        c.nearestCluster = c2


    def determineClusters(self, pointli, opt):
        """ Partition pointli into clusters; each cluster contains all points
            that are less than opt.within_cluster_dist from at least one
            other point in the cluster
        """
        clusterli = []
        for p1 in pointli:
            if opt.stop_requested:
                return []
            if p1.cluster:
                continue
            for p2 in pointli:
                if p1 != p2 and p1.dist(p2) <= geometry.toPixelUnits(
                                                    opt.within_cluster_dist,
                                                    self.pixelwidth):
                    if p2.cluster != None:
                        p1.cluster = p2.cluster
                        clusterli[p1.cluster].append(p1)
                        break
            else:
                p1.cluster = len(clusterli)
                clusterli.append(ClusterData([p1]))
        return clusterli

            
    def __parse(self, opt):
        """ Parse profile data from input file
        """
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = readFile(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        self.psdli = []
        self.holeli = []
        while li:
            s = li.pop(0).replace("\n","").strip()
            if s.split(" ")[0].upper() == "IMAGE":
                self.src_img = s.split(" ")[1]
            elif s.split(" ")[0].upper() == "SYNAPSE_ID":
                try:
                    self.ID = s.split(" ")[1]
                except IndexError, ValueError:
                    ProfileWarning(self, "Profile ID not defined or invalid")
            elif s.split(" ")[0].upper() == "COMMENT":
                try:
                    self.comment = s.split(" ", 1)[1]
                except IndexError:
                    self.comment = ''
            elif s.split(" ")[0].upper() == "PIXELWIDTH":
                try: 
                    self.pixelwidth = float(s.split(" ")[1])
                    self.metric_unit = s.split(" ")[2]               
                except ValueError:
                    raise ProfileError(self,
                                       "PIXELWIDTH is not a valid number")                   
            elif s.split(" ")[0].upper() == "POSTSYNAPTIC_PROFILE":
                try:
                    self.postsynProfile = s.split(" ", 1)[1]
                except IndexError:
                    pass
            elif s.split(" ")[0].upper() == "PRESYNAPTIC_PROFILE":
                try:
                    self.presynProfile = s.split(" ", 1)[1]
                except IndexError:
                    pass
            elif s.upper() == "POSTSYNAPTIC_ELEMENT":
                self.posel = geometry.SegmentedPath(self.__getCoords(li,
                                                "postsynaptic element"))
            elif s.upper() == "PRESYNAPTIC_ELEMENT":
                self.prsel = geometry.SegmentedPath(self.__getCoords(li,
                                                "presynaptic element"))
            elif s.upper() == "POSTSYNAPTIC_DENSITY":
                self.psdli.append(PSD(self.__getCoords(li, "PSD")))
            elif s.upper() == "HOLE":
                self.holeli.append(geometry.SegmentedPath(self.__getCoords(li,
                                                          "hole")))
            elif s.upper() == "PARTICLES":
                self.pli = PointList(self.__getCoords(li, "particle"),
                                        "particle")
            elif s.upper() == "GRID":
                self.gridli = PointList(self.__getCoords(li, "grid"),
                                           "grid")
            elif s.upper() == "RANDOM_POINTS":
                self.randomli = PointList(self.__getCoords(li, "random"),
                                             "random")
            elif s[0] != "#":          # unless specifically commented out           
                ProfileWarning(self, "Unrecognized string '" + s +
                                    "' in input file")
        # Now, let's see if everything was found
        self.__checkParsedData(opt)

    def __checkParsedData(self, opt):
        """ See if the profile data was parsed correctly, and print info on the
            parsed data to standard output.            
        """
        self.checkVarDefault('src_img', "Source image", "N/A")        
        self.checkVarDefault('ID', "Profile ID", "N/A")
        self.checkVarDefault('comment', "Comment", "")
        self.checkVarVal('metric_unit', "Metric unit", 'metric_unit', opt)        
        self.checkRequiredVar('pixelwidth', "Pixel width", self.metric_unit)
        self.checkVarDefault('postsynProfile', "Postsynaptic profile", "N/A")
        self.checkVarDefault('presynProfile', "Presynaptic profile", "N/A")
        self.checkListVar('posel', 'Postsynaptic element', 'nodes', 2)
        self.checkListVar('prsel', 'Presynaptic element', 'nodes', 2)
        self.checkTableVar('psdli', "Postsynaptic density", 
                             "Postsynaptic densities", 1, 2)
        self.checkListVar('pli', 'Particles', '', 0)        
        self.checkTableVar('holeli', "Hole", "Holes", 0, 2)
        self.checkVarExists('gridli', "Grid", 'useGrid', opt)
        self.checkVarExists('randomli', "Random points", 'useRandom', opt)
        for n, h in enumerate(self.holeli):
            if not h.isSimplePolygon():
                raise ProfileError(self,
                                   "Profile hole %d is not a simple polygon" 
                                    % (n+1))
            for n2, h2 in enumerate(self.holeli[n+1:]):
                if h.overlapsPolygon(h2):
                    raise ProfileError(self,
                                       "Profile hole %d overlaps with hole %d "
                                       % (n+1, n+n2+2))                                    

    def checkRequiredVar(self, var_to_check, var_str, post_str):
        """ Confirm that a required variable exists; else, raise ProfileError.
        """
        if not hasattr(self, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        else:            
            sys.stdout.write("  %s: %s %s\n" 
                             % (var_str, self.__dict__[var_to_check], post_str))

    def checkListLen(self, var, min_len):
        """ Returns True if var is a list and has at least min_len elements,
            else False
        """        
        return isinstance(var, list) and len(var) >= min_len
        
    def checkListVar(self, var_to_check, var_str, post_str, min_len):
        """ Confirms that var_to_check exists, is a list and has at least 
            min_len elements; if var_to_check does not exist and min_len <= 0, 
            assigns an empty list to var_to_check. Else, raise a ProfileError.
        """
        if not hasattr(self, var_to_check):
            if min_len > 0:
                raise ProfileError(self, "%s not found in input file"
                                         % var_str_1)
            else:
                self.__dict__[var_to_check] = []  
        elif not self.checkListLen(self.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = " " + post_str            
        sys.stdout.write("  %s%s: %d\n" 
                         % (var_str, post_str, 
                            len(self.__dict__[var_to_check])))
                             
    def checkTableVar(self, var_to_check, var_str_singular, var_str_plural, 
                        min_len_1, min_len_2):
        """ Confirms that var_to_check exists, is a list and has at least 
            min_len_1 elements, and that each of these has at least min_len_2 
            subelements; if var_to_check does not exist and min_len_1 <= 0, 
            assigns an empty list to var_to_check. Else, raise ProfileError.
        """    
        if not hasattr(self, var_to_check):
            if min_len_1 > 0:
                raise ProfileError(self, "%s not found in input file"
                                         % var_str_plural)
            else:
                self.__dict__[var_to_check] = []        
        elif not self.checkListLen(self.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file"
                               % var_str_plural.lower())
        else:
            for element in self.__dict__[var_to_check]:
                if not self.checkListLen(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates"
                                       % var_str_singular)                    
        sys.stdout.write("  %s: %d\n" % (var_str_plural, 
                                         len(self.__dict__[var_to_check])))
                             

    def checkVarDefault(self, var_to_check, var_str, default=""):
        """ Checks if var_to_check exists; if not, assign the default value
        to var_to_check. Never raises a ProfileError.
        """
        if not hasattr(self, var_to_check):
            self.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str, self.__dict__[var_to_check]))   
        
    
    def checkVarExists(self, var_to_check, var_str, optflag, opt):
        """ Checks for consistency between profiles with respect to the
            existence of var_to_check (i.e., var_to_check must be present 
            either in all profiles or in none).  
            
            If optflag is not set (i.e., this is the first profile), then
            set optflag to True or False depending on the existence of
            var_to_check. If optflag is already set (for consequent profiles),
            var_to_check must (if optflag is True) or must not (if optflag is 
            False) exist. If not so, raise ProfileError.
        """
        if not hasattr(opt, optflag):
            if hasattr(self, var_to_check):
                opt.__dict__[optflag] = True
            else:
                opt.__dict__[optflag] = False
        if opt.__dict__[optflag]:
            if hasattr(self, var_to_check):
                sys.stdout.write("  %s: yes\n" % var_str)
            else:
                raise ProfileError(self, "%s not found in input file" % var_str)
        elif hasattr(self, var_to_check):
            raise ProfileError(self, "%s found but not expected" % var_str)
        else:
            sys.stdout.write("  %s: no\n" % var_str)

    def checkVarVal(self, var_to_check, var_str, optvar, opt):
        """ Checks for consistency between profiles with respect to the
            value of var_to_check (i.e., var_to_check must be present and 
            have equal value in all profiles).  
            
            If optvar is not set (i.e., this is the first profile), then
            set optflag to the value of var_to_check. If optvar is already set 
            (for consequent profiles), the value of var_to_check must be equal 
            to that of optvar. If not so, raise ProfileError.
        """
        if not hasattr(self, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        if not hasattr(opt, optvar):
            opt.__dict__[optvar] = self.__dict__[var_to_check]
        elif self.__dict__[var_to_check] == opt.__dict__[optvar]:
            sys.stdout.write("  %s: %s\n" 
                             % (var_str, self.__dict__[var_to_check]))
        else:
            raise ProfileError(self, "%s value '%s'  differs from the value "
                                     "specified ('%s') in the first input file" 
                               % (var_str, self.__dict__[var_to_check],
                                  opt.__dict__[optvar]))

        
    def __checkPaths(self):
        """ Make sure that posel, prsel or psd do not intersect with themselves
        """
        def checkPath(path, s):
            for n1 in range(0, len(path)-3):
                for n2 in range(0, len(path)-1):
                    if n1 not in (n2, n2+1) and n1+1 not in (n2, n2+1):
                        if geometry.segmentIntersection(path[n1], path[n1+1],
                                                        path[n2], path[n2+1]):
                            raise ProfileError(self, "%s invalid (crosses itself)" % s)
            return True
                    
        for path in [(self.posel, "Postsynaptic element"), 
                   (self.prsel, "Presynaptic element")]:
            checkPath(*path)
        for path in self.psdli:
            checkPath(path, "PSD")
        for path in self.holeli:
            checkPath(path, "Hole")
        sys.stdout.write("  Paths are ok.\n")


    def __orientPrsel(self):
        """ The posel and prsel should be oriented the same way, so that 
            the line between posel[0] and prsel[0] does not cross the line
            between posel[-1] and prsel[-1]. If the lines do cross, reverse 
            the order of the nodes in prsel.
        """
        if geometry.segmentIntersection(self.posel[0], self.prsel[0],
                                        self.posel[-1], self.prsel[-1]):
            self.prsel.reverse()
                           

    def __getTotalPosm(self):
        """ Construct a path comprising the postsynaptic membrane of all PSDs 
            and perforations.
            Assume well-behaved PSDs.
        """

        def distToLeftEndnode(p):
            path = geometry.SegmentedPath()
            project, seg_project = p.projectOnPathOrEndnode(self.posel)
            path.extend([self.posel[0], project])
            for n in range(1, seg_project):
                path.insert(len(path)-1, self.posel[n])
            return path.length()

        left_mindist = right_mindist = self.posel.length()
        for p in itertools.chain(*[[psd[0], psd[-1]] for psd in self.psdli]):
            left_d = distToLeftEndnode(p)
            right_d = self.posel.length() - left_d
            if left_d <= left_mindist:
                left_mindist = left_d
                left_p = p
            elif right_d <= right_mindist:
                right_mindist = right_d
                right_p = p
        pseudoPSD = PSD(geometry.SegmentedPath([left_p, right_p]))
        return pseudoPSD.getPosm(self)


            
    def __getCoords(self, strli, coordType=""):
        """ Pop point coordinates from list strli.
            When an element of strli is not a valid point,
            a warning is issued.
        """
        pointli = []
        s = strli.pop(0).replace("\n","").replace(" ","").strip()
        while s != "END":
            try:
                p = geometry.Point(float(s.split(",")[0]),
                                   float(s.split(",")[1]))
                if pointli and (p == pointli[-1] or 
                                (coordType == 'particle' and p in pointli)):
                    sys.stdout.write("Duplicate %s coordinates %s: skipping "
                                     "2nd instance\n" % (coordType, p))
                else:
                    pointli.append(Point(p.x, p.y, ptype=coordType))
            except ValueError:
                if s[0] != "#":
                    ProfileWarning(self, "'%s' not valid %s coordinates"
                                   % (s, coordType))
                else:
                    pass 
            s = strli.pop(0).replace("\n","").strip()
        # For some reason, sometimes the endnodes have the same coordinates;
        # in that case, delete the last endnode to avoid division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]): 
            del pointli[-1]                                            
        return pointli        

    def __saveResults(self, opt):
        """ Output results from a single profile to file
        """ 
        
        def m(x):
            try:
                return toMetricUnits(x, self.pixelwidth)
            except ZeroDivisionError:
                return None
                
        def m2(x): 
            try:
                return toMetricUnits(x, self.pixelwidth**2) # for area units...
            except ZeroDivisionError:
                return None
       
        def fwrite(*args):
            f.writerow(args)
                    
        
        try:
            self.outputfn = os.path.join(opt.output_dir,
                                         os.path.basename(self.inputfn)
                                         + opt.output_filename_suffix
                                         + opt.output_filename_ext)
            
            if (os.path.exists(self.outputfn) and
                opt.action_if_output_file_exists == 'enumerate'):
                    self.outputfn = enumFilename(self.outputfn, 2)
            sys.stdout.write("Writing to '%s'...\n" % self.outputfn)
            if opt.output_file_format == "csv":
                csv_format = { 'dialect' : 'excel', 'lineterminator' : '\n'}
                if opt.csv_delimiter == 'tab':
                    csv_format['delimiter'] = '\t'
                f = unicode_csv.writer(file(self.outputfn, "w"),
                                        **opt.csv_format)
            elif opt.output_file_format == 'excel':
                import xls
                f = xls.writer(self.outputfn)            
            fwrite("Table 1. Profile-centric data")
            fwrite("Source image:", self.src_img)
            fwrite("Comment:",self.comment) 
            fwrite("Pixel width:", tostr(float(self.pixelwidth), 2), 
                                   self.metric_unit)
            fwrite("Postsynaptic structure:", self.postsynProfile)
            fwrite("Presynaptic structure:", self.presynProfile)
            fwrite("Length of postsynaptic element:", m(self.posel.length()))
            fwrite("Length of presynaptic element:", m(self.posel.length()))
            fwrite("Number of PSDs:", len(self.psdli))            
            fwrite("Total postsynaptic membrane length incl perforations:", 
                    m(self.totalPosm.length()))
            fwrite("Total postsynaptic membrane length excl perforations:", 
                    sum([m(psd.posm.length()) for psd in self.psdli]))
            fwrite("Total PSD area:", sum([m2(psd.psdposm.area()) 
                                           for psd in self.psdli]))
            fwrite("Number of particles (total):", len(self.pli))
            fwrite("Number of particles in PSD:",
                   len([p for p in self.pli if p.isWithinPSD]))
            fwrite("Number of particles within %s metric units of PSD:" 
                    % opt.spatial_resolution,
                   len([p for p in self.pli if p.isAssociatedWithPSD]))
            fwrite("Particles in PSD / PSD area (*1e6):", 
                   1e6 * len([p for p in self.pli if p.isWithinPSD])
                   / sum([m2(psd.psdposm.area()) for psd in self.psdli]))
           
            fwrite("Table 2. Particle-centric data")                                   
            columnheadings = ["Particle number (as appearing in input file)",
                              "Particle coordinates (in pixels)",
                              "Distance to membrane of postsynaptic element",
                              "Distance to membrane of presynaptic element",
                              "Lateral location",
                              "Lateral distance to nearest PSD",
                              "Normalized lateral distance to PSD",
                              "Particle within PSD",
                              "Particle within %(opt.spatial_resolution) "
                                  "metric units of PSD"]
            fwrite(*columnheadings) 
            f.writerows([[n+1, 
                          str(p), 
                          tostr(m(p.distToPosel), 2),
                          tostr(m(p.distToPrsel), 2),
                          p.lateralLocation,
                          m(p.lateralDistPSD), 
                          p.normLateralDistPSD,
                          yes_or_no(p.isWithinPSD),
                          yes_or_no(p.isAssociatedWithPSD)]
                          for n, p in enumerate(self.pli)])
            fwrite("Table 3. Random point-centric data")                                   
            columnheadings = ["Point number (as appearing in input file)",
                              "Point coordinates (in pixels)",
                              "Distance to membrane of postsynaptic element",
                              "Distance to membrane of presynaptic element",
                              "Lateral location",
                              "Lateral distance to nearest PSD",
                              "Normalized lateral distance to PSD",
                              "Point within PSD",
                              "Point within %(opt.spatial_resolution) "
                                  "metric units of PSD"]
            fwrite(*columnheadings) 
            f.writerows([[n+1, 
                          str(p), 
                          tostr(m(p.distToPosel), 2),
                          tostr(m(p.distToPrsel), 2),
                          p.lateralLocation,
                          m(p.lateralDistPSD), 
                          p.normLateralDistPSD,
                          yes_or_no(p.isWithinPSD),
                          yes_or_no(p.isAssociatedWithPSD)]
                          for n, p in enumerate(self.randomli)])      
            fwrite("Table 4. Grid-centric data")                                   
            columnheadings = ["Point number (as appearing in input file)",
                              "Point coordinates (in pixels)",
                              "Distance to membrane of postsynaptic element",
                              "Distance to membrane of presynaptic element",
                              "Lateral location",
                              "Lateral distance to nearest PSD",
                              "Normalized lateral distance to PSD",
                              "Point within PSD",
                              "Point within %(opt.spatial_resolution) "
                                  "metric units of PSD"]
            fwrite(*columnheadings) 
            f.writerows([[n+1, 
                          str(p), 
                          tostr(m(p.distToPosel), 2),
                          tostr(m(p.distToPrsel), 2),
                          p.lateralLocation,
                          m(p.lateralDistPSD), 
                          p.normLateralDistPSD,
                          yes_or_no(p.isWithinPSD),
                          yes_or_no(p.isAssociatedWithPSD)]
                          for n, p in enumerate(self.gridli)])      
            f.close()            
        except IOError:
            raise ProfileError(self, "Unable to write to file '%s'"
                               % self.outputfn)
        sys.stdout.write("Done.\n")
        return 1

    #def warn(profile, msg):
    #    sys.stdout.write("Warning: %s\n" % msg)
    #    profile.warnflag = True
                                                     
    def _getCleftSimplistic(self):
        """ Return the polygon outline of the synaptic cleft, by concatenating 
            posm and prsm (by joining them by two vertices across the cleft).
        """
        pol = geometry.SegmentedPath()
        if len(self.prsel) == 0:
            return pol          
        pol.extend(self.posm[:])
        if (pol[0].dist(self.prsm[0]) < pol[0].dist(self.prsm[-1])):
            pol.reverse()
        pol.extend(self.prsm)
        return pol

            
    def cleftWidthCenter(self):
        try:
            pt, foo = self.posm.centerPoint().projectOnPathOrEndnode(self.prsel)
            return abs(Point(pt.x, pt.y).perpendDist(self.posm))
        except (TypeError, IndexError):
            return None            
# end of class ProfileData
        
class OptionData:
    def __init__(self):
        self.input_file_list = []
        self.spatial_resolution = 25
        self.shell_width = 200   # Skip points farther than this from the
                                # postsynaptic element
        self.outputs = {'profile summary': True, 'point summary': True,
                        'random summary': True, 'session summary': True,
                        'individual profiles': False}
        self.output_file_format = "excel"
        self.output_filename_ext = ".xls"
        self.input_filename_ext = ".syn"
        self.output_filename_suffix = ''
        self.output_filename_other_suffix = ''
        self.output_filename_date_suffix = True
        self.output_filename_use_other_suffix = False
        self.csv_delimiter = 'comma'
        self.action_if_output_file_exists = 'overwrite'
        self.output_dir = ''
        self.determine_clusters = False
        self.within_cluster_dist = 50
        self.run_monte_carlo = False
        self.monte_carlo_runs = 99
        self.monte_carlo_simulation_window = 'whole profile'
        self.monte_carlo_strict_location = False
        self.determine_interpoint_dists = False
        self.interpoint_dist_mode = 'nearest neighbour'
        self.interpoint_relations = {'particle - particle': True,
                                    'random - particle': True,
                                    'particle - simulated': False,
                                    'simulated - particle': False,
                                    'simulated - simulated': False}
        self.interpoint_shortest_dist = True
        self.interpoint_lateral_dist = False
  
# end of class OptionData

class ProfileError(exceptions.Exception):
    def __init__(self, profile, msg):
        self.args = (profile, msg + ".")


def ProfileWarning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True

def ProfileMessage(profile, msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)
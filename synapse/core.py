import sys
import exceptions
import random
import itertools
import unicode_csv
import geometry
import file_io
import stringconv


# Convenience functions
def dot_progress(x, linelength=80, char='.'):
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
        if isinstance(x, geometry.Point):
            geometry.Point.__init__(self, x.x, x.y)
        else:
            geometry.Point.__init__(self, x, y)
        self.profile = profile
        if self.profile is not None:
            self.opt = self.profile.opt
        else:
            self.opt = None
        self.skipped = False
        self.ptype = ptype
        self.dist_to_posel = None
        self.dist_to_prsel = None
        self.lateral_location = None
        self.strict_lateral_location = None
        self.lateral_dist_psd = None
        self.norm_lateral_dist_psd = None
        self.closest_psd = None
        self.axodendritic_location = None
        self.is_within_psd = None
        self.is_associated_with_psd = None
        self.associated_psd = None
        self.is_within_postsynaptic_membrane_shell = None
        self.is_postsynaptic_membrane_associated = None
        self.is_presynaptic_membrane_associated = None
        self.cluster = None

    def get_dist_to_posel(self):
        """ Determines the perpendicular distance to the outline of the
            postsynptic element
        """
        dist_to_posel = self.perpend_dist(self.profile.posel,
                                          posloc=self.profile.postsynloc)
        if dist_to_posel is None and self.ptype == "particle":
            self.skipped = True
            if self.ptype == "particle":
                profile_warning(self.profile,
                                "Unable to project on postsynaptic element\n"
                                "  => skipping particle at %s" % self)
        return dist_to_posel

    def get_dist_to_prsel(self):
        """ Determines the perpendicular distance to the outline of the
            presynaptic element
        """
        if len(self.profile.prsel) > 0:
            return self.perpend_dist(self.profile.prsel, 
                                     posloc=self.profile.postsynloc)
        else:
            return None

    def get_lateral_location(self):
        """ Determines the lateral location, defined as the location
            along the mediolateral synaptic axis of the projection of
            the point on the postsynaptic membrane.
        """
        # A point is classified as synaptic if its projection on the
        # postsynaptic membrane is within the spatial resolution limit 
        # of the postsynaptic membrane
        for psd in self.profile.psdli:
            if (self.lateral_dist_syn(self.profile.posel, psd.posm) -
                    self.profile.spatial_resolution_in_pixels <=
                    psd.posm.length() / 2.0):
                return "synaptic"
        # If not synaptic but still within the extent of the synaptic
        # membrane
        if (self.lateral_dist_syn(self.profile.posel, self.profile.total_posm) <
                (self.profile.total_posm.length() / 2.0)):
            return "within perforation"
        # Otherwise it may be perisynaptic, defined here as being
        # outside the synapse but within half a diameter of the nearest
        # PSD... or perhaps it should be a fixed distance, e g 200 nm?
        # I don't want to determine which PSD is closest and which edge
        # of that PSD faces the extrasynapse, so I just determine the
        # distance to the nearest edge; if that is an internal edge,
        # the point will have been identified as within perforation
        # above.
        for psd in self.profile.psdli:
            if (self.lateral_dist_syn(self.profile.posel, psd.posm) <=
                    psd.posm.length()):
                return "perisynaptic"
        # If none of the above
        return "extrasynaptic"

    def get_strict_lateral_location(self):
        """ Determines the lateral location, defined as the location
            along the mediolateral synaptic axis of the projection of
            the particle on the postsynaptic membrane. Does not take
            spatial resolution into account, such that a particle is
            considered synaptic only when strictly projectable on the
            postsynaptic membrane.
        """
        # A point is classified as synaptic if its projection on the
        # postsynaptic membrane is within the postsynaptic membrane
        for psd in self.profile.psdli:
            if (self.lateral_dist_syn(self.profile.posel, psd.posm) <=
                    psd.posm.length() / 2.0):
                return "synaptic"
        # If not synaptic but still within the extent of the synaptic
        # membrane
        if (self.lateral_dist_syn(self.profile.posel, self.profile.total_posm) <
                (self.profile.total_posm.length() / 2.0)):
            return "within perforation"
        # Otherwise it may be perisynaptic, defined here as being
        # outside the synapse but within half a diameter of the nearest
        # PSD... or perhaps it should be a fixed distance, e g 200 nm?
        # I don't want to determine which PSD is closest and which edge
        # of that PSD faces the extrasynapse, so I just determine the
        # distance to the nearest edge; if that is an internal edge,
        # the point will have been identified as within perforation
        # above.
        for psd in self.profile.psdli:
            if (self.lateral_dist_syn(self.profile.posel, psd.posm) <=
                    psd.posm.length()):
                return "perisynaptic"
        # If none of the above
        return "extrasynaptic"

    def _axodendritic_location(self):
        """ Determines the particle's location along the axodendritic axis.
        """
        if self.dist_to_posel >= 0:
            return "postsynaptic"
        elif self.dist_to_prsel is None:  # if there is no presynaptic membrane,
            return "not postsynaptic"  # that's about all we can say
        elif self.dist_to_prsel <= 0:
            return "presynaptic"
        else:
            return "neither pre- or postsynaptic"

    def _lateral_dist_psd(self):
        mindist = sys.maxint
        closest_psd = None
        for psd in self.profile.psdli:
            d = (self.lateral_dist_syn(self.profile.posel, psd.posm) -
                 (psd.posm.length() / 2))
            if d < mindist:
                mindist = d
                closest_psd = psd
        if not closest_psd:
            raise ProfileError(self.profile,
                               "could not determine lateral distance "
                               "to a PSD of particle at %s" % self)
        mindist = self.lateral_dist_syn(self.profile.posel, closest_psd.posm)
        normdist = mindist / (closest_psd.posm.length() / 2)
        return mindist, normdist, closest_psd

    def _psd_association(self):
        is_within_psd = False
        is_associated_with_psd = False
        associated_psd = None
        mindist = sys.maxint
        for psd in self.profile.psdli:
            if self.is_within_polygon(psd.psdposm):
                is_within_psd = True
                is_associated_with_psd = True
                associated_psd = psd
                break
            dist = self.perpend_dist_closed_path(psd.psdposm, 
                                                 dont_care_if_on_or_off_seg=
                                                 True)
            if dist <= self.profile.spatial_resolution_in_pixels:
                is_associated_with_psd = True
                if dist < mindist:
                    associated_psd = psd
                    mindist = dist
        return is_within_psd, is_associated_with_psd, associated_psd

    def _is_within_postsynaptic_membrane_shell(self):
        return (self.dist_to_posel is not None
                and abs(self.dist_to_posel) <=
                geometry.to_pixel_units(self.opt.shell_width,
                                        self.profile.pixelwidth))

    def _membrane_association(self):
        posm_associated = False
        presm_associated = False
        if (self.dist_to_posel is not None and abs(self.dist_to_posel) <=
                self.profile.spatial_resolution_in_pixels):
            posm_associated = True
        if (self.dist_to_prsel is not None and abs(self.dist_to_prsel) <=
                self.profile.spatial_resolution_in_pixels):
            presm_associated = True
        return posm_associated, presm_associated

    def is_within_hole(self):
        """  Determine whether self is inside a profile hole
        """
        for h in self.profile.holeli:
            if self.is_within_polygon(h):
                return True
        return False

    def determine_stuff(self):
        if self.is_within_hole():
            # no warning for random and grid points
            if self.ptype == "particle":
                profile_warning(self.profile, "Discarding point at %s: Located "
                                              "within a profile hole" % self)
            self.skipped = True
            return
        self.dist_to_posel = self.get_dist_to_posel()
        if self.dist_to_posel > geometry.to_pixel_units(
                self.opt.shell_width, self.profile.pixelwidth):
            if self.ptype == "particle":
                profile_message("Discarding particle at %s: Located outside "
                                "shell" % self)
            self.skipped = True
            return
        self.dist_to_prsel = self.get_dist_to_prsel()
        self.lateral_location = self.get_lateral_location()
        self.strict_lateral_location = self.get_strict_lateral_location()
        self.lateral_dist_psd, self.norm_lateral_dist_psd, \
            self.closest_psd = self._lateral_dist_psd()
        self.axodendritic_location = self._axodendritic_location()
        self.is_within_psd, self.is_associated_with_psd, \
            self.associated_psd = self._psd_association()
        self.is_within_postsynaptic_membrane_shell = \
            self._is_within_postsynaptic_membrane_shell()
        self.is_postsynaptic_membrane_associated, \
            self.is_presynaptic_membrane_associated = \
            self._membrane_association()

    def determine_nearest_neighbour(self, pointli):
        # Assumes that only valid (projectable, within shell etc) points 
        # are in pointli
        mindist = float(sys.maxint)
        for p in pointli:
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
        if not mindist < float(sys.maxint):
            return None
        else:
            nearest_neighbour_dist = mindist
            return nearest_neighbour_dist

    def determine_nearest_lateral_neighbour(self, pointli):
        # Assumes that only valid (projectable, within shell etc) points 
        # are in pointli
        mindist = float(sys.maxint)
        for p in pointli:
            if p is not self:
                d = self.lateral_dist_to_point(p, self.profile.posel)
                if d < mindist:
                    mindist = d
        if not mindist < float(sys.maxint):
            return None
        else:
            nearest_lateral_neighbour_dist = mindist
            return nearest_lateral_neighbour_dist

    def lateral_dist_to_point(self, p2, sel):
        """ Determine lateral distance to a point p2 along sel.
            Overrides function in geometry.Point, which only works
            with a closed path.
        """
        path = geometry.OpenPath()
        p2_project, p2_seg_project = p2.project_on_path_or_endnode(sel)
        project, seg_project = self.project_on_path_or_endnode(sel)
        path.extend([project, p2_project])
        if p2_seg_project < seg_project:
            path.reverse()
        for n in range(min(p2_seg_project, seg_project) + 1,
                       max(p2_seg_project, seg_project)):
            path.insert(len(path) - 1, sel[n])
        return path.length()

    def lateral_dist_syn(self, sel, sm):
        """ Determine lateral distance to center of (post- or pre-)
            synaptic membrane. 
        """
        return self.lateral_dist_to_point(sm.center_point(), sel)


class PointList(list):
    def __init__(self, pointli, ptype):
        try:
            self.extend([Point(p.x, p.y, ptype) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError('not a list of Point elements')
        

class PSD(geometry.ClosedPath):
    def __init__(self, pointli, profile):
        super(geometry.ClosedPath, self).__init__(pointli)
        self.profile = profile
        self.posm = geometry.OpenPath()
        self.prsm = geometry.OpenPath()
        self.psdposm = geometry.ClosedPath()
        self.cleft = geometry.ClosedPath()
        self.cleft_width = None

    def adjust_psd(self):
        """ Adjust PSD coordinates so that the PSD ends exactly on 
            the postsynaptic element. 
        """
        # Partition psd into paths defined by the intersections with posel,
        # beginning and ending with the projections of the end nodes of psd.
        # Each path will then be completely on one side of posel.
        pathli = [geometry.OpenPath()]
        pathli[0].append(self[0].project_on_path_or_endnode(
            self.profile.posel)[0])
        p = geometry.Point()
        for n in range(0, len(self) - 1):
            pathli[-1].append(self[n])
            for d in range(len(self.profile.posel) - 1):
                p = geometry.segment_intersection(self[n], self[n + 1],
                                                  self.profile.posel[d],
                                                  self.profile.posel[d + 1])
                if p:
                    break
            if p:
                pathli[-1].append(p)
                pathli.append(geometry.OpenPath())
                pathli[-1].append(p)
        pathli[-1].append(self[-1])
        pathli[-1].append(self[-1].project_on_path_or_endnode(
            self.profile.posel)[0])
        # Now, look for the longest path. This is assumed to be the intended
        # psd. (Perhaps area is more relevant. However, this is messier because
        # we need to determine the part of dm enclosed by path for each path.)
        maxlength = 0
        longest_path = geometry.OpenPath()
        for path in pathli:
            length = path.length()
            if length > maxlength:
                longest_path = path
                maxlength = length
        del self[:]
        self.extend(longest_path)
        return self

    def get_posm(self):
        """ Return postsynaptic membrane (i e the part of the 
            postsynaptic element which is bounded by the PSD). Assume
            that the PSD is adjusted so as to end exactly on 
            postsynaptic element, and that the PSD is oriented in the
            same direction as posm (ie node1 <= node2).
        """
        posm = geometry.OpenPath()
        p1, node1 = self[0].project_on_path_or_endnode(self.profile.posel)
        p2, node2 = self[-1].project_on_path_or_endnode(self.profile.posel)
        if not None in (node1, node2):
            posm.extend(self.profile.posel[node1 + 1:node2 + 1])
            posm.insert(0, p1)
            posm.append(p2)
        if len(posm) == 0 or posm.length() == 0:
            raise ProfileError(self.profile, "Could not determine postsynaptic "
                                             "membrane")
        # It appears that sometimes the first or last nodes get 
        # duplicated, resulting in a zero vector and division by zero 
        # when projecting on posm. Simply checking for the duplicated
        # node and removing it if present obviously works, but perhaps
        # I should delve deeper into the problem. The node should not 
        # be duplicated in the first place. 
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

    def get_prsm(self):
        prsm = geometry.OpenPath()
        p1, node1 = self[0].project_on_path_or_endnode(self.profile.prsel)
        p2, node2 = self[-1].project_on_path_or_endnode(self.profile.prsel)
        if not None in (node1, node2):
            # Now we have the projections of the psd endnodes on prsel.
            # But if the prsel is far from a psd endnode, and the psd 
            # is very curved, the projection might not coincide with 
            # the prsel endnode, even if the prsel endnode is apposed
            # to the psd. Therefore, see if the prsel endnode 
            # projections are on the posm, and if these are closer than
            # the projections of psd endnodes on prsel. If so, use them
            # instead. Phew! (04/03/2003)
            p0, foo = self.profile.prsel[0].project_on_path(
                self.posm)  # should be self.posm? 10.11.07 / ML
            p_n, foo = self.profile.prsel[-1].project_on_path(self.posm)
            if p0 and p0.dist(p1) < self[0].dist(p1):
                p1 = self.profile.prsel[0]
                node1 = 0
            if p_n and p_n.dist(p2) < self[-1].dist(p2):
                p2 = self.profile.prsel[-1]
                node2 = len(self.profile.prsel)
            prsm.extend(self.profile.prsel[node1 + 1:node2 + 1])
            prsm.insert(0, p1)
            prsm.append(p2)
        if len(prsm) == 0 or prsm.length() == 0:
            raise ProfileError(self.profile, "Could not determine presynaptic "
                                             "membrane")
        return prsm

    def get_psd_posm(self):
        """ Make a closed-path PSD by concatenating the PSD and the 
            posm. Assume that PSD and posm share end nodes and that 
            they are oriented in the same direction. 
        """
        pol = geometry.ClosedPath()
        pol.extend(self.posm)
        pol.reverse()
        pol.extend(self[1:-1])
        return pol

    def get_cleft(self):
        """ Return the polygon outline of the synaptic cleft, defined 
            as the polygon outlining the area bounded by the prsm, the 
            segments formed by projecting the end nodes of prsm onto 
            posm, and the part of posm between these projections.  
        """
        pol = geometry.ClosedPath()
        if len(self.profile.prsel) == 0:
            return pol
        proj0, seg0 = self.prsm[0].project_on_path_or_endnode(self.posm)
        projn, segn = self.prsm[-1].project_on_path_or_endnode(self.posm)
        pol.extend(self.prsm)
        if seg0 < segn:
            pol.reverse()
            pol.append(proj0)
            pol.extend(self.posm[seg0 + 1:segn + 1])
            pol.append(projn)
        else:
            pol.append(projn)
            pol.extend(self.posm[segn + 1:seg0 + 1])
            pol.append(proj0)
        return pol

    def cleft_width_average(self):
        """ Return synaptic cleft width, calculated by dividing the area 
            of the cleft outline by the length of posm. 
        """
        try:
            return self.cleft.area() / self.posm.length()
        except (TypeError, IndexError):
            return None


# End of class PSD

class ClusterData(list):
    def __init__(self, pli=()):
        try:
            self.extend([Point(p.x, p.y) for p in pli])
        except (AttributeError, IndexError):
            raise TypeError('not a point list')
        self.convex_hull = geometry.SegmentedPath()

    def lateral_dist_syn(self, c2, mem):
        """ Determine lateral distance to a cluster c2 along postsynaptic
            element membrane.
        """
        centroid = Point(self.convex_hull.centroid())
        centroid2 = Point(c2.convex_hull.centroid())
        return centroid.lateral_dist_to_point(centroid2, mem)


class ProfileData:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.opt = opt
        self.outputfn = ""
        self.warnflag = False
        self.errflag = False
        self.spatial_resolution_in_pixels = None
        self.total_posm = geometry.OpenPath()
        self.postsynloc = geometry.Point()
        self.pli = []
        self.gridli = []
        self.randomli = []

    def process(self):
        """ Parse profile data from a file and determine distances
        """
        try:
            self._parse()
            self._check_paths()
            sys.stdout.write("Processing profile...\n")
            self.spatial_resolution_in_pixels = geometry.to_pixel_units(
                self.opt.spatial_resolution,
                self.pixelwidth)
            self.prsel.orient_to_path(self.posel)
            for p in self.psdli:
                p.adjust_psd()
                p.orient_to_path(self.posel)
                p.psdloc = p.center_point()
                p.posm = p.get_posm()
                p.prsm = p.get_prsm()
                p.psdposm = p.get_psd_posm()
                p.cleft = p.get_cleft()
                p.cleft_width = p.cleft_width_average()
            self.total_posm = self._get_total_posm()
            self.postsynloc = self.psdli[0].psdloc
            for p in self.pli:
                p.determine_stuff(self.pli)
            self.pli = [p for p in self.pli if not p.skipped]
            if self.opt.useGrid:
                for g in self.gridli:
                    g.determine_stuff(self.gridli)
                self.gridli = [g for g in self.gridli if not g.skipped]
            if self.opt.useRandom:
                for r in self.randomli:
                    r.determine_stuff(self.randomli)
                self.randomli = [r for r in self.randomli if not r.skipped]
            self._get_inter_distlis()
            self._get_clusters()
            self._get_monte_carlo()
            if self.opt.stop_requested:
                return
            sys.stdout.write("Done.\n")
            if self.opt.outputs['individual profiles']:
                self._save_results()
        except ProfileError, (self, msg):
            sys.stdout.write("Error: %s\n" % msg)
            self.errflag = True

    def _get_monte_carlo(self):
        if self.opt.run_monte_carlo is False:
            self.mcli = []
        else:
            sys.stdout.write("Running Monte Carlo simulations...\n")
            self.mcli = self._run_monte_carlo()

    def _get_inter_distlis(self):
        if not self.opt.determine_interpoint_dists:
            return
        if not True in [val for key, val in
                        self.opt.interpoint_relations.items()
                        if "simulated" not in key]:
            return
        sys.stdout.write("Determining interpoint distances...\n")
        if self.opt.interpoint_relations["particle - particle"]:
            self.pp_distli, self.pp_latdistli = \
                self._get_same_interpoint_distances(self.pli)
        if self.opt.useRandom and \
                self.opt.interpoint_relations["random - particle"]:
            self.rp_distli, self.rp_latdistli = \
                self._get_interpoint_distances2(self.randomli, self.pli)

    def _get_clusters(self):
        if not (self.opt.determine_clusters and
                self.opt.within_cluster_dist > 0):
            return
        sys.stdout.write("Determining clusters...\n")
        self.clusterli = self._determine_clusters(self.pli)
        self._process_clusters(self.clusterli)

    def _get_same_interpoint_distances(self, pointli):
        dli = []
        latdli = []
        for i in range(0, len(pointli)):
            if self.opt.stop_requested:
                return [], []
            if self.opt.interpoint_dist_mode == 'all':
                for j in range(i + 1, len(pointli)):
                    if self.opt.interpoint_shortest_dist:
                        dli.append(pointli[i].dist(pointli[j]))
                    if self.opt.interpoint_lateral_dist:
                        latdli.append(pointli[i].lateral_dist_to_point(
                            pointli[j],
                            self.posel))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(pointli[i].determine_nearest_neighbour(pointli))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(
                        pointli[i].determine_nearest_lateral_neighbour(pointli))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def _get_interpoint_distances2(self, pointli, pointli2=()):
        dli = []
        latdli = []
        for p in pointli:
            if self.opt.stop_requested:
                return [], []
            if self.opt.interpoint_dist_mode == 'all':
                for p2 in pointli2:
                    if self.opt.interpoint_shortest_dist:
                        dli.append(p.dist(p2))
                    if self.opt.interpoint_lateral_dist:
                        latdli.append(p.lateral_dist_to_point(p2, self.posel))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(p.determine_nearest_neighbour(pointli2))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(
                        p.determine_nearest_lateral_neighbour(pointli2))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def _run_monte_carlo(self):
        
        def is_valid(rand_p):
            d = rand_p.perpend_dist(self.posel)
            if (d is None or abs(d) > shell_width or rand_p.is_within_hole(self)
                    or rand_p in mcli[n]["pli"]):
                return False
            if self.opt.monte_carlo_simulation_window == "whole profile":
                return True
            if self.opt.monte_carlo_strict_location:
                latloc = rand_p.get_strict_lateral_location(self)
            else:
                latloc = rand_p.get_lateral_location(self)
            if (self.opt.monte_carlo_simulation_window == "synapse" and
                    latloc in ("synaptic", "within perforation")):
                return True
            if (self.opt.monte_carlo_simulation_window ==
                    "synapse - perforations" and latloc == "synaptic"):
                return True
            if (self.opt.monte_carlo_simulation_window ==
                    "synapse + perisynapse" and latloc in
                    ("synaptic", "within perforation", "perisynaptic")):
                return True
            # TODO : include PSD as a window option
            return False

        pli = self.pli
        if self.opt.monte_carlo_strict_location:
            locvar = "strict_lateral_location"
        else:
            locvar = "lateral_location"
        if self.opt.monte_carlo_simulation_window == "whole profile":
            # particles outside shell have been discarded
            numpoints = len(pli)
        elif self.opt.monte_carlo_simulation_window == "synapse":
            numpoints = len([p for p in pli if p.__dict__[locvar] in
                             ("synaptic", "within perforation")])
        elif self.opt.monte_carlo_simulation_window == "synapse - perforations":
            numpoints = len([p for p in pli
                             if p.__dict__[locvar] == "synaptic"])
        elif self.opt.monte_carlo_simulation_window == "synapse + perisynapse":
            numpoints = len([p for p in pli if p.__dict__[locvar] in
                             ("synaptic", "within perforation",
                              "perisynaptic")])
        else:
            return []
        allpaths = geometry.SegmentedPath(
            [p for path in (self.posel, self.prsel,
                            [q for psd in self.psdli for q in psd],
                            [r for h in self.holeli for r in h])
             for p in path])
        box = allpaths.bounding_box()
        shell_width = geometry.to_pixel_units(self.opt.shell_width,
                                              self.pixelwidth)
        mcli = []
        for n in range(0, self.opt.monte_carlo_runs):
            if self.opt.stop_requested:
                return []
            dot_progress(n)
            mcli.append({"pli": [],
                         "simulated - simulated": {"dist": [], "latdist": []},
                         "simulated - particle": {"dist": [], "latdist": []},
                         "particle - simulated": {"dist": [], "latdist": []},
                         "clusterli": []})
            p = Point()
            for _ in range(0, numpoints):
                while True:
                    x = random.randint(int(box[0].x - shell_width),
                                       int(box[1].x + shell_width) + 1)
                    y = random.randint(int(box[0].y - shell_width),
                                       int(box[2].y + shell_width) + 1)
                    p = Point(x, y)
                    if is_valid(p):
                        break
                # escape the while loop when a valid simulated point 
                # is found
                mcli[n]["pli"].append(p)
            for p in mcli[n]["pli"]:
                p.determine_stuff(mcli[n]["pli"])
            if self.opt.interpoint_relations["simulated - simulated"]:
                distlis = self._get_same_interpoint_distances(mcli[n]["pli"])
                mcli[n]["simulated - simulated"]["dist"].append(distlis[0])
                mcli[n]["simulated - simulated"]["latdist"].append(distlis[1])
            if self.opt.interpoint_relations["simulated - particle"]:
                distlis = self._get_interpoint_distances2(mcli[n]["pli"], pli)
                mcli[n]["simulated - particle"]["dist"].append(distlis[0])
                mcli[n]["simulated - particle"]["latdist"].append(distlis[1])
            if self.opt.interpoint_relations["particle - simulated"]:
                distlis = self._get_interpoint_distances2(pli, mcli[n]["pli"])
                mcli[n]["particle - simulated"]["dist"].append(distlis[0])
                mcli[n]["particle - simulated"]["latdist"].append(distlis[1])
        if self.opt.determine_clusters:
            for n, li in enumerate(mcli):
                dot_progress(n)
                mcli[n]["clusterli"] = self._determine_clusters(li["pli"])
                self._process_clusters(mcli[n]["clusterli"])
        sys.stdout.write("\n")
        return mcli

    def _process_clusters(self, clusterli):
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.convex_hull = geometry.convex_hull(c)
            c.distToPath = c.convex_hull.centroid().perpend_dist(self.posel)
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.nearestCluster = ClusterData()
            if len(clusterli) == 1:
                c.distToNearestCluster = -1
                return
            c.distToNearestCluster = sys.maxint
            for c2 in clusterli:
                if c2 != c:
                    d = c.lateral_dist_syn(c2, self.posel)
                    if d < c.distToNearestCluster:
                        c.distToNearestCluster = d
                        c.nearestCluster = c2

    def _determine_clusters(self, pointli):
        """ Partition pointli into clusters; each cluster contains all
            points that are less than opt.within_cluster_dist from at 
            least one other point in the cluster
        """
        clusterli = []
        for p1 in pointli:
            if self.opt.stop_requested:
                return []
            if p1.cluster:
                continue
            for p2 in pointli:
                if p1 != p2 and p1.dist(p2) <= geometry.to_pixel_units(
                        self.opt.within_cluster_dist, self.pixelwidth):
                    if p2.cluster is not None:
                        p1.cluster = p2.cluster
                        clusterli[p1.cluster].append(p1)
                        break
            else:
                p1.cluster = len(clusterli)
                clusterli.append(ClusterData([p1]))
        return clusterli

    def _parse(self):
        """ Parse profile data from input file
        """
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = file_io.read_file(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        self.psdli = []
        self.holeli = []
        while li:
            s = li.pop(0).replace("\n", "").strip()
            if s.split(" ")[0].upper() == "IMAGE":
                self.src_img = s.split(" ")[1]
            elif s.split(" ")[0].upper() == "SYNAPSE_ID":
                try:
                    self.ID = s.split(" ")[1]
                except (IndexError, ValueError):
                    profile_warning(self, "Profile ID not defined or invalid")
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
                    raise ProfileError(self, "PIXELWIDTH is not a valid number")
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
                self.posel = geometry.OpenPath(
                    self._get_coords(li, "postsynaptic element"))
            elif s.upper() == "PRESYNAPTIC_ELEMENT":
                self.prsel = geometry.OpenPath(
                    self._get_coords(li, "presynaptic element"))
            elif s.upper() == "POSTSYNAPTIC_DENSITY":
                self.psdli.append(PSD(self._get_coords(li, "PSD"), self))
            elif s.upper() == "HOLE":
                self.holeli.append(geometry.ClosedPath(
                    self._get_coords(li, "hole")))
            elif s.upper() == "PARTICLES":
                self.pli = PointList(
                    self._get_coords(li, "particle"), "particle")
            elif s.upper() == "GRID":
                self.gridli = PointList(
                    self._get_coords(li, "grid"), "grid")
            elif s.upper() == "RANDOM_POINTS":
                self.randomli = PointList(
                    self._get_coords(li, "random"), "random")
            elif s[0] != "#":  # unless specifically commented out           
                profile_warning(self, "Unrecognized string '" + s +
                                "' in input file")
        # Now, let's see if everything was found
        self._check_parsed_data()

    def _check_parsed_data(self):
        """ See if the profile data was parsed correctly, and print 
            info on the parsed data to standard output.            
        """
        self._check_var_default('src_img', "Source image", "N/A")
        self._check_var_default('ID', "Profile ID", "N/A")
        self._check_var_default('comment', "Comment", "")
        self._check_var_val('metric_unit', "Metric unit", 'metric_unit')
        self._check_required_var('pixelwidth', "Pixel width", self.metric_unit)
        self._check_var_default('postsynProfile', "Postsynaptic profile", "N/A")
        self._check_var_default('presynProfile', "Presynaptic profile", "N/A")
        self._check_list_var('posel', 'Postsynaptic element', 'nodes', 2)
        self._check_list_var('prsel', 'Presynaptic element', 'nodes', 2)
        self._check_table_var('psdli', "Postsynaptic density",
                              "Postsynaptic densities", 1, 2)
        self._check_list_var('pli', 'Particles', '', 0)
        self._check_table_var('holeli', "Hole", "Holes", 0, 2)
        self._check_var_exists('gridli', "Grid", 'useGrid')
        self._check_var_exists('randomli', "Random points", 'useRandom')
        for n, h in enumerate(self.holeli):
            if not h.isSimplePolygon():
                raise ProfileError(self,
                                   "Profile hole %d is not a simple polygon"
                                   % (n + 1))
            for n2, h2 in enumerate(self.holeli[n + 1:]):
                if h.overlapsPolygon(h2):
                    raise ProfileError(self,
                                       "Profile hole %d overlaps with hole %d "
                                       % (n + 1, n + n2 + 2))

    def _check_required_var(self, var_to_check, var_str, post_str):
        """ Confirm that a required variable exists; else, raise 
            ProfileError.
        """
        if not hasattr(self, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        else:
            sys.stdout.write("  %s: %s %s\n" % (var_str, 
                                                self.__dict__[var_to_check], 
                                                post_str))

    @staticmethod
    def _check_list_len(var, min_len):
        """ Returns True if var is a list and has at least min_len
            elements, else False
        """
        return isinstance(var, list) and len(var) >= min_len

    def _check_list_var(self, var_to_check, var_str, post_str, min_len):
        """ Confirms that var_to_check exists, is a list and has at 
            least min_len elements; if var_to_check does not exist and
            min_len <= 0, assigns an empty list to var_to_check. Else, 
            raise a ProfileError.
        """
        if not hasattr(self, var_to_check):
            if min_len > 0:
                raise ProfileError(self, "%s not found in input file"
                                   % var_str)
            else:
                self.__dict__[var_to_check] = []
        elif not self._check_list_len(self.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = " " + post_str
        sys.stdout.write("  %s%s: %d\n" % (var_str, post_str,
                                           len(self.__dict__[var_to_check])))

    def _check_table_var(self, var_to_check, var_str_singular, var_str_plural,
                         min_len_1, min_len_2):
        """ Confirms that var_to_check exists, is a list and has at
            least min_len_1 elements, and that each of these has at 
            least min_len_2 subelements; if var_to_check does not exist
            and min_len_1 <= 0, assigns an empty list to var_to_check. 
            Else, raise ProfileError.
        """
        if not hasattr(self, var_to_check):
            if min_len_1 > 0:
                raise ProfileError(self, "%s not found in input file"
                                   % var_str_plural)
            else:
                self.__dict__[var_to_check] = []
        elif not self._check_list_len(self.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file"
                               % var_str_plural.lower())
        else:
            for element in self.__dict__[var_to_check]:
                if not self._check_list_len(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates"
                                       % var_str_singular)
        sys.stdout.write("  %s: %d\n" % (var_str_plural, 
                                         len(self.__dict__[var_to_check])))

    def _check_var_default(self, var_to_check, var_str, default=""):
        """ Checks if var_to_check exists; if not, assign the default value
        to var_to_check. Never raises a ProfileError.
        """
        if not hasattr(self, var_to_check):
            self.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str, self.__dict__[var_to_check]))

    def _check_var_exists(self, var_to_check, var_str, optflag):
        """ Checks for consistency between profiles with respect to the
            existence of var_to_check (i.e., var_to_check must be present 
            either in all profiles or in none).  
            
            If optflag is not set (i.e., this is the first profile), then
            set optflag to True or False depending on the existence of
            var_to_check. If optflag is already set (for consequent profiles),
            var_to_check must (if optflag is True) or must not (if optflag is 
            False) exist. If not so, raise ProfileError.
        """
        if not hasattr(self.opt, optflag):
            if hasattr(self, var_to_check):
                self.opt.__dict__[optflag] = True
            else:
                self.opt.__dict__[optflag] = False
        if self.opt.__dict__[optflag]:
            if hasattr(self, var_to_check):
                sys.stdout.write("  %s: yes\n" % var_str)
            else:
                raise ProfileError(self, "%s not found in input file" % var_str)
        elif hasattr(self, var_to_check):
            raise ProfileError(self, "%s found but not expected" % var_str)
        else:
            sys.stdout.write("  %s: no\n" % var_str)

    def _check_var_val(self, var_to_check, var_str, optvar):
        """ Checks for consistency between profiles with respect to the
            value of var_to_check (i.e., var_to_check must be present 
            and have equal value in all profiles).  
            
            If optvar is not set (i.e., this is the first profile), 
            then set optflag to the value of var_to_check. If optvar is
            already set (for consequent profiles), the value of 
            var_to_check must be equal to that of optvar. If not so, 
            raise ProfileError.
        """
        if not hasattr(self, var_to_check):
            raise ProfileError(self, "%s not found in input file" % var_str)
        if not hasattr(self.opt, optvar):
            self.opt.__dict__[optvar] = self.__dict__[var_to_check]
        elif self.__dict__[var_to_check] == self.opt.__dict__[optvar]:
            sys.stdout.write("  %s: %s\n" % (var_str, 
                                             self.__dict__[var_to_check]))
        else:
            raise ProfileError(self, "%s value '%s'  differs from the value "
                                     "specified ('%s') in the first input file"
                               % (var_str, self.__dict__[var_to_check],
                                  self.opt.__dict__[optvar]))

    def _check_paths(self):
        """ Make sure that posel, prsel or psd do not intersect with
            themselves
        """

        def check_path(_path, s):
            for n1 in range(0, len(_path) - 3):
                for n2 in range(0, len(_path) - 1):
                    if n1 not in (n2, n2 + 1) and n1 + 1 not in (n2, n2 + 1):
                        if geometry.segment_intersection(_path[n1],
                                                         _path[n1 + 1],
                                                         _path[n2],
                                                         _path[n2 + 1]):
                            raise ProfileError(
                                self, "%s invalid (crosses itself)" % s)
            return True

        for path in [(self.posel, "Postsynaptic element"),
                     (self.prsel, "Presynaptic element")]:
            check_path(*path)
        for path in self.psdli:
            check_path(path, "PSD")
        for path in self.holeli:
            check_path(path, "Hole")
        sys.stdout.write("  Paths are ok.\n")

    def __orient_prsel(self):
        """ The posel and prsel should be oriented the same way, so 
            that the line between posel[0] and prsel[0] does not cross 
            the line between posel[-1] and prsel[-1]. If the lines do 
            cross, reverse the order of the nodes in prsel.
        """
        if geometry.segment_intersection(self.posel[0], self.prsel[0],
                                         self.posel[-1], self.prsel[-1]):
            self.prsel.reverse()

    def _get_total_posm(self):
        """ Construct a path comprising the postsynaptic membrane of all 
            PSDs and perforations.

            Assume well-behaved PSDs.
        """

        def dist_to_left_endnode(_p):
            path = geometry.OpenPath()
            project, seg_project = _p.project_on_path_or_endnode(self.posel)
            path.extend([self.posel[0], project])
            for n in range(1, seg_project):
                path.insert(len(path) - 1, self.posel[n])
            return path.length()

        left_mindist = right_mindist = self.posel.length()
        left_p = Point()
        right_p = Point()
        for p in itertools.chain(*[[psd[0], psd[-1]] for psd in self.psdli]):
            left_d = dist_to_left_endnode(p)
            right_d = self.posel.length() - left_d
            if left_d <= left_mindist:
                left_mindist = left_d
                left_p = p
            elif right_d <= right_mindist:
                right_mindist = right_d
                right_p = p
        pseudo_psd = PSD(geometry.SegmentedPath([left_p, right_p]), self)
        return pseudo_psd.get_posm()

    def _get_coords(self, strli, coord_type=""):
        """ Pop point coordinates from list strli.
            When an element of strli is not a valid point,
            a warning is issued.
        """
        pointli = []
        s = strli.pop(0).replace("\n", "").replace(" ", "").strip()
        while s != "END":
            try:
                p = geometry.Point(float(s.split(",")[0]),
                                   float(s.split(",")[1]))
                if pointli and (p == pointli[-1] or
                                (coord_type == 'particle' and p in pointli)):
                    sys.stdout.write(
                        "Duplicate %s coordinates %s: skipping "
                        "2nd instance\n" % (coord_type, p))
                else:
                    pointli.append(Point(p.x, p.y, ptype=coord_type))
            except ValueError:
                if s[0] != "#":
                    profile_warning(self, "'%s' not valid %s coordinates"
                                    % (s, coord_type))
                else:
                    pass
            s = strli.pop(0).replace("\n", "").strip()
        # For some reason, sometimes the endnodes have the same
        # coordinates; in that case, delete the last endnode to avoid
        # division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]):
            del pointli[-1]
        return pointli

    def _save_results(self):
        """ Output results from a single profile to file
        """

        def m(x):
            try:
                return geometry.to_metric_units(x, self.pixelwidth)
            except ZeroDivisionError:
                return None

        def m2(x):
            # for area units...
            try:
                return geometry.to_metric_units(x, self.pixelwidth ** 2)
            except ZeroDivisionError:
                return None

        def fwrite(*args):
            f.writerow(args)

        try:
            self.outputfn = \
                file_io.os.path.join(self.opt.output_dir,
                                     file_io.os.path.basename(self.inputfn) +
                                     self.opt.output_filename_suffix +
                                     self.opt.output_filename_ext)

            if (file_io.os.path.exists(self.outputfn) and
                    self.opt.action_if_output_file_exists == 'enumerate'):
                self.outputfn = file_io.enum_filename(self.outputfn, 2)
            sys.stdout.write("Writing to '%s'...\n" % self.outputfn)
            if self.opt.output_file_format == "csv":
                csv_format = {'dialect': 'excel', 'lineterminator': '\n'}
                if self.opt.csv_delimiter == 'tab':
                    csv_format['delimiter'] = '\t'
                f = unicode_csv.writer(file(self.outputfn, "w"),
                                       **self.opt.csv_format)
            elif self.opt.output_file_format == 'excel':
                import xls

                f = xls.writer(self.outputfn)
            fwrite("Table 1. Profile-centric data")
            fwrite("Source image:", self.src_img)
            fwrite("Comment:", self.comment)
            fwrite("Pixel width:", stringconv.tostr(float(self.pixelwidth), 2),
                   self.metric_unit)
            fwrite("Postsynaptic structure:", self.postsynProfile)
            fwrite("Presynaptic structure:", self.presynProfile)
            fwrite("Length of postsynaptic element:", m(self.posel.length()))
            fwrite("Length of presynaptic element:", m(self.posel.length()))
            fwrite("Number of PSDs:", len(self.psdli))
            fwrite("Total postsynaptic membrane length incl perforations:",
                   m(self.total_posm.length()))
            fwrite("Total postsynaptic membrane length excl perforations:",
                   sum([m(psd.posm.length()) for psd in self.psdli]))
            fwrite("Total PSD area:", sum([m2(psd.psdposm.area())
                                           for psd in self.psdli]))
            fwrite("Number of particles (total):", len(self.pli))
            fwrite("Number of particles in PSD:",
                   len([p for p in self.pli if p.isWithinPSD]))
            fwrite("Number of particles within %s metric units of PSD:"
                   % self.opt.spatial_resolution,
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
            f.writerows([[n + 1,
                          str(p),
                          stringconv.tostr(m(p.distToPosel), 2),
                          stringconv.tostr(m(p.distToPrsel), 2),
                          p.lateralLocation,
                          m(p.lateralDistPSD),
                          p.normLateralDistPSD,
                          stringconv.yes_or_no(p.isWithinPSD),
                          stringconv.yes_or_no(p.isAssociatedWithPSD)]
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
            f.writerows([[n + 1,
                          str(p),
                          stringconv.tostr(m(p.distToPosel), 2),
                          stringconv.tostr(m(p.distToPrsel), 2),
                          p.lateralLocation,
                          m(p.lateralDistPSD),
                          p.normLateralDistPSD,
                          stringconv.yes_or_no(p.isWithinPSD),
                          stringconv.yes_or_no(p.isAssociatedWithPSD)]
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
            f.writerows([[n + 1,
                          str(p),
                          stringconv.tostr(m(p.distToPosel), 2),
                          stringconv.tostr(m(p.distToPrsel), 2),
                          p.lateralLocation,
                          m(p.lateralDistPSD),
                          p.normLateralDistPSD,
                          stringconv.yes_or_no(p.isWithinPSD),
                          stringconv.yes_or_no(p.isAssociatedWithPSD)]
                         for n, p in enumerate(self.gridli)])
            f.close()
        except IOError:
            raise ProfileError(self, "Unable to write to file '%s'"
                               % self.outputfn)
        sys.stdout.write("Done.\n")
        return 1
# end of class ProfileData


class OptionData:
    def __init__(self):
        self.input_file_list = []
        self.spatial_resolution = 25
        self.shell_width = 200  # Skip points farther than this from the
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
        self.stop_requested = False

    def reset(self):
        """ Resets all options to default, and removes those that are not
            set in __init__().
        """
        self.__dict__ = {}
        self.__init__()
# end of class OptionData


class ProfileError(exceptions.Exception):
    def __init__(self, profile, msg):
        self.args = (profile, msg + ".")


def profile_warning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True


def profile_message(msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)
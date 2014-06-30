import sys
import geometry
from err_warn import ProfileError, profile_warning, profile_message


class PointList(list):
    def __init__(self, pointli, ptype, profile):
        try:
            self.extend([Point(p.x, p.y, ptype, profile) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError('not a list of Point elements')


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

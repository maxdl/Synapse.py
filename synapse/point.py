import sys
from . import geometry
from .err_warn import ProfileError, profile_warning, profile_message


def lazy_property(fn):
    """Decorator that makes a property lazily evaluated.
       From https://stevenloria.com/lazy-properties/.
    """
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


class PointList(list):
    def __init__(self, pointli, ptype, profile):
        super().__init__()
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
        self.discard = False
        self.ptype = ptype
        self.cluster = None
        self.is_within_psd = None
        self.lateral_dist_psd = None
        self.norm_lateral_dist_psd = None
        self.nearest_psd = None
        self.is_associated_with_psd = None
        self.associated_psd = None
        self.nearest_neighbour_dist = None
        self.nearest_neighbour_point = geometry.Point()
        self.nearest_lateral_neighbour_dist = None
        self.nearest_lateral_neighbour_point = geometry.Point()
        self.nearest_neighbour = geometry.Point()

    def determine_stuff(self):
        def mark_to_discard(msg):
            if self.ptype == 'particle':  # don't warn if random points
                profile_message("Discarding particle at %s: %s" % (self, msg))
            self.discard = True
            self.profile.n_discarded[self.ptype] += 1
            return

        if self.dist_to_posel is None:
            mark_to_discard("Could not project on postsynaptic element")
            return
        if self.dist_to_prsel is None:
            mark_to_discard("Could not project on presynaptic element")
            return
        if not self.is_within_shell:
            mark_to_discard("Located outside the shell")
            return
        if self.is_within_hole:
            mark_to_discard("Located within a profile hole")
            return
        __ = self.lateral_location
        __ = self.strict_lateral_location
        __ = self.axodendritic_location
        __ = self.is_within_postsynaptic_membrane_shell
        __ = self.is_postsynaptic_membrane_associated
        __ = self.is_presynaptic_membrane_associated
        self.get_lateral_dist_psd()
        self.get_psd_association()

    @lazy_property
    def dist_to_posel(self):
        """Return perpendicular distance to postsynaptic element membrane"""
        return self.perpend_dist(self.profile.posel, posloc=self.profile.posloc)

    @lazy_property
    def dist_to_prsel(self):
        """Return perpendicular distance to postsynaptic element membrane"""
        if len(self.profile.prsel) > 0:
            return self.perpend_dist(self.profile.prsel, posloc=self.profile.posloc)
        else:
            return None

    @lazy_property
    def is_within_hole(self):
        """  Determine whether self is inside a profile hole
        """
        for h in self.profile.holeli:
            if self.is_within_polygon(h):
                return True
        return False

    @lazy_property
    def is_within_shell(self):
        """Determine whether self is within shell"""
        return (self.dist_to_posel is not None and
                abs(self.dist_to_posel) <= geometry.to_pixel_units(self.opt.shell_width,
                                                                   self.profile.pixelwidth))

    @lazy_property
    def is_within_postsynaptic_membrane_shell(self):
        return (self.dist_to_posel is not None
                and abs(self.dist_to_posel) <=
                geometry.to_pixel_units(self.opt.shell_width,
                                        self.profile.pixelwidth))

    @lazy_property
    def lateral_location(self):
        """ Determines the lateral location, defined as the location
            along the mediolateral synaptic axis of the projection of
            the point on the postsynaptic membrane.
            
            Assumes that PSDs and the postsynaptic membrane have been determined.
        """
        # A point is classified as synaptic if its projection on the
        # postsynaptic membrane is within the spatial resolution limit
        # of the postsynaptic membrane
        for psd in self.profile.psdli:
            if (self.lateral_dist_syn(self.profile.posel, psd.posm) -
                    geometry.to_pixel_units(self.opt.spatial_resolution,
                                            self.profile.pixelwidth) <= psd.posm.length() / 2):
                return "synaptic"
        # If not synaptic but still within the extent of the synaptic
        # membrane
        if (self.lateral_dist_syn(self.profile.posel, self.profile.total_posm) <
                self.profile.total_posm.length() / 2):
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
            if (self.lateral_dist_syn(self.profile.posel, psd.posm) -
                    geometry.to_pixel_units(self.opt.spatial_resolution,
                                            self.profile.pixelwidth) <= psd.posm.length()):
                return "perisynaptic"
        # If none of the above
        return "extrasynaptic"

    @lazy_property
    def strict_lateral_location(self):
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
            if self.lateral_dist_syn(self.profile.posel, psd.posm) <= psd.posm.length() / 2.0:
                return "synaptic"
        # If not synaptic but still within the extent of the synaptic
        # membrane
        if (self.lateral_dist_syn(self.profile.posel, self.profile.total_posm) <
                self.profile.total_posm.length() / 2.0):
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
            if self.lateral_dist_syn(self.profile.posel, psd.posm) <= psd.posm.length():
                return "perisynaptic"
        # If none of the above
        return "extrasynaptic"

    @lazy_property
    def axodendritic_location(self):
        """ Determines the particle's location along the axodendritic axis.
        """
        if self.dist_to_posel >= 0:
            return "postsynaptic"
        elif self.dist_to_prsel is None:   # if there is no presynaptic membrane,
            return "not postsynaptic"      # that's about all we can say
        elif self.dist_to_prsel <= 0:
            return "presynaptic"
        else:
            return "neither pre- or postsynaptic"

    def get_lateral_dist_psd(self):
        mindist = sys.maxsize
        nearest_psd = None
        for psd in self.profile.psdli:
            d = self.lateral_dist_syn(self.profile.posel, psd.posm) - psd.posm.length() / 2
            if d < mindist:
                mindist = d
                nearest_psd = psd
        if not nearest_psd:
            raise ProfileError(self.profile,
                               "could not determine lateral distance to a PSD of particle at %s"
                               % self)
        mindist = self.lateral_dist_syn(self.profile.posel, nearest_psd.posm)
        normdist = mindist / (nearest_psd.posm.length() / 2)
        self.lateral_dist_psd = mindist
        self.norm_lateral_dist_psd = normdist
        self.nearest_psd = nearest_psd

    def get_psd_association(self):
        if self.is_within_psd is not None:
            return
        is_within_psd = False
        is_associated_with_psd = False
        associated_psd = None
        mindist = sys.maxsize
        for psd in self.profile.psdli:
            if self.is_within_polygon(psd.psdposm):
                is_within_psd = True
                is_associated_with_psd = True
                associated_psd = psd
                break
            dist = self.perpend_dist_closed_path(psd.psdposm, dont_care_if_on_or_off_seg=True)
            if dist <= geometry.to_pixel_units(self.opt.spatial_resolution,
                                               self.profile.pixelwidth):
                is_associated_with_psd = True
                if dist < mindist:
                    associated_psd = psd
                    mindist = dist
        self.is_within_psd = is_within_psd
        self.is_associated_with_psd = is_associated_with_psd
        self.associated_psd = associated_psd

    @lazy_property
    def is_postsynaptic_membrane_associated(self):
        if (self.dist_to_posel is not None and
                abs(self.dist_to_posel) <= geometry.to_pixel_units(self.opt.spatial_resolution,
                                                                   self.profile.pixelwidth)):
            return True
        else:
            return False

    @lazy_property
    def is_presynaptic_membrane_associated(self):
        if (self.dist_to_prsel is not None and
                abs(self.dist_to_prsel) <= geometry.to_pixel_units(self.opt.spatial_resolution,
                                                                   self.profile.pixelwidth)):
            return True
        else:
            return False

    def get_nearest_neighbour(self, pointli):
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxsize)
        for p in pointli:
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
        if not mindist < float(sys.maxsize):
            return None
        else:
            nearest_neighbour_dist = mindist
            return nearest_neighbour_dist

    def get_nearest_lateral_neighbour(self, pointli):
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxsize)
        for p in pointli:
            if p is not self:
                d = self.lateral_dist_to_point(p, self.profile.posel)
                if d < mindist:
                    mindist = d
        if not mindist < float(sys.maxsize):
            return None
        else:
            nearest_lateral_neighbour_dist = mindist
            return nearest_lateral_neighbour_dist

    def lateral_dist_to_point(self, p2, sel):
        """ Determine lateral distance to a point p2 along sel.
            Overrides function in geometry.Point, which only works
            with a closed path.
        """
        path = geometry.SegmentedPath()
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


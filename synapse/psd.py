from . import geometry
from .err_warn import ProfileError


class PSD(geometry.SegmentedPath):
    def __init__(self, pointli, profile):
        super().__init__(pointli)
        self.profile = profile
        self.posm = geometry.SegmentedPath()
        self.prsm = geometry.SegmentedPath()
        self.psdposm = geometry.SegmentedPath()
        self.cleft = geometry.SegmentedPath()
        self.cleft_width = None

    def adjust_psd(self):
        """ Adjust PSD coordinates so that the PSD ends exactly on
            the postsynaptic element.
        """
        # Partition psd into paths defined by the intersections with posel,
        # beginning and ending with the projections of the end nodes of psd.
        # Each path will then be completely on one side of posel.
        pathli = [geometry.SegmentedPath()]
        pathli[0].append(self[0].project_on_path_or_endnode(self.profile.posel)[0])
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
                pathli.append(geometry.SegmentedPath())
                pathli[-1].append(p)
        pathli[-1].append(self[-1])
        pathli[-1].append(self[-1].project_on_path_or_endnode(self.profile.posel)[0])
        # Now, look for the longest path. This is assumed to be the intended
        # psd. (Perhaps area is more relevant. However, this is messier because
        # we need to determine the part of dm enclosed by path for each path.)
        maxlength = 0
        longest_path = geometry.SegmentedPath()
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
        posm = geometry.SegmentedPath()
        p1, node1 = self[0].project_on_path_or_endnode(self.profile.posel)
        p2, node2 = self[-1].project_on_path_or_endnode(self.profile.posel)
        if None not in (node1, node2):
            posm.extend(self.profile.posel[node1 + 1:node2 + 1])
            posm.insert(0, p1)
            posm.append(p2)
        if len(posm) == 0 or posm.length() == 0:
            raise ProfileError(self.profile, "Could not determine postsynaptic membrane")
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
        prsm = geometry.SegmentedPath()
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
            raise ProfileError(self.profile, "Could not determine presynaptic membrane")
        return prsm

    def get_psd_posm(self):
        """ Make a closed-path PSD by concatenating the PSD and the
            posm. Assume that PSD and posm share end nodes and that
            they are oriented in the same direction.
        """
        pol = geometry.SegmentedPath()
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
        pol = geometry.SegmentedPath()
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



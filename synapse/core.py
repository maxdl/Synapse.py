import sys
import random
import itertools
from . import geometry
from . import file_io
from . import point
from . import psd
from .err_warn import ProfileError, profile_warning


# Convenience functions

def dot_progress(line_length=80, char='.', reset=False):
    """Simple progress indicator on sys.stdout"""
    if not hasattr(dot_progress, 'counter'):
        dot_progress.counter = 0
    if reset:
        dot_progress.counter = 0
        sys.stdout.write('\n')
    dot_progress.counter += 1
    sys.stdout.write(char)
    if dot_progress.counter == line_length:
        dot_progress.counter = 0
        sys.stdout.write('\n')


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


#
# Classes
#

class ClusterData(list):
    def __init__(self, pointli=None):
        super().__init__()
        if pointli is None:
            pointli = []
        try:
            self.extend([point.Point(p.x, p.y) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError("not a point list")
        self.convex_hull = geometry.SegmentedPath()

    def lateral_dist_syn(self, c2, mem):
        """ Determine lateral distance to a cluster c2 along postsynaptic
            element membrane.
        """
        centroid = point.Point(self.convex_hull.centroid())
        centroid2 = point.Point(c2.convex_hull.centroid())
        return centroid.lateral_dist_to_point(centroid2, mem)


class ProfileData:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.src_img = None
        self.opt = opt
        self.presyn_profile = None
        self.postsyn_profile = None
        self.holeli = []
        self.psdli = []
        self.pli = []
        self.posel = geometry.SegmentedPath()
        self.prsel = geometry.SegmentedPath()
        self.randomli = []
        self.mcli = []
        self.clusterli = []
        self.pp_distli, self.pp_latdistli = [], []
        self.rp_distli, self.rp_latdistli = [], []
        self.n_discarded = {'particle': 0, 'random': 0}
        self.comment = ''
        self.pixelwidth = None
        self.metric_unit = ''
        self.posloc = geometry.Point()
        self.negloc = geometry.Point()
        self.total_synm = geometry.SegmentedPath()
        self.perimeter = None
        self.feret = None
        self.warnflag = False
        self.errflag = False
        self.spatial_resolution_in_pixels = None
        self.total_posm = geometry.SegmentedPath()
        self.postsynloc = geometry.Point()
        self.pli = []
        self.gridli = None
        self.randomli = None

    def process(self):
        """ Parse profile data from a file and determine distances
        """
        try:
            sys.stdout.write("Processing profile...\n")
            self.__parse()
            self.__check_paths()
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
            self.total_posm = self.__get_total_posm()
            self.postsynloc = self.psdli[0].psdloc
            self.__compute_stuff()
            if self.opt.determine_interpoint_dists:
                sys.stdout.write("Determining interparticle distances...\n")
                self.__determine_interdistlis()
            if self.opt.determine_clusters:
                sys.stdout.write("Determining clusters...\n")
                self.clusterli = self.__determine_clusters(self.pli)
            if self.opt.run_monte_carlo:
                sys.stdout.write("Running Monte Carlo simulations...\n")
                self.__run_monte_carlo()
            if self.opt.stop_requested:
                return
            sys.stdout.write("Done.\n")
        except ProfileError as err:
            sys.stdout.write("Error: %s\n" % err.msg)
            self.errflag = True

    def __compute_stuff(self):
        for p in self.pli:
            p.determine_stuff()
        self.pli = [p for p in self.pli if not p.discard]
        for p in self.randomli:
            p.determine_stuff()
        self.randomli = [p for p in self.randomli if not p.discard]
        for ptype in ('particle', 'random'):
            if ptype == 'random' and not self.opt.use_random:
                continue
            ptypestr = 'particles' if ptype == 'particle' else ptype + ' points'
            sys.stdout.write("  Number of %s discarded: %d\n"
                             % (ptypestr, self.n_discarded[ptype]))

    def __determine_interdistlis(self):
        if True not in [val for key, val in self.opt.interpoint_relations.items()
                        if 'simulated' not in key]:
            return
        if self.opt.interpoint_relations['particle - particle']:
            self.pp_distli, self.pp_latdistli = self.__get_same_interpoint_distances(self.pli)
        if self.opt.use_random and self.opt.interpoint_relations['random - particle']:
            self.rp_distli, self.rp_latdistli = self.__get_interpoint_distances2(self.randomli,
                                                                                 self.pli)

    def __get_same_interpoint_distances(self, pointli):
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
                            pointli[j], self.posel))
            elif self.opt.interpoint_dist_mode == 'nearest neighbour':
                if self.opt.interpoint_shortest_dist:
                    dli.append(pointli[i].get_nearest_neighbour(pointli))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(pointli[i].get_nearest_lateral_neighbour(
                        pointli))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def __get_interpoint_distances2(self, pointli, pointli2=None):
        if pointli2 is None:
            pointli2 = []
        dli = []
        latdli = []
        for i, p in enumerate(pointli):
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
                    dli.append(p.get_nearest_neighbour(pointli2))
                if self.opt.interpoint_lateral_dist:
                    latdli.append(p.get_nearest_lateral_neighbour(pointli2))
        dli = [d for d in dli if d is not None]
        latdli = [d for d in latdli if d is not None]
        return dli, latdli

    def __run_monte_carlo(self):
        
        def is_valid(p_candidate):
            d = p_candidate.perpend_dist(self.posel)
            if (d is None or abs(d) > border or p_candidate.is_within_hole()
                    or p_candidate in mcli[n]['pli']):
                return False
            if self.opt.monte_carlo_simulation_window == 'profile':
                return True
            if self.opt.monte_carlo_strict_location:
                latloc = p_candidate.get_strict_lateral_location()
            else:
                latloc = p_candidate.get_lateral_location()
            if (self.opt.monte_carlo_simulation_window == 'synapse' and
                    latloc in ('synaptic', 'within perforation')):
                return True
            if (self.opt.monte_carlo_simulation_window ==
                    'synapse - perforations' and latloc == 'synaptic'):
                return True
            if (self.opt.monte_carlo_simulation_window ==
                    'synapse + perisynapse' and latloc in
                    ('synaptic', 'within perforation', 'perisynaptic')):
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
            numpoints = len([p for p in pli
                             if getattr(p, locvar) in ("synaptic", "within perforation")])
        elif self.opt.monte_carlo_simulation_window == "synapse - perforations":
            numpoints = len([p for p in pli if getattr(p, locvar) == "synaptic"])
        elif self.opt.monte_carlo_simulation_window == "synapse + perisynapse":
            numpoints = len([p for p in pli
                             if getattr(p, locvar) in
                             ("synaptic", "within perforation", "perisynaptic")])
        else:
            return []
        allpaths = geometry.SegmentedPath(
            [p for path in (self.posel, self.prsel,
                            [q for _psd in self.psdli for q in _psd],
                            [r for h in self.holeli for r in h])
             for p in path])
        box = allpaths.bounding_box()
        border = geometry.to_pixel_units(min(self.opt.shell_width, self.opt.spatial_resolution))
        mcli = []
        for n in range(0, self.opt.monte_carlo_runs):
            if self.opt.stop_requested:
                return []
            dot_progress(n)
            mcli.append({'pli': [],
                         'simulated - simulated': {'dist': [], 'latdist': []},
                         'simulated - particle': {'dist': [], 'latdist': []},
                         'particle - simulated': {'dist': [], 'latdist': []},
                         'clusterli': []})
            p = point.Point()
            for __ in range(0, numpoints):
                while True:
                    x = random.randint(int(box[0].x - border), int(box[1].x + border) + 1)
                    y = random.randint(int(box[0].y - border), int(box[2].y + border) + 1)
                    p = point.Point(x, y, ptype='simulated', profile=self)
                    if p not in mcli[n]['pli'] and is_valid(p):
                        break
                # escape the while loop when a valid simulated point is found
                mcli[n]['pli'].append(p)
            for p in mcli[n]['pli']:
                p.determine_stuff()
        if self.opt.determine_interpoint_dists:
            for n in range(0, self.opt.monte_carlo_runs):
                if self.opt.stop_requested:
                    return []
                dot_progress(n)
                if self.opt.interpoint_relations['simulated - simulated']:
                        distlis = self.__get_same_interpoint_distances(mcli[n]['pli'])
                        mcli[n]['simulated - simulated']['dist'].append(distlis[0])
                        mcli[n]['simulated - simulated']['latdist'].append(distlis[1])
                if self.opt.interpoint_relations['simulated - particle']:
                    distlis = self.__get_interpoint_distances2(mcli[n]['pli'], pli)
                    mcli[n]['simulated - particle']['dist'].append(distlis[0])
                    mcli[n]['simulated - particle']['latdist'].append(distlis[1])
                if self.opt.interpoint_relations['particle - simulated']:
                    distlis = self.__get_interpoint_distances2(pli, mcli[n]['pli'])
                    mcli[n]['particle - simulated']['dist'].append(distlis[0])
                    mcli[n]['particle - simulated']['latdist'].append(distlis[1])
        if self.opt.determine_clusters:
            sys.stdout.write("\n-> Determining simulated clusters:\n   ")
            for n, li in enumerate(mcli):
                if self.opt.stop_requested:
                    return []
                dot_progress(n)
                mcli[n]['clusterli'] = self.__determine_clusters(li['pli'])
                self.__process_clusters(mcli[n]['clusterli'])
        self.mcli = mcli
        sys.stdout.write("\n")
        return mcli

    def __process_clusters(self, clusterli):
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.convex_hull = geometry.convex_hull(c)
            c.dist_to_posel = c.convex_hull.centroid().perpend_dist(self.posel)
        for c in clusterli:
            if self.opt.stop_requested:
                return
            c.nearest_cluster = ClusterData()
            if len(clusterli) == 1:
                c.dist_to_nearest_cluster = -1
                return
            c.dist_to_nearest_cluster = sys.maxsize
            for c2 in clusterli:
                if c2 != c:
                    d = c.lateral_dist_to_cluster(c2, self.posel)
                    if d < c.dist_to_nearest_cluster:
                        c.dist_to_nearest_cluster = d
                        c.nearest_cluster = c2

    def __determine_clusters(self, pointli):
        """ Partition pointli into clusters; each cluster contains all points
            that are less than opt.within_cluster_dist from at least one
            other point in the cluster
        """
        if self.opt.within_cluster_dist < 0:
            return
        clusterli = []
        for p1 in pointli:
            if self.opt.stop_requested:
                return []
            if p1.cluster:
                continue
            for p2 in pointli:
                if p1 != p2 and p1.dist(p2) <= geometry.to_pixel_units(self.opt.within_cluster_dist,
                                                                       self.pixelwidth):
                    if p2.cluster is not None:
                        p1.cluster = p2.cluster
                        clusterli[p1.cluster].append(p1)
                        break
            else:
                p1.cluster = len(clusterli)
                clusterli.append(ClusterData([p1]))
        self.__process_clusters(clusterli)
        return clusterli

    def __parse(self):
        """ Parse profile data from input file 
        """
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = file_io.read_file(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        self.psdli = []
        self.holeli = []
        while li:
            s = li.pop(0).replace('\n', '').strip()
            if s.split(' ')[0].upper() == 'IMAGE':
                self.src_img = s.split(' ')[1]
            elif s.split(' ')[0].upper() in ('PROFILE_ID', 'SYNAPSE_ID'):
                try:
                    self.id = int(s.split(' ')[1])
                except (IndexError, ValueError):
                    profile_warning(self, "Profile id not defined or invalid")
            elif s.split(' ')[0].upper() == 'COMMENT':
                try:
                    self.comment = s.split(' ', 1)[1]
                except IndexError:
                    self.comment = ''
            elif s.split(' ')[0].upper() == 'PIXELWIDTH':
                try:
                    self.pixelwidth = float(s.split(' ')[1])
                    self.metric_unit = s.split(' ')[2]
                except (IndexError, ValueError):
                    raise ProfileError(self, "PIXELWIDTH is not a valid number")
            elif s.split(" ")[0].upper() == "POSTSYNAPTIC_PROFILE":
                try:
                    self.postsyn_profile = s.split(" ", 1)[1]
                except IndexError:
                    pass
            elif s.split(" ")[0].upper() == "PRESYNAPTIC_PROFILE":
                try:
                    self.presyn_profile = s.split(" ", 1)[1]
                except IndexError:
                    pass
            elif s.upper() == "POSTSYNAPTIC_ELEMENT":
                self.posel = geometry.SegmentedPath(self.__get_coords(li, "postsynaptic element"))
            elif s.upper() == "PRESYNAPTIC_ELEMENT":
                self.prsel = geometry.SegmentedPath(self.__get_coords(li, "presynaptic element"))
            elif s.upper() in ("POSTSYNAPTIC_DENSITY", "PSD_OR_ACTIVE_ZONE", "PSD"):
                self.psdli.append(psd.PSD(self.__get_coords(li, "PSD"), self))
            elif s.upper() in ("PROFILE_HOLE", "HOLE"):
                self.holeli.append(geometry.SegmentedPath(self.__get_coords(li, "hole")))
            elif s.upper() in ("POINTS", "PARTICLES"):
                self.pli = point.PointList(self.__get_coords(li, "particle"), "particle", self)
            elif s.upper() == "RANDOM_POINTS":
                self.randomli = point.PointList(self.__get_coords(li, "random"), "random", self)
            elif s.upper() == 'GRID':
                # Retrieve coordinates to dummy variable as they will not be used
                __ = point.PointList(self.__get_coords(li, 'grid'), 'grid', self)
                profile_warning(self, "Grid found; however, as grids are no longer supported " 
                                      "it will be discarded")
            elif s[0] != "#":  # unless specifically commented out
                profile_warning(self, "Unrecognized string '" + s +
                                "' in input file")
        # Now, let's see if everything was found
        self.__check_parsed_data()

    def __check_parsed_data(self):
        """See if the profile data was parsed correctly, and print info
        on the parsed data to stdout.
        """
        self.__check_var_default('src_img', "Source image", "N/A")
        self.__check_var_default('id', "Profile ID", "N/A")
        self.__check_var_default('comment', "Comment", "")
        self.__check_var_val('metric_unit', "Metric unit", 'metric_unit')
        self.__check_required_var('pixelwidth', "Pixel width", self.metric_unit)
        self.__check_var_default('postsyn_profile', "Postsynaptic profile","N/A")
        self.__check_var_default('presyn_profile', "Presynaptic profile", "N/A")
        self.__check_list_var('posel', 'Postsynaptic element', 'nodes', 2)
        self.__check_list_var('prsel', 'Presynaptic element', 'nodes', 2)
        self.__check_table_var('psdli', "Postsynaptic density", "Postsynaptic densities", 1, 2)
        self.__check_list_var('pli', 'Particles', '', 0)
        self.__check_table_var('holeli', "Hole", "Holes", 0, 2)
        self.__check_var_exists('randomli', "Random points", 'use_random')

    def __check_required_var(self, var_to_check, var_str, post_str):
        """Confirm that self has a required variable; else, raise
        ProfileError.
        """
        if not self.__dict__[var_to_check]:
            raise ProfileError(self, "%s not found in input file" % var_str)
        else:
            sys.stdout.write("  %s: %s %s\n" % (var_str, self.__dict__[var_to_check], post_str))

    @staticmethod
    def __check_list_len(var, min_len):
        """Return True if var is a list and has at least min_len
        elements, else False.
        """
        return isinstance(var, list) and len(var) >= min_len

    def __check_list_var(self, var_to_check, var_str, post_str, min_len):
        """Confirms that self has a var_to_check that is a list and
        has at least min_len elements; if var_to_check does not exist
        and min_len <= 0, assigns an empty list to var_to_check. Else,
        raise a ProfileError.
        """
        if not self.__dict__[var_to_check]:
            if min_len > 0:
                raise ProfileError(self, "%s not found in input file" % var_str)
            else:
                self.__dict__[var_to_check] = []
        elif not self.__check_list_len(self.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = " " + post_str
        sys.stdout.write("  %s%s: %d\n" % (var_str, post_str, len(self.__dict__[var_to_check])))

    def __check_table_var(self, var_to_check, var_str_singular,
                          var_str_plural, min_len_1, min_len_2):
        """Confirms that var_to_check exists, is a list and has at
        least min_len_1 elements, and that each of these has at least
        min_len_2 subelements; if var_to_check does not exist and
        min_len_1 <= 0, assigns an empty list to var_to_check. Else,
        raise ProfileError.
        """
        if not self.__dict__[var_to_check]:
            if min_len_1 > 0:
                raise ProfileError(self, "%s not found in input file" % var_str_plural)
            else:
                self.__dict__[var_to_check] = []
        elif not self.__check_list_len(self.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file" % var_str_plural.lower())
        else:
            for element in self.__dict__[var_to_check]:
                if not self.__check_list_len(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates" % var_str_singular)
        sys.stdout.write("  %s: %d\n" % (var_str_plural, len(self.__dict__[var_to_check])))

    def __check_var_default(self, var_to_check, var_str, default=""):
        """Checks if var_to_check exists; if not, assign the default
        value to var_to_check. Never raises a ProfileError.
        """
        if not self.__dict__[var_to_check]:
            self.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str, self.__dict__[var_to_check]))

    def __check_var_exists(self, var_to_check, var_str, optflag):
        """Checks for consistency between profiles with respect to the
        existence of var_to_check (i.e., var_to_check must be present
        either in all profiles or in none).

        If optflag is not set (i.e., this is the first profile), then
        set optflag to True or False depending on the existence of
        var_to_check. If optflag is already set (for consequent
        profiles), var_to_check must (if optflag is True) or must not
        (if optflag is False) exist. If not so, raise ProfileError.
        """
        if not hasattr(self.opt, optflag):
            if self.__dict__[var_to_check]:
                self.opt.__dict__[optflag] = True
            else:
                self.opt.__dict__[optflag] = False
        if self.opt.__dict__[optflag]:
            if self.__dict__[var_to_check]:
                sys.stdout.write("  %s: yes\n" % var_str)
            else:
                raise ProfileError(self, "%s expected but not found in input file" % var_str)
        elif self.__dict__[var_to_check]:
            raise ProfileError(self, "%s found but not expected" % var_str)
        else:
            sys.stdout.write("  %s: no\n" % var_str)

    def __check_var_val(self, var_to_check, var_str, optvar):
        """Checks for consistency between profiles with respect to the
        value of var_to_check (i.e., var_to_check must be present and
        have equal value in all profiles).

        If optvar is not set (i.e., this is the first profile), then
        set optflag to the value of var_to_check. If optvar is already
        set (for consequent profiles), the value of var_to_check must
        be equal to that of optvar. If not so, raise ProfileError.
        """
        if not self.__dict__[var_to_check]:
            raise ProfileError(self, "%s not found in input file" % var_str)
        if not hasattr(self.opt, optvar):
            self.opt.__dict__[optvar] = self.__dict__[var_to_check]
        elif self.__dict__[var_to_check] == self.opt.__dict__[optvar]:
            pass  # really no point in pointing out that it's ok
            # sys.stdout.write("  %s: %s\n"
            #                  % (var_str, parent.__dict__[var_to_check]))
        else:
            raise ProfileError(self, "%s value '%s'  differs from the value "
                                     "specified ('%s') in the first input file"
                               % (var_str, self.__dict__[var_to_check],
                                  self.opt.__dict__[optvar]))

    def __check_paths(self):
        """Check if profile border and holes intersect with themselves."""

        def check_path(_path, s):
            for p in range(0, len(_path) - 3):
                for q in range(0, len(_path) - 1):
                    if p not in (q, q + 1) and p + 1 not in (q, q + 1):
                        if geometry.segment_intersection(_path[p],
                                                         _path[p + 1],
                                                         _path[q],
                                                         _path[q + 1]):
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
        for n, h in enumerate(self.holeli):
            if not h.is_simple_polygon():
                raise ProfileError(self, "Profile hole %d is not a simple polygon" % (n + 1))
            for n2, h2 in enumerate(self.holeli[n + 1:]):
                if h.overlaps_polygon(h2):
                    raise ProfileError(self, "Profile hole %d overlaps with hole %d "
                                       % (n + 1, n + n2 + 2))
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

    def __get_total_posm(self):
        """ Construct a path comprising the postsynaptic membrane of all 
            PSDs and perforations.

            Assume well-behaved PSDs.
        """

        def dist_to_left_endnode(_p):
            path = geometry.SegmentedPath()
            project, seg_project = _p.project_on_path_or_endnode(self.posel)
            path.extend([self.posel[0], project])
            for n in range(1, seg_project):
                path.insert(len(path) - 1, self.posel[n])
            return path.length()

        left_mindist = right_mindist = self.posel.length()
        left_p = point.Point()
        right_p = point.Point()
        for p in itertools.chain(*[[_psd[0], _psd[-1]] for _psd in self.psdli]):
            left_d = dist_to_left_endnode(p)
            right_d = self.posel.length() - left_d
            if left_d <= left_mindist:
                left_mindist = left_d
                left_p = p
            elif right_d <= right_mindist:
                right_mindist = right_d
                right_p = p
        pseudo_psd = psd.PSD(geometry.SegmentedPath([left_p, right_p]), self)
        return pseudo_psd.get_posm()

    def __get_coords(self, strli, coord_type=""):
        """Pop point coordinates from list strli.

        When an element of strli is not a valid point, a warning is
        issued.
        """
        pointli = []
        s = strli.pop(0).replace('\n', '').replace(' ', '').strip()
        while s != 'END':
            try:
                p = geometry.Point(float(s.split(',')[0]), float(s.split(',')[1]))
                if pointli and (p == pointli[-1] or (coord_type == 'particle' and p in pointli)):
                    sys.stdout.write("Duplicate %s coordinates %s: skipping "
                                     "2nd instance\n" % (coord_type, p))
                else:
                    pointli.append(point.Point(p.x, p.y, ptype=coord_type))
            except ValueError:
                if s[0] != '#':
                    profile_warning(self, "'%s' not valid %s coordinates" % (s, coord_type))
                else:
                    pass
            s = strli.pop(0).replace('\n', '').strip()
        # For some reason, sometimes the endnodes have the same coordinates;
        # in that case, delete the last endnode to avoid division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]):
            del pointli[-1]
        return pointli
# end of class ProfileData

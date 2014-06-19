# -*- coding: utf-8 -*-

import sys
import os.path
import time
from core import *
import geometry
from file_io import *
import version
import stringconv as sc
    
#
# Functions
#


def save_output(profileli, opt):
    """ Save a summary of results of evaluated profiles
    """
    def m(x, pixelwidth):
        return geometry.to_metric_units(x, pixelwidth)

    def m2(x, pixelwidth):
        # for area units...
        return geometry.to_metric_units(x, pixelwidth**2)  # for area units...

    def na(x):
        if x in (None, -1):
            return "N/A"
        else:
            return x

    def write_session_summary():
        with file_io.FileWriter("session.summary", opt) as f:
            f.writerow(["%s version:" % version.title,
                       "%s (Last modified %s %s, %s)"
                       % ((version.version,) + version.date)])
            f.writerow(["Number of evaluated profiles:", len(eval_proli)])
            if err_fli:
                f.writerow(["Number of non-evaluated profiles:",
                            len(err_fli)])
            f.writerow(["Metric unit:", eval_proli[0].metric_unit])
            f.writerow(["Spatial resolution:", opt.spatial_resolution,
                        eval_proli[0].metric_unit])
            f.writerow(["Shell width:", opt.shell_width,
                        eval_proli[0].metric_unit])
            f.writerow(["Interpoint distances calculated:",
                        sc.yes_or_no(opt.determine_interpoint_dists)])
            if opt.determine_interpoint_dists:
                f.writerow(["Interpoint distance mode:",
                            opt.interpoint_dist_mode])
                f.writerow(["Shortest interpoint distances:",
                            sc.yes_or_no(opt.interpoint_shortest_dist)])
                f.writerow(["Lateral interpoint distances:",
                            sc.yes_or_no(opt.interpoint_lateral_dist)])
            f.writerow(["Monte Carlo simulations performed:",
                        sc.yes_or_no(opt.run_monte_carlo)])
            if opt.run_monte_carlo:
                f.writerow(["Number of Monte Carlo runs:",
                            opt.monte_carlo_runs])
                f.writerow(["Monte Carlo simulation window:",
                            opt.monte_carlo_simulation_window])
                f.writerow(["Strict localization in simulation window:",
                            sc.yes_or_no(opt.monte_carlo_strict_location)])
            f.writerow(["Clusters determined:",
                        sc.yes_or_no(opt.determine_clusters)])
            if opt.determine_clusters:
                f.writerow(["Within-cluster distance:",
                            opt.within_cluster_dist,
                            eval_proli[0].metric_unit])
            if clean_fli:
                f.writerow(["Input files processed cleanly:"])
                f.writerows([[fn] for fn in clean_fli])
            if nop_fli:
                f.writerow(["Input files processed but which generated no "
                            "point distances:"])
                f.writerows([[fn] for fn in nop_fli])
            if warn_fli:
                f.writerow(["Input files processed but which generated "
                            "warnings (see log for details):"])
                f.writerows([[fn] for fn in warn_fli])
            if err_fli:
                f.writerow(["Input files not processed or not included in "
                            "summary (see log for details):"])
                f.writerows([[fn] for fn in err_fli])

    def write_profile_summary():
        with file_io.FileWriter("profile.summary", opt) as f:
            f.writerow(["Postsynaptic element length",
                        "Presynaptic element length",
                        "Number of PSDs:",             
                        "Total postsynaptic membrane length incl perforations:", 
                        "Total postsynaptic membrane length excl perforations:", 
                        "Total PSD area:",
                        "Particles (total)",
                        "Particles in PSD",
                        "Particles within %s %s of PSD" 
                        % (opt.spatial_resolution, eval_proli[0].metric_unit),
                        "Shell particles strictly synaptic and postsynaptic",
                        "Shell particles strictly synaptic and postsynaptic "
                        "or associated with postsynaptic membrane",
                        "Synaptic particles associated w/ postsynaptic "
                        "membrane",
                        "Synaptic particles associated w/ presynaptic membrane",
                        "Perisynaptic particles associated w/ postsynaptic "
                        "membrane",
                        "Perisynaptic particles associated w/ presynaptic "
                        "membrane",
                        "Within-perforation particles associated w/ "
                        "postsynaptic membrane",
                        "Within-perforation particles associated w/ "
                        "presynaptic membrane",
                        "Presynaptic profile",
                        "Postsynaptic profile",
                        "Input file"])
            f.writerows([[m(pro.posel.length(), pro.pixelwidth),
                          m(pro.prsel.length(), pro.pixelwidth),
                          len(pro.psdli),
                          m(pro.totalPosm.length(), pro.pixelwidth),
                          sum([m(psd.posm.length(), pro.pixelwidth)
                               for psd in pro.psdli]),
                          sum([m2(psd.psdposm.area(), pro.pixelwidth)
                               for psd in pro.psdli]),
                          len(pro.pli),
                          len([p for p in pro.pli if p.isWithinPSD]),
                          len([p for p in pro.pli if p.isAssociatedWithPSD]),
                          len([p for p in pro.pli
                               if p.strictLateralLocation == "synaptic" and
                               p.axodendriticLocation == "postsynaptic" and
                               p.isWithinPostsynapticMembraneShell]),
                          len([p for p in pro.pli
                               if p.strictLateralLocation == "synaptic" and
                               (p.axodendriticLocation == "postsynaptic" and
                                p.isWithinPostsynapticMembraneShell) or
                               p.isPostsynapticMembraneAssociated]),
                          len([p for p in pro.pli
                               if p.lateralLocation == "synaptic" and
                               p.isPostsynapticMembraneAssociated]),
                          len([p for p in pro.pli
                               if p.lateralLocation == "synaptic" and
                               p.isPresynapticMembraneAssociated]),
                          len([p for p in pro.pli
                               if p.lateralLocation == "perisynaptic" and
                               p.isPostsynapticMembraneAssociated]),
                          len([p for p in pro.pli
                               if p.lateralLocation == "perisynaptic" and
                               p.isPresynapticMembraneAssociated]),
                          len([p for p in pro.pli
                               if p.lateralLocation == "within perforation"
                               and p.isPostsynapticMembraneAssociated]),
                          len([p for p in pro.pli
                               if p.lateralLocation == "within perforation"
                               and p.isPresynapticMembraneAssociated]),
                          pro.presynProfile,
                          pro.postsynProfile,
                          os.path.basename(pro.inputfn)] for pro in eval_proli])

    def write_point_summary(ptype):
        if ptype == "particle":
            pli = "pli"
            pstr = "particle"
        elif ptype == "random":
            if not opt.useRandom:
                return
            else:
                pli = "randomli"
                pstr = "point"
        elif ptype == "grid":
            if not opt.useGrid:
                return
            else:
                pli = "gridli"
                pstr = "point"
        else:
            return
        with file_io.FileWriter("%s.summary" % ptype, opt) as f:
            f.writerow(["%s number (as appearing in input file)"
                        % pstr.capitalize(),
                        "Axodendritic location",                        
                        "Distance to postsynaptic element membrane",
                        "Distance to presynaptic element membrane",                        
                        "Lateral location",
                        "Strict lateral location",
                        "Lateral distance to nearest PSD center", 
                        "Normalized lateral distance to nearest PSD center",
                        "Within PSD",
                        "Within %s %s of PSD"
                        % (opt.spatial_resolution, eval_proli[0].metric_unit),
                        "Total postsynaptic membrane length incl perforations",
                        "Total postsynaptic membrane length excl perforations",
                        "Length of laterally closest PSD",
                        "Presynaptic profile",
                        "Postsynaptic profile",
                        "Input file",
                        "Comment"])
            f.writerows([[n+1,
                          p.axodendriticLocation,                           
                          m(p.distToPosel, pro.pixelwidth),
                          m(p.distToPrsel, pro.pixelwidth),
                          p.lateralLocation,
                          p.strictLateralLocation,
                          m(p.lateralDistPSD, pro.pixelwidth),
                          p.normLateralDistPSD,                          
                          sc.yes_or_no(p.isWithinPSD),
                          sc.yes_or_no(p.isAssociatedWithPSD),
                          m(pro.totalPosm.length(), pro.pixelwidth),
                          m(sum([psd.posm.length() for psd in pro.psdli]),
                            pro.pixelwidth),
                          m(p.closestPSD.posm.length(), pro.pixelwidth),
                          pro.presynProfile,
                          pro.postsynProfile,
                          os.path.basename(pro.inputfn),
                          pro.comment] for pro in eval_proli for n, p in
                         enumerate(pro.__dict__[pli])])

    def write_cluster_summary():
        if not opt.determine_clusters:
            return
        with FileWriter("cluster.summary", opt) as f:
            f.writerow(["Cluster number",
                       "Number of particles in cluster",
                       "Distance to postsynaptic membrane of centroid",
                       "Distance to nearest cluster along postsynaptic element "
                       "membrane",
                       "Profile ID",
                       "Input file",
                       "Comment"])
            f.writerows([[n + 1,
                        len(c),
                        m(c.distToPath, pro.pixelwidth),
                        m(na(c.distToNearestCluster), pro.pixelwidth),
                        pro.ID,
                        os.path.basename(pro.inputfn),
                        pro.comment]for pro in eval_proli for n, c in
                        enumerate(pro.clusterli)])

    def write_interpoint_summaries():

        def _m(x):
            return m(x, pro.pixelwidth)

        if not opt.determine_interpoint_dists:
            return
        ip_rels = dict([(key, val)
                        for key, val in opt.interpoint_relations.items()
                        if val and "simulated" not in key])
        if not opt.useRandom:
            for key, val in opt.interpoint_relations.items():
                if "random" in key and val:
                    del ip_rels[key]
        if (len(ip_rels) == 0 or not
           (opt.interpoint_shortest_dist or opt.interpoint_lateral_dist)):
            return
        table = []
        if opt.interpoint_dist_mode == 'all':
            s = "all distances"
        else:
            s = "nearest neighbour distances"
        table.append(["Mode: " + s])
        headerli = ip_rels.keys()
        prefixli = []
        for key, val in ip_rels.items():
            prefix = key[0] + key[key.index("- ") + 2] + "_"
            #if prefix[0] == prefix[1]:
            #    prefix = prefix.replace(prefix[0], "", 1)
            prefixli.append(prefix)
        if opt.interpoint_shortest_dist and opt.interpoint_lateral_dist:
            headerli.extend(headerli)
            prefixli.extend(map(lambda t: t + "lat", prefixli))
        topheaderli = []
        if opt.interpoint_shortest_dist:
            topheaderli.append("Shortest distances")
            if opt.interpoint_lateral_dist:
                topheaderli.extend([""] * (len(ip_rels) - 1))
        if opt.interpoint_lateral_dist:
            topheaderli.append("Lateral distances along postsynaptic element "
                               "membrane")
        table.extend([topheaderli, headerli])
        cols = [[] for c in prefixli]
        for pro in eval_proli:
            for n, li in enumerate([pro.__dict__[prefix + "distli"]
                                    for prefix in prefixli]):
                cols[n].extend(map(_m, li))
        # transpose cols and append to table
        table.extend(map(lambda *col: [e if e is not None else "" for e in col],
                         *cols))
        with file_io.FileWriter("interpoint.distances", opt) as f:
            f.writerows(table)

    def write_mc_dist_to_psd(ptype):

        def m_li(*_li):
            return [m(x, pro.pixelwidth) for x in _li]

        if not opt.run_monte_carlo:
            return
        table = []
        if ptype == "metric":
            table.append(["Lateral distances in %s to the nearest PSD"
                          % eval_proli[0].metric_unit])
        elif ptype == "normalized":
            table.append(["Normalized lateral distances to the nearest PSD"])
        table.append(["Run %d" % (n + 1)
                     for n in range(0, opt.monte_carlo_runs)])
        for pro in eval_proli:
            if ptype == "metric":
                table.extend(map(m_li, *[[p.lateralDistPSD for p in li["pli"]]
                                         for li in pro.mcli]))
            elif ptype == "normalized":
                table.extend(map(None, *[[p.normLateralDistPSD
                                          for p in li["pli"]]
                                         for li in pro.mcli]))
        with FileWriter("simulated.PSD.%s.lateral.distances" % ptype, opt) as f:
            f.writerows(table)

    def write_mc_dist_to_posel():

        def m_li(*_li):
            return [m(x, pro.pixelwidth) for x in _li]

        if not opt.run_monte_carlo:
            return
        table = [["Run %d" % (n + 1)
                  for n in range(0, opt.monte_carlo_runs)]]
        for pro in eval_proli:
            table.extend(map(m_li, *[[p.distToPosel for p in li["pli"]]
                                     for li in pro.mcli]))
        with FileWriter("simulated.postsynaptic.element.membrane.distances",
                        opt) as f:
            f.writerows(table)

    def write_mc_ip_dists(dist_type):

        def m_li(*_li):
            return [m(x, pro.pixelwidth) for x in _li]

        if not opt.run_monte_carlo:
            return
        for ip_type in [key for key, val in opt.interpoint_relations.items()
                        if "simulated" in key and val]:
            if ((dist_type == "shortest" and not opt.interpoint_shortest_dist)
                or
               (dist_type == "lateral" and not opt.interpoint_lateral_dist)):
                return
            if dist_type == "lateral":
                short_dist_type = "lat"
            else:
                short_dist_type = ""
            table = [["Run %d" % (n + 1)
                      for n in range(0, opt.monte_carlo_runs)]]
            for pro in eval_proli:
                table.extend(map(m_li,
                                 *[p for li in pro.mcli
                                 for p in li[ip_type]
                                 ["%sdist" % short_dist_type]]))
            with file_io.FileWriter("%s.interpoint.%s.distance.summary"
                                    % (ip_type.replace(" ", ""),
                                       dist_type), opt) as f:
                f.writerows(table)

    def write_mc_cluster_summary():
        if not (opt.determine_clusters and opt.run_monte_carlo):
            return
        table = [["N particles in cluster", "Run",
                  "Distance to postsynaptic element membrane from centroid",
                  "Distance to nearest cluster",
                  "Profile ID",
                  "Input file",
                  "Comment"]]
        for pro in eval_proli:
            for n in range(0, opt.monte_carlo_runs):
                for c in pro.mcli[n]["clusterli"]:
                    table.append([len(c), n + 1,
                                 m(c.dist_to_posel, pro.pixelwidth),
                                 m(na(c.dist_to_nearest_cluster),
                                   pro.pixelwidth),
                                 pro.ID,
                                 os.path.basename(pro.inputfn),
                                 pro.comment])
        with file_io.FileWriter("simulated.cluster.summary", opt) as f:
            f.writerows(table)

    sys.stdout.write("\nSaving summaries...\n")
    opt.save_result = {'any_saved': False, 'any_err': False}
    eval_proli = [profile for profile in profileli if not profile.errflag]
    clean_fli = [profile.inputfn for profile in profileli
                 if not (profile.errflag or profile.warnflag)]
    warn_fli = [profile.inputfn for profile in profileli if profile.warnflag]
    err_fli = [profile.inputfn for profile in profileli if profile.errflag]
    nop_fli = [profile.inputfn for profile in profileli if not profile.pli]
    write_session_summary()
    write_profile_summary()
    write_point_summary("particle")
    write_point_summary("random")
    write_point_summary("grid")
    write_interpoint_summaries()
    write_cluster_summary()
    write_mc_dist_to_posel()
    write_mc_dist_to_psd("metric")
    write_mc_dist_to_psd("normalized")
    write_mc_ip_dists("shortest")
    write_mc_ip_dists("lateral")
    write_mc_cluster_summary()
    if opt.save_result['any_err']:
        sys.stdout.write("Note: One or more summaries could not be saved.\n")
    if opt.save_result['any_saved']:
        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("No summaries saved.\n")


def reset_options(opt):
    """ Deletes certain options that should always be set anew for each run
        (each time the "Start" button is pressed)
    """
    if hasattr(opt, "metric_unit"):
        delattr(opt, "metric_unit")
    if hasattr(opt, "use_grid"):
        delattr(opt, "use_grid")
    if hasattr(opt, "use_random"):
        delattr(opt, "use_random")


def show_options(opt):
    sys.stdout.write("{} version: {} (Last modified {} {}, {})\n".format(
                     version.title, version.version, *version.date))
    sys.stdout.write("Output file format: %s\n" % opt.output_file_format)
    sys.stdout.write("Suffix of output files: %s\n"
                     % opt.output_filename_suffix)
    sys.stdout.write("Output directory: %s\n" % opt.output_dir)
    sys.stdout.write("Spatial resolution: %d\n" % opt.spatial_resolution)
    sys.stdout.write("Shell width: %d metric units\n" % opt.shell_width)
    sys.stdout.write("Interpoint distances calculated: %s\n"
                     % sc.yes_or_no(opt.determine_interpoint_dists))
    if opt.determine_interpoint_dists:
        sys.stdout.write("Interpoint distance mode: %s\n"
                         % opt.interpoint_dist_mode.capitalize())
        sys.stdout.write("Shortest interpoint distances: %s\n"
                         % sc.yes_or_no(opt.interpoint_shortest_dist))
        sys.stdout.write("Lateral interpoint distances: %s\n"
                         % sc.yes_or_no(opt.interpoint_lateral_dist))
    sys.stdout.write("Monte Carlo simulations performed: %s\n"
                     % sc.yes_or_no(opt.run_monte_carlo))
    if opt.run_monte_carlo:
        sys.stdout.write("Number of Monte Carlo runs: %d\n"
                         % opt.monte_carlo_runs)
        sys.stdout.write("Monte Carlo simulation window: %s\n"
                         % opt.monte_carlo_simulation_window)
        sys.stdout.write("Strict localization in simulation window: %s\n"
                         % sc.yes_or_no(opt.monte_carlo_strict_location))
    sys.stdout.write("Clusters determined: %s\n" %
                     sc.yes_or_no(opt.determine_clusters))
    if opt.determine_clusters:
        sys.stdout.write("Within-cluster distance: %d\n"
                         % opt.within_cluster_dist)


def get_output_format(opt):
    if opt.output_file_format == 'excel':
        import imp
        try:
            imp.find_module("pyExcelerator")
        except ImportError:
            sys.stdout.write("Unable to write Excel files: resorting to csv "
                             "format.\n")
            opt.output_file_format = "csv"
    if opt.output_file_format == 'csv':
        opt.output_filename_ext = ".csv"
        opt.csv_format = {'dialect': 'excel', 'lineterminator': '\n',
                          'encoding': sys.getfilesystemencoding()}
        if opt.csv_delimiter == 'tab':
            opt.csv_format['delimiter'] = '\t'
    if opt.output_filename_date_suffix:
        import datetime
        opt.outputFilenameSuffix = "." + datetime.date.today().isoformat()
    if opt.output_filename_other_suffix != '':
        opt.outputFilenameSuffix += "." + opt.output_filename_other_suffix
        
      
def main_proc(parent):
    """ Process profile data files
    """
    opt = parent.opt
    if not opt.input_file_list:
        sys.stdout.write("No input files.\n")
        return 0
    i, n = 0, 0
    profileli = []
    sys.stdout.write("--- Session started %s local time ---\n" % time.ctime())
    # Remove duplicate filenames
    for f in opt.input_file_list:
        if opt.input_file_list.count(f) > 1:
            sys.stdout.write("Duplicate input filename %s:\n   => "
                             "removing first occurrence in list\n" % f)
            opt.input_file_list.remove(f)
    get_output_format(opt)
    reset_options(opt)
    show_options(opt)
    while True:
        if i < len(opt.input_file_list):
            inputfn = opt.input_file_list[i]
            i += 1
        else: 
            sys.stdout.write("\nNo more input files...\n")
            break
        parent.process_queue.put(("new_file", inputfn))       
        profileli.append(ProfileData(inputfn, opt))
        profileli[-1].process(opt)
        if opt.stop_requested:
            sys.stdout.write("\n--- Session aborted by user %s local time ---\n" 
                             % time.ctime())
            return 3            
        if not profileli[-1].errflag:
            n += 1
            if profileli[-1].warnflag:
                sys.stdout.write("Warning(s) found while processing "
                                 "input file.\n")
                continue
        else:
            sys.stdout.write("Error(s) found while processing input file =>\n"
                             "  => No distances could be determined.\n")
            continue
    # no more input files
    errfli = [pro.inputfn for pro in profileli if pro.errflag]
    warnfli = [pro.inputfn for pro in profileli if pro.warnflag]
    if errfli:
        sys.stdout.write("\n%s input %s generated one or more errors:\n"
                         % (stringconv.plurality("This", len(errfli)),
                            stringconv.plurality("file", len(errfli))))        
        sys.stdout.write("%s\n" % "\n".join([fn for fn in errfli]))
    if warnfli:
        sys.stdout.write("\n%s input %s generated one or more warnings:\n" 
                         % (stringconv.plurality("This", len(warnfli)),
                            stringconv.plurality("file", len(warnfli))))        
        sys.stdout.write("%s\n" % "\n".join([fn for fn in warnfli]))
    if n > 0:
        parent.process_queue.put(("saving_summaries", ""))
        save_output(profileli, opt)
    else:
        sys.stdout.write("\nNo files processed.\n")
    sys.stdout.write("--- Session ended %s local time ---\n" % time.ctime())
    parent.process_queue.put(("done", ""))
    if errfli: 
        return 0
    elif warnfli:
        return 2
    else:
        return 1
# End of main.py

#
#    Module      : main.py
#    Description : Processes input data and generates output
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

from __future__ import with_statement
import sys
import os.path
import time
import datetime
from classes import *
import geometry
from fileIO import *
import version
import stringconv
    
#
# Functions
#

def saveOutput(profileli, opt):
    """ Save a summary of results of evaluated profiles
    """
    def m(x, pixelwidth):
        return geometry.toMetricUnits(x, pixelwidth)

    def m2(x, pixelwidth): 
        return geometry.toMetricUnits(x, pixelwidth**2)  # for area units...
    
    def m_inv(x):
        try:
            return 1 / m(1 / x)
        except (TypeError, ZeroDivisionError):
            return None

    def na(x):
        if x in (None, -1):
            return "N/A"
        else:
            return x

    def writeSessionSummary():
        with FileWriter("session.summary", opt) as f:
            f.writerow(["%s version:" % version.title, 
                       "%s (Last modified %s %s, %s)" 
                       % ((version.version,) + version.date)])
            f.writerow(["Number of evaluated profiles:", len(eval_proli)])
            if err_fli:
                f.writerow(["Number of non-evaluated profiles:", len(err_fli)])
            f.writerow(["Metric unit:", eval_proli[0].metric_unit])
            f.writerow(["Spatial resolution:", opt.spatial_resolution,
                        eval_proli[0].metric_unit])
            f.writerow(["Shell width:", opt.shell_width,
                        eval_proli[0].metric_unit])
            f.writerow(["Interpoint distances calculated:", 
                        yes_or_no(opt.determine_interpoint_dists)])
            if opt.determine_interpoint_dists:
                f.writerow(["Interpoint distance mode:", 
                            opt.interpoint_dist_mode])
                f.writerow(["Shortest interpoint distances:", 
                            yes_or_no(opt.interpoint_shortest_dist)])
                f.writerow(["Lateral interpoint distances:", 
                            yes_or_no(opt.interpoint_lateral_dist)])
            f.writerow(["Monte Carlo simulations performed:", 
                        yes_or_no(opt.run_monte_carlo)])
            if opt.run_monte_carlo:
                f.writerow(["Number of Monte Carlo runs:",
                            opt.monte_carlo_runs])
                f.writerow(["Monte Carlo simulation window:", 
                            opt.monte_carlo_simulation_window])
                f.writerow(["Strict localization in simulation window:", 
                            yes_or_no(opt.monte_carlo_strict_location)])
            f.writerow(["Clusters determined:", 
                        yes_or_no(opt.determine_clusters)])
            if opt.determine_clusters:
               f.writerow(["Within-cluster distance:", 
                           opt.within_cluster_dist, eval_proli[0].metric_unit])
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
     

    def writeProfileSummary():
        with FileWriter("profile.summary", opt) as f:
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
                        "Shell particles strictly synaptic and postsynaptic or"\
                             " associated with postsynaptic membrane",                        
                        "Synaptic particles associated w/ postsynaptic membrane",
                        "Synaptic particles associated w/ presynaptic membrane",
                        "Perisynaptic particles associated w/ postsynaptic membrane",
                        "Perisynaptic particles associated w/ presynaptic membrane",
                        "Within-perforation particles associated w/ "\
                            "postsynaptic membrane",
                        "Within-perforation particles associated w/ "\
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
                                    if p.strictLateralLocation == "synaptic" 
                                       and p.axodendriticLocation == "postsynaptic"
                                       and p.isWithinPostsynapticMembraneShell]),                                                    
                          len([p for p in pro.pli
                                    if p.strictLateralLocation == "synaptic" 
                                       and 
                                       (p.axodendriticLocation == "postsynaptic"
                                        and p.isWithinPostsynapticMembraneShell)
                                       or p.isPostsynapticMembraneAssociated]),                          
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


    def writePointSummary(pType):
        if pType == "particle":
            pli = "pli"
            pstr = "particle"
        elif pType == "random":
            if not opt.useRandom:
                return
            else:
                pli = "randomli"
                pstr = "point"
        elif pType == "grid":
            if not opt.useGrid:
                return
            else:
                pli = "gridli"
                pstr = "point"
        else:
            return
        with FileWriter("%s.summary" % pType, opt) as f:
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
                          yes_or_no(p.isWithinPSD),
                          yes_or_no(p.isAssociatedWithPSD),
                          m(pro.totalPosm.length(), pro.pixelwidth),
                          m(sum([psd.posm.length() for psd in pro.psdli]),
                            pro.pixelwidth),
                          m(p.closestPSD.posm.length(), pro.pixelwidth),
                          pro.presynProfile,
                          pro.postsynProfile,
                          os.path.basename(pro.inputfn),
                          pro.comment]
                          for pro in eval_proli for n, p in 
                                                enumerate(pro.__dict__[pli])])


    def writeClusterSummary():
        if not opt.determine_clusters:
            return
        with FileWriter("cluster.summary", opt) as f:
            f.writerow(["Cluster number",
                       "Number of particles in cluster",
                       "Distance to postsynaptic membrane of centroid",
                       "Distance to nearest cluster along postsynaptic element membrane",
                       "Profile ID",
                       "Input file",
                       "Comment"])
            f.writerows([[n + 1,
                        len(c),
                        m(c.distToPath, pro.pixelwidth),
                        m(na(c.distToNearestCluster), pro.pixelwidth),
                        pro.ID,
                        os.path.basename(pro.inputfn),
                        pro.comment]
                        for pro in eval_proli
                        for n, c in enumerate(pro.clusterli)])


    def writeInterpointSummaries():

        def _m(x):
            return m(x, pro.pixelwidth)

        if not opt.determine_interpoint_dists:
            return
        ipRels = dict([(key, val)
                        for key, val in opt.interpoint_relations.items()
                        if val and "simulated" not in key])
        if not opt.useRandom:
            for key, val in opt.interpoint_relations.items():
                if "random" in key and val:
                    del ipRels[key]
        if (len(ipRels) == 0 or not
            (opt.interpoint_shortest_dist or opt.interpoint_lateral_dist)):
            return
        table = []
        if opt.interpoint_dist_mode == 'all':
            s = "all distances"
        else:
            s = "nearest neighbour distances"
        table.append(["Mode: " + s])
        headerli = ipRels.keys()
        prefixli = []
        for key, val in ipRels.items():
            prefix = key[0] + key[key.index("- ") + 2] + "_"
            #if prefix[0] == prefix[1]:
            #    prefix = prefix.replace(prefix[0], "", 1)
            prefixli.append(prefix)
        if opt.interpoint_shortest_dist and opt.interpoint_lateral_dist:
            headerli.extend(headerli)
            prefixli.extend(map(lambda s: s + "lat", prefixli))
        topheaderli = []
        if opt.interpoint_shortest_dist:
            topheaderli.append("Shortest distances")
            if opt.interpoint_lateral_dist:
                topheaderli.extend([""] * (len(ipRels)-1))
        if opt.interpoint_lateral_dist:
            topheaderli.append("Lateral distances along postsynaptic element membrane")
        table.extend([topheaderli, headerli])
        cols = [[] for c in prefixli]
        for pro in eval_proli:
            for n, li in enumerate([pro.__dict__[prefix + "distli"]
                                    for prefix in prefixli]):
                cols[n].extend(map(_m, li))
        # transpose cols and append to table
        table.extend(map(lambda *col:[e if e != None else "" for e in col],
                                  *cols))
        with FileWriter("interpoint.distances", opt) as f:
            f.writerows(table)


    def writeMonteCarloDistToPSD(type):

        def m_li(*li):
            return [m(x, pro.pixelwidth) for x in li]

        if not opt.run_monte_carlo:
            return
        table = []
        if type == "metric":
            table.append(["Lateral distances in %s to the nearest PSD"
                          % eval_proli[0].metric_unit])
        elif type == "normalized":
            table.append(["Normalized lateral distances to the nearest PSD"])
        table.append(["Run %d" % (n + 1)
                     for n in range(0, opt.monte_carlo_runs)])
        for pro in eval_proli:
            if type == "metric":
                table.extend(map(m_li, *[[p.lateralDistPSD
                                            for p in li["pli"]]
                                            for li in pro.mcli]))
            elif type == "normalized":
                table.extend(map(None, *[[p.normLateralDistPSD
                                            for p in li["pli"]]
                                            for li in pro.mcli]))
        with FileWriter("simulated.PSD.%s.lateral.distances" % type, opt) as f:
            f.writerows(table)


    def writeMonteCarloDistToPosel():

        def m_li(*li):
            return [m(x, pro.pixelwidth) for x in li]

        if not opt.run_monte_carlo:
            return
        table = []
        table.append(["Run %d" % (n + 1)
                     for n in range(0, opt.monte_carlo_runs)])
        for pro in eval_proli:
            table.extend(map(m_li, *[[p.distToPosel for p in li["pli"]]
                                        for li in pro.mcli]))
        with FileWriter("simulated.postsynaptic.element.membrane.distances",
                        opt) as f:
            f.writerows(table)


    def writeMonteCarloIPDists(dist_type):

        def m_li(*li):
            return [m(x, pro.pixelwidth) for x in li]

        if not opt.run_monte_carlo:
            return
        for ip_type in [key for key, val in opt.interpoint_relations.items()
                        if "simulated" in key and val]:
            if ((dist_type == "shortest" and not opt.interpoint_shortest_dist) or
                (dist_type == "lateral" and not opt.interpoint_lateral_dist)):
                return
            if dist_type == "lateral":
               short_dist_type = "lat"
            else:
               short_dist_type = ""
            table = []
            table.append(["Run %d" % (n + 1)
                         for n in range(0, opt.monte_carlo_runs)])
            for pro in eval_proli:
                table.extend(map(m_li,
                                 *[p for li in pro.mcli
                                    for p in li[ip_type]
                                               ["%sdist" % short_dist_type]]))
            with FileWriter("%s.interpoint.%s.distance.summary"
                            % (ip_type.replace(" ", ""), dist_type), opt) as f:
                f.writerows(table)


    def writeMonteCarloClusterSummary():
        if not (opt.determine_clusters and opt.run_monte_carlo):
            return
        table = []
        table.append(["N particles in cluster", "Run",
                     "Distance to postsynaptic element membrane from centroid",
                     "Distance to nearest cluster",
                     "Profile ID",
                     "Input file",
                     "Comment"])
        for pro in eval_proli:
            for n in range (0, opt.monte_carlo_runs):
                for c in pro.mcli[n]["clusterli"]:
                    table.append([len(c), n + 1,
                                 m(c.distToPath, pro.pixelwidth),
                                 m(na(c.distToNearestCluster), pro.pixelwidth),
                                 pro.ID,
                                 os.path.basename(pro.inputfn),
                                 pro.comment])
        with FileWriter("simulated.cluster.summary", opt) as f:
            f.writerows(table)


    sys.stdout.write("\nSaving summaries...\n")
    opt.save_result = {'any_saved': False, 'any_err': False}
    eval_proli = [pro for pro in profileli if not pro.errflag]
    clean_fli = [pro.inputfn for pro in profileli if not (pro.errflag or pro.warnflag)]
    warn_fli = [pro.inputfn for pro in profileli if pro.warnflag]
    err_fli = [pro.inputfn for pro in profileli if pro.errflag]
    nop_fli = [pro.inputfn for pro in eval_proli if not pro.pli]
    if opt.output_file_format == 'excel':
        import xls
    elif opt.output_file_format == 'csv':
        csv_format = { 'dialect' : 'excel', 'lineterminator' : '\n'}
        if opt.csv_delimiter == 'tab':
            csv_format['delimiter'] = '\t'  
    writeSessionSummary()
    writeProfileSummary()
    writePointSummary("particle")
    writePointSummary("random")
    writePointSummary("grid")
    writeInterpointSummaries()
    writeClusterSummary()
    writeMonteCarloDistToPosel()
    writeMonteCarloDistToPSD("metric")
    writeMonteCarloDistToPSD("normalized")
    writeMonteCarloIPDists("shortest")
    writeMonteCarloIPDists("lateral")
    writeMonteCarloClusterSummary()
    if opt.save_result['any_err'] == True:
        sys.stdout.write("Note: One or more summaries could not be saved.\n")
    if opt.save_result['any_saved'] == True:
        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("No summaries saved.\n")


def showOptions(opt):   
    sys.stdout.write("%s version: %s (Last modified %s %s, %s)\n" 
                      % ((version.title, version.version) + version.date))                           
    sys.stdout.write("Output file format: %s\n" % opt.output_file_format)
    sys.stdout.write("Suffix of output files: %s\n"
                     % opt.output_filename_suffix)
    sys.stdout.write("Output directory: %s\n" % opt.output_dir)
    sys.stdout.write("Spatial resolution: %d\n" % opt.spatial_resolution)
    sys.stdout.write("Shell width: %d metric units\n" % opt.shell_width)
    sys.stdout.write("Interpoint distances calculated: %s\n" 
                     % yes_or_no(opt.determine_interpoint_dists))
    if opt.determine_interpoint_dists:
        sys.stdout.write("Interpoint distance mode: %s\n" 
                         % opt.interpoint_dist_mode.capitalize())
        sys.stdout.write("Shortest interpoint distances: %s\n" 
                         % yes_or_no(opt.interpoint_shortest_dist))
        sys.stdout.write("Lateral interpoint distances: %s\n" 
                         % yes_or_no(opt.interpoint_lateral_dist))
    sys.stdout.write("Monte Carlo simulations performed: %s\n" 
                     % yes_or_no(opt.run_monte_carlo))
    if opt.run_monte_carlo:
        sys.stdout.write("Number of Monte Carlo runs: %d\n" 
                         % opt.monte_carlo_runs)
        sys.stdout.write("Monte Carlo simulation window: %s\n" 
                         % opt.monte_carlo_simulation_window)
        sys.stdout.write("Strict localization in simulation window: %s\n" 
                         % yes_or_no(opt.monte_carlo_strict_location))
    sys.stdout.write("Clusters determined: %s\n" % 
                     yes_or_no(opt.determine_clusters))
    if opt.determine_clusters:
       sys.stdout.write("Within-cluster distance: %d\n" 
                        % opt.within_cluster_dist)
   
def getOutputFormat(opt):
    if opt.output_file_format == 'excel':
        try:
            import xls
        except ImportError:
            sys.stdout.write("Unable to write Excel files: resorting to csv "
                             "format.\n")
            opt.output_file_format = "csv"
    if opt.output_file_format == 'csv':
        opt.output_filename_ext = ".csv"
        opt.csv_format = { 'dialect' : 'excel', 'lineterminator' : '\n',
                       'encoding': sys.getfilesystemencoding() }
        if opt.csv_delimiter == 'tab':
            opt.csv_format['delimiter'] = '\t'
    if opt.output_filename_date_suffix:
        import datetime
        opt.outputFilenameSuffix = "." + datetime.date.today().isoformat()
    if opt.output_filename_other_suffix != '':
        opt.outputFilenameSuffix += "." + opt.output_filename_other_suffix
        
      
def mainProc(parent, opt):
    """ Process profile data files
    """
    
    def removeDuplicateFilenames(fli):
        """ Remove duplicate filenames in input file list
        """
        for f in fli:
            if fli.count(f) > 1:
                sys.stdout.write("Duplicate input filename %s:\n   => " 
                                 "removing first occurrence in list\n" % f)
                fli.remove(f)    
                

    if not opt.input_file_list:
        sys.stdout.write("No input files.\n")
        return 0                 
    i, n = 0, 0
    profileli = []
    sys.stdout.write("--- Session started %s local time ---\n" 
                      % time.ctime())
    removeDuplicateFilenames(opt.input_file_list)
    getOutputFormat(opt)
    showOptions(opt)
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
        sys.stdout.write("\n%s input %s generated one or more "
                        "errors:\n" 
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
        saveOutput(profileli, opt)
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

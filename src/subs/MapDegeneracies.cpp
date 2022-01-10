//

#include "MapDegeneracies.h"

template <class S>
MapDegeneracies<S>::MapDegeneracies(const vector<unique_ptr<typename S::Extension>> &extensions_,
                                    const typename S::Collection &mapCollection, const set<string> &mapTypes,
                                    bool defaulted)
        : extensions(extensions_) {
    // Get all map names meeting criteria
    for (auto mapname : mapCollection.allMapNames()) {
        // Skip fixed maps
        if (mapCollection.getFixed(mapname)) continue;
        // Skip non-defaulted if we only want defaulted maps:
        if (defaulted && !mapCollection.getDefaulted(mapname)) continue;

        auto m = mapCollection.cloneMap(mapname);
        if (mapTypes.count(m->getType())) maps[mapname];
        delete m;
    }

    // Find all extensions using a relevant map
    for (int iextn = 0; iextn < extensions.size(); iextn++) {
        auto extptr = extensions[iextn].get();
        if (!extptr) continue;
        for (auto dependent : mapCollection.dependencies(extptr->mapName)) {
            if (maps.count(dependent)) {
                // This extension uses this dependent map
                maps[dependent].insert(iextn);
                extns[iextn].insert(dependent);
            }
        }
    }
    return;
}

template <class S>
void MapDegeneracies<S>::eraseMap(string mapname) {
    // Change all extensions this map touches
    for (auto iextn : maps[mapname]) {
        extns[iextn].erase(mapname);
        // Get rid of the extension if it's now empty
        if (extns[iextn].empty()) extns.erase(iextn);
    }
    // And kill the map
    maps.erase(mapname);
    return;
}

template <class S>
set<int> MapDegeneracies<S>::findNondegenerate() const {
    set<int> out;
    for (auto &m : extns)
        if (m.second.size() == 1) out.insert(m.first);
    return out;
}

template <class S>
list<string> MapDegeneracies<S>::replaceWithIdentity(const set<string> &candidates) {
    list<string> mapsToReplace;

    while (!maps.empty()) {
        set<string> goodmaps;
        auto goodextns = findNondegenerate();
        if (goodextns.empty()) {
            // We have maps with no clear degeneracy breaking path.
            // Find the candidate map that we could fix...
            // Choose the one that will "unlock" the most other extensions by
            // reducing their map counts from 2 to 1.
            int maxFix = -1;
            string maxName;
            for (auto mapname : candidates) {
                if (maps.count(mapname) == 0) continue;  // Skip if candidate is no longer degenerate
                int nFix = 0;
                for (auto iextn : maps[mapname])
                    if (extns[iextn].size() == 2) nFix++;
                if (nFix > maxFix) {
                    maxFix = nFix;
                    maxName = mapname;
                }
            }
            //**/cerr << "--maxFix " << maxFix << "--maxName " << maxName << endl;
            // If there is no way to un-degenerate any other maps by fixing one,
            // we are screwed (??? unless we want to start looking for pairs to
            // fix simultaneously, by finding all exposures with all but one map
            // that are candidates to fix)
            if (maxFix <= 0) {
                cerr << "WARNING: no clear path for breaking degeneracy of these maps:" << endl;
                for (auto m : maps) cerr << "  " << m.first << endl;
                return mapsToReplace;
            }

            // Otherwise we will try to turn this map into Identity and continue
            mapsToReplace.push_back(maxName);
            eraseMap(maxName);
        } else {
            // The map that is in each nondegen extension is no longer degenerate.
            for (auto iextn : goodextns) {
                // ?? check for over-constrained ??
                goodmaps.insert(extns[iextn].begin(), extns[iextn].end());
            }
            for (auto mapname : goodmaps) eraseMap(mapname);
        }
    }
    return mapsToReplace;
}

template <class S>
list<set<int>> MapDegeneracies<S>::initializationOrder() {
    list<set<int>> out;
    while (!maps.empty()) {
        auto goodextns = findNondegenerate();
        cerr << ".." << to_string(maps.size()) << " maps, " << to_string(goodextns.size()) << " clean extns"
             << endl;
        if (goodextns.empty()) {
            // We have maps with no clear degeneracy breaking path.
            cerr << "ERROR: no path to initialize these maps without degeneracies:" << endl;
            for (auto &m : maps) cerr << "  " << m.first << endl;
            exit(1);
        }

        // Do everything that can be done - for a given choose all relevant exposures in
        // a given exposure.  Exposure with the most data to contribute now?
        map<string, map<int, set<int>>> mapCounts;
        // First index in map name, 2nd is which exposure it's from, the set is
        // the relevant exposures.

        for (auto iextn : goodextns) {
            string itsMap = *extns[iextn].begin();
            int mapUses = maps[itsMap].size();
            if (mapUses == 1) {
                // This extn's map is not used anywhere else.  Go ahead and add to init list
                // and erase the map
                set<int> tmp = {iextn};
                out.push_back(tmp);
                eraseMap(itsMap);
                //**/cerr << "..singleton " << itsMap << endl;
            } else {
                int mapExpo = extensions[iextn]->exposure;
                mapCounts[itsMap][mapExpo].insert(iextn);
            }
        }
        for (auto &m : mapCounts) {
            // For each unique mapname, initialize using exposure with the most
            // available extensions
            int maxCount = 0;
            int maxExpo = 0;
            for (auto &m2 : m.second) {
                if (m2.second.size() > maxCount) {
                    maxCount = m2.second.size();
                    maxExpo = m2.first;
                }
            }
            cerr << "...Do map " << m.first << " from exposure " << to_string(maxExpo) << " with "
                 << to_string(maxCount) << " exposures" << endl; /**/
            out.push_back(m.second[maxExpo]);
            eraseMap(m.first);
        }
    }
    return out;
}

template class MapDegeneracies<Astro>;

template class MapDegeneracies<Photo>;

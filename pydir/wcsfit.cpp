#include <pybind11/pybind11.h>
#include "pybind11/stl.h"
#include "pybind11/eigen.h"
#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h>

#include "FitSubroutines.h"
#include "Instrument.h"
#include "TPVMap.h"
#include "WCSFoFRoutine.h"
#include "WCSFitRoutine.h"


using namespace pybind11::literals;
namespace py = pybind11;

void spaceReplace2(string &s) {
    stripWhite(s);
    cerr << "stripWhite done" << endl;
    std::string regex_ = "[[:space:]]+";
    std::regex e2(regex_, std::regex::extended);
    cerr << "check done" << endl;
    s = regexReplace("[[:space:]]+", "_", s);
}

void floatPrint(double &i, int &j) {
    cerr << "double is: " << i << endl;
    cerr << "int is: " << j << endl;
    cerr << "done" << endl;
}

// Make lambda for Instrument and PolyMap so this is unnecessary?
template <class T = float>
void declareBounds(py::module &m) {
    using Class = Bounds<T>;

    py::class_<Class>(m, "Bounds").def(py::init<const T, const T, const T, const T>());
}


PYBIND11_MODULE(wcsfit, m) {
    m.def("readWCSs", &readWCSs);
    m.attr("DEGREE") = DEGREE;
    m.attr("ARCSEC") = ARCSEC;
    m.def("spaceReplace2", &spaceReplace2);
    m.def("floatPrint", &floatPrint);
    
    m.doc() = "pybind11 example plugin"; // optional module docstring


    declareBounds<double>(m);
    // Not needed?
    py::class_<img::FTable>(m, "FTable")
            .def(py::init<long>())
            .def("addColumnDouble", &img::FTable::addColumn<double>, py::arg("values"), py::arg("columnName"),
                 py::arg("repeat") = -1, py::arg("stringLength") = -1)
            .def("addColumnStr", &img::FTable::addColumn<std::string>, py::arg("values"),
                 py::arg("columnName"), py::arg("repeat") = -1, py::arg("stringLength") = -1);

    //////// Classes needed for setting up the WCS //////////
    // Not needed?
    py::class_<astrometry::PixelMap>(m, "PixelMap");
    // Not needed?
    py::class_<astrometry::SphericalCoords>(m, "SphericalCoords");
    // Make lambda to combine the next three?
    py::class_<astrometry::SphericalICRS, astrometry::SphericalCoords>(m, "SphericalICRS")
            .def(py::init<>())
            .def(py::init<double, double>());

    py::class_<astrometry::Orientation>(m, "Orientation")
            .def(py::init<>())
            .def(py::init<astrometry::SphericalCoords &, double>(), py::arg("pole_"), py::arg("pa") = 0)
            .def("set", &astrometry::Orientation::set)
            .def("getPole", &astrometry::Orientation::getPole);

    py::class_<astrometry::Gnomonic, astrometry::SphericalCoords>(m, "Gnomonic")
    //py::class_<astrometry::Gnomonic>(m, "Gnomonic")
            //.def(py::init<astrometry::Orientation const &, bool>(), py::arg("o"),
            //     py::arg("shareOrient") = false)
            .def(py::init([](double ra, double dec){
                astrometry::SphericalICRS pole(ra, dec);
                astrometry::Orientation orientIn(pole);
                return astrometry::Gnomonic(orientIn.getPole(), orientIn);
            }));
            //.def(py::init<astrometry::SphericalCoords const &, astrometry::Orientation &>())
            //.def(py::init<>());

    py::class_<astrometry::LinearMap, astrometry::PixelMap>(m, "LinearMap")
            .def(py::init<astrometry::DVector::Base const &>());

    py::class_<poly2d::Poly2d>(m, "Poly2d")
            .def(py::init<int>())
            .def(py::init<astrometry::DMatrix::Base const &>())
            .def("nCoeffs", &poly2d::Poly2d::nCoeffs)
            .def("getC", &poly2d::Poly2d::getC)
            .def("setC", &poly2d::Poly2d::setC)
            .def("vectorIndex", &poly2d::Poly2d::vectorIndex);

    py::class_<astrometry::PolyMap, astrometry::PixelMap>(m, "PolyMap")
            .def(py::init<poly2d::Poly2d, poly2d::Poly2d, std::string, Bounds<double>, double>());

    py::class_<astrometry::SubMap, astrometry::PixelMap>(m, "SubMap")
            .def(py::init<list<astrometry::PixelMap *> const &, std::string, bool>());

    py::class_<astrometry::Wcs, std::shared_ptr<astrometry::Wcs>>(m, "Wcs")
            .def(py::init<astrometry::PixelMap *, astrometry::SphericalCoords const &, std::string, double,
                          bool>(),
                 py::arg("pm_"), py::arg("nativeCoords_"), py::arg("name") = "", py::arg("wScale_") = DEGREE,
                 py::arg("shareMap_") = false)
            .def("reprojectTo",
            // Don't need this scoped_estream_redirect thing I think
                 [](astrometry::Wcs &self, const astrometry::SphericalCoords &targetCoords_) {
                     py::scoped_estream_redirect stream(
                             std::cerr,                                 // std::ostream&
                             py::module_::import("sys").attr("stderr")  // Python output
                     );
                     self.reprojectTo(targetCoords_);
                 })
            .def("getTargetCoords", &astrometry::Wcs::getTargetCoords)
            .def("toWorld", [](astrometry::Wcs &self, double x, double y) {
                double ra;
                double dec;
                self.toWorld(x, y, ra, dec);
                std::vector<double> radec = {ra, dec};
                return radec;
            });

    /////////// Data Classes /////////////////
    py::class_<Fields>(m, "Fields")
            .def(py::init<std::vector<std::string>, std::vector<double>, std::vector<double>,
                          std::vector<double>>());
    // don't need?
    py::class_<FieldsHelper>(m, "FieldsHelper")
            .def(py::init<std::vector<std::string>, std::vector<double>, std::vector<double>,
                          std::vector<double>>());

    py::class_<Instrument, std::shared_ptr<Instrument>>(m, "Instrument")
            .def(py::init<std::string>(), py::arg("name_") = "")
            .def_readwrite("name", &Instrument::name)
            //.def_readwrite("nDevices", &Instrument::nDevices)
            .def_readwrite("band", &Instrument::band)
            //.def_readwrite("deviceNames", &Instrument::deviceNames)
            .def("addDevice", &Instrument::addDevice)
            .def("hasDevice", [](Instrument &self, std::string name) { return self.deviceNames.has(name); });

    py::class_<ExposuresHelper>(m, "ExposuresHelper")
            .def(py::init<std::vector<std::string>, std::vector<int>, std::vector<int>, std::vector<double>,
                          std::vector<double>, std::vector<double>, std::vector<double>,
                          std::vector<double>>())
            .def_readwrite("instrumentNumbers", &ExposuresHelper::instrumentNumbers)
            .def_readwrite("fieldNumbers", &ExposuresHelper::fieldNumbers);

    //declareBounds<double>(m);

    py::class_<astrometry::YAMLCollector>(m, "YAMLCollector")
            .def(py::init<std::string, std::string>())
            .def(
                    "addInput",
                    [](astrometry::YAMLCollector &self, std::string is, std::string filter, bool prepend) {
                        std::istringstream iss(is);
                        self.addInput(iss, filter, prepend);
                    },
                    py::arg("is"), py::arg("filter") = "", py::arg("prepend") = false);
    //.def("addMap", &YAMLCollector::addMap, py::arg("name"), py::arg("dict"),
    //     py::arg("criticalTime")=nullptr);
    //.def("addMap", [](YAMLCollector & self, string name, py::dict pyDict) {
    //                  astrometry::YAMLCollector::Dictionary d;
    //                  for (auto e : pyDict) {
    //                      d[e.first] = e.second;
    //                  }
    //                  self.addMap(name, d);});

    py::class_<astrometry::IdentityMap, astrometry::PixelMap>(m, "IdentityMap")
            .def(py::init<>())
            .def("getName", &astrometry::IdentityMap::getName);
    //don't need?
    py::class_<ExtensionObjectSet>(m, "ExtensionObjectSet").def(py::init<std::string>());

    ////////////// WCS-fitting class ////////////////////////
    py::class_<WCSFit>(m, "WCSFit")
            .def(py::init<Fields &, std::vector<std::shared_ptr<Instrument>>, ExposuresHelper,
                          std::vector<int>, std::vector<int>, astrometry::YAMLCollector,
                          std::vector<std::shared_ptr<astrometry::Wcs>>, std::vector<int>,
                          std::vector<LONGLONG>, std::vector<LONGLONG>, std::vector<int>, double, double, int,
                          std::string, std::string, bool, int>(),
                 py::arg("fields"), py::arg("instruments"), py::arg("exposures"),
                 py::arg("extensionExposureNumbers"), py::arg("extensionDevices"), py::arg("inputYAML"),
                 py::arg("wcss"), py::arg("sequence"), py::arg("extns"), py::arg("objects"),
                 py::arg("exposureColorPriorities") = std::vector<int>(), py::arg("sysErr") = 2.0,
                 py::arg("refSysErr") = 2.0, py::arg("minMatches") = 2, py::arg("skipObjectsFile") = "",
                 py::arg("fixMaps") = "", py::arg("usePM") = true, py::arg("verbose") = 0)
            .def("setObjects", &WCSFit::setObjects, py::arg("i"), py::arg("tableMap"), py::arg("xKey"),
                 py::arg("yKey"), py::arg("xyErrKeys"), py::arg("idKey") = "", py::arg("pmCovKey") = "",
                 py::arg("magKey") = "", py::arg("magKeyElement") = 0, py::arg("magErrKey") = "",
                 py::arg("magErrKeyElement") = 0, py::arg("pmRaKey") = "", py::arg("pmDecKey") = "",
                 py::arg("parallaxKey") = "")

            .def_readwrite("verbose", &WCSFit::verbose)
            .def("setRefWCSNames", &WCSFit::setRefWCSNames)  // to delete
            .def("setupMaps", &WCSFit::setupMaps)            // to delete
            .def("setMatches", &WCSFit::setMatches)          // to delete

            .def("fit", &WCSFit::fit, py::arg("maxError") = 100., py::arg("minFitExposures") = 200,
                 py::arg("reserveFraction") = 0.2, py::arg("randomNumberSeed") = 1234,
                 py::arg("minimumImprovement") = 0.02, py::arg("clipThresh") = 5.0,
                 py::arg("chisqTolerance") = 0.001, py::arg("clipEntireMatch") = false,
                 py::arg("divideInPlace") = false, py::arg("purgeOutput") = false,
                 py::arg("minColor") = -10.0, py::arg("maxColor") = 10.0)
            .def("saveResults", &WCSFit::saveResults);

    /// The following are WCSFoF-specific classes:
    ///////////// Friends of Friends Fitting Class //////////////////
    py::class_<FoFClass>(m, "FoFClass")
            .def(py::init<Fields &, std::vector<std::shared_ptr<Instrument>>, ExposuresHelper,
                          std::vector<double>, double>(),
                 py::arg("fields"), py::arg("instruments"), py::arg("exposures"), py::arg("fieldExtents"),
                 py::arg("matchRadius"))
            .def("reprojectWCS", &FoFClass::reprojectWCS)
            //.def("addCatalog", py::overload_cast<std::shared_ptr<astrometry::Wcs> &, std::string, int, int,
            .def("addCatalog",
                 py::overload_cast<astrometry::Wcs const &, std::string, int, int, int, int, long,
                                   std::vector<bool>, std::vector<double>, std::vector<double>,
                                   std::vector<long>>(&FoFClass::addCatalog))
            .def("sortMatches", &FoFClass::sortMatches, py::arg("fieldNumber"), py::arg("minMatches") = 2,
                 py::arg("allowSelfMatches") = false)
            .def_readwrite("sequence", &FoFClass::sequence)
            .def_readwrite("extn", &FoFClass::extn)
            .def_readwrite("obj", &FoFClass::obj);

    // The rest of WCSFit class:
    //.def(py::init<string>())
    /*
    .def_readwrite("fields", &WCSFit::fields)
    .def_readwrite("instruments", &WCSFit::instruments)
    .def("setExposures", [](WCSFit & self, vector<std::shared_ptr<Exposure>> expos, double sysErr, double
    refSysErr) { py::scoped_estream_redirect stream( std::cerr,                               // std::ostream&
            py::module_::import("sys").attr("stderr") // Python output
        );
        self.setExposures(expos, sysErr, refSysErr);
    })

    .def("addMap", &WCSFit::addMap)




    .def("getMatchLength", &WCSFit::getMatchLength)
    //.def("fit", [](WCSFit & self) {
    //    py::scoped_estream_redirect stream(
    //        std::cerr,                               // std::ostream&
    //        py::module_::import("sys").attr("stdout") // Python output
    //    );
    //    self.fit();
    //})

    .def_readwrite("randomNumberSeed", &WCSFit::randomNumberSeed)
    .def_readwrite("clipThresh", &WCSFit::clipThresh)
    .def_readwrite("maxError", &WCSFit::maxError)
    .def_readwrite("minMatches", &WCSFit::minMatches)
    .def_readwrite("minFitExposures", &WCSFit::minFitExposures)
    .def_readwrite("clipEntireMatch", &WCSFit::clipEntireMatch)
    .def_readwrite("chisqTolerance", &WCSFit::chisqTolerance)
    .def_readwrite("divideInPlace", &WCSFit::divideInPlace)
    .def_readwrite("purgeOutput", &WCSFit::purgeOutput)
    .def_readwrite("minColor", &WCSFit::minColor)
    .def_readwrite("maxColor", &WCSFit::maxColor)
    .def_readwrite("exposures", &WCSFit::exposures)

    .def_readwrite("extensions", &WCSFit::extensions)
    .def_readwrite("fieldNames", &WCSFit::fieldNames)
    .def_readwrite("fieldProjections", &WCSFit::fieldProjections)
    .def_readwrite("fieldEpochs", &WCSFit::fieldEpochs)
    .def_readwrite("matches", &WCSFit::matches);
    //.def_readwrite("inputYAML", &WCSFit::inputYAML)*/

    /*
    m.doc() = "tmp docstring"; // optional module docstring

    m.attr("REF_INSTRUMENT") = REF_INSTRUMENT;
    m.attr("PM_INSTRUMENT") = PM_INSTRUMENT;
    m.attr("ARCSEC") = ARCSEC;
    m.attr("MILLIARCSEC") = MILLIARCSEC;

    //m.attr("RESIDUAL_UNIT") = RESIDUAL_UNIT;
    //m.attr("WCS_UNIT") = WCS_UNIT;
    m.def("spaceReplace", &spaceReplace);
    m.def("loadPixelMapParser", &loadPixelMapParser);
    //declareReadMatches<Astro>(m);
    m.def("readMatches", py::overload_cast<vector<int>&, vector<LONGLONG>&, vector<LONGLONG>&,
                                         Astro::MCat&, vector<Astro::Extension*>&,
                                         vector<Astro::ColorExtension*>&, ExtensionObjectSet const&,
                                         int, bool>
                                         (&readMatches<Astro>));
    m.def("readObjects_oneExtension", &readObjects_oneExtension<Astro>);
    m.def("reportStatistics", &Astro::reportStatistics);
    m.def("readFields", &readFields);


    py::class_<astrometry::SphericalCoords>(m, "SphericalCoords");

    // TODO: check which of these are actually needed


    py::class_<Exposure, std::shared_ptr<Exposure>>(m, "Exposure")
        .def(py::init<string, astrometry::SphericalCoords const&>())
        .def_readwrite("name", &Exposure::name)
        .def_readwrite("projection", &Exposure::projection)
        .def_readwrite("field", &Exposure::field)
        .def_readwrite("instrument", &Exposure::instrument)
        .def_readwrite("airmass", &Exposure::airmass)
        .def_readwrite("exptime", &Exposure::exptime)
        .def_readwrite("mjd", &Exposure::mjd)
        .def_readwrite("apcorr", &Exposure::apcorr)
        .def_readwrite("epoch", &Exposure::epoch)
        .def_readwrite("projection", &Exposure::projection)
        .def("setAstrometricCovariance", [](Exposure & self, astrometry::Matrix22::Base const&
    matrix){self.astrometricCovariance = matrix;}) .def("getAstrometricCovariance", [](Exposure const&
    self){return static_cast <astrometry::Matrix22::Base const&>(self.astrometricCovariance);})
        .def("addToAstrometricCovariance", [](Exposure & self, astrometry::Matrix22::Base const&
    matrix){self.astrometricCovariance += matrix;});

    py::class_<SPExtension, std::shared_ptr<SPExtension>>(m, "SPExtension")
        .def(py::init<>())
        //.def("addWcs", &Class::addWcs)
        .def_readwrite("exposure", &SPExtension::exposure)
        .def_readwrite("device", &SPExtension::device)
        //.def_readwrite("map", &SPExtension::map)
        .def_readwrite("airmass", &SPExtension::airmass)
        .def_readwrite("apcorr", &SPExtension::apcorr)
        .def_readwrite("magshift", &SPExtension::magshift)
        .def_readwrite("startWcs", &SPExtension::startWcs)
        //.def_readwrite("wcs", &Class::wcs)
        .def_readwrite("wcsName", &SPExtension::wcsName)
        .def_readwrite("mapName", &SPExtension::mapName);
        //.def_readwrite("keepers", &Class::keepers);


    py::class_<astrometry::Detection>(m, "Dectection");

    py::class_<astrometry::PixelMap>(m, "PixelMap");


    declareColorExtension<astrometry::Match>(m);








    py::class_<astrometry::Match>(m, "Match");



    py::class_<NameIndex>(m, "NameIndex")
        .def(py::init<>())
        .def("append", &NameIndex::append)
        .def("nameOf", &NameIndex::nameOf);


    //declareDictionary<string, string>(m);
    //py::class_<YAMLCollector::Dictionary, std::map<string, string>>(m, "Dictionary");

    // This may not be necessary long-term
    // the "read" function still doesn't work
    py::class_<PixelMapCollection>(m, "PixelMapCollection")
        .def(py::init<>())
        //.def("read", &PixelMapCollection::read, py::arg("is"), py::arg("namePrefix")="")
        .def("read", [](PixelMapCollection & self, string is, string NamePrefix){
                        istringstream iss(is);
                        self.read(iss, NamePrefix);},
                        py::arg("is"), py::arg("namePrefix")="")
        .def("allWcsNames", &PixelMapCollection::allWcsNames)
        .def("cloneWcs", &PixelMapCollection::cloneWcs);




    // The following are WCSFoF-specific classes:
    py::class_<Point>(m, "Point")
        .def(py::init<double, double, long int, long int, long int>())
        .def_readwrite("extensionNumber", &Point::extensionNumber)
        .def_readwrite("objectNumber", &Point::objectNumber)
        .def_readwrite("exposureNumber", &Point::exposureNumber);

    py::class_<Field>(m, "Field")
        .def(py::init<>())
        .def_readwrite("name", &Field::name)
        .def_readwrite("matchRadius", &Field::matchRadius)
        .def_readwrite("projection", &Field::projection)
        .def_readwrite("extent", &Field::extent)
        .def_readwrite("catalogs", &Field::catalogs);

    declareList<Point const*>(m, "PointList");

    declareMatch<Point, 2>(m);

    declareMatchSet<Point, 2>(m);

    declareCatalog<Point, 2>(m, "PointCat");

    py::bind_map<map<string, PointCat>>(m, "PointCatDict");

    py::class_<Device, Bounds<double>>(m, "Device")
        .def(py::init<>())
        //.def(py::init<double const, double const, double const, double const>())
        .def("setXMin", &Device::setXMin)
        .def("setXMax", &Device::setXMax)
        .def("setYMin", &Device::setYMin)
        .def("setYMax", &Device::setYMax)
        .def_readwrite("name", &Device::name);

    declareVector<Device>(m);

    py::bind_vector<Instr, vector<Device>>(m, "Instr")
        .def(py::init<string>())
        .def_readwrite("name", &Instr::name);

    py::class_<Expo>(m, "Expo")
        .def(py::init<>())
        .def_readwrite("name", &Expo::name)
        .def_readwrite("field", &Expo::field)
        .def_readwrite("instrument", &Expo::instrument)
        .def_readwrite("pointing", &Expo::pointing);

    py::class_<FoFClass>(m, "FoFClass")
        .def(py::init<>())
        .def_readwrite("useAffinities", &FoFClass::useAffinities)
        .def_readwrite("minMatches", &FoFClass::minMatches)
        .def_readwrite("allowSelfMatches", &FoFClass::allowSelfMatches)
        .def_readwrite("fields", &FoFClass::fields)
        .def_readwrite("instruments", &FoFClass::instruments)
        .def_readwrite("exposures", &FoFClass::exposures)
        .def_readwrite("allPoints", &FoFClass::allPoints)
        .def_readwrite("sequence", &FoFClass::sequence)
        .def_readwrite("extn", &FoFClass::extn)
        .def_readwrite("obj", &FoFClass::obj)
        .def_readwrite("matchRadius", &FoFClass::matchRadius)
        .def_readwrite("a", &FoFClass::a)
        .def_readwrite("b", &FoFClass::b)
        .def_readwrite("c", &FoFClass::c)
        .def("addTest", &FoFClass::addTest)
        .def("addCatalog", &FoFClass::addCatalog)
        .def("sortMatches", &FoFClass::sortMatches,
             py::call_guard<py::scoped_ostream_redirect,
                            py::scoped_estream_redirect>())
        .def("writeMatches", &FoFClass::writeMatches);
    */
}

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"
#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h>

#include "FitSubroutines.h"
#include "Instrument.h"
#include "TPVMap.h"
//#include "WCSFoF_match.h"
#include "WCSFit_fit.h"

using namespace pybind11::literals;
namespace py = pybind11;


template<class T = float>
void declareBounds(py::module &m) {
    using Class = Bounds<T>;

    py::class_<Class>(m, "Bounds")
        .def(py::init<const T, const T, const T, const T>());

}

using InstrumentVector = std::vector<std::unique_ptr<Instrument>>;

PYBIND11_MAKE_OPAQUE(InstrumentVector);

void declareInstrumentVector(py::module & mod) {
    py::class_<InstrumentVector> cls(mod, "InstrumentVector");
    cls.def(py::init<>());
    cls.def(
        "add",
        [](InstrumentVector & self, std::string name, std::string band, py::dict devices) {
            std::unique_ptr<Instrument> instrument(new Instrument(name));
            instrument->band = band;
            for (auto item : devices) {
                auto device_name = py::cast<std::string>(item.first);
                auto domain = py::cast<Bounds<double>>(item.second);
                instrument->addDevice(device_name, domain);
            }
        },
        "name"_a, "band"_a, "devices"_a,
        "Add a new instrument to the vector.\n"
        "\n"
        "'devices' is a dict mapping detector name to the bounding box of \n"
        "that detector.\n"
    );
}

/*template<class T1, class T2> 
void declareExtension(py::module &m) {
    using Class = ExtensionBase<T1, T2>;

    py::class_<Class, shared_ptr<Class>>(m, "Extension")
        .def(py::init<>())
        .def("addWcs", &Class::addWcs)
        .def_readwrite("exposure", &Class::exposure)
        .def_readwrite("device", &Class::device)
        .def_readwrite("map", &Class::map)
        .def_readwrite("airmass", &Class::airmass)
        .def_readwrite("apcorr", &Class::apcorr)
        .def_readwrite("magshift", &Class::magshift)
        .def_readwrite("startWcs", &Class::startWcs)
        .def_readwrite("wcs", &Class::wcs)
        .def_readwrite("wcsName", &Class::wcsName)
        .def_readwrite("mapName", &Class::mapName)
        .def_readwrite("keepers", &Class::keepers);
}*/

/*

void declareFields(py::module & mod) {
    py::class_<Fields> cls(mod, "Fields");
    cls.def(
        py::init<std::vector<std::string>, std::vector<double>, std::vector<double>, std::vector<double>>(),
        "names"_a,
        "ra"_a,
        "dec"_a,
        "epochs"_a
    );
}



template<class T>
void declareColorExtension(py::module &m) {
    using Class = ColorExtensionBase<T>;

    py::class_<Class>(m, "ColorExtension")
        .def(py::init<>());
}


template<class T>
void declareVector(py::module &m) {
    using Class = std::vector<T>;

    py::class_<Class>(m, "deviceVector");
}

template<class T>
void declareList(py::module &m, string ListName) {
    using Class = std::list<T>;

    py::class_<Class>(m, ListName.c_str());
}

template<class P, int DIM = 2>
void declareMatch(py::module &m) {
    using Class = fof::Match<P, DIM>;
    //using Class2 = std::list<P const*>;
    //declareList<P>(m);

    py::class_<Class>(m, "FoFMatch")
        //.def(py::init<>())
        //.def_static("fromVector", [](Class & self, std::vector<P> v){
        //    for 
        //    self{ std::begin(v), std::end(v)};
        //})
        .def("toVector", [](Class & self){
            //std::vector<const P*> v{ std::begin(self), std::end(self) };
            std::vector<const P*> v{ std::begin(self), std::end(self) };
            return v;
        });
    //py::bind_vector<Class, Class2>(m, "FoFMatch");
}

template<class P, int DIM>
void declareMatchSet(py::module &m) {
    using Class = set<fof::Match<P, DIM>*>;

    py::class_<Class>(m, "MatchSet");
}

template<class P, int DIM>
void declareCatalog(py::module &m, string CatalogType) {
    using Class = fof::Catalog<P, DIM>;

    py::class_<Class>(m, CatalogType.c_str())
        .def("toVector", [](Class & self){
            std::vector<fof::Match<P, DIM>*> v{ std::begin(self), std::end(self) };
            return v;
        });
    //py::class_<Class>(m, CatalogType.c_str());
}

template<class S>
void declareReadMatches(py::module &m) {
    m.def("readMatches", readMatches<S>);
}

//template<typename T1, typename T2>
//void declareDictionary(py::module &m) {
//    py::class_<YAMLCollector::Dictionary, std::map<T1, T2>>();
//}
*/
PYBIND11_MODULE(wcsfit, m) {

    m.def("readWCSs", &readWCSs);
    m.attr("DEGREE") = DEGREE;

    py::class_<img::FTable>(m, "FTable")
        .def(py::init<long>())
        .def("addColumnDouble", &FTable::addColumn<double>, py::arg("values"), py::arg("columnName"),
             py::arg("repeat")=-1, py::arg("stringLength")=-1)
        .def("addColumnStr", &FTable::addColumn<string>, py::arg("values"), py::arg("columnName"),
             py::arg("repeat")=-1, py::arg("stringLength")=-1);

    py::class_<astrometry::PixelMap>(m, "PixelMap");

    py::class_<astrometry::SphericalCoords>(m, "SphericalCoords");

    py::class_<astrometry::SphericalICRS, astrometry::SphericalCoords>(m, "SphericalICRS")
        .def(py::init<>())
        .def(py::init<double, double>());

    py::class_<astrometry::Orientation>(m, "Orientation")
        .def(py::init<>())
        .def(py::init<astrometry::SphericalCoords&, double>(), py::arg("pole_"), py::arg("pa")=0)
        .def("set", &astrometry::Orientation::set)
        .def("getPole", &astrometry::Orientation::getPole);

    py::class_<astrometry::Gnomonic, astrometry::SphericalCoords>(m, "Gnomonic")
        .def(py::init<astrometry::Orientation const&, bool>(), py::arg("o"), py::arg("shareOrient")=false)
        .def(py::init<astrometry::SphericalCoords const&, astrometry::Orientation&>())
        .def(py::init<>());

    py::class_<astrometry::LinearMap, astrometry::PixelMap>(m, "LinearMap")
        //.def("setVector", [](astrometry::LinearMap & self, astrometry::DVector::Base const& vector){self.v = vector;});
        .def(py::init<astrometry::DVector::Base const&>());

    py::class_<poly2d::Poly2d>(m, "Poly2d")
        .def(py::init<int>())
        .def(py::init<astrometry::DMatrix::Base const&>())
        .def("nCoeffs", &poly2d::Poly2d::nCoeffs)
        .def("getC", &poly2d::Poly2d::getC)
        .def("setC", &poly2d::Poly2d::setC)
        .def("vectorIndex", &poly2d::Poly2d::vectorIndex);

    py::class_<astrometry::PolyMap, astrometry::PixelMap>(m, "PolyMap")
        .def(py::init<poly2d::Poly2d, poly2d::Poly2d, string, Bounds<double>, double>());

    py::class_<astrometry::SubMap, astrometry::PixelMap>(m, "SubMap")
        .def(py::init<list<astrometry::PixelMap*> const&, string, bool>());

    py::class_<astrometry::Wcs, shared_ptr<astrometry::Wcs>>(m, "Wcs")
        .def(py::init<astrometry::PixelMap*, astrometry::SphericalCoords const&, string, double, bool>(),
            py::arg("pm_"), py::arg("nativeCoords_"), py::arg("name")="", py::arg("wScale_")=DEGREE,
            py::arg("shareMap_")=false)
        .def("reprojectTo", [](astrometry::Wcs & self, const SphericalCoords& targetCoords_) {
            py::scoped_estream_redirect stream(
                std::cerr,                               // std::ostream&
                py::module_::import("sys").attr("stderr") // Python output
            );
            self.reprojectTo(targetCoords_);
        })
        .def("getTargetCoords", &astrometry::Wcs::getTargetCoords)
        .def("toWorld", [](astrometry::Wcs & self, double x, double y) {
            double ra;
            double dec;
            self.toWorld(x, y, ra, dec);
            vector<double> radec = {ra, dec};
            return radec;
        });

    declareInstrumentVector(m);
    py::class_<Fields>(m, "Fields")
        .def(py::init<vector<string>, vector<double>, vector<double>, vector<double>>());

    py::class_<FieldsHelper>(m, "FieldsHelper")
        .def(py::init<vector<string>, vector<double>, vector<double>, vector<double>>());

    py::class_<Instrument, shared_ptr<Instrument>>(m, "Instrument")
        .def(py::init<string>(), py::arg("name_")="")
        .def_readwrite("name", &Instrument::name)
        //.def_readwrite("nDevices", &Instrument::nDevices)
        .def_readwrite("band", &Instrument::band)
        //.def_readwrite("deviceNames", &Instrument::deviceNames)
        .def("addDevice", &Instrument::addDevice)
        .def("hasDevice", [](Instrument & self, string name) {
                          return self.deviceNames.has(name);
        });

    py::class_<ExposuresHelper>(m, "ExposuresHelper")
        .def(py::init<vector<string>, vector<int>, vector<int>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>>());

    declareBounds<double>(m);

    py::class_<YAMLCollector>(m, "YAMLCollector")
        .def(py::init<string, string>())
        .def("addInput", [](YAMLCollector & self, string is, string filter, bool prepend){
                        istringstream iss(is);
                        self.addInput(iss, filter, prepend);},
                        py::arg("is"), py::arg("filter")="", py::arg("prepend")=false);
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

    py::class_<ExtensionObjectSet>(m, "ExtensionObjectSet")
        .def(py::init<string>());

    //declareExtension<astrometry::SubMap, astrometry::Detection>(m);

    py::class_<FitClass>(m, "FitClass")
        .def(py::init<Fields &, vector<shared_ptr<Instrument>>, ExposuresHelper, vector<int>, vector<int>, 
                      YAMLCollector, vector<shared_ptr<astrometry::Wcs>>, vector<int>, vector<LONGLONG>, vector<LONGLONG>, double, double, int, ExtensionObjectSet, string, bool, int>(),
             py::arg("fields_"), py::arg("instruments_"), py::arg("exposures_"),
             py::arg("extensionExposureNumbers"), py::arg("extensionDevices"), py::arg("inputYAML"),
             py::arg("wcss"), py::arg("sequence"), py::arg("extns"), py::arg("objects"),
             py::arg("sysErr")=2.0, py::arg("refSysErr")=2.0, py::arg("minMatches")=2, py::arg("matchSkipSet")=ExtensionObjectSet(""),
             py::arg("fixMaps")="", py::arg("usePM")=true, py::arg("verbose")=0)
        .def("setObjects", &FitClass::setObjects) 
        //.def_readwrite("reserveFraction", &FitClass::reserveFraction)
        .def_readwrite("verbose", &FitClass::verbose)
        .def("setRefWCSNames", &FitClass::setRefWCSNames) // to delete
        .def("setupMaps", &FitClass::setupMaps) // to delete
        .def("setMatches", &FitClass::setMatches) // to delete
        .def("setExtensions", &FitClass::setExtensions) // to delete
        .def("fit", &FitClass::fit, py::arg("maxError")=100., py::arg("minFitExposures")=200,
             py::arg("reserveFraction")=0.2, py::arg("randomNumberSeed")=1234,
             py::arg("minimumImprovement")=0.02, py::arg("clipThresh")=5.0, py::arg("chisqTolerance")=0.001,
             py::arg("clipEntireMatch")=false, py::arg("divideInPlace")=false, py::arg("purgeOutput")=false,
             py::arg("minColor")=-10.0, py::arg("maxColor")=10.0)
        ;
        //.def(py::init<string>())
        /*
        .def_readwrite("fields", &FitClass::fields) 
        .def_readwrite("instruments", &FitClass::instruments) 
        .def("setExposures", [](FitClass & self, vector<shared_ptr<Exposure>> expos, double sysErr, double refSysErr) {
            py::scoped_estream_redirect stream(
                std::cerr,                               // std::ostream&
                py::module_::import("sys").attr("stderr") // Python output
            );
            self.setExposures(expos, sysErr, refSysErr);
        })
        
        .def("addMap", &FitClass::addMap)
        
        
        
        
        .def("getMatchLength", &FitClass::getMatchLength)
        //.def("fit", [](FitClass & self) {
        //    py::scoped_estream_redirect stream(
        //        std::cerr,                               // std::ostream&
        //        py::module_::import("sys").attr("stdout") // Python output
        //    );
        //    self.fit();
        //})
        
        .def_readwrite("randomNumberSeed", &FitClass::randomNumberSeed)
        .def_readwrite("clipThresh", &FitClass::clipThresh)
        .def_readwrite("maxError", &FitClass::maxError)
        .def_readwrite("minMatches", &FitClass::minMatches)
        .def_readwrite("minFitExposures", &FitClass::minFitExposures)
        .def_readwrite("clipEntireMatch", &FitClass::clipEntireMatch)
        .def_readwrite("chisqTolerance", &FitClass::chisqTolerance)
        .def_readwrite("divideInPlace", &FitClass::divideInPlace)
        .def_readwrite("purgeOutput", &FitClass::purgeOutput)
        .def_readwrite("minColor", &FitClass::minColor)
        .def_readwrite("maxColor", &FitClass::maxColor)
        .def_readwrite("exposures", &FitClass::exposures)
        
        .def_readwrite("extensions", &FitClass::extensions)
        .def_readwrite("fieldNames", &FitClass::fieldNames)
        .def_readwrite("fieldProjections", &FitClass::fieldProjections)
        .def_readwrite("fieldEpochs", &FitClass::fieldEpochs)
        .def_readwrite("matches", &FitClass::matches);
        //.def_readwrite("inputYAML", &FitClass::inputYAML)*/

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
    
    
    py::class_<Exposure, shared_ptr<Exposure>>(m, "Exposure")
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
        .def("setAstrometricCovariance", [](Exposure & self, astrometry::Matrix22::Base const& matrix){self.astrometricCovariance = matrix;})
        .def("getAstrometricCovariance", [](Exposure const& self){return static_cast <astrometry::Matrix22::Base const&>(self.astrometricCovariance);})
        .def("addToAstrometricCovariance", [](Exposure & self, astrometry::Matrix22::Base const& matrix){self.astrometricCovariance += matrix;});

    py::class_<SPExtension, shared_ptr<SPExtension>>(m, "SPExtension")
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

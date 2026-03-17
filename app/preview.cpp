#include "preview_common.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>

namespace {

void printUsage(const char* exe)
{
    std::cerr
        << "Usage: " << exe << " --in <out_dir> [options]\n"
        << "\n"
        << "Options:\n"
        << "  --out <file.png|file.ppm>    Output image path (default: preview.png)\n"
        << "  --width <int>                Image width  (default: 1024)\n"
        << "  --height <int>               Image height (default: 512)\n"
        << "  --sun-zenith-deg <float>     Sun zenith angle in degrees (default: 30)\n"
        << "  --sun-azimuth-deg <float>    Sun azimuth angle in degrees (default: 0)\n"
        << "  --exposure <float>           Exposure multiplier before tonemap (default: 10)\n"
        << "  --gamma <float>              Output gamma (default: 2.2)\n"
        << "  --include-rayleigh <0|1>     Include single Rayleigh (default: 1)\n"
        << "  --include-mie <0|1>          Include single Mie (default: 1)\n"
        << "  --include-multiple <0|1>     Include multiple scattering (default: 1)\n"
        << "  --black-below-horizon <0|1>  Output black for mu <= 0 (default: 1)\n"
        << "\n"
        << "Optional consistency checks / overrides:\n"
        << "  --sky-mu <int>               Must match LUT header if specified\n"
        << "  --sky-mu-s <int>             Must match LUT header if specified\n"
        << "  --sky-nu <int>               Must match LUT header if specified\n"
        << "  --mie-g <float>\n"
        << "  --mu-s-min <float>\n"
        << "  --observer-altitude-m <float>\n";
}

bool parseBool01(const std::string& s, bool& out)
{
    if(s == "0") {
        out = false;
        return true;
    }
    if(s == "1") {
        out = true;
        return true;
    }
    return false;
}

} // namespace

int main(int argc, char** argv)
{
    if(argc <= 1) {
        printUsage(argv[0]);
        return 1;
    }

    std::filesystem::path inDir;
    std::filesystem::path outPath = "preview.png";

    preview::RenderOptions options;
    std::optional<std::uint32_t> checkSkyMu;
    std::optional<std::uint32_t> checkSkyMuS;
    std::optional<std::uint32_t> checkSkyNu;
    std::optional<float> overrideMieG;
    std::optional<float> overrideMuSMin;
    std::optional<float> overrideObserverAltitudeM;

    for(int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto requireValue = [&](const char* name) -> const char* {
            if(i + 1 >= argc) {
                std::cerr << "Missing value after " << name << "\n";
                std::exit(1);
            }
            return argv[++i];
        };

        if(arg == "--in") {
            inDir = requireValue("--in");
        }
        else if(arg == "--out") {
            outPath = requireValue("--out");
        }
        else if(arg == "--width") {
            options.width = static_cast<std::uint32_t>(std::stoul(requireValue("--width")));
        }
        else if(arg == "--height") {
            options.height = static_cast<std::uint32_t>(std::stoul(requireValue("--height")));
        }
        else if(arg == "--sun-zenith-deg") {
            options.sunZenithDeg = std::stof(requireValue("--sun-zenith-deg"));
        }
        else if(arg == "--sun-azimuth-deg") {
            options.sunAzimuthDeg = std::stof(requireValue("--sun-azimuth-deg"));
        }
        else if(arg == "--exposure") {
            options.exposure = std::stof(requireValue("--exposure"));
        }
        else if(arg == "--gamma") {
            options.gamma = std::stof(requireValue("--gamma"));
        }
        else if(arg == "--sky-mu") {
            checkSkyMu = static_cast<std::uint32_t>(std::stoul(requireValue("--sky-mu")));
        }
        else if(arg == "--sky-mu-s") {
            checkSkyMuS = static_cast<std::uint32_t>(std::stoul(requireValue("--sky-mu-s")));
        }
        else if(arg == "--sky-nu") {
            checkSkyNu = static_cast<std::uint32_t>(std::stoul(requireValue("--sky-nu")));
        }
        else if(arg == "--mie-g") {
            overrideMieG = std::stof(requireValue("--mie-g"));
        }
        else if(arg == "--mu-s-min") {
            overrideMuSMin = std::stof(requireValue("--mu-s-min"));
        }
        else if(arg == "--observer-altitude-m") {
            overrideObserverAltitudeM = std::stof(requireValue("--observer-altitude-m"));
        }
        else if(arg == "--include-rayleigh") {
            if(!parseBool01(requireValue("--include-rayleigh"), options.includeRayleigh)) {
                std::cerr << "--include-rayleigh must be 0 or 1\n";
                return 1;
            }
        }
        else if(arg == "--include-mie") {
            if(!parseBool01(requireValue("--include-mie"), options.includeMie)) {
                std::cerr << "--include-mie must be 0 or 1\n";
                return 1;
            }
        }
        else if(arg == "--include-multiple") {
            if(!parseBool01(requireValue("--include-multiple"), options.includeMultiple)) {
                std::cerr << "--include-multiple must be 0 or 1\n";
                return 1;
            }
        }
        else if(arg == "--black-below-horizon") {
            if(!parseBool01(requireValue("--black-below-horizon"), options.blackBelowHorizon)) {
                std::cerr << "--black-below-horizon must be 0 or 1\n";
                return 1;
            }
        }
        else if(arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        }
        else {
            std::cerr << "Unknown argument: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    if(inDir.empty()) {
        std::cerr << "--in is required.\n";
        return 1;
    }

    std::string error;
    preview::PreviewMetadata meta;
    if(!preview::loadMetadataJson(inDir / "metadata.json", meta, error)) {
        std::cerr << error << "\n";
        return 1;
    }

    if(overrideMieG.has_value()) {
        meta.miePhaseFunctionG = *overrideMieG;
    }
    if(overrideMuSMin.has_value()) {
        meta.muSMin = *overrideMuSMin;
    }
    if(overrideObserverAltitudeM.has_value()) {
        meta.observerAltitudeM = *overrideObserverAltitudeM;
    }
    if(checkSkyMu.has_value()) {
        meta.skyMu = *checkSkyMu;
    }
    if(checkSkyMuS.has_value()) {
        meta.skyMuS = *checkSkyMuS;
    }
    if(checkSkyNu.has_value()) {
        meta.skyNu = *checkSkyNu;
    }

    preview::FinalLuts luts;
    if(!preview::loadFinalLuts(inDir, meta, luts, error)) {
        std::cerr << error << "\n";
        return 1;
    }

    preview::ImageRgb32f image;
    preview::renderEquirectangularSky(meta, luts, options, image);

    if(!preview::saveImage8(outPath, image, options.exposure, options.gamma, error)) {
        std::cerr << error << "\n";
        return 1;
    }


    std::cout << "Preview image written to: " << outPath.string() << "\n";
    std::cout << "Sky dims  : nu=" << meta.skyNu
              << " mu=" << meta.skyMu
              << " mu_s=" << meta.skyMuS
              << " lambda=" << meta.wavelengthsNm.size() << "\n";
    std::cout << "Table type: sky headers are interpreted as [nu][mu][mu_s][lambda], lambda fastest.\n";
    return 0;
}

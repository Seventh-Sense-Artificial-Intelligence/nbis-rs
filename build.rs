use std::env;

fn android_abi_from_target(target: &str) -> Option<&'static str> {
    if target.contains("aarch64") {
        Some("arm64-v8a")
    } else if target.contains("armv7") {
        Some("armeabi-v7a")
    } else if target.contains("x86_64") {
        Some("x86_64")
    } else if target.contains("i686") {
        Some("x86")
    } else {
        None
    }
}

// build.rs
fn main() {
    println!("cargo:rerun-if-env-changed=CLIPPY");

    let target = env::var("TARGET").unwrap_or_default();
    let is_android = target.contains("android");

    // ---- CMake for OpenCV ----
    let mut cmake = cmake::Config::new("ext/opencv-2.4.13.6");

    if is_android {
        let ndk = env::var("ANDROID_NDK_HOME").expect("ANDROID_NDK_HOME not set");
        let abi = if target.contains("aarch64") {
            "arm64-v8a"
        } else if target.contains("armv7") {
            "armeabi-v7a"
        } else {
            panic!("Unsupported Android ABI: {}", target);
        };

        cmake
            .define("CMAKE_SYSTEM_NAME", "Android")
            .define("CMAKE_SYSTEM_VERSION", "21") // minSdkVersion
            .define("CMAKE_ANDROID_ARCH_ABI", abi)
            .define("CMAKE_ANDROID_NDK", &ndk)
            .define("ANDROID_NATIVE_API_LEVEL", "21")
            .define("ANDROID_ABI", abi)
            .define("ANDROID_STL", "c++_static")
            .define("CMAKE_CXX_STANDARD", "14")
            .define("CMAKE_TOOLCHAIN_FILE", format!("{}/build/cmake/android.toolchain.cmake", ndk));
    }

    cc::Build::new()
        .file("ext/nbis/nfiq/src/lib/nfiq/nfiq.c")
        .file("ext/nbis/nfiq/src/lib/nfiq/nfiqgbls.c")
        .file("ext/nbis/nfiq/src/lib/nfiq/nfiqread.c")
        .file("ext/nbis/nfiq/src/lib/nfiq/znorm.c")
        .file("ext/nbis/commonbis/src/lib/util/memalloc.c")
        .file("ext/nbis/commonbis/src/lib/util/ssxstats.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/dataio.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/fileexst.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/filehead.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/fileroot.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/filesize.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/filetail.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/findfile.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/newext.c")
        .file("ext/nbis/commonbis/src/lib/ioutil/readutil.c")
        .file("ext/nbis/commonbis/src/lib/cblas/sgemv.c")
        .file("ext/nbis/commonbis/src/lib/cblas/xerbla.c")
        .file("ext/nbis/commonbis/src/lib/cblas/lsame.c")
        .file("ext/nbis/pcasys/src/lib/mlp/runmlp.c")
        .file("ext/nbis/pcasys/src/lib/mlp/acs.c")
        .file("ext/nbis/pcasys/src/lib/mlp/mlpcla.c")
        .include("ext/nbis/pcasys/include")
        .include("ext/nbis/nfiq/include")
        .include("ext/nbis/mindtct/include")
        .include("ext/nbis/commonbis/include")
        .include("ext/nbis/imgtools/include")
        .define("NOVERBOSE", None) // you probably don’t want stdout spam
        .flag_if_supported("-w") // for GCC/Clang: suppress *all* warnings
        .compile("nfiq");

    cc::Build::new()
        .file("ext/nbis/bozorth/src/lib/bozorth3/bozorth3.c")
        .file("ext/nbis/bozorth/src/lib/bozorth3/bz_alloc.c")
        .file("ext/nbis/bozorth/src/lib/bozorth3/bz_drvrs.c")
        .file("ext/nbis/bozorth/src/lib/bozorth3/bz_gbls.c")
        .file("ext/nbis/bozorth/src/lib/bozorth3/bz_io.c")
        .file("ext/nbis/bozorth/src/lib/bozorth3/bz_sort.c")
        .include("ext/nbis/commonbis/include")
        .file("ext/nbis/bozorth/src/lib/bozorth3/bozorth_glue.c")
        .include("ext/nbis/bozorth/include") // to find bozorth.h
        .define("NOVERBOSE", None) // you probably don’t want stdout spam
        .flag_if_supported("-w") // for GCC/Clang: suppress *all* warnings
        .compile("bozorth");

    cc::Build::new()
        .file("ext/nbis/mindtct/src/lib/mindtct/log.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/line.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/contour.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/imgutil.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/quality.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/block.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/loop.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/mytime.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/minutia.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/link.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/matchpat.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/binar.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/morph.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/chaincod.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/detect.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/dft.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/free.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/globals.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/init.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/isempty.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/remove.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/ridges.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/shape.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/sort.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/util.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/maps.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/xytreps.c")
        .file("ext/nbis/mindtct/src/lib/mindtct/getmin.c")
        .include("ext/nbis/mindtct/include") // to find bozorth.h
        .define("NOVERBOSE", None) // you probably don’t want stdout spam
        .flag_if_supported("-w") // for GCC/Clang: suppress *all* warnings
        .compile("mindtct");

    let dst = cmake
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("BUILD_PNG", "OFF")
        .define("BUILD_JPEG", "OFF")
        .define("BUILD_TIFF", "OFF")
        .define("BUILD_WEBP", "OFF")
        .define("BUILD_OPENJPEG", "OFF")
        .define("WITH_FFMPEG", "OFF")
        .define("BUILD_opencv_videoio", "OFF")
        .define(
            "BUILD_LIST",
            "core,imgproc,video,features2d,flann,calib3d,objdetect,legacy,highgui",
        ) // only build needed modules
        .define("BUILD_EXAMPLES", "OFF")
        .define("BUILD_TESTS", "OFF")
        .define("BUILD_ZLIB", "OFF")
        .define("BUILD_PERF_TESTS", "OFF")
        .build();

    if is_android {
        use std::fs;

        let abi = android_abi_from_target(&target)
            .expect("Unsupported Android target");

        let out_dir = dst.display().to_string();
        let opencv_lib_dir = dst.join("build").join("lib").join(abi);

        // List of expected OpenCV static libs
        for lib in &["opencv_core", "opencv_imgproc", "opencv_highgui", "opencv_legacy"] {
            let src = opencv_lib_dir.join(format!("lib{}.a", lib));
            let dst = std::path::Path::new(&out_dir).join(format!("lib{}.a", lib));
            if src.exists() {
                fs::copy(&src, &dst).unwrap_or_else(|_| panic!("Failed to copy {}", lib));
            } else {
                panic!("Missing static lib: {}", src.display());
            }
        }

        // Set Android-specific link search path
        println!("cargo:rustc-link-search=native={}", out_dir);
    } else {
        // macOS path
        println!("cargo:rustc-link-search=native={}/lib", dst.display());
    }



    println!("cargo:rustc-link-lib=static=opencv_core");
    println!("cargo:rustc-link-lib=static=opencv_imgproc");
    println!("cargo:rustc-link-lib=static=opencv_highgui");
    println!("cargo:rustc-link-lib=static=opencv_legacy"); // needed by highgui in 2.4.x
    println!("cargo:rustc-link-lib=z");

    cc::Build::new()
        .cpp(true)
        .file("ext/nbis/misc/sivv/src/SIVVCore.cpp")
        .file("ext/nbis/misc/sivv/src/SIVVGraph.cpp")
        .file("ext/nbis/misc/sivv/src/sivv_wrapper.cpp")
        .include("ext/nbis/misc/sivv/include")
        .include("ext/opencv-2.4.13.6/include")
        .include("ext/opencv-2.4.13.6/include/opencv")
        .include("ext/opencv-2.4.13.6/modules/core/include")
        .include("ext/opencv-2.4.13.6/modules/imgproc/include")
        .include("ext/opencv-2.4.13.6/modules/video/include")
        .include("ext/opencv-2.4.13.6/modules/features2d/include")
        .include("ext/opencv-2.4.13.6/modules/flann/include")
        .include("ext/opencv-2.4.13.6/modules/calib3d/include")
        .include("ext/opencv-2.4.13.6/modules/objdetect/include")
        .include("ext/opencv-2.4.13.6/modules/legacy/include")
        .include("ext/opencv-2.4.13.6/modules/highgui/include")
        .define("NOVERBOSE", None) // you probably don’t want stdout spam
        .flag_if_supported("-w") // for GCC/Clang: suppress *all* warnings
        .flag_if_supported("-Wno-everything") // extra if using
        .compile("sivv");

    // Automatically re-run build.rs if these files change
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/src/lib/bozorth3/bozorth3.c");
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/include");
}

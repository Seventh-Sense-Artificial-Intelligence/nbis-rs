use std::{env, path::Path};

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
    let is_linux = target.contains("linux") && !target.contains("android");

    // ---- CMake for OpenCV ----
    let mut cmake = cmake::Config::new("ext/opencv-4.10.0");

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
            .define("INSTALL_CREATE_DISTRIB", "ON")
            .define("CMAKE_INSTALL_PREFIX", "opencv_install")
            .define("CMAKE_INSTALL_INCLUDEDIR", "opencv_install/sdk/native/jni/include")
            .define("CMAKE_INSTALL_LIBDIR", "opencv_install/sdk/native/libs")
            .define("BUILD_ANDROID_PROJECTS", "OFF")
            .define("BUILD_ANDROID_EXAMPLES", "OFF")
            .define("BUILD_opencv_java", "OFF")
            .build_target("install")
            .define(
                "CMAKE_TOOLCHAIN_FILE",
                format!("{}/build/cmake/android.toolchain.cmake", ndk),
            );
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
        // Disbale image codecs we don't need
        .define("BUILD_PNG", "OFF")
        .define("BUILD_JPEG", "OFF")
        .define("BUILD_TIFF", "OFF")
        .define("BUILD_WEBP", "OFF")
        .define("BUILD_OPENJPEG", "OFF")
        // For Mac
        .define("WITH_TEGRA", "OFF")
        .define("WITH_CAROTENE", "OFF") // ← stop building the carotene_o4t HAL
        .define("WITH_LAPACK", "OFF")
        .define("WITH_OPENCL", "OFF")
        // disable unnecessary modules
        .define("WITH_FFMPEG", "OFF")
        .define("WITH_CUDA", "OFF")
        .define("WITH_CUDNN", "OFF")
        .define("WITH_GSTREAMER", "OFF")
        .define("WITH_V4L", "OFF")
        .define("WITH_V4L2", "OFF")
        .define("WITH_LIBV4L", "OFF")
        .define("WITH_IPP", "OFF")
        .define("BUILD_IPP_IW", "OFF")
        .define("WITH_ITT", "OFF")
        .define("BUILD_opencv_hal", "ON")
        .define("BUILD_opencv_python2", "OFF")
        .define("BUILD_opencv_python3", "OFF")
        .define("BUILD_opencv_python_bindings_generator", "OFF")
        .define("BUILD_opencv_videoio", "OFF")
        .define(
            "BUILD_LIST",
            "core,imgproc",
        ) // only build needed modules
        .define("BUILD_EXAMPLES", "OFF")
        .define("BUILD_TESTS", "OFF")
        .define("BUILD_ZLIB", "OFF")
        .define("BUILD_PERF_TESTS", "OFF")
        .define("CMAKE_CXX_STANDARD", "14")
        .build();

    // Print dst as a warning for the user
    //eprintln!("OpenCV build directory: {}", dst.display());

    cc::Build::new().cpp(true)
        .file("ext/nbis/misc/sivv/src/SIVVCore.cpp")
        .file("ext/nbis/misc/sivv/src/SIVVGraph.cpp")
        .file("ext/nbis/misc/sivv/src/sivv_wrapper.cpp")
        .include("ext/nbis/misc/sivv/include")
        .include(dst.join("include/opencv4"))
        // Additional includes for Android
        .include(dst.join("build/opencv_install/sdk/native/jni/include"))
        .define("NOVERBOSE", None) // you probably don’t want stdout spam
        .flag_if_supported("-w") // for GCC/Clang: suppress *all* warnings
        .flag_if_supported("-Wno-everything") // extra if using
        .compile("sivv");

    if is_android || is_linux {
        use std::fs;

        let out_dir = dst.display().to_string();
        let opencv_lib_dir = if is_android {
            let abi = android_abi_from_target(&target).expect("Unsupported Android target");
            dst.join("build").join("lib").join(abi)
        } else {
            dst.join("build").join("lib")
        };

        // List of expected OpenCV static libs
        for entry in fs::read_dir(&opencv_lib_dir).expect("Failed to read OpenCV lib dir") {
            let entry = entry.expect("Failed to read dir entry");
            let path = entry.path();

            if path.extension().is_some_and(|ext| ext == "a") {
                let filename = path.file_name().unwrap();
                let dst_path = Path::new(&out_dir).join(filename);
                fs::copy(&path, &dst_path)
                    .unwrap_or_else(|_| panic!("Failed to copy static lib {}", path.display()));
            }
        }

        // Set Android-specific link search path
        println!("cargo:rustc-link-search=native={}", out_dir);
    } else {
        // macOS path
        println!("cargo:rustc-link-search=native={}/lib", dst.display());
    }

    println!("cargo:rustc-link-lib=static=opencv_imgproc");
    println!("cargo:rustc-link-lib=static=opencv_core");
    println!("cargo:rustc-link-lib=z");

    // Automatically re-run build.rs if these files change
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/src/lib/bozorth3/bozorth3.c");
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/include");
}

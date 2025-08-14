use std::{
    env, fs, io,
    path::{Path, PathBuf},
};

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

fn copy_nfiq2_dirs() {
    let source = Path::new("ext/opencv-4.10.0");
    let dst = Path::new("ext/NFIQ2-2.3.0/opencv");
    copy_dir_recursive(source, dst).expect("failed to copy OpenCV dir");

    let source = Path::new("ext/FingerJetFXOSE");
    let dst = Path::new("ext/NFIQ2-2.3.0/fingerjetfxose");
    copy_dir_recursive(source, dst).expect("failed to copy FingerJetFXOSE dir");

    let source = Path::new("ext/digestpp");
    let dst = Path::new("ext/NFIQ2-2.3.0/digestpp");
    copy_dir_recursive(source, dst).expect("failed to copy digestpp dir");

    let source = Path::new("ext/libbiomeval-10.0");
    let dst = Path::new("ext/NFIQ2-2.3.0/libbiomeval");
    copy_dir_recursive(source, dst).expect("failed to copy libbiomeval-10.0 dir");
}

fn build_nfiq2() -> PathBuf {
    let target = env::var("TARGET").unwrap_or_default();
    let is_android = target.contains("android");
    let is_linux = target.contains("linux") && !target.contains("android");
    //let is_windows = target.contains("windows");
    let is_macos = target.contains("apple") || target.contains("darwin");
    copy_nfiq2_dirs();

    // ---- CMake for NFIQ2 ----
    let mut cmake = cmake::Config::new("ext/NFIQ2-2.3.0");
    cmake
        .define("CMAKE_BUILD_TYPE", "Release")
        .define("CMAKE_INSTALL_PREFIX", "NFIQ2-2.3.0/install")
        .define("EMBED_RANDOM_FOREST_PARAMETERS", "ON")
        .define("EMBEDDED_RANDOM_FOREST_PARAMETER_FCT", "3")
        .define("BUILD_NFIQ2_CLI", "OFF");

    if is_android {
        let ndk = env::var("ANDROID_NDK_ROOT").expect("ANDROID_NDK_ROOT not set");
        let abi = android_abi_from_target(&target)
            .expect("Unsupported Android ABI. Supported ABIs: arm64-v8a, armeabi-v7a, x86_64, x86");
        cmake.define("ANDROID_ABI", abi);
        cmake.define(
            "CMAKE_TOOLCHAIN_FILE",
            format!("{ndk}/build/cmake/android.toolchain.cmake"),
        );
    }

    let dst = cmake.build();

    // Define the include and library paths for NFIQ2
    let nfiq2_include_path = dst.join("build/install_staging/nfiq2/include");
    let nfiq2_lib_path = dst.join("build/install_staging/nfiq2/lib");

    // On Android, OpenCV libraries are in a different location
    let opencv_android_lib_path = if is_android {
        let abi = android_abi_from_target(&target).expect("Unsupported Android ABI");
        Some(dst.join(format!(
            "build/install_staging/nfiq2/sdk/native/staticlibs/{abi}",
        )))
    } else {
        None
    };

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let opencv_lib_path = out_dir.join("build/install_staging/nfiq2/lib");

    // 1) Compile the C++ FFI wrapper
    cc::Build::new()
        .cpp(true) // switch to a C++ compiler
        .flag_if_supported("-std=c++14") // or c++11/17, whichever you need
        .include(&nfiq2_include_path) // where nfiq2.hpp lives
        .include("src/cwrapper")
        .file("src/cwrapper/nfiq_wrapper.cpp") // your FFI source
        .define("NOVERBOSE", None) // you probably don’t want stdout spam
        .flag_if_supported("-w") // for GCC/Clang: suppress *all* warnings
        .compile("nfiq2_ffi"); // emits libnfiq2_ffi.a

    // 2) Link against both the wrapper and the NFIQ2 / OpenCV libs
    println!("cargo:rustc-link-lib=static=nfiq2_ffi");
    println!(
        "cargo:rustc-link-search=native={}",
        nfiq2_lib_path.display()
    );
    println!(
        "cargo:rustc-link-search=native={}",
        opencv_lib_path.display()
    );

    // Include the OpenCV libs path for Android
    if let Some(opencv_android_lib_path) = opencv_android_lib_path {
        println!(
            "cargo:rustc-link-search=native={}",
            opencv_android_lib_path.display()
        );
    }

    println!("cargo:rustc-link-lib=static=nfiq2");

    if is_linux {
        println!("cargo:rustc-link-lib=dylib=stdc++");
        println!("cargo:rustc-link-lib=dylib=z");
    }

    if is_android {
        println!("cargo:rustc-link-lib=z");
        println!("cargo:rustc-link-lib=android");
        println!("cargo:rustc-link-lib=c++_shared");
    }

    if is_macos {
        println!("cargo:rustc-link-lib=framework=Accelerate");
        println!("cargo:rustc-link-lib=framework=OpenCL");
        //println!("cargo:rustc-link-lib=static=zlib");
    }

    dst
}

// build.rs
fn main() {
    println!("cargo:rerun-if-env-changed=CLIPPY");

    let target = env::var("TARGET").unwrap_or_default();
    let is_android = target.contains("android");
    let is_linux = target.contains("linux") && !target.contains("android");
    let is_windows = target.contains("windows");

    let dst = build_nfiq2();
    // dst: /home/coder/nbis-rs/target/release/build/nbis-rs-4685d910ce23e274/out

    let mut bozorth_cc = cc::Build::new();
    bozorth_cc
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
        ;

    if is_windows {
        bozorth_cc
            .file("ext/sys_time/time.cpp")
            .include("ext/sys_time");
    }

    bozorth_cc.compile("bozorth");

    let mut mindtct_cc = cc::Build::new();
    mindtct_cc
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
        ;

    if is_windows {
        mindtct_cc
            .file("ext/sys_time/time.cpp")
            .include("ext/sys_time");
    }

    mindtct_cc.compile("mindtct");

    let mut sivv_cpp = cc::Build::new();
    sivv_cpp
        .cpp(true)
        .flag("-std=c++11")
        .file("ext/nbis/misc/sivv/src/SIVVCore.cpp")
        .file("ext/nbis/misc/sivv/src/sivv_wrapper.cpp")
        .include("ext/nbis/misc/sivv/include")
        .include(dst.join("build/install_staging/nfiq2/include"))
        // Additional includes for Android
        .include(dst.join("build/opencv_install/sdk/native/jni/include"))
        .define("NOVERBOSE", None) // you probably don’t want stdout spam
        .flag_if_supported("-w") // for GCC/Clang: suppress *all* warnings
        .flag_if_supported("-Wno-everything") // extra if using
        ;

    if is_windows {
        sivv_cpp
            .file("ext/sys_time/time.cpp")
            .include("ext/sys_time");
    }

    if is_android {
        // dst: /home/coder/nbis-rs/target/aarch64-linux-android/release/build/nbis-rs-0311d84ca63bc87e/out
        let ocv_header_path = dst.join("build/install_staging/nfiq2/sdk/native/jni/include");
        sivv_cpp.include(ocv_header_path);
    }

    sivv_cpp.compile("sivv");

    if is_android || is_linux {
        let opencv_lib_dir = if is_android {
            let abi = android_abi_from_target(&target).expect("Unsupported Android target");
            dst.join("build/install_staging/nfiq2/lib").join(abi)
        } else {
            dst.join("build/install_staging/nfiq2/lib")
        };

        // Set lib path
        println!(
            "cargo:rustc-link-search=native={}",
            opencv_lib_dir.display()
        );
    } else {
        // macOS path
        println!("cargo:rustc-link-search=native={}/lib", dst.display());
    }

    if !is_windows {
        println!("cargo:rustc-link-lib=static=opencv_imgproc");
        println!("cargo:rustc-link-lib=static=opencv_ml");
        println!("cargo:rustc-link-lib=static=opencv_imgcodecs");
        println!("cargo:rustc-link-lib=static=opencv_imgproc");
        println!("cargo:rustc-link-lib=static=opencv_core");
        println!("cargo:rustc-link-lib=static=FRFXLL_static");
    } else {
        // dst.join("build/install_staging/nfiq2/lib");
        // let lib_src_dir_str = format!("{}/build/lib", &dst.display());
        let lib_src_dir_str = format!("{}/build/install_staging/nfiq2", &dst.display());
        let lib_src_dir = Path::new(&lib_src_dir_str);
        let lib_dst_dir = Path::new("ext/nfiq2_libs");

        // Copy the directory
        if let Err(e) = copy_dir_recursive(lib_src_dir, lib_dst_dir) {
            panic!("Failed to copy directory: {e}");
        }

        println!("cargo:rustc-link-search=native=C:/msys64/mingw64/lib");
        println!("cargo:rustc-link-search=native={}/lib", &lib_src_dir_str);
        println!("cargo:rustc-link-search=native={}/x64/mingw/staticlib", &lib_src_dir_str);
        println!("cargo:rustc-link-lib=static=opencv_imgproc4100");
        println!("cargo:rustc-link-lib=static=opencv_ml4100");
        println!("cargo:rustc-link-lib=static=opencv_imgcodecs4100");
        println!("cargo:rustc-link-lib=static=opencv_core4100");
        println!("cargo:rustc-link-lib=static=FRFXLL_static");
        println!("cargo:rustc-link-lib=static=openblas");
        println!("cargo:rustc-link-lib=static=gomp");
        println!("cargo:rustc-link-lib=static=stdc++");
    }

    println!("cargo:rustc-link-lib=z");

    // Automatically re-run build.rs if these files change
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/src/lib/bozorth3/bozorth3.c");
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/include");
}

fn copy_dir_recursive(src: &Path, dst: &Path) -> io::Result<()> {
    // If dst exists, delete it
    if dst.exists() {
        fs::remove_dir_all(dst)?;
    }

    // Create the destination directory if it doesn't exist
    if !dst.exists() {
        fs::create_dir_all(dst)?;
    }

    // Read the source directory
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let src_path = entry.path();
        let dst_path = dst.join(entry.file_name());

        if src_path.is_dir() {
            // Recursively copy subdirectories
            copy_dir_recursive(&src_path, &dst_path)?;
        } else {
            // Copy files
            fs::copy(&src_path, &dst_path)?;
        }
    }

    Ok(())
}

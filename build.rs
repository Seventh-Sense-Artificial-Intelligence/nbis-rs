// build.rs
fn main() {
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

    // Automatically re-run build.rs if these files change
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/src/lib/bozorth3/bozorth3.c");
    println!("cargo:rerun-if-changed=ext/nbis/bozorth/include");
}

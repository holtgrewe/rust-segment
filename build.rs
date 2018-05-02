extern crate bindgen;
extern crate fs_utils;

use fs_utils::copy::copy_directory;

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    if !out.join("haarseg").exists() {
        copy_directory("haarseg", &out).unwrap();
    }

    if Command::new("make")
        .current_dir(out.join("haarseg"))
        .arg("CFLAGS=-g -Wall -O2 -fPIC")
        .arg("lib-static")
        .arg("-B")
        .status()
        .unwrap()
        .success() != true
    {
        panic!("failed to build libhaarseg");
    }

    let bindings = bindgen::Builder::default()
        // https://github.com/servo/rust-bindgen/issues/687
        .blacklist_type("FP_NAN")
        .blacklist_type("FP_INFINITE")
        .blacklist_type("FP_ZERO")
        .blacklist_type("FP_SUBNORMAL")
        .blacklist_type("FP_NORMAL")
        // https://github.com/servo/rust-bindgen/issues/550
        .blacklist_type("max_align_t")
        .header("wrapper.h")
        .generate()
        .expect("Unable to generate bindings.");
    bindings
        .write_to_file(out.join("bindings.rs"))
        .expect("Could not write bindings.");

    println!(
        "cargo:rustc-flags=-L {out}/haarseg -L {out} -l static=haarseg",
        out = out.to_str().unwrap()
    );
}

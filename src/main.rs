use clap::{App, Arg, ArgMatches};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};

#[allow(non_snake_case)]
#[derive(Debug)]
struct Datf {
    A: u32,
    T: u32,
    G: u32,
    C: u32,
    N: u32,
    GC: f64,
    LEN: u32,
}

fn cmd() -> ArgMatches<'static> {
    let opt = App::new("Count_ATGC")
        .version("0.1.0")
        .author("sharkLoc")
        .about("a program to count ATGCN from fasta file.")
        .arg(
            Arg::with_name("fasta")
                .short("f")
                .long("fasta")
                .takes_value(true)
                .required(true)
                .help("path to fasta file name"),
        )
        .arg(
            Arg::with_name("out")
                .short("o")
                .long("out")
                .takes_value(true)
                .default_value("summary.txt")
                .help("output file name"),
        )
        .get_matches();

    opt
}

fn main() {
    let args = cmd();
    let file = args.value_of("fasta").expect("fasta file arg error");
    let out = args.value_of("out").expect("output file arg error");

    let mut fp = File::open(file).expect("open failed!");
    let mut op = File::create(out).expect("create failed!");

    let mut content = String::new();
    fp.read_to_string(&mut content).expect("read error");
    let mut block: Vec<&str> = content.split('>').collect();
    block.remove(0); // skip first ">"

    let mut hash = HashMap::new();
    let mut chrs = Vec::new();
    let mut total = Datf {
        A: 0,
        T: 0,
        G: 0,
        C: 0,
        N: 0,
        GC: 0.0,
        LEN: 0,
    };

    for each in block {
        let mut chr_seq: Vec<&str> = each.split('\n').collect();
        let head: Vec<&str> = chr_seq.remove(0).split(' ').collect();
        chrs.push(head[0].clone());
        eprintln!("start run chromosome: {}", &head[0]);

        let mut df = Datf {
            A: 0,
            T: 0,
            G: 0,
            C: 0,
            N: 0,
            GC: 0.0,
            LEN: 0,
        };

        for line in chr_seq {
            for base in line.chars() {
                df.LEN += 1;

                if base == 'A' || base == 'a' {
                    df.A += 1;
                }
                if base == 'T' || base == 't' {
                    df.T += 1;
                }
                if base == 'G' || base == 'g' {
                    df.G += 1;
                }
                if base == 'C' || base == 'c' {
                    df.C += 1;
                }
                if base == 'N' || base == 'n' {
                    df.N += 1;
                }
            }
        }
        df.GC = (df.G + df.C) as f64 / df.LEN as f64;
        hash.insert(head[0], df);
    }
    

    for (_k, v) in &hash {
        total.A += v.A;
        total.T += v.T;
        total.G += v.G;
        total.C += v.C;
        total.N += v.N;
        total.LEN += v.LEN;
    }
    total.GC = (total.G + total.C) as f64 / total.LEN as f64;
    hash.insert("ALL", total);
    chrs.push("ALL");

    op.write(
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "chr_Name", "base_A", "base_T", "base_G", "base_C", "base_N", "GC_Rate", "seq_Len"
        )
        .as_bytes(),
    )
    .expect("head error");
    
    for n in &chrs {
        let val = hash.get(n).unwrap();
        //println!("{}\t{:?}", n, val);
        op.write(
            format!(
                "{}\t{:?}\t{:?}\t{:?}\t{:?}\t{:?}\t{:.2}%\t{:?}\n",
                n,
                val.A,
                val.T,
                val.G,
                val.C,
                val.N,
                val.GC * 100.0,
                val.LEN
            )
            .as_bytes(),
        )
        .expect("write failed!");
    }
}

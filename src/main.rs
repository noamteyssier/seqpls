use std::io::Write;
use std::sync::Arc;

use anyhow::{Result, bail};
use clap::Parser;
use memchr::memmem::Finder;
use paraseq::fastq::Reader;
use paraseq::fastx::Record;
use paraseq::parallel::{
    PairedParallelProcessor, PairedParallelReader, ParallelProcessor, ParallelReader,
};
use parking_lot::Mutex;

type Patterns = Vec<Finder<'static>>;
type Expressions = Vec<regex::bytes::Regex>;
type BoxedWriter = Box<dyn Write + Send>;
const DEFAULT_BUFFER_SIZE: usize = 64 * 1024;

#[derive(Clone)]
struct GrepProcessor {
    /// Patterns to search for
    mp1: Patterns, // in primary
    mp2: Patterns, // in secondary
    pat: Patterns, // in either

    /// Regex expressions to match on
    re1: Expressions, // in primary
    re2: Expressions, // in secondary
    re: Expressions,  // in either

    /// Invert the pattern selection
    invert: bool,

    /// Local write buffers
    local_buffer: Vec<u8>,
    local_counter: usize,

    /// Global values
    global_writer: Arc<Mutex<BoxedWriter>>,
    global_counter: Arc<Mutex<usize>>,
}
impl GrepProcessor {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mp1: Patterns,
        mp2: Patterns,
        pat: Patterns,
        re1: Expressions,
        re2: Expressions,
        re: Expressions,
        invert: bool,
        output: BoxedWriter,
    ) -> Self {
        Self {
            mp1,
            mp2,
            pat,
            re1,
            re2,
            re,
            invert,
            global_writer: Arc::new(Mutex::new(output)),
            local_buffer: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
            local_counter: 0,
            global_counter: Arc::new(Mutex::new(0)),
        }
    }

    fn search_primary(&self, seq: &[u8]) -> bool {
        if self.mp1.is_empty() || seq.is_empty() {
            return true;
        }
        self.mp1.iter().all(|pat| pat.find(seq).is_some())
    }

    fn search_secondary(&self, seq: &[u8]) -> bool {
        if self.mp2.is_empty() || seq.is_empty() {
            return true;
        }
        self.mp2.iter().any(|pat| pat.find(seq).is_some())
    }

    fn search_either(&self, r1: &[u8], r2: &[u8]) -> bool {
        if self.pat.is_empty() {
            return true;
        }
        self.pat
            .iter()
            .any(|pat| pat.find(r1).is_some() || pat.find(r2).is_some())
    }

    fn regex_primary(&self, seq: &[u8]) -> bool {
        if self.re1.is_empty() || seq.is_empty() {
            return true;
        }
        self.re1.iter().any(|re| re.find(seq).is_some())
    }

    fn regex_secondary(&self, seq: &[u8]) -> bool {
        if self.re2.is_empty() || seq.is_empty() {
            return true;
        }
        self.re2.iter().any(|re| re.find(seq).is_some())
    }

    fn regex_either(&self, r1: &[u8], r2: &[u8]) -> bool {
        if self.re.is_empty() {
            return true;
        }
        self.re
            .iter()
            .any(|re| re.find(r1).is_some() || re.find(r2).is_some())
    }

    pub fn pattern_match(&self, primary: &[u8], secondary: &[u8]) -> bool {
        let pred = self.search_primary(primary)
            && self.search_secondary(secondary)
            && self.search_either(primary, secondary)
            && self.regex_primary(primary)
            && self.regex_secondary(secondary)
            && self.regex_either(primary, secondary);
        if self.invert { !pred } else { pred }
    }

    pub fn write_record<Rf: Record>(&mut self, record: Rf) -> Result<()> {
        self.local_buffer.write(b"@")?;
        self.local_buffer.extend_from_slice(record.id());
        self.local_buffer.write(b"\n")?;
        self.local_buffer.extend_from_slice(record.seq());
        self.local_buffer.write(b"\n+\n")?;
        self.local_buffer.extend_from_slice(record.qual().unwrap());
        self.local_buffer.write(b"\n")?;
        Ok(())
    }
}
impl ParallelProcessor for GrepProcessor {
    fn process_record<Rf: paraseq::fastx::Record>(
        &mut self,
        record: Rf,
    ) -> paraseq::parallel::Result<()> {
        if self.pattern_match(record.seq(), &[]) {
            self.write_record(record)?;
            self.local_counter += 1;
        }
        Ok(())
    }
    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        {
            let mut global_writer = self.global_writer.lock();
            global_writer.write_all(&self.local_buffer)?;
            global_writer.flush()?;
        }
        self.local_buffer.clear();

        *self.global_counter.lock() += self.local_counter;
        self.local_counter = 0;
        Ok(())
    }
}
impl PairedParallelProcessor for GrepProcessor {
    fn process_record_pair<Rf: Record>(
        &mut self,
        record1: Rf,
        record2: Rf,
    ) -> paraseq::parallel::Result<()> {
        if self.pattern_match(record1.seq(), record2.seq()) {
            self.write_record(record1)?;
            self.write_record(record2)?;
            self.local_counter += 1;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        {
            let mut global_writer = self.global_writer.lock();
            global_writer.write_all(&self.local_buffer)?;
            global_writer.flush()?;
        }
        self.local_buffer.clear();

        *self.global_counter.lock() += self.local_counter;
        self.local_counter = 0;
        Ok(())
    }
}

fn match_output(path: Option<String>) -> Result<BoxedWriter> {
    if let Some(path) = path {
        Ok(Box::new(std::fs::File::create(path)?))
    } else {
        Ok(Box::new(std::io::stdout()))
    }
}

fn grep_paired(
    r1_path: &str,
    r2_path: &str,
    output: Option<String>,
    num_threads: usize,
    args: &GrepArgs,
) -> Result<()> {
    let (r1_handle, _comp) = niffler::send::from_path(r1_path)?;
    let (r2_handle, _comp) = niffler::send::from_path(r2_path)?;

    let r1_reader = Reader::new(r1_handle);
    let r2_reader = Reader::new(r2_handle);
    let output = match_output(output)?;

    let processor = GrepProcessor::new(
        args.bytes_mp1(),
        args.bytes_mp2(),
        args.bytes_pat(),
        args.bytes_reg1(),
        args.bytes_reg2(),
        args.bytes_reg(),
        args.invert,
        output,
    );

    r1_reader.process_parallel_paired(r2_reader, processor, num_threads)?;

    Ok(())
}

fn grep_single(
    r1_path: &str,
    output: Option<String>,
    num_threads: usize,
    args: &GrepArgs,
) -> Result<()> {
    let (r1_handle, _comp) = niffler::send::from_path(r1_path)?;

    let r1_reader = Reader::new(r1_handle);
    let output = match_output(output)?;

    let processor = GrepProcessor::new(
        args.bytes_mp1(),
        args.bytes_mp2(),
        args.bytes_pat(),
        args.bytes_reg1(),
        args.bytes_reg2(),
        args.bytes_reg(),
        args.invert,
        output,
    );

    r1_reader.process_parallel(processor, num_threads)?;

    Ok(())
}

fn main() -> Result<()> {
    let args = GrepCommand::parse();
    args.grep.validate()?;

    if args.inputs.len() == 2 {
        grep_paired(
            &args.inputs[0],
            &args.inputs[1],
            args.output,
            args.threads,
            &args.grep,
        )
    } else if args.inputs.len() == 1 {
        grep_single(&args.inputs[0], args.output, args.threads, &args.grep)
    } else {
        bail!("Must provide either 1 or 2 input files")
    }
}

/// Grep a BINSEQ file and output to FASTQ or FASTA.
#[derive(Parser, Debug)]
pub struct GrepCommand {
    /// Input file paths (two files are assumed paired)
    #[clap(num_args= 1..=2, required = true)]
    pub inputs: Vec<String>,

    /// Output file path [default: stdout]
    #[clap(short = 'o', long)]
    pub output: Option<String>,

    /// Number of threads to use [default: 1]
    #[clap(short = 'T', long, default_value_t = 1)]
    pub threads: usize,

    #[clap(flatten)]
    pub grep: GrepArgs,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "SEARCH OPTIONS")]
pub struct GrepArgs {
    /// Fixed string pattern to search for in primary sequence
    #[clap(short = 'e', long)]
    pub pat1: Vec<String>,

    /// Fixed string pattern to search for in secondary sequence
    #[clap(short = 'E', long)]
    pub pat2: Vec<String>,

    /// Pattern to search for in either sequence
    #[clap(short = 'F', long)]
    pub pat: Vec<String>,

    /// Regex expression to search for in primary sequence
    #[clap(short = 'r', long)]
    pub reg1: Vec<String>,

    /// Regex expression to search for in secondary sequence
    #[clap(short = 'R', long)]
    pub reg2: Vec<String>,

    /// Regex expression to search for in either sequence
    #[clap(short = 'P', long)]
    pub reg: Vec<String>,

    /// Invert pattern criteria (like grep -v)
    #[clap(short = 'v', long)]
    pub invert: bool,
}
impl GrepArgs {
    pub fn validate(&self) -> Result<()> {
        if self.pat1.is_empty()
            && self.pat2.is_empty()
            && self.pat.is_empty()
            && self.reg1.is_empty()
            && self.reg2.is_empty()
            && self.reg.is_empty()
        {
            anyhow::bail!("At least one pattern must be specified");
        }
        Ok(())
    }
    pub fn bytes_mp1(&self) -> Vec<Finder<'static>> {
        self.pat1
            .iter()
            .map(|s| Finder::new(s.as_bytes()))
            .map(|f| f.into_owned())
            .collect()
    }
    pub fn bytes_mp2(&self) -> Vec<Finder<'static>> {
        self.pat2
            .iter()
            .map(|s| Finder::new(s.as_bytes()))
            .map(|f| f.into_owned())
            .collect()
    }
    pub fn bytes_pat(&self) -> Vec<Finder<'static>> {
        self.pat2
            .iter()
            .map(|s| Finder::new(s.as_bytes()))
            .map(|f| f.into_owned())
            .collect()
    }
    pub fn bytes_reg1(&self) -> Vec<regex::bytes::Regex> {
        self.reg1
            .iter()
            .map(|s| regex::bytes::Regex::new(s).expect("Could not build regex from pattern: {s}"))
            .collect()
    }
    pub fn bytes_reg2(&self) -> Vec<regex::bytes::Regex> {
        self.reg2
            .iter()
            .map(|s| regex::bytes::Regex::new(s).expect("Could not build regex from pattern: {s}"))
            .collect()
    }
    pub fn bytes_reg(&self) -> Vec<regex::bytes::Regex> {
        self.reg
            .iter()
            .map(|s| regex::bytes::Regex::new(s).expect("Could not build regex from pattern: {s}"))
            .collect()
    }
}

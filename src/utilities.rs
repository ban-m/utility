use bio;
use dtw;
use histogram_minimizer::pc;
use knn_predictor::knn::KNN;
use rayon::prelude::*;
use squiggler;
use std;
use std::collections::HashMap;
use std::fs::{read_dir, File};
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

/// Get sam file from specified file. Each entry can be accessed via its id.
/// Each entry consists of two element, flag and location.
/// For semantics of flag, see samfile specification elsewhere. Same for location.
pub fn get_sam(path: &str) -> std::io::Result<HashMap<String, (u32, usize)>> {
    Ok(BufReader::new(File::open(&Path::new(path))?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let contents: Vec<_> = e.split('\t').collect();
            let id = match contents[0].split('_').nth(0) {
                Some(id) => id.to_string(),
                None => return None,
            };
            let flag: u32 = match contents[1].parse() {
                Ok(res) if res == 0 || res == 16 => res,
                _ => return None,
            };
            let location: usize = match contents[3].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            let quality: usize = match contents[4].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            if quality > 10 {
                Some((id, (flag, location)))
            } else {
                None
            }
        })
        .collect())
}

/// Set up template complement strand from specified fasta file.
/// If you use fastq file, please call setup_tempate_complement_fastq instead.
/// (template,reverse)
pub fn setup_template_complement(path: &Path) -> std::result::Result<(String, String), ()> {
    let seq = bio::io::fasta::Reader::from_file(path).map_err(|_| ())?;
    let reference: Vec<u8> = seq
        .records()
        .filter_map(|e| e.ok())
        .fold(Vec::new(), |mut acc, e| {
            acc.extend_from_slice(&mut e.seq());
            acc
        });
    let reverse = bio::alphabets::dna::revcomp(&reference);
    let temp = String::from_utf8(reference).map_err(|_| ())?;
    let rev = String::from_utf8(reverse).map_err(|_| ())?;
    Ok((temp, rev))
}

/// Set up template complement strand from specified fasta file.
/// If you use fastq file, please call setup_tempate_complement_fastq instead.
/// (template,reverse)
pub fn setup_template_complement_autofit(
    path: &Path,
    refsize: usize,
) -> std::result::Result<(String, String), ()> {
    let seq = bio::io::fasta::Reader::from_file(path).map_err(|_| ())?;
    let reference: Vec<u8> = seq
        .records()
        .filter_map(|e| e.ok())
        .fold(Vec::new(), |mut acc, e| {
            acc.extend_from_slice(&mut e.seq());
            acc
        });
    let reference = extend_vector(reference, refsize / 2);
    let reverse = bio::alphabets::dna::revcomp(&reference);
    let temp = String::from_utf8(reference).map_err(|_| ())?;
    let rev = String::from_utf8(reverse).map_err(|_| ())?;
    Ok((temp, rev))
}

fn extend_vector(text: Vec<u8>, length: usize) -> Vec<u8> {
    let mut result = vec![];
    loop {
        if result.len() < length {
            result.append(&mut text.clone());
        } else {
            break result.into_iter().take(length).collect();
        }
    }
}

/// deplicate and concatinate string until it reaches specified length
pub fn extend_strand(strand: &str, len: usize) -> String {
    let mut result = String::new();
    loop {
        if result.len() < len {
            result += strand;
        } else {
            return result;
        }
    }
}

/// Convert given string to signal by using specified model and
/// normalize it.
pub fn convert_to_squiggle(strand: &str, model: &squiggler::Squiggler) -> Vec<f32> {
    let sig: Vec<_> = model
        .get_signal_from_fasta(&strand)
        .into_iter()
        .map(|e| e.2)
        .collect();
    dtw::normalize(&sig, dtw::NormalizeType::Z)
}

/// Parse the given mode. Currently for only Sub, SakoeChiba and scouting.
pub fn get_mode(mode: &str) -> std::result::Result<dtw::Mode, String> {
    if mode.starts_with("Sub") {
        Ok(dtw::Mode::Sub)
    } else if mode.starts_with("Chiba") {
        let bandwidth: usize = match mode.split(',').nth(1).and_then(|e| e.parse().ok()) {
            Some(b) => b,
            None => return Err("Given chiba but couldn't parse bandwidth correctly".to_string()),
        };
        Ok(dtw::Mode::SakoeChiba(bandwidth))
    } else if mode.starts_with("Scouting") {
        let contents: Vec<usize> = mode
            .split(',')
            .skip(1)
            .filter_map(|e| e.parse().ok())
            .collect();
        if contents.len() == 2 {
            Ok(dtw::Mode::Scouting(contents[0], contents[1]))
        } else {
            return Err("Given scouting but couldn't parse argument.".to_string());
        }
    } else if mode.starts_with("FastSub") {
        match mode.split(',').skip(1).nth(0).and_then(|e| e.parse().ok()) {
            Some(res) => Ok(dtw::Mode::FastSub(res)),
            None => return Err("Not Valid radious".to_string()),
        }
    } else {
        return Err(format!("invalid mode name:{}", mode));
    }
}

/// Extract events sequence from specified directry and
/// split from 50 to 50 + querysize and its id.
/// Note that this function doesn't normalize these events.
#[cfg(feature = "python")]
pub fn get_queries(
    querypath: &str,
    querysize: usize,
    takenum: usize,
) -> Result<Vec<(String, Vec<f32>)>, std::io::Error> {
    use fast5wrapper;
    Ok(read_dir(&Path::new(querypath))?
        .filter_map(|e| e.ok())
        .filter_map(|e| e.path().to_str().map(|e| e.to_string()))
        .filter_map(|e| match fast5wrapper::get_read_id(&e) {
            Ok(id) => match fast5wrapper::get_event_for(&e, 50, querysize)
                .map(|e| e.into_iter().map(|e| e[2]).collect())
            {
                Ok(res) => Some((id, res)),
                Err(why) => {
                    eprintln!("{:?}", why);
                    return None;
                }
            },
            Err(why) => {
                eprintln!("{:?}", why);
                return None;
            }
        })
        .take(takenum)
        .collect())
}

/// Merge given sma file with given query file.
pub fn merge_queries_and_sam<T>(
    queries: &Vec<(String, Vec<T>)>,
    sam: &HashMap<String, (u32, usize)>,
) -> Vec<(Vec<T>, u32, usize)>
where
    T: Copy,
{
    queries
        .into_iter()
        .filter_map(|&(ref id, ref query)| {
            sam.get(id)
                .map(|&(flag, location)| (query.clone(), flag, location))
        })
        .collect()
}

/// Get dataset from specified queries and template, reverse strand ,mode, refsize, power, metric and reference size.
pub fn get_dataset(
    queries: &Vec<(Vec<f32>, u32, usize)>,
    temp: &Vec<f32>,
    rev: &Vec<f32>,
    refsize: usize,
    querysize: usize,
    mode: &dtw::Mode,
    power: f32,
    metric: &str,
) -> Vec<(f64, f64)> {
    queries
        .par_iter()
        .map(|&(ref query, flag, location)| {
            if flag == 0 {
                (query, flag, location)
            } else {
                (
                    query,
                    flag,
                    if rev.len() > location + refsize {
                        rev.len() - location - refsize
                    } else {
                        0
                    },
                )
            }
        })
        .filter_map(|(query, flag, location)| {
            let pos_ref = if flag == 0 { temp } else { rev };
            let neg_ref = if flag == 0 { rev } else { temp };
            if let Some(pos) = get_score(
                &query, location, pos_ref, refsize, querysize, mode, power, metric,
            ) {
                if let Some(neg) = get_score(
                    &query, location, neg_ref, refsize, querysize, mode, power, metric,
                ) {
                    Some((pos, neg))
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect()
}

/// Inner function to call from get_datasets, but
/// I put it in public scope because
/// sometimes we just want to
/// compute only for very small amount of query, difference setting.
#[inline]
pub fn get_score(
    query: &Vec<f32>,
    location: usize,
    reference: &Vec<f32>,
    refsize: usize,
    querysize: usize,
    mode: &dtw::Mode,
    power: f32,
    metric: &str,
) -> Option<f64> {
    let query = dtw::normalize(&query, dtw::NormalizeType::Z);
    // important
    let query = squiggler::dedup(&query, power);
    if query.len() < querysize {
        return None;
    };
    let query: Vec<_> = query.into_iter().take(querysize).collect();
    let offset = 500;
    let subref = if location < offset {
        &reference[0..refsize]
    } else if location + refsize - offset > reference.len() {
        &reference[reference.len() - offset - refsize..]
    } else {
        &reference[location - offset..location - offset + refsize]
    };
    dtw::utils::dtw_wrapper(&query, subref, mode, metric, &None, &None).map(|e| e as f64)
}

/// function to compute cross validataion
/// for specified K packs.
/// queries should consists of (event,flag,location) without normalization.
/// Returns: true positive,false positive,positive num,total num
pub fn cross_validation(
    queries: &Vec<(Vec<f32>, u32, usize)>,
    temp: &Vec<f32>,
    rev: &Vec<f32>,
    refsize: usize,
    querysize: usize,
    mode: &dtw::Mode,
    power: f32,
    metric: &str,
    model: &str,
    k: usize,
) -> (u32, u32, u32, u32) {
    let (temp, rev) = {
        let mut current_length = rev.len();
        let mut elonged_temp = temp.clone();
        let mut elonged_rev = rev.clone();
        loop {
            if current_length > refsize {
                break;
            } else {
                current_length += rev.len();
                elonged_temp.extend(&temp.clone());
                elonged_rev.extend(&rev.clone());
            };
        }
        (elonged_temp, elonged_rev)
    };
    let dataset = get_dataset(
        queries, &temp, &rev, refsize, querysize, mode, power, metric,
    );
    validate(&dataset, model, k)
}

/// k-cross-varidation by specified model:KNN or PC.
/// Returns: true positive,false positive,positive num,total num
pub fn validate(dataset: &Vec<(f64, f64)>, model: &str, k: usize) -> (u32, u32, u32, u32) {
    compute_k_folds(dataset, k)
        .iter()
        .map(|&(ref train, ref test)| {
            let train = train.iter().fold(vec![], |mut acc, &(pos, neg)| {
                acc.push((pos, true));
                acc.push((neg, false));
                acc
            });
            let test = test.iter().fold(vec![], |mut acc, &(pos, neg)| {
                acc.push((pos, true));
                acc.push((neg, false));
                acc
            });
            validate_for_single_pack(&train, &test, model)
        })
        .fold((0, 0, 0, 0), |acc, x| {
            (acc.0 + x.0, acc.1 + x.1, acc.2 + x.2, acc.3 + x.3)
        })
}

/// validation for single pack
/// Returns: true positive,false positive,positive num,total num
pub fn validate_for_single_pack(
    train: &Vec<(f64, bool)>,
    test: &Vec<(f64, bool)>,
    model: &str,
) -> (u32, u32, u32, u32) {
    let result = if model == "KNN" {
        let predictor = KNN::new(&train, 15); //15-NN
        test.iter()
            .map(|&(score, is_true)| (predictor.predict(score), is_true))
            .collect::<Vec<_>>()
    } else {
        let negdata = train
            .iter()
            .filter(|&&(_, is_pos)| !is_pos)
            .map(|&(score, _)| score)
            .collect();
        let posdata = train
            .iter()
            .filter(|&&(_, is_pos)| is_pos)
            .map(|&(score, _)| score)
            .collect();
        let predictor = pc::PC::new(negdata, posdata);
        test.iter()
            .map(|&(score, is_true)| (predictor.classify(score), is_true))
            .collect::<Vec<_>>()
    };
    result
        .into_iter()
        .fold((0, 0, 0, 0), |acc, (predict, answer)| {
            if predict && answer {
                (acc.0 + 1, acc.1, acc.2 + 1, acc.3 + 1)
            } else if predict && !answer {
                (acc.0, acc.1 + 1, acc.2, acc.3 + 1)
            } else if !predict && answer {
                (acc.0, acc.1, acc.2 + 1, acc.3 + 1)
            } else {
                (acc.0, acc.1, acc.2, acc.3 + 1)
            }
        })
}

/// compute k folds
pub fn compute_k_folds(
    dataset: &Vec<(f64, f64)>,
    k: usize,
) -> Vec<(Vec<(f64, f64)>, Vec<(f64, f64)>)> {
    let window_size = dataset.len() / k;
    let mut result = std::vec::Vec::with_capacity(k);
    for i in 0..(k - 1) {
        let testset: Vec<_> = dataset[i * window_size..(i + 1) * window_size]
            .iter()
            .map(|e| e.clone())
            .collect();
        let trainset: Vec<_> = dataset
            .iter()
            .enumerate()
            .filter(|&(idx, _)| idx < i * window_size || idx >= (i + 1) * window_size)
            .map(|(_, &e)| e.clone())
            .collect();
        result.push((trainset, testset));
    }
    let testset: Vec<_> = dataset[(k - 1) * window_size..]
        .iter()
        .map(|e| e.clone())
        .collect();
    let trainset: Vec<_> = dataset[..(k - 1) * window_size]
        .iter()
        .map(|e| e.clone())
        .collect();
    result.push((trainset, testset));
    result
}

/// Parse given file into events streams. Each event stream should be separated by "---" charactor.
pub fn parse_events(file: &Path) -> std::io::Result<Vec<Vec<(f32, f32, usize, usize)>>> {
    let mut result = vec![];
    let mut temp = vec![];
    for line in BufReader::new(File::open(file)?)
        .lines()
        .filter_map(|e| e.ok())
    {
        let line = if line.contains("---") {
            if !temp.is_empty() {
                result.push(temp);
            }
            temp = vec![];
            let line = line.replace("---", "");
            line
        } else {
            line
        };
        let contents: Vec<_> = line.split(',').collect();
        if contents.len() < 4 {
            continue;
        }
        let mean: f32 = match contents[0].parse() {
            Ok(res) => res,
            Err(_) => continue,
        };
        let sd: f32 = match contents[1].parse() {
            Ok(res) => res,
            Err(_) => continue,
        };
        let index: f32 = match contents[2].parse() {
            Ok(res) => res,
            Err(_) => continue,
        };
        let length: f32 = match contents[3].parse() {
            Ok(res) => res,
            Err(_) => continue,
        };
        temp.push((mean, sd, index as usize, length as usize));
    }
    result.push(temp);
    Ok(result)
}

/// Parse given file into events stream.
pub fn parse_event(file: &Path) -> std::io::Result<Vec<(f32, f32, usize, usize)>> {
    Ok(BufReader::new(File::open(file)?)
        .lines()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let contents = e.split(',').collect::<Vec<_>>();
            if contents.len() < 4 {
                return None;
            }
            let mean: f32 = match contents[0].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            let sd: f32 = match contents[1].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            let index: f32 = match contents[2].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            let length: f32 = match contents[3].parse() {
                Ok(res) => res,
                Err(_) => return None,
            };
            Some((mean, sd, index as usize, length as usize))
        })
        .collect())
}

/// Asess if the given query is actual event stream
pub fn is_actual_event(data: &[f32], thr: f32) -> bool {
    let max = data.iter().fold(0., |a, &b| if a > b { a } else { b });
    let min = data.iter().fold(10000., |a, &b| if a > b { b } else { a });
    (max - min) > thr
}

/// Get signals from the given directly.
pub fn get_signals(
    querypath: &str,
    querysize: usize,
    takenum: usize,
) -> Result<Vec<(String, Vec<u32>)>, std::io::Error> {
    eprintln!("{}", querypath);
    Ok(read_dir(&Path::new(querypath))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter_map(|e| {
            let f = match File::open(e) {
                Ok(res) => res,
                Err(why) => {
                    eprintln!("{:?}", why);
                    return None;
                }
            };
            let contents: Vec<_> = BufReader::new(&f).lines().filter_map(|e| e.ok()).collect();
            if contents.len() < 2500 {
                return None;
            }
            let id = match contents[0].split(':').nth(0) {
                Some(res) => res,
                None => "mock",
            }
            .to_string();
            let signals: Vec<_> = contents[1..]
                .into_iter()
                .take(querysize)
                .filter_map(|e| e.parse().ok())
                .collect();
            Some((id, signals))
        })
        .take(takenum)
        .collect())
}

/// Get signals from the given directly.
pub fn get_signals_with_filename(
    querypath: &str,
    querysize: usize,
    takenum: usize,
) -> Result<Vec<(std::path::PathBuf, Vec<u32>)>, std::io::Error> {
    eprintln!("{}", querypath);
    Ok(read_dir(&Path::new(querypath))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|e| e.is_file())
        .filter_map(|e| {
            if let Some(ext) = e.extension() {
                if ext != "dat" {
                    return None;
                }
            }
            let f = match File::open(e.clone()) {
                Ok(res) => res,
                Err(why) => {
                    eprintln!("{:?}", why);
                    return None;
                }
            };
            let contents: Vec<_> = BufReader::new(&f).lines().filter_map(|e| e.ok()).collect();
            if contents.len() < 2500 {
                return None;
            }
            let signals: Vec<_> = contents[1..]
                .into_iter()
                .take(querysize)
                .filter_map(|e| e.parse().ok())
                .collect();
            Some((e, signals))
        })
        .take(takenum)
        .collect())
}

/// Get GC content of the given reference
pub fn get_gc_content(genome: &str) -> f32 {
    let gc = genome
        .matches(|e| e == 'G' || e == 'C' || e == 'g' || e == 'c')
        .count();
    gc as f32 / genome.len() as f32
}

#[derive(Clone, Copy, Debug)]
/// A struct to express a event after segmentation for raw signal.
/// It is almost the same as the inner expression of MinKNOW's event struct.
pub struct Event {
    /// The mean value of this events.
    pub mean: f64,
    /// The standard deviation
    pub stdv: f64,
    /// The start index from the start of the raw signal.
    pub start: usize,
    /// The end index from the start of the raw signal.
    pub length: usize,
}
use std::fmt;
impl fmt::Display for Event {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{},{},{},{}",
            self.mean, self.stdv, self.start, self.length
        )
    }
}
fn compute_sum_sumsq(data: &[u32]) -> (Vec<u64>, Vec<u64>) {
    // compute accumulate sum and squared sum.
    let (mut sum, mut sumsq) = (
        Vec::with_capacity(data.len()),
        Vec::with_capacity(data.len()),
    );
    let (mut current_sum, mut current_sumsq) = (0, 0);
    for &x in data.iter() {
        let x = x as u64;
        current_sum += x;
        current_sumsq += x * x;
        sum.push(current_sum);
        sumsq.push(current_sumsq);
    }
    (sum, sumsq)
}

const EPSILON: f64 = 0.000000000000000000000000001;
fn compute_tstat(sums: &[u64], sumsqs: &[u64], window: usize) -> Vec<f64> {
    // compute t-statistics value for each successive window.
    (0..sums.len())
        .map(|i| {
            if i < window {
                0.
            } else if i < sums.len() - window {
                let lsums = if i == window { 0 } else { sums[i - window - 1] };
                let lsumsq = if i == window {
                    0
                } else {
                    sumsqs[i - window - 1]
                };
                let sum1 = sums[i - 1] - lsums;
                let sum2 = sums[i + window - 1] - sums[i - 1];
                let sumsq1 = sumsqs[i - 1] - lsumsq;
                let sumsq2 = sumsqs[i + window - 1] - sumsqs[i - 1];
                let deltasq = (sum2 - sum1).pow(2);
                let denom = sumsq1 + sumsq2 - (sum1.pow(2) + sum2.pow(2)) / window as u64;
                let denom = if denom as f64 > EPSILON {
                    denom as f64
                } else {
                    EPSILON
                };
                (deltasq as f64 / denom).sqrt()
            } else {
                0.
            }
        })
        .collect()
}

fn short_long_peak_detector(
    tstat1: &[f64],
    tstat2: &[f64],
    windows: &[usize; 2],
    thresholds: &[f64; 2],
    peak_hight: f64,
) -> Vec<usize> {
    let len = tstat1.len();
    let def_peak_pos = -1;
    let def_peak_val = 10000.;
    let mut peak_values = [def_peak_val, def_peak_val];
    let mut peak_pos = [def_peak_pos, def_peak_pos];
    let mut maskedto = -1;
    let mut validpeak = [false, false];
    let mut times = vec![];
    times.push(0);
    for i in 0..len {
        for k in 0..2 {
            if k == 1 && maskedto >= i as i32 {
                continue;
            }

            let current_value = if k == 0 { tstat1[i] } else { tstat2[i] };
            if peak_pos[k] == def_peak_pos {
                if current_value < peak_values[k] {
                    peak_values[k] = current_value;
                } else if current_value - peak_values[k] > peak_hight {
                    peak_values[k] = current_value;
                    peak_pos[k] = i as i32;
                }
            } else {
                if current_value > peak_values[k] {
                    peak_values[k] = current_value;
                    peak_pos[k] = i as i32;
                }
                if k == 0 {
                    if peak_values[0] > thresholds[0] {
                        maskedto = peak_pos[0] + windows[0] as i32;
                        peak_pos[1] = def_peak_pos;
                        peak_values[1] = def_peak_val;
                        validpeak[1] = false;
                    }
                }
                if peak_values[k] - current_value > peak_hight && peak_values[k] > thresholds[k] {
                    validpeak[k] = true;
                }
                if validpeak[k] && (i as i32 - peak_pos[k] > windows[k] as i32 / 2) {
                    times.push(peak_pos[k]);
                    peak_pos[k] = def_peak_pos;
                    peak_values[k] = current_value;
                    validpeak[k] = false;
                }
            }
        }
    }
    times
        .into_iter()
        .filter(|&e| e >= 0)
        .map(|e| e as usize)
        .collect()
}
/// Trim signal assumed as a noise front of the query.
/// To do that, the difference of signals is first calculate,
/// then, skip the signal while the differences are few.
pub fn trim_front(signal: &Vec<u32>, thr: u32, count: usize) -> usize {
    let mut good_so_far = 0;
    let location = signal
        .windows(2)
        .into_iter()
        .map(|e| {
            if e[1] > e[0] {
                e[1] - e[0]
            } else {
                e[0] - e[1]
            }
        })
        .take_while(|diff| {
            good_so_far += if diff > &thr { 1 } else { 0 };
            good_so_far < count
        })
        .count();
    if location > signal.len() {
        return signal.len();
    } else {
        location
    }
}

/// Trim signal assumed as a noise front of the query.
/// To do that, the difference of signals is first calculate,
/// then, skip the signal while the differences are few.
pub fn trim_front_vec(signal: &[i32], thr: u32, count: usize) -> Vec<u32> {
    let thr = thr as i32;
    let mut good_so_far = 0;
    let position = signal
        .windows(2)
        .into_iter()
        .map(|e| (e[1] - e[0]).abs())
        .take_while(|diff| {
            good_so_far += if diff > &thr { 1 } else { 0 };
            good_so_far < count
        })
        .count();
    if position >= signal.len() {
        vec![]
    } else {
        let (median, median_of_dev) = get_median_median_of_dev(&signal[position..]);
        signal[position..]
            .iter()
            .filter(|&e| (e - median) < 7 * median_of_dev)
            .map(|&e| e as u32)
            .collect()
    }
}

fn get_median_median_of_dev(signal: &[i32]) -> (i32, i32) {
    let mut signal: Vec<_> = signal.iter().collect();
    let len = signal.len();
    signal.sort();
    let median = if len % 2 == 0 {
        (signal[len / 2 - 1] + signal[len / 2]) / 2
    } else {
        *signal[len / 2]
    };
    let mut devs: Vec<_> = signal.iter().map(|&e| (e - median).abs()).collect();
    devs.sort();
    if len % 2 == 0 {
        (median, (devs[len / 2 - 1] + devs[len / 2]) / 2)
    } else {
        (median, devs[len / 2])
    }
}
/// This is a clone from ONT MinKNOW Core cpp source code.
pub fn event_detect_mult_ttest(
    data: &[u32],
    windows: &[usize; 2],
    thresholds: &[f64; 2],
    peak_hight: f64,
) -> Vec<Event> {
    let (sums, sumsqs) = compute_sum_sumsq(&data);
    let len = data.len();
    if len == 0 {
        vec![]
    } else if len <= 2 * windows[1].max(windows[0]) as usize {
        let len = data.len() as f64;
        let mean = sums[0] as f64 / len;
        let stdv = (sumsqs[0] as f64 / len - mean.powi(2)).sqrt();
        let start = 0;
        let length = data.len();
        let event = Event {
            mean,
            stdv,
            start,
            length,
        };
        vec![event]
    } else {
        let (tstat1, tstat2) = (
            compute_tstat(&sums, &sumsqs, windows[0]),
            compute_tstat(&sums, &sumsqs, windows[1]),
        );
        let starts = short_long_peak_detector(&tstat1, &tstat2, windows, thresholds, peak_hight);
        (0..starts.len())
            .map(|i| {
                let start = starts[i] - 1;
                let lsums = if i == 0 { 0 } else { sums[start] };
                let lsumsq = if i == 0 { 0 } else { sumsqs[start] };
                let end = if i == starts.len() - 1 {
                    data.len() - 1
                } else {
                    starts[i + 1] - 1
                };
                let length = end - start;
                let mean = (sums[end] - lsums) as f64 / length as f64;
                let stdv = ((sumsqs[end] - lsumsq) as f64 / length as f64 - mean.powi(2)).sqrt();
                Event {
                    mean,
                    stdv,
                    start,
                    length,
                }
            })
            .collect()
    }
}

// use serde::Serialize;
use std::fs::{OpenOptions};
use std::io::{Write, BufWriter};
use anyhow::Result;


const QE: f64 = 1.602176565e-19;  // C, electron charge
const EPS0: f64 = 8.85418782e-12; // C/V/m vac permiability
const ME: f64 = 9.10938215e-31;   // kg electron mass

#[derive(serde::Serialize)]
struct Sample {
       x0: f64,
       dx: f64,
       phi: Vec<f64>,
       rho: Vec<f64>,
       ef: Vec<f64>
}

#[derive(serde::Serialize)]
struct SeriesElement {
       x_i: f64,
       phi_i: f64,
       rho_i: f64,
       ef_i: f64
}

// Will delete the contents of the file if it already exists
fn write_json(sample :&Sample, filename: &str) -> Result<()> {
   let file = OpenOptions::new()
       .write(true)
       .create(true)
       .truncate(true)
       .open(filename)?;
       
   let mut writer = BufWriter::new(file);
   let line = serde_json::to_string(sample)?;
   writeln!(writer, "{}", line)?;

   writer.flush()?;
   Ok(())
}


// Will delete the contents of the file if it already exists

// Potential improvements here: conversion to struct in order to encode json
// might not be necessary...
// some extra vars are created as well
// could assert that the vectors are all the same length
fn write_series_json(sample :&Sample, filename: &str) -> Result<()> {
   let file = OpenOptions::new()
       .write(true)
       .create(true)
       .truncate(true)
       .open(filename)?;
       
   let mut writer = BufWriter::new(file);

   for (i, ((phi_i, rho_i), ef_i)) in sample.phi.iter()
       	   	    	    	      .zip(sample.rho.iter())
				      .zip(sample.ef.iter())
				      .enumerate() {
       let s = SeriesElement { x_i: sample.x0 + (i as f64) * sample.dx,
       	       		       phi_i: *phi_i,
			       rho_i: *rho_i,
			       ef_i: *ef_i };
       let line = serde_json::to_string(&s)?;
       writeln!(writer, "{}", line)?;
   }
   writer.flush()?;
   Ok(())
}

// could assert that the vectors are all the same length
fn solvePotentialDirect(dx: f64, phi: &mut Vec<f64>, rho: &Vec<f64>) -> Result<()> {
   let ni: usize = phi.len();
   let mut a: Vec<f64> = vec![0.0; ni];  // coef phi[i-1]
   let mut b: Vec<f64> = vec![0.0; ni];  // coef phi[i]
   let mut c: Vec<f64> = vec![0.0; ni];  // coef phi[i+1]
   let mut d: Vec<f64> = vec![0.0; ni];  // rhs

   // precompute some common values
   let inv_sq = 1.0 / (dx * dx);
   let two_inv_sq = -2.0 / (dx * dx);
   
   for i in 0..ni {
      if i == 0 || i == ni-1 {
      	 b[i] = 1.0;
	 d[i] = 0.0;
      } else {
         a[i] = inv_sq;
	 b[i] = two_inv_sq;
	 c[i] = inv_sq;
	 d[i] = -rho[i] / EPS0;
      }
   }

   c[0] = c[0] / b[0];
   d[0] = d[0] / b[0];

   for i in 1..ni {
      if i < (ni - 1) {
         c[i] = c[i] / (b[i] - a[i] * c[i-1]);
      }
      d[i] = (d[i] - a[i] * d[i-1])/(b[i] - a[i] * c[i-1]);
   }

   phi[ni-1] = d[ni-1];
   for i in (0..ni-1).rev() {
      phi[i] = d[i] - c[i] * phi[i+1];
   }
   
   Ok(())   
}

fn main() -> Result<()> {
    println!("Starting execution");
    let ni = 21;
    let x0: f64 = 0.0;
    let xm: f64 = 0.1;
    let dc: f64 = (xm - x0) / (ni - 1) as f64;

    let mut phi: Vec<f64> = vec![0.0; ni];
    let mut rho: Vec<f64> = vec![QE*1e12; ni];
    let mut ef: Vec<f64> = vec![0.0; ni];

    solvePotentialDirect(dc, &mut phi, &rho)?;
    
    let s1 = Sample{x0: x0, dx: dc, phi: phi, rho: rho, ef: ef};

    write_json(&s1, "output.json")?;
    write_series_json(&s1, "output_series.json")?;
    println!("Done writing files");
    Ok(())
    
    
}

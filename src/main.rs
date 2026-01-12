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

fn print_series(sample :&Sample) -> Result<()> {
   for (i, ((phi_i, rho_i), ef_i)) in sample.phi.iter()
       	   	    	    	      .zip(sample.rho.iter())
				      .zip(sample.ef.iter())
				      .enumerate() {
       let xi:f64 = sample.x0 + (i as f64) * sample.dx; 
       println!("X_i {xi:.6}, phi_i {:.3}, rho_i {:.3}, ef_i {:.3}",
                *phi_i, *rho_i, *ef_i);
   }
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

// want this to return a value indicating successful convergence
fn solvePotentialGsSOR(dx: f64, phi: &mut Vec<f64>, rho: &Vec<f64>, max_iter: i32) -> Result<i32, String> {

   let mut L2: f64 = 1e12;
   let L2_conv: f64 = 1e-6;
   let dx2: f64 = dx * dx;
   let w: f64 = 1.4;  //make this a param?
   let ni: usize = phi.len();

   let found = {
      let mut result = None;
      
      for iter in 0..max_iter {
         phi[0] = 0.0;
      	 phi[ni-1] = 0.0;

         for i in 1..(ni-1) {
            let g: f64 = 0.5 * (phi[i-1] + phi[i+1] + dx2 * rho[i] / EPS0);
	    phi[i] = phi[i] + w * (g - phi[i]);
         }

         if iter % 50 == 0 {
      	    let mut sum: f64 = 0.0;
	    for i in 1..(ni-1) {
	       let res: f64 = -rho[i]/EPS0 - (phi[i-1] - 2.0 * phi[i] + phi[i+1])/dx2;
	       sum += res * res;
	    }
	    let L2: f64 = (sum/(ni as f64)).sqrt();
	    if L2 < L2_conv {
	       result = Some(iter);
	       break;
	    }
	 }
      }
      result
   };
   match found {
      Some(i) => {
         println!("GS SOR solver converged after {i} iterations");
	 Ok(i)
      }
      None => {
         println!("GS SOR solver didn't converge, L2 residual {L2:.6}");
	 Err("GS SOR didn't converge".into())
      }
   }
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
    let s1 = Sample{x0: x0, dx: dc, phi: phi.clone(), rho: rho.clone(), ef: ef.clone()};
    print_series(&s1);

    // reset values for next sovler
    phi = vec![0.0; ni];
    rho = vec![QE*1e12; ni];
    ef = vec![0.0; ni];
    match solvePotentialGsSOR(dc, &mut phi, &rho, 1000) {
       Ok(i) => {
          let s1 = Sample{x0: x0, dx: dc, phi: phi.clone(), rho: rho.clone(), ef: ef.clone()};
          print_series(&s1);
       }
       Err(s) => println!("{s}")
    }

    let s1 = Sample{x0: x0, dx: dc, phi: phi.clone(), rho: rho.clone(), ef: ef.clone()};
    write_json(&s1, "output.json")?;
    write_series_json(&s1, "output_series.json")?;
    println!("Done writing files");
    Ok(())
    
    
}

// external
use anyhow::Result;
use std::fs::{OpenOptions};
use std::io::{Write, BufWriter};

// internal
use crate::constants::*;
use super::traits::OneDimFieldInterp;


#[derive(serde::Serialize)]
struct LinearOneDimCell {
    x: f64,
    phi: f64,
    rho: f64,
    ef: f64
}

#[derive(serde::Serialize)]
pub struct LinearOneDimField {
    x0: f64,
    dx: f64,
    cells: Vec<LinearOneDimCell>
}

// Qs
// Can you have optional arguments?
impl OneDimFieldInterp for LinearOneDimField {
    fn new(n: usize, x0: f64, dx: f64, rho0: f64) -> Result<Self> {
        let cells: Vec<LinearOneDimCell> = (0..n)
            .map(|i| {
                let x = x0 + i as f64 * dx;
                LinearOneDimCell {
                    x,
                    phi: 0.0,
                    rho: rho0,
                    ef: 0.0
                }
            }).collect();
        let mut field = Self {x0:x0, dx:dx, cells:cells };
        field.solve_potential()?;
        field.compute_ef()?;
        Ok(field)
    }

    fn len(&self) -> usize {self.cells.len() }

    // Assumes that x values remain unmodified
    fn reset(&mut self, rho0: f64) -> Result<()> {
        for c in &mut self.cells {
            c.phi = 0.0;
            c.rho = rho0;
            c.ef = 0.0;
        }
        self.solve_potential()?;
        self.compute_ef()?;
        Ok(())
    }

    fn phi_at(&self, x: f64) -> f64 {
        let li = self.x_to_l(x);
        self.gather_phi(li)
    }

    fn rho_at(&self, x: f64) -> f64 {
        let li = self.x_to_l(x);
        self.gather_rho(li)
    }

    fn ef_at(&self, x: f64) -> f64 {
        let li = self.x_to_l(x);
        self.gather_ef(li)
    }


    fn get_phi_max(&self) -> f64 {
        self.cells
            .iter().map(|c| c.phi)
            .reduce(f64::max).unwrap()
    }

    fn x0(&self) -> f64 { self.x0 }
    fn dx(&self) -> f64 { self.dx }

    // Will delete the contents of the file if it already exists
    fn write_field_to_json(&self, filename : &str) -> Result<()> {
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(filename)?;
       
        let mut writer = BufWriter::new(file);
        let line = serde_json::to_string(self)?; //??? is self a LinearOneDimField here?
        writeln!(writer, "{}", line)?;

        writer.flush()?;
        Ok(())
    }


    fn print(&self) {
        for cell in self.cells.iter() {
            println!("X_i {:.6}, phi_i {:.3}, rho_i {:.3e}, ef_i {:.3}",
                     cell.x, cell.phi, cell.rho, cell.ef);
        }
    }

}

impl LinearOneDimField {
    // private methods?  Also, missing lots of error handling right now -
    // need to figure out the right fix for that
    fn x_to_l(&self, x: f64) -> f64 {
        (x - self.x0) / self.dx
    }

    fn gather_ef(&self, li: f64) -> f64 {
        let i : usize = li as usize; // should we check that li is non-negative?
        let di : f64 = li.fract();
        self.cells[i].ef * (1.0 - di) + self.cells[i+1].ef * di
    }

    fn gather_phi(&self, li: f64) -> f64 {
        let i : usize = li as usize; // should we check that li is non-negative?
        let di : f64 = li.fract();
        self.cells[i].phi * (1.0 - di) + self.cells[i+1].phi * di
    }

    fn gather_rho(&self, li: f64) -> f64 {
        let i : usize = li as usize; // should we check that li is non-negative?
        let di : f64 = li.fract();
        self.cells[i].rho * (1.0 - di) + self.cells[i+1].rho * di
    }

    fn solve_potential(&mut self) -> Result<()> {
        // Hardcode which solver to use for now...
        // Direct solver
        //self.solve_potential_direct()?;

        // Iterative solver - Gauss Seidel
        match self.solve_potential_gs_sor(1000) {
            Ok(i) => {
                println!("GS SOR solver converged after {i} iterations");
                Ok(())
            }
            Err(s) => {
                println!("{s}");
                Err(anyhow::anyhow!(s))
            }
        }
    }

    fn solve_potential_direct(&mut self) -> Result<()> {
        let ni: usize = self.len();
        let mut a: Vec<f64> = vec![0.0; ni];  // coef phi[i-1]
        let mut b: Vec<f64> = vec![0.0; ni];  // coef phi[i]
        let mut c: Vec<f64> = vec![0.0; ni];  // coef phi[i+1]
        let mut d: Vec<f64> = vec![0.0; ni];  // rhs

        // precompute some common values
        let inv_sq = 1.0 / (self.dx * self.dx);
        let two_inv_sq = -2.0 / (self.dx * self.dx);
        
        for i in 0..ni {
            if i == 0 || i == ni-1 {
      	        b[i] = 1.0;
	        d[i] = 0.0;
            } else {
                a[i] = inv_sq;
	        b[i] = two_inv_sq;
	        c[i] = inv_sq;
	        d[i] = -self.cells[i].rho / EPS0;
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

        self.cells[ni-1].phi = d[ni-1];
        for i in (0..ni-1).rev() {
            self.cells[i].phi = d[i] - c[i] * self.cells[i+1].phi;
        }
   
        Ok(())   
    }

    // want this to return a value indicating successful convergence
    fn solve_potential_gs_sor(&mut self, max_iter: i32) -> Result<i32, String> {
        let mut l2: f64 = 1e12;
        let l2_conv: f64 = 1e-6;
        let dx2: f64 = self.dx * self.dx;
        let w: f64 = 1.4;  //make this a param?
        let ni: usize = self.len();
        
        let found = {
            let mut result = None;
            
            for iter in 0..max_iter {
                self.cells[0].phi = 0.0;
      	        self.cells[ni-1].phi = 0.0;

                for i in 1..(ni-1) {
                    let g: f64 = 0.5 * (self.cells[i-1].phi + self.cells[i+1].phi + dx2 * self.cells[i].rho / EPS0);
	            self.cells[i].phi = self.cells[i].phi + w * (g - self.cells[i].phi);
                }
                
                if iter % 50 == 0 {
      	            let mut sum: f64 = 0.0;
	            for i in 1..(ni-1) {
	                let res: f64 = -self.cells[i].rho/EPS0 - (self.cells[i-1].phi - 2.0 * self.cells[i].phi + self.cells[i+1].phi)/dx2;
	                sum += res * res;
	            }
	            l2 = (sum/(ni as f64)).sqrt();
	            if l2 < l2_conv {
	                result = Some(iter);
	                break;
	            }
	        }
            }
            result
        };
        match found {
            Some(i) => { Ok(i) }
            None => {
	        Err(format!("GS SOR didn't converge.  L2 residual {l2:.6}"))
            }
        }
    }

    fn compute_ef(&mut self) -> Result <()> {
        let ni = self.len();
        let denom = 2.0 * self.dx;
        
        for i in 1..(ni - 1) {
            self.cells[i].ef = -(self.cells[i+1].phi - self.cells[i-1].phi) / denom;
        }
        // do 2nd order calc at ends as well
        self.cells[0].ef = (3.0 * self.cells[0].phi
                            - 4.0 * self.cells[1].phi
                            + self.cells[2].phi)
            / denom;
        self.cells[ni - 1].ef = (-self.cells[ni - 3].phi
                                 + 4.0 * self.cells[ni - 2].phi
                                 - 3.0 * self.cells[ni - 1].phi)
            / denom;
        Ok(())
    }
}

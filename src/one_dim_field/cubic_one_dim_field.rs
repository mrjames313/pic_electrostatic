// external
use anyhow::Result;
use std::fs::{OpenOptions};
use std::io::{Write, BufWriter};

// internal
use crate::constants::*;
use super::traits::OneDimFieldInterp;


#[derive(serde::Serialize)]
struct CubicOneDimCell {
    x: f64,
    phi: f64,
    rho: f64,
    ef: f64,
    // now for the cubic approx params
    phi_a: f64,
    phi_b: f64,
    phi_c: f64,
    phi_d: f64,
    rho_a: f64,
    rho_b: f64,
    rho_c: f64,
    rho_d: f64,
    ef_a: f64,
    ef_b: f64,
    ef_c: f64,
    ef_d: f64,
}

#[derive(serde::Serialize)]
pub struct CubicOneDimField {
    x0: f64,
    dx: f64,
    cells: Vec<CubicOneDimCell>
}

// Qs
// Can you have optional arguments?
impl OneDimFieldInterp for CubicOneDimField {
    fn new(n: usize, x0: f64, dx: f64, rho0: f64) -> Result<Self> {
        let cells: Vec<CubicOneDimCell> = (0..n)
            .map(|i| {
                let x = x0 + i as f64 * dx;
                CubicOneDimCell {
                    x,
                    phi: 0.0,
                    rho: rho0,
                    ef: 0.0,
                    phi_a: 0.0,
                    phi_b: 0.0,
                    phi_c: 0.0,
                    phi_d: 0.0,
                    rho_a: 0.0,
                    rho_b: 0.0,
                    rho_c: 0.0,
                    rho_d: 0.0,
                    ef_a: 0.0,
                    ef_b: 0.0,
                    ef_c: 0.0,
                    ef_d: 0.0,
                }
            }).collect();
        let mut field = Self {x0:x0, dx:dx, cells:cells };
        field.solve_potential()?;
        field.compute_ef()?;
        field.compute_cubic_params_natural()?;
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
        self.compute_cubic_params_natural()?;
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

    // This may be just an approximation - possible that some
    // cubic interpolation goes above any knot value
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
        let line = serde_json::to_string(self)?; //??? is self a CubicOneDimField here?
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

impl CubicOneDimField {
    // private methods?  Also, missing lots of error handling right now -
    // need to figure out the right fix for that
    fn x_to_l(&self, x: f64) -> f64 {
        (x - self.x0) / self.dx
    }

    fn gather_ef(&self, li: f64) -> f64 {
        let i : usize = li as usize; // should we check that li is non-negative?
        let dx : f64 = li.fract() * self.dx;
        self.cells[i].ef_a + dx * (self.cells[i].ef_b +
                                   dx * (self.cells[i].ef_c +
                                         dx * self.cells[i].ef_d))
            
    }

    fn gather_phi(&self, li: f64) -> f64 {
        let i : usize = li as usize; // should we check that li is non-negative?
        let dx : f64 = li.fract() * self.dx;
        self.cells[i].phi_a + dx * (self.cells[i].phi_b +
                                   dx * (self.cells[i].phi_c +
                                         dx * self.cells[i].phi_d))
    }

    fn gather_rho(&self, li: f64) -> f64 {
        let i : usize = li as usize; // should we check that li is non-negative?
        let dx : f64 = li.fract() * self.dx;
        self.cells[i].rho_a + dx * (self.cells[i].rho_b +
                                   dx * (self.cells[i].rho_c +
                                         dx * self.cells[i].rho_d))
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

    fn compute_cubic_params_natural(&mut self) -> Result <()> {
        let n = self.len() - 1;
        let mut m_phi = vec![0.0_f64; n+1];
        let mut m_rho = vec![0.0_f64; n+1];
        let mut m_ef = vec![0.0_f64; n+1];

        // assumes n >= 2
        let inv_dx2 = 1.0 / (self.dx * self.dx);
        let interior = n-1;
        let mut cprime_phi = vec![0.0_f64; interior];
        let mut cprime_rho = vec![0.0_f64; interior];
        let mut cprime_ef = vec![0.0_f64; interior];
        let mut dprime_phi = vec![0.0_f64; interior];
        let mut dprime_rho = vec![0.0_f64; interior];
        let mut dprime_ef = vec![0.0_f64; interior];

        let rhs1_phi = 6.0 * inv_dx2 * (self.cells[2].phi - 2.0 * self.cells[1].phi + self.cells[0].phi);
        let rhs1_rho = 6.0 * inv_dx2 * (self.cells[2].rho - 2.0 * self.cells[1].rho + self.cells[0].rho);
        let rhs1_ef = 6.0 * inv_dx2 * (self.cells[2].ef - 2.0 * self.cells[1].ef + self.cells[0].ef);
        let denom1 = 4.0;
        // Are all the cprimes the same???
        cprime_phi[0] = 1.0 / denom1;
        dprime_phi[0] = rhs1_phi / denom1;
        cprime_rho[0] = 1.0 / denom1;
        dprime_rho[0] = rhs1_rho / denom1;
        cprime_ef[0] = 1.0 / denom1;
        dprime_ef[0] = rhs1_ef / denom1;

        for k in 1..interior {
            let i = k+1;
            let rhs_phi = 6.0 * inv_dx2 * (self.cells[i+1].phi - 2.0 * self.cells[i].phi + self.cells[i-1].phi);
            let rhs_rho = 6.0 * inv_dx2 * (self.cells[i+1].rho - 2.0 * self.cells[i].rho + self.cells[i-1].rho);
            let rhs_ef = 6.0 * inv_dx2 * (self.cells[i+1].ef - 2.0 * self.cells[i].ef + self.cells[i-1].ef);

            let denom_phi = 4.0 - 1.0 * cprime_phi[k-1];
            let denom_rho = 4.0 - 1.0 * cprime_rho[k-1];
            let denom_ef = 4.0 - 1.0 * cprime_ef[k-1];

            cprime_phi[k] = 1.0 / denom_phi;
            cprime_rho[k] = 1.0 / denom_rho;
            cprime_ef[k] = 1.0 / denom_ef;

            dprime_phi[k] = (rhs_phi - 1.0 * dprime_phi[k-1]) / denom_phi;
            dprime_rho[k] = (rhs_rho - 1.0 * dprime_rho[k-1]) / denom_rho;
            dprime_ef[k] = (rhs_ef - 1.0 * dprime_ef[k-1]) / denom_ef;
        }

        m_phi[n-1] = dprime_phi[interior - 1];
        m_rho[n-1] = dprime_rho[interior - 1];
        m_ef[n-1] = dprime_ef[interior - 1];

        for k in (0..interior - 1).rev() {
            m_phi[k+1] = dprime_phi[k] - cprime_phi[k] * m_phi[k+2];
            m_rho[k+1] = dprime_rho[k] - cprime_rho[k] * m_rho[k+2];
            m_ef[k+1] = dprime_ef[k] - cprime_ef[k] * m_ef[k+2];
            
        }

        // The last point doesn't get cubic params, if I'm understanding correctly
        for i in 0..n {
            self.cells[i].phi_a = self.cells[i].phi; // so... why have both if the first param is the same as the value
            self.cells[i].phi_b = (self.cells[i+1].phi - self.cells[i].phi) / self.dx -
                (self.dx / 6.0) * (2.0 * m_phi[i] + m_phi[i + 1]);
            self.cells[i].phi_c = 0.5 * m_phi[i];
            self.cells[i].phi_d = (m_phi[i+1] - m_phi[i]) / (6.0 * self.dx);
            self.cells[i].rho_a = self.cells[i].rho; // so... why have both if the first param is the same as the value
            self.cells[i].rho_b = (self.cells[i+1].rho - self.cells[i].rho) / self.dx -
                (self.dx / 6.0) * (2.0 * m_rho[i] + m_rho[i + 1]);
            self.cells[i].rho_c = 0.5 * m_rho[i];
            self.cells[i].rho_d = (m_rho[i+1] - m_rho[i]) / (6.0 * self.dx);
            self.cells[i].ef_a = self.cells[i].ef; // so... why have both if the first param is the same as the value
            self.cells[i].ef_b = (self.cells[i+1].ef - self.cells[i].ef) / self.dx -
                (self.dx / 6.0) * (2.0 * m_ef[i] + m_ef[i + 1]);
            self.cells[i].ef_c = 0.5 * m_ef[i];
            self.cells[i].ef_d = (m_ef[i+1] - m_ef[i]) / (6.0 * self.dx);
        }
        Ok(())
    }
}

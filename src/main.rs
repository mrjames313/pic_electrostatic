// use serde::Serialize;
use std::fs::{OpenOptions};
use std::io::{Write, BufWriter};
use anyhow::Result;


const QE: f64 = 1.602176565e-19;  // C, electron charge
const EPS0: f64 = 8.85418782e-12; // C/V/m vac permiability
const ME: f64 = 9.10938215e-31;   // kg electron mass


#[derive(serde::Serialize)]
struct Cell {
    x: f64,
    phi: f64,
    rho: f64,
    ef: f64
}

#[derive(serde::Serialize)]
struct Field {
    x0: f64,
    dx: f64,
    cells: Vec<Cell>
}

// Qs
// Can you have optional arguments?
impl Field {
    fn new(n: usize, x0: f64, dx: f64, rho0: f64) -> Self {
        let cells: Vec<Cell> = (0..n)
            .map(|i| {
                let x = x0 + i as f64 * dx;
                Cell {
                    x,
                    phi: 0.0,
                    rho: rho0,
                    ef: 0.0
                }
            }).collect();
        Self {x0:x0, dx:dx, cells:cells }
    }

    fn len(&self) -> usize {self.cells.len() }

    // Assumes that x values remain unmodified
    fn reset(&mut self, rho0: f64) {
        for c in &mut self.cells {
            c.phi = 0.0;
            c.rho = rho0;
            c.ef = 0.0;
        }
    }
}

// Will delete the contents of the file if it already exists
fn write_field_to_json(field : &Field, filename : &str) -> Result<()> {
   let file = OpenOptions::new()
       .write(true)
       .create(true)
       .truncate(true)
       .open(filename)?;
       
   let mut writer = BufWriter::new(file);
   let line = serde_json::to_string(field)?;
   writeln!(writer, "{}", line)?;

   writer.flush()?;
   Ok(())

}


fn print_field(field :&Field) -> Result<()> {
    for cell in field.cells.iter() {
        println!("X_i {:.6}, phi_i {:.3}, rho_i {:.3}, ef_i {:.3}",
                  cell.x, cell.phi, cell.rho, cell.ef);
    }
    Ok(())
}


fn solve_potential_direct(field : &mut Field) -> Result<()> {
   let ni: usize = field.len();
   let mut a: Vec<f64> = vec![0.0; ni];  // coef phi[i-1]
   let mut b: Vec<f64> = vec![0.0; ni];  // coef phi[i]
   let mut c: Vec<f64> = vec![0.0; ni];  // coef phi[i+1]
   let mut d: Vec<f64> = vec![0.0; ni];  // rhs

   // precompute some common values
   let inv_sq = 1.0 / (field.dx * field.dx);
   let two_inv_sq = -2.0 / (field.dx * field.dx);
   
   for i in 0..ni {
      if i == 0 || i == ni-1 {
      	 b[i] = 1.0;
	 d[i] = 0.0;
      } else {
         a[i] = inv_sq;
	 b[i] = two_inv_sq;
	 c[i] = inv_sq;
	 d[i] = -field.cells[i].rho / EPS0;
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

   field.cells[ni-1].phi = d[ni-1];
   for i in (0..ni-1).rev() {
      field.cells[i].phi = d[i] - c[i] * field.cells[i+1].phi;
   }
   
   Ok(())   
}

// want this to return a value indicating successful convergence
fn solve_potential_gs_sor(field : &mut Field, max_iter: i32) -> Result<i32, String> {
   let mut L2: f64 = 1e12;
   let L2_conv: f64 = 1e-6;
   let dx2: f64 = field.dx * field.dx;
   let w: f64 = 1.4;  //make this a param?
   let ni: usize = field.len();

   let found = {
      let mut result = None;
      
      for iter in 0..max_iter {
         field.cells[0].phi = 0.0;
      	 field.cells[ni-1].phi = 0.0;

         for i in 1..(ni-1) {
            let g: f64 = 0.5 * (field.cells[i-1].phi + field.cells[i+1].phi + dx2 * field.cells[i].rho / EPS0);
	    field.cells[i].phi = field.cells[i].phi + w * (g - field.cells[i].phi);
         }

         if iter % 50 == 0 {
      	    let mut sum: f64 = 0.0;
	    for i in 1..(ni-1) {
	       let res: f64 = -field.cells[i].rho/EPS0 - (field.cells[i-1].phi - 2.0 * field.cells[i].phi + field.cells[i+1].phi)/dx2;
	       sum += res * res;
	    }
	    L2 = (sum/(ni as f64)).sqrt();
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

fn compute_EF(field : &mut Field) -> Result <()> {
    let ni = field.len();
    let denom = 2.0 * field.dx;
    
    for i in (1..ni - 1) {
        field.cells[i].ef = -(field.cells[i+1].phi - field.cells[i-1].phi) / denom;
    }
    // do 2nd order calc at ends as well
    field.cells[0].ef = (3.0 * field.cells[0].phi
                         - 4.0 * field.cells[1].phi
                         + field.cells[2].phi)
        / denom;
    field.cells[ni - 1].ef = (-field.cells[ni - 3].phi
                              + 4.0 * field.cells[ni - 2].phi
                              - 3.0 * field.cells[ni - 1].phi)
        / denom;
    Ok(())
}
    

fn main() -> Result<()> {
    println!("Starting execution");
    let ni = 21;
    let x0: f64 = 0.0;
    let xm: f64 = 0.1;
    let dx: f64 = (xm - x0) / (ni - 1) as f64;

    let mut field = Field::new(ni, x0, dx, QE*1e12); 

    solve_potential_direct(&mut field)?;
    print_field(&field);

    field.reset(QE*1e12);
    
    match solve_potential_gs_sor(&mut field, 1000) {
       Ok(i) => {
           print_field(&field);
       }
       Err(s) => println!("{s}")
    }

    compute_EF(&mut field)?;
    print_field(&field);
    
    write_field_to_json(&field, "output.json")?;

    println!("Done writing files");
    Ok(())
    
    
}

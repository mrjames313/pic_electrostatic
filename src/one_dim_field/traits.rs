use anyhow::Result;

pub trait OneDimFieldInterp {

    fn new(n: usize, x0: f64, dx: f64, rho0: f64) -> Result<Self>
        where Self: Sized;
    
    fn reset(&mut self, rho0: f64)-> Result<()>;
    
    fn phi_at(&self, x: f64) -> f64;
    fn rho_at(&self, x: f64) -> f64;
    fn ef_at(&self, x: f64) -> f64;

    fn get_phi_max(&self) -> f64;
    
    fn len(&self) -> usize;
    
    fn x0(&self) -> f64;
    fn dx(&self) -> f64;

    fn write_field_to_json(&self, filename : &str) -> Result<()>;
    fn print(&self);
    
}


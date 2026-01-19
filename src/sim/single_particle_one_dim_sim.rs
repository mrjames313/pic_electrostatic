use crate::constants::*;
use crate::one_dim_field::OneDimFieldInterp;

pub struct SingleParticleOneDimTrace {
    pub ts_vals: Vec<usize>,
    pub x_vals: Vec<f64>,
    pub v_vals: Vec<f64>,
    pub phi_vals: Vec<f64>,
    pub ke_vals: Vec<f64>,
    pub pe_vals: Vec<f64>,
    pub te_vals: Vec<f64>,
    pub min_x: f64,
    pub max_x: f64,
    pub min_e: f64,
    pub max_e: f64,    
}

impl SingleParticleOneDimTrace {
    pub fn with_capacity(steps: usize) -> Self {
        Self {
            ts_vals: Vec::with_capacity(steps),
            x_vals: Vec::with_capacity(steps),
            v_vals: Vec::with_capacity(steps),
            phi_vals: Vec::with_capacity(steps),
            ke_vals: Vec::with_capacity(steps),
            pe_vals: Vec::with_capacity(steps),
            te_vals: Vec::with_capacity(steps),
            min_x: f64::INFINITY,
            max_x: f64::NEG_INFINITY,
            min_e: f64::INFINITY,
            max_e: f64::NEG_INFINITY,
        }
    }

    pub fn update_trace(&mut self, ts: usize, x: f64, v: f64, phi: f64,
                        ke: f64, pe: f64, te: f64)
    {
        if x < self.min_x {self.min_x = x; }
        if x > self.max_x {self.max_x = x; }
        if te < self.min_e {self.min_e = te; }
        if te > self.max_e {self.max_e = te; }

        self.ts_vals.push(ts);
        self.x_vals.push(x);
        self.v_vals.push(v);
        self.phi_vals.push(phi);
        self.ke_vals.push(ke);
        self.pe_vals.push(pe);
        self.te_vals.push(te);
    }

}

pub fn simulate<F>(field: &F, x_init: f64, v_init: f64,
                   dt: f64, steps: usize) -> SingleParticleOneDimTrace
where F: crate::one_dim_field::OneDimFieldInterp
{

    let mut x = x_init;
    let mut v = v_init;
    let mut x_old;
    let phi_max = field.get_phi_max();
    let mut trace = SingleParticleOneDimTrace::with_capacity(steps);

    // leapfrog half-step
    let mut ef_p = field.ef_at(x);
    v -= 0.5 * (-QE / ME) * ef_p * dt;

    for ts in 0..steps {
        ef_p = field.ef_at(x);
        x_old = x;
        v += (-QE / ME) * ef_p * dt;
        x += v * dt;

        let x_inter : f64 = 0.5 * (x + x_old); // interpolate x to same time offset as v

        // now get some other physics data
        let phi_p = field.phi_at(x_inter);

        let ke : f64 = 0.5 * ME * v * v / QE;
        let pe : f64 = (-QE) * (phi_p - phi_max) / QE; // This is weird, QE and /QE?
        let te : f64 = ke + pe;

        trace.update_trace(ts, x, v, phi_p, ke, pe, te);
        if ts % 50 == 0 {
            println!("{ts}, {x_inter:.6}, {v:.4}, {phi_p:.4}, {ke:.5}, {pe:.5}, {te:.5}");
        }
    }
    trace
}


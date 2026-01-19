// external
use anyhow::Result;
use plotly::{Plot, Scatter};
use plotly::common::Mode;

// internal
use pic_electrostatic::constants::*;
use pic_electrostatic::one_dim_field::{OneDimFieldInterp,LinearOneDimField};

    
fn main() -> Result<()> {
    println!("Starting execution");
    let ni = 21;
    let x0: f64 = 0.0;
    let xm: f64 = 0.1;
    let dx: f64 = (xm - x0) / (ni - 1) as f64;

    // ****************   Create the field from potentials *************
    
    let mut field = match LinearOneDimField::new(ni, x0, dx, QE*1e12) {
        Ok(s) => s,
        Err(_) => {
            println!("Failed to create a viable field");
            return Err(anyhow::anyhow!("badness"));
        }
    };

    // Just demonstrating a reset of Rho field - not necessary
//    match field.reset(QE*1e12) {
//        Ok(()) => {},
//        Err(_) => {
//            println!("Field became invalid after resetting values");
//            return Err(anyhow::anyhow!("badness"));
//        }
//    }

    field.print();
    field.write_field_to_json("output.json")?;

    
    // ********  Now simulate a particle (electron) in the field ********
    let mut x : f64 = 4.0 * dx; // start 4 cells in
    let mut v : f64 = 0.0;
    let mut x_old : f64;
    let dt : f64 = 1e-10;

    // Get phi_max
    let phi_max : f64 = field.get_phi_max();
    
    // offset velocity by half timestep
    let mut ef_p : f64 = field.ef_at(x);
    v -= 0.5 * (-QE / ME) * ef_p * dt;

    // keep some stats
    let mut min_x : f64 = f64::INFINITY;
    let mut max_x : f64 = f64::NEG_INFINITY;
    let mut min_energy : f64 = f64::INFINITY;
    let mut max_energy : f64 = f64::NEG_INFINITY;

    let loop_iters : usize = 1201;
    //let loop_iters : usize = 201;

    let mut ts_vals = Vec::with_capacity(loop_iters);
    let mut x_vals = Vec::with_capacity(loop_iters);
    let mut v_vals = Vec::with_capacity(loop_iters);
    let mut phi_vals = Vec::with_capacity(loop_iters);
    let mut ke_vals = Vec::with_capacity(loop_iters);
    let mut pe_vals = Vec::with_capacity(loop_iters);
    let mut te_vals = Vec::with_capacity(loop_iters);
    
    for ts in 0..loop_iters {
        ef_p = field.ef_at(x);

        x_old = x;
        v += (-QE / ME) * ef_p * dt;
        x += v * dt;

        let x_inter : f64 = 0.5 * (x + x_old); // interpolate x to same time offset as v

        // now get some other physics data
        let phi_p = field.phi_at(x_inter);

        let ke : f64 = 0.5 * ME * v * v / QE;
        let pe : f64 = (-QE) * (phi_p - phi_max) / QE; // This is weird, QE and /QE?
        let total_e : f64 = ke + pe;
        
        if x < min_x {min_x = x; }
        if x > max_x {max_x = x; }
        if total_e < min_energy {min_energy = total_e; }
        if total_e > max_energy {max_energy = total_e; }

        ts_vals.push(ts);
        x_vals.push(x_inter);
        v_vals.push(v);
        phi_vals.push(phi_p);
        ke_vals.push(ke);
        pe_vals.push(pe);
        te_vals.push(total_e);
        
        if ts % 50 == 0 {
            println!("{ts}, {x_inter:.4}, {v:.2}, {phi_p:.2}, {ke:.4}, {pe:.4}, {total_e:.4}");
        }
    }

    // make some plots to debug energy
    let trace = Scatter::new(x_vals, phi_vals)
        .mode(Mode::Lines)
        .name("total_e(x)");
    let mut plot = Plot::new();
    plot.add_trace(trace);
    //plot.show();
    let path = "/home/ubuntu/plots/plotly_tmp.html";
    plot.write_html(path);
    
    println!("Min x {min_x}, max x {max_x}, min total energy {min_energy}, max total energy {max_energy}");
    
    println!("Done writing files");
    Ok(())
    
    
}

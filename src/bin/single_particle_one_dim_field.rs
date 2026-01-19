// external
use anyhow::Result;

// internal
use pic_electrostatic::constants::*;
use pic_electrostatic::one_dim_field::{OneDimFieldInterp,LinearOneDimField,CubicOneDimField};
use pic_electrostatic::sim::{SingleParticleOneDimTrace, simulate};
use pic_electrostatic::plot::plot_line_html;

fn main() -> Result<()> {
    println!("Starting execution");

    // ************ Set up the parameters for the field ******************
    let ni = 21;
    let x0: f64 = 0.0;
    let xm: f64 = 0.1;
    let dx: f64 = (xm - x0) / (ni - 1) as f64;

    // ****************   Create the field from potentials *************
    //let mut field = match LinearOneDimField::new(ni, x0, dx, QE*1e12) {
    let mut field = match CubicOneDimField::new(ni, x0, dx, QE*1e12) {
        Ok(s) => s,
        Err(_) => {
            println!("Failed to create a viable field");
            return Err(anyhow::anyhow!("badness"));
        }
    };

    // Demonstrating a reset of Rho field - not necessary
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
    let x_init : f64 = 4.0 * dx; // start 4 cells in
    let v_init : f64 = 0.0;
    let dt : f64 = 1e-10;
    let loop_iters : usize = 120001;
    //let loop_iters : usize = 401;

    // looks like the type gets automatically inferred...
    let sim_trace : SingleParticleOneDimTrace = simulate(
        &field, x_init, v_init, dt, loop_iters);

    // *********** make plots and print stats ******************
    let path = "/home/ubuntu/plots/plotly_tmp.html";
    let name = "total_e(x)";
    plot_line_html(sim_trace.x_vals, sim_trace.te_vals, name, path);
    
    println!("Min x {}, max x {}, min total energy {}, max total energy {}",
             sim_trace.min_x, sim_trace.max_x,
             sim_trace.min_e, sim_trace.max_e);
    
    println!("Done writing files");
    Ok(())
    
    
}

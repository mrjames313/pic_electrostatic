use anyhow::Result;

fn main() -> Result <()> {

    let m: f64 = 9.10938215e-31;
    let q: f64 = 1.602176565e-19;
    let mut x: f64 = 0.0;
    let mut v: f64 = 0.0;

    let dt: f64 = 1.0e-9;
    let E: f64 = -100.0; // electric field


    // leapfrog method - back v up by 0.5 dt
    v -= 0.5 * (q / m) * E * dt;

    for i in (0..10) {
        // to be more correct, would interpolate v to next value
        // before printing
        println!("{:.9}: {x:.6},  {v:.3}", i as f64 * dt);
        v += (q/m) * E * dt;
        x += v * dt;
    }
    Ok(())
}

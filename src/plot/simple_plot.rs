use plotly::{common::Mode, Plot, Scatter};

pub fn plot_line_html(x: Vec<f64>, y: Vec<f64>, name: &str, path: &str) {
    let trace = Scatter::new(x, y)
        .mode(Mode::Lines)
        .name(name);
    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.write_html(path);
}

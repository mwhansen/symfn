use symfn::{
    clear_caches,
    skew_lr::{expand_skew, take_peak_layer_states},
    Partition,
};
fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}
fn main() {
    let which = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "0".into())
        .parse::<usize>()
        .unwrap();
    let shapes: [&[u32]; 4] = [
        &[10, 8, 6, 4],
        &[12, 10, 8, 6],
        &[8, 7, 6, 5, 4, 3],
        &[16, 13, 10, 7],
    ];
    let mu = p(shapes[which]);
    let w = mu.part(0);
    let outer: Vec<u32> = mu
        .parts()
        .iter()
        .map(|x| x + w)
        .chain(mu.parts().iter().copied())
        .collect();
    let inner: Vec<u32> = mu
        .parts()
        .iter()
        .map(|_| w)
        .chain(std::iter::repeat_n(0, mu.len()))
        .collect();
    clear_caches();
    let _ = take_peak_layer_states();
    let r = expand_skew(&p(&outer), &Partition::new(inner));
    eprintln!(
        "{mu}^2\t{} terms\tpeak_states {}",
        r.len(),
        take_peak_layer_states()
    );
}

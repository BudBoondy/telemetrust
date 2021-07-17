
pub trait MeanFilter: Statistics {
    type Output;
    type Inner;
    fn hampel_filter(&mut self) -> Self::Output;
}

pub trait Statistics {
    type Inner;
    fn average(&self) -> Self::Inner;
    fn median(&self) -> Self::Inner;
    fn variance(&self) -> Self::Inner;
}
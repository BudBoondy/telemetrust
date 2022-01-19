pub use crate::stats::Statistics;
pub use crate::stats::MeanFilter;

impl Statistics for Vec<f64> {
    type Inner = f64;
    fn average(&self) -> f64{
        if self.is_empty(){
            return 0.0
        }
        let self_clone = self.clone();
        self_clone.into_iter().sum::<f64>() / (self.len() as f64)
    }

    fn median(&self) -> f64{
        match self.len(){
            0 => return 0.0,
            1 => return self[0],
            _ => ()
        }
        let mut self_clone = self.clone();
        self_clone.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if self_clone.len() % 2 != 0{
            self_clone[self_clone.len()/2 + 1]
        } else {
            (self_clone[(self_clone.len()/2)-1] + self_clone[(self_clone.len()/2)])/2.0
        }
    }

    fn variance(&self) -> f64 {
        if self.is_empty(){
            return 0.0
        }
        let avg = self.average();
        let variance = self.iter().map(|x| x * x).sum::<f64>() / self.len() as f64 - (avg * avg);
        variance
    }
}

impl MeanFilter for Vec<f64>{
    type Output = Vec<f64>;
    type Inner = f64;
    fn hampel_filter(&mut self) -> Vec<f64>{
        let gaussian_scale = 1.4826;
        let median = self.median();
        let mad = gaussian_scale * self.clone().into_iter().map(|x| (x - median).abs()).collect::<Vec<f64>>().median();
        
        for i in 0..self.len(){
            if (self[i] - median) > 3.0*mad{
                self[i] = median;
            }
        }
        self.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_hampel_vec_f64(){
        let vec_f64_odd = vec![2.0, 5.5, 7.4, 6.9, 10.2, 3.3, 8.8];
        assert_eq!((vec_f64_odd.average()*10.0).round()/10.0, 6.3);
        assert_eq!(vec_f64_odd.median(), 7.4);
        let vec_f64_even = vec![2.0, 5.5, 7.4, 6.9, 10.2, 3.3, 8.8, 1.2];

        assert_eq!((vec_f64_even.average()*10.0).round()/10.0, 5.7);
        assert_eq!(vec_f64_even.median(), 6.2);

        let vec_no_elements: Vec<f64> = vec![];
        assert_eq!((vec_no_elements.average()*10.0).round()/10.0, 0.0);
        assert_eq!(vec_no_elements.median(), 0.0);

        let vec_one_element = vec![2.5];
        assert_eq!((vec_one_element.average()*10.0).round()/10.0, 2.5);
        assert_eq!(vec_one_element.median(), 2.5);
    }

    #[test]
    fn test_variance(){
        let vec_f64 = vec![1.0, 3.0, 5.0, 9.0, 12.0];
        assert_eq!(16.0, vec_f64.variance())
    }

    #[test]
    fn test_hampel_filter(){
        let mut vec_hampel = vec![1.0, 3.0, 5.0, 9.0, 12.0, 5.0, 11.0, 50.0, 20.0, 17.5, 70.0, 12.5];
        // assert median
        let median = vec_hampel.median();
        // assert two elements (50.0, and 70.0)
        vec_hampel.hampel_filter();
        assert_eq!(vec_hampel[7], median);
        assert_eq!(vec_hampel[10], median);
    }
}
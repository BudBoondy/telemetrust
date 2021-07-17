// MLE method
pub struct MLE;

use super::*;
use rulinalg::matrix;
use rulinalg::matrix::Matrix;
use rulinalg::vector::Vector;

impl TriangulationFunction for MLE {
    fn triangulate(b: Vec<Bearing>) -> Result<Location, EqualAngleError> {
        let y = b.iter().map(|y| y.loc.lon).collect::<Vec<f64>>();
        let x = b.iter().map(|x| x.loc.lat).collect::<Vec<f64>>();
        let bearings  = b.iter().map(|x| x.angle).collect::<Vec<f64>>();
        let theta: Vec<f64> = bearings.iter().map(|br| std::f64::consts::PI/180.0*(90.0-br)).collect();
        let si: Vec<f64> = theta.iter().map(|t| t.sin()).collect();
        let ci: Vec<f64> = theta.iter().map(|t| t.cos()).collect();
        let xi = x;
        let yi = y;
        let zi = subtract_vec(&dot_vec(&si, &xi), &dot_vec(&ci, &yi));
        let first: f64 = si.iter().map(|i| i * i).sum();
        let second_first: f64 = -1_f64 * si.iter().zip(ci.iter()).map(|(s,c)| s * c).sum::<f64>();
        let third_second: f64 = ci.iter().map(|i| i * i).sum();  
        let a =  Matrix::new(2, 2, vec![
            first, second_first,
            second_first, third_second
        ]);
        let b = Vector::new(vec![si.iter().zip(zi.iter()).map(|(s, z)| s * z).sum(), -1_f64 * ci.iter().zip(zi.iter()).map(|(c, z)| c * z).sum::<f64>()]); 
        let xy = match a.solve(b){
            Ok(x) => x,
            Err(e) => {println!("{:?}", e); return Err(EqualAngleError)} 
        };
        let mut x0 = xy[0];
        let mut y0 = xy[1];
        let mut x: f64 = 0.0; 
        let mut y: f64 = 0.0; 
        let mut counter = 0;

        while(x.abs()-x0.abs()).abs() > 0.0001 && (y.abs()-y0.abs()).abs() > 0.0001 && counter <= 999{
            if x != 0.0 && y != 0.0{
                x0 = x;
                y0 = y;
            } else { 
                x = x0;
                y = y0;
            }
            let di: Vec<f64> = (add_vec(&sub(&xi,x).iter().map(|val: &f64| val * val).collect::<Vec<f64>>(), &sub(&yi, y).iter().map(|val: &f64| val * val).collect::<Vec<f64>>())).iter().map(|val| val.sqrt()).collect::<Vec<f64>>();
            let di_pow3 = di.iter().map(|v| v * v * v).collect::<Vec<f64>>();
            let s_star = div_vec(&sub(&yi, y), &di_pow3);
            let c_star = div_vec(&sub(&xi, x), &di_pow3);
            
            let ffirst: f64 = si.iter().zip(s_star.iter()).map(|(i,j)| i * j).sum();
            let fsecond: f64 = -1_f64 * si.iter().zip(c_star.iter()).map(|(i,j)| i * j).sum::<f64>();

            let sfirst: f64 = -1_f64 * ci.iter().zip(s_star.iter()).map(|(s,c)| s * c).sum::<f64>();
            let ssecond: f64 = ci.iter().zip(c_star.iter()).map(|(s,c)| s * c).sum::<f64>();
            let q =  matrix![
                round_to_one(ffirst), round_to_one(sfirst);
                round_to_one(fsecond), round_to_one(ssecond)
            ];
            
            let p = Vector::new(vec![
                round_to_one(s_star.iter().zip(zi.iter()).map(|(s, z)| s * z).sum::<f64>()), 
                round_to_one(-1_f64 * c_star.iter().zip(zi.iter()).map(|(c, z)| c * z).sum::<f64>())]); 
            
            let xyn = match q.solve(p){
                Ok(n) => n,
                Err(e) => {println!("{:?}", e); return Err(EqualAngleError)} 
            };
            
            x = xyn[0];
            y = xyn[1];      
            
            counter += 1;
            
            if x.is_nan() || y.is_nan() || counter == 999{ 
                return Err(EqualAngleError)
            }
        }
        Ok(Location{lat: x, lon: y})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_triangulation(){
        // first test
        let base_sample_one: Vec<Bearing> = vec![
            Bearing{loc: Location{lat: -103.099, lon: 33.774}, name: "".to_owned(), angle: 322.0},
            Bearing{loc: Location{lat: -103.098, lon: 33.774}, name: "".to_owned(), angle: 296.0},
            Bearing{loc: Location{lat: -103.098, lon: 33.775}, name: "".to_owned(), angle: 290.0},
        ];
        assert_eq!(Location{lat: -103.10058, lon: 33.77577}, triangulate::<MLE>(base_sample_one).unwrap());

        // second test
        let base_sample_two: Vec<Bearing> = vec![
            Bearing{loc: Location{lat: -103.080, lon: 33.795}, name: "".to_owned(), angle: 69.0},
            Bearing{loc: Location{lat: -103.080, lon: 33.796}, name: "".to_owned(), angle: 95.0},
            Bearing{loc: Location{lat: -103.080, lon: 33.795}, name: "".to_owned(), angle: 51.0},
        ];
        assert_eq!(Location{lat: -103.07850, lon: 33.79587}, triangulate::<MLE>(base_sample_two).unwrap());

        // third test
        let base_sample_three: Vec<Bearing> = vec![
            Bearing{loc: Location{lat: -99.874, lon: 36.506}, name: "".to_owned(), angle: 137.0},
            Bearing{loc: Location{lat: -99.874, lon: 36.500}, name: "".to_owned(), angle: 51.0},
            Bearing{loc: Location{lat: -99.867, lon: 36.5077}, name: "".to_owned(), angle: 215.0},
        ];
        assert_eq!(Location{lat: -99.87076, lon: 36.50256}, triangulate::<MLE>(base_sample_three).unwrap());

        // fourth test
        let base_sample_four: Vec<Bearing> = vec![
            Bearing{loc: Location{lat: -99.873, lon: 36.5077}, name: "".to_owned(), angle: 97.0},
            Bearing{loc: Location{lat: -99.867, lon: 36.5077}, name: "".to_owned(), angle: 220.0},
            Bearing{loc: Location{lat: -99.874, lon: 36.5042}, name: "".to_owned(), angle: 50.0},
        ];
        assert_eq!(Location{lat: -99.86714, lon: 36.50754}, triangulate::<MLE>(base_sample_four).unwrap());
    }
}
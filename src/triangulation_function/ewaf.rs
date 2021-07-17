// algorithm based on 
//  https://www.movable-type.co.uk/scripts/latlong.html
//  http://www.edwilliams.org/avform147.htm#Intersection

use itertools::Itertools;
pub struct EWAF;

use super::*;

pub fn great_circle(b: Bearing) -> Vec<f64> {
   let lat = b.loc.lat.to_radians(); 
   let lng = b.loc.lon.to_radians();
   let angle = b.angle.to_radians();
   let great_circle = vec![ 
        lng.sin() * angle.cos() - lat.sin() * lng.cos() * angle.sin(),
        -1.0 * lng.cos() * angle.cos() - lat.sin() * lng.sin() * angle.sin(),
        lat.cos() * angle.sin()
    ];
   great_circle
}

fn triangulate_bearings_pair(b1: &Bearing, b2: &Bearing) -> Result<Location, EqualAngleError> {
    if (b1.angle - b2.angle).abs() < 0.05{
        return Err(super::EqualAngleError)
    }
    let n1 = super::to_normalized_vector(b1.loc); 
    let angle1 = b1.angle.to_radians();
    let n2 = super::to_normalized_vector(b2.loc); 
    let angle2 = b2.angle.to_radians();
    // north pole
    let n = vec!(0.0, 0.0, 1.0);                    

    let de1 = super::unit(super::cross_prod(&n, &n1));            
    let dn1 = super::cross_prod(&n1, &de1);
    let de1_sin0 = de1.iter().map(|x| x * angle1.sin()).collect::<Vec<f64>>();
    let dn1_cos0 = dn1.iter().map(|x| x * angle1.cos()).collect::<Vec<f64>>();
    let d1 = dn1_cos0.iter()
        .zip(de1_sin0.iter())
        .map(|(a, b)| a + b)
        .collect::<Vec<f64>>();

    let c1 = cross_prod(&n1, &d1);                  
    let de2 = unit(cross_prod(&n, &n2));
    let dn2 = cross_prod(&n2, &de2);    
    
    
    let de2_sin0 = de2.iter().map(|x| x * angle2.sin()).collect::<Vec<f64>>();
    let dn2_cos0 = dn2.iter().map(|x| x * angle2.cos()).collect::<Vec<f64>>();
    
    let d2 = dn2_cos0.iter()
        .zip(de2_sin0.iter())
        .map(|(a, b)| a + b)
        .collect::<Vec<f64>>();   
    let c2 = cross_prod(&n2, &d2);
    let ni = cross_prod(&c1, &c2);
    Ok(to_lat_lng(vec!(ni[0], ni[1], ni[2])))
}

impl TriangulationFunction for EWAF{
    fn triangulate(bearings: Vec<Bearing>) -> Result<Location, EqualAngleError> {
        let location_vec: Vec<Location> = bearings.into_iter()
        .combinations(2)
        .into_iter()
        .filter_map(|x| {
            let ret_val;
            let val_first = || -> Option<Location> {triangulate_bearings_pair(&x[0], &x[1]).ok()};
            let val_second = || -> Option<Location> {triangulate_bearings_pair(&x[1], &x[0]).ok()};

            ret_val = match val_first(){
                Some(l) if l.distance_to(&x[0].loc) < 1000.0 => l,
                _ => match val_second(){
                    Some(s) if s.distance_to(&x[0].loc) < 1000.0 => s,
                    _ => return None,
                },
            };
            Some(ret_val)
        })
        .collect();

        if location_vec.is_empty(){
            return Err(EqualAngleError)
        }
        let location = fold_location_vec(location_vec);
        Ok(location)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_triangulation(){
        assert_eq!(Location{lat: 51.1297000, lon: 1.3214000}.to_string_precision(4), 
            triangulate_bearings_pair(&Bearing{loc: Location{lat: 50.7175, lon: 1.65139}, angle: 333.3508, name: "irrelevant".to_owned()}, &Bearing{loc: Location{lat: 50.9250, lon: 1.7094}, angle: 310.1414, name: "irrelevant".to_owned()})
            .unwrap()
            .to_string_precision(4));
    }
    
    #[test]
    fn test_real_triangulation(){
        assert_eq!(Location{lat: 50.824255556259786, lon: 8.784189999999997}.to_string_precision(6), 
            triangulate_bearings_pair(&Bearing{loc: Location{lat: 50.814781, lon: 8.769190}, angle: 45.0, name: "irrelevant".to_owned()}, &Bearing{loc: Location{lat: 50.814781, lon: 8.799190}, angle: 315.0, name: "irrelevant".to_owned()})
        .unwrap()
        .to_string_precision(6))
    }
    #[test]
    fn test_real_triangulation_inf(){
        let result = EWAF::triangulate(vec![Bearing{loc: Location{lat: 50.814781, lon: 8.769190}, angle: 45.0, name: "irrelevant".to_owned()}, Bearing{loc: Location{lat: 50.814781, lon: 8.799190}, angle: 45.0, name: "irrelevant".to_owned()}]);
        match result{
            Ok(_) => panic!("did not expect OK result"),
            Err(_) => assert_eq!(true, true)
        }
        let expected = Err(EqualAngleError);
        assert_eq!(expected, result);
    }

    #[test]
    fn test_great_circle(){
        assert_eq!(vec!("-0.794".to_owned(),"0.129".to_owned(),"0.594".to_owned()), great_circle(Bearing{loc: Location{lat: 53.3206, lon: -1.7297}, angle: 96.0, name: "irrelevant".to_owned()})
            .iter()
            .map(|x| format!("{:.3}", x))
            .collect::<Vec<String>>());
    }
    
    #[test]
    fn test_distance_real(){
        assert_eq!(format!("{:.1}", 1490.12), format!("{:.1}", Location{lat: 50.814781, lon: 8.769190}.distance_to(&Location{lat: 50.824255556259786, lon: 8.784189999999997})))
    }

    #[test]
    fn test_distance(){
        assert_eq!(format!("{:.1}", 662700.0/1000.0), format!("{:.1}", Location{lat: 52.205, lon: 0.119}.distance_to(&Location{lat: 58.0257, lon: 2.351})/1000.0))
    }

    #[test]
    fn test_to_lat_lang(){
        assert_eq!(Location{lat: 45.0000000,lon: 045.0000000}.to_string() ,to_lat_lng(vec!(0.5000, 0.5000, 0.7071)).to_string())
    }

    #[test]
    fn test_vector_folding(){
        assert_eq!(Location{lat: 5.0, lon: 4.0}, fold_location_vec(vec![Location{lat: 8.52, lon: 7.1}, Location{lat: 1.48, lon: 0.9}]));
    }
}
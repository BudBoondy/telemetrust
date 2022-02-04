// Based on https://github.com/henry-dang/triangulation
pub struct Huber;

use super::*;
use rulinalg::matrix::Matrix;
use rulinalg::vector::Vector;

impl TriangulationFunction for Huber {
    fn triangulate(b: Vec<Bearing>) -> Result<Location, TriangulationError> {
        let y = b.iter().map(|y| y.loc.lon).collect::<Vec<f64>>();
        let x = b.iter().map(|x| x.loc.lat).collect::<Vec<f64>>();
        let bearings = b.iter().map(|x| x.angle).collect::<Vec<f64>>();
        let theta: Vec<f64> = bearings
            .iter()
            .map(|br| std::f64::consts::PI / 180.0 * (90.0 - br))
            .collect();
        let si: Vec<f64> = theta.iter().map(|t| t.sin()).collect();
        let ci: Vec<f64> = theta.iter().map(|t| t.cos()).collect();
        let xi = x;
        let yi = y;
        let zi = subtract_vec(&dot_vec(&si, &xi), &dot_vec(&ci, &yi));
        let mut wi = vec![1.0, 1.0, 1.0];
        let first = wi.iter().zip(si.iter()).map(|(w, s)| w * s.powi(2)).sum();
        let second = -1_f64
            * wi.iter()
                .zip(si.iter().zip(ci.iter()))
                .map(|(w, (s, c))| w * s * c)
                .sum::<f64>();
        let fourth = wi.iter().zip(ci.iter()).map(|(w, c)| w * c.powi(2)).sum(); // clone
        let a = Matrix::new(
            2,
            2,
            vec![
                round_to_one(first),
                round_to_one(second),
                round_to_one(second),
                round_to_one(fourth),
            ],
        );
        let b_first = round_to_one(
            wi.iter()
                .zip(si.iter().zip(zi.iter()))
                .map(|(w, (s, z))| w * s * z)
                .sum(),
        );
        let b_second = round_to_one(
            -wi.iter()
                .zip(ci.iter().zip(zi.iter()))
                .map(|(w, (c, z))| w * c * z)
                .sum::<f64>(),
        );
        let b = Vector::new(vec![b_first, b_second]);
        let xy = match a.solve(b) {
            Ok(x) => x,
            Err(_) => return Err(TriangulationError::UnsolvableEquation),
        };

        let mut x0 = xy[0];
        let mut y0 = xy[1];
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut counter = 0;

        while (x.abs() - x0.abs()).abs() > 0.0001
            && (y.abs() - y0.abs()).abs() > 0.0001
            && counter <= 999
        {
            if x != 0.0 && y != 0.0 {
                x0 = x;
                y0 = y;
            } else {
                x = x0;
                y = y0;
            }

            let mui = sub(&yi, y)
                .iter()
                .zip(sub(&xi, x).iter())
                .map(|(y, x)| y.atan2(*x))
                .collect::<Vec<f64>>();
            let vector_sub = subtract_vec(&theta, &mui);
            let vector_cos = vector_sub.into_iter().map(|v| v.cos());
            let cw: f64 = wi.iter().zip(vector_cos).map(|(w, v)| w * v).sum::<f64>();
            let kappa_inv = 2.0 * (1.0 - cw)
                + (1.0 - cw).powi(2) * (0.48794 - 0.82905 * cw - 1.3915 * cw.powi(2)) / cw;
            let kappa = 1.0 / kappa_inv;
            let phi = subtract_vec(&theta, &mui);
            let t: Vec<f64> = (sub(&phi.iter().map(|p| p.cos()).collect::<Vec<f64>>(), 1.0))
                .iter()
                .map(|x| (x * 2.0 * kappa).abs().sqrt())
                .zip(phi.iter().map(|x| x.rem_euclid(2.0 * std::f64::consts::PI)))
                .map(|(tv, rv)| match rv {
                    x if (0.0..=std::f64::consts::PI).contains(&x) => tv,
                    _ => -tv,
                })
                .collect();
            let psih: Vec<f64> = t
                .iter()
                .map(|te| {
                    te.signum()
                        * match te.abs() < 1.5 {
                            true => te.abs(),
                            false => 1.5,
                        }
                })
                .collect();
            wi = div_vec(&psih, &t);
            let di: Vec<f64> = (add_vec(
                &sub(&xi, x)
                    .iter()
                    .map(|val: &f64| val * val)
                    .collect::<Vec<f64>>(),
                &sub(&yi, y)
                    .iter()
                    .map(|val: &f64| val * val)
                    .collect::<Vec<f64>>(),
            ))
            .iter()
            .map(|val| val.sqrt())
            .collect::<Vec<f64>>();
            let di_pow3 = di.iter().map(|v| v * v * v).collect::<Vec<f64>>();
            let s_star = div_vec(&sub(&yi, y), &di_pow3);
            let c_star = div_vec(&sub(&xi, x), &di_pow3);

            let first = wi
                .iter()
                .zip(si.iter().zip(s_star.iter()))
                .map(|(w, (s, ss))| w * s * ss)
                .sum();
            let second = -1_f64
                * wi.iter()
                    .zip(si.iter().zip(c_star.iter()))
                    .map(|(w, (s, cs))| w * s * cs)
                    .sum::<f64>();
            let third = -1_f64
                * wi.iter()
                    .zip(ci.iter().zip(s_star.iter()))
                    .map(|(w, (c, ss))| w * c * ss)
                    .sum::<f64>();
            let fourth = wi
                .iter()
                .zip(ci.iter().zip(c_star.iter()))
                .map(|(w, (c, cs))| w * c * cs)
                .sum();
            let q = Matrix::new(2, 2, vec![first, third, second, fourth]);
            let p = Vector::new(vec![
                wi.iter()
                    .zip(s_star.iter().zip(zi.iter()))
                    .map(|(w, (s, z))| w * s * z)
                    .sum::<f64>(),
                -1_f64
                    * wi.iter()
                        .zip(c_star.iter().zip(zi.iter()))
                        .map(|(w, (c, z))| w * c * z)
                        .sum::<f64>(),
            ]);

            let xyn = match q.solve(p) {
                Ok(n) => n,
                Err(_) => return Ok(Location { lat: x0, lon: y0 }), // there is a solution, might not be optimal
            };
            x = xyn[0];
            y = xyn[1];

            counter += 1;

            if x.is_nan() || y.is_nan() || counter == 999 {
                return Err(TriangulationError::NoTriangulation);
            }
        }
        Ok(Location { lat: x, lon: y })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_triangulation() {
        // first test
        let base_sample_one: Vec<Bearing> = vec![
            Bearing {
                loc: Location {
                    lat: -103.099,
                    lon: 33.774,
                },
                angle: 322.0,
            },
            Bearing {
                loc: Location {
                    lat: -103.098,
                    lon: 33.774,
                },
                angle: 296.0,
            },
            Bearing {
                loc: Location {
                    lat: -103.098,
                    lon: 33.775,
                },
                angle: 290.0,
            },
        ];

        assert_eq!(
            Location {
                lat: -103.10058,
                lon: 33.77577
            },
            triangulate::<Huber>(base_sample_one).unwrap()
        );

        // second test
        let base_sample_two: Vec<Bearing> = vec![
            Bearing {
                loc: Location {
                    lat: -103.080,
                    lon: 33.795,
                },
                angle: 69.0,
            },
            Bearing {
                loc: Location {
                    lat: -103.080,
                    lon: 33.796,
                },
                angle: 95.0,
            },
            Bearing {
                loc: Location {
                    lat: -103.080,
                    lon: 33.795,
                },
                angle: 51.0,
            },
        ];
        assert_eq!(
            Location {
                lat: -103.07850,
                lon: 33.79587
            },
            triangulate::<Huber>(base_sample_two).unwrap()
        );

        // third test
        let base_sample_three: Vec<Bearing> = vec![
            Bearing {
                loc: Location {
                    lat: -99.874,
                    lon: 36.506,
                },
                angle: 137.0,
            },
            Bearing {
                loc: Location {
                    lat: -99.874,
                    lon: 36.500,
                },
                angle: 51.0,
            },
            Bearing {
                loc: Location {
                    lat: -99.867,
                    lon: 36.5077,
                },
                angle: 215.0,
            },
        ];
        assert_eq!(
            Location {
                lat: -99.87076,
                lon: 36.50256
            },
            triangulate::<Huber>(base_sample_three).unwrap()
        );

        // fourth test
        let base_sample_four: Vec<Bearing> = vec![
            Bearing {
                loc: Location {
                    lat: -99.873,
                    lon: 36.5077,
                },
                angle: 97.0,
            },
            Bearing {
                loc: Location {
                    lat: -99.867,
                    lon: 36.5077,
                },
                angle: 220.0,
            },
            Bearing {
                loc: Location {
                    lat: -99.874,
                    lon: 36.5042,
                },
                angle: 50.0,
            },
        ];
        assert_eq!(
            Location {
                lat: -99.86714,
                lon: 36.50754
            },
            triangulate::<Huber>(base_sample_four).unwrap()
        );
    }
}

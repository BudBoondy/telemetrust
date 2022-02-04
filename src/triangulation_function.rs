use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::fmt;

pub mod andrews;
pub mod ewaf;
pub mod huber;
pub mod mle;

const RADIUS: f64 = 6371000.0;

pub trait TriangulationFunction {
    fn triangulate(s: Vec<Bearing>) -> Result<Location, TriangulationError>;
}

pub fn triangulate<F>(b: Vec<Bearing>) -> Result<Location, TriangulationError>
where
    F: TriangulationFunction,
{
    F::triangulate(b)
}

impl fmt::Display for TriangulationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "No triangulation.")
    }
}

impl fmt::Debug for TriangulationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{ file: {}, line: {} }}", file!(), line!())
    }
}

pub enum TriangulationError {
    EqualAngle,
    UnsolvableEquation,
    UnsolvableState,
    NoTriangulation,
}

#[derive(Debug, Clone)]
pub struct Bearing {
    pub loc: Location,
    pub angle: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct Location {
    pub lat: f64,
    pub lon: f64,
}

impl std::ops::Add<Location> for Location {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            lat: self.lat + other.lat,
            lon: self.lon + other.lon,
        }
    }
}

impl std::ops::Div<f64> for Location {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        if rhs == 0.0 {
            panic!("cannot divide by 0");
        }
        Self {
            lat: self.lat / rhs,
            lon: self.lon / rhs,
        }
    }
}

/*
impl Statistics for Vec<Location>{
    type Inner = Location;

    fn average(&self) -> Location {}
    fn mean(&self) -> Location {}
    fn variance(&self) -> Location {}

}
*/

impl Location {
    pub fn to_string_precision(self, i: usize) -> String {
        format!("{:.prec$}:{:.prec$}", self.lat, self.lon, prec = i)
    }

    // faster algorithm?
    pub fn distance_to_spherical(self, other: Location) -> f64 {
        let n1 = to_normalized_vector(self);
        let n2 = to_normalized_vector(other);

        let sin0: f64 = length(cross_prod(&n1, &n2));
        let cos0 = dot(n1, n2);
        let a = sin0.atan2(cos0); // tanδ = |n₁×n₂| / n₁⋅n₂

        a * RADIUS
    }

    pub fn distance_to(&self, other: &Location) -> f64 {
        let phi1 = self.lat * PI / 180.0;
        let phi2 = other.lat * PI / 180.0;
        let delta_phi = (other.lat - self.lat) * PI / 180.0;
        let delta_lambda = (other.lon - self.lon) * PI / 180.0;
        let a = (delta_phi / 2.0).sin() * (delta_phi / 2.0)
            + phi1.cos() * phi2.cos() * (delta_lambda / 2.0).sin() * (delta_lambda / 2.0).sin();
        let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
        RADIUS * c
    }
}

// adapt to support round errors after the fifth number
impl PartialEq for Location {
    fn eq(&self, other: &Self) -> bool {
        (self.lat - other.lat).abs() < 0.01 && (self.lon - other.lon).abs() < 0.01
    }
}

impl ToString for Location {
    fn to_string(&self) -> String {
        format!("{:.1}:{:.1}", self.lat, self.lon)
    }
}

pub fn interpolate_location(x: f64, x1: f64, x2: f64, fx1: Location, fx2: Location) -> Location {
    let lat = fx1.lat + ((x - x1) / (x2 - x1)) * (fx2.lat - fx1.lat);
    let lon = fx1.lon + ((x - x1) / (x2 - x1)) * (fx2.lon - fx1.lon);
    Location { lat, lon }
}

pub fn to_lat_lng(v: Vec<f64>) -> Location {
    let x = v[0];
    let y = v[1];
    let z = v[2];
    let lat = z.atan2((x * x + y * y).sqrt());
    let lng = y.atan2(x);
    Location {
        lat: lat.to_degrees(),
        lon: lng.to_degrees(),
    }
}

fn to_normalized_vector(l: Location) -> Vec<f64> {
    let lat_n = l.lat.to_radians();
    let lng_n = l.lon.to_radians();
    let sin_lat = lat_n.sin();
    let cos_lat = lat_n.cos();
    let sin_lng = lng_n.sin();
    let cos_lng = lng_n.cos();

    // right-handed vector: x -> 0°E,0°N; y -> 90°E,0°N, z -> 90°N
    let x = cos_lat * cos_lng;
    let y = cos_lat * sin_lng;
    let z = sin_lat;
    vec![x, y, z]
}

fn unit(v: Vec<f64>) -> Vec<f64> {
    let norm = v.len() as f64;
    v.iter().map(|x| x / norm).collect::<Vec<f64>>()
}

fn dot(v1: Vec<f64>, v2: Vec<f64>) -> f64 {
    let x = v1[0] * v2[0];
    let y = v1[1] * v2[1];
    let z = v1[2] * v2[2];
    x + y + z
}

fn dot_vec(v1: &[f64], v2: &[f64]) -> Vec<f64> {
    v1.iter().zip(v2.iter()).map(|(a, b)| a * b).collect()
}

fn subtract_vec(v1: &[f64], v2: &[f64]) -> Vec<f64> {
    v1.iter().zip(v2.iter()).map(|(a, b)| a - b).collect()
}

pub fn fold_location_vec(loc: Vec<Location>) -> Location {
    let mut normalized_loc = Location { lat: 0.0, lon: 0.0 };
    let mut normalized_loc = loc.iter().fold(normalized_loc, |_a, b| {
        normalized_loc.lat += b.lat;
        normalized_loc.lon += b.lon;
        normalized_loc
    });
    normalized_loc.lat /= loc.len() as f64;
    normalized_loc.lon /= loc.len() as f64;
    normalized_loc
}

fn cross_prod(v1: &[f64], v2: &[f64]) -> Vec<f64> {
    let x = v1[1] * v2[2] - v1[2] * v2[1];
    let y = v1[2] * v2[0] - v1[0] * v2[2];
    let z = v1[0] * v2[1] - v1[1] * v2[0];
    vec![x, y, z]
}

fn add_vec(v1: &[f64], v2: &[f64]) -> Vec<f64> {
    v1.iter().zip(v2.iter()).map(|(a, b)| a + b).collect()
}

fn div_vec(v1: &[f64], v2: &[f64]) -> Vec<f64> {
    v1.iter().zip(v2.iter()).map(|(a, b)| a / b).collect()
}

fn sub(v1: &[f64], f: f64) -> Vec<f64> {
    v1.iter().map(|a| f - a).collect::<Vec<f64>>()
}

fn round_to_one(v: f64) -> f64 {
    (v * 10.0).round() / 10.0
}

trait DegRad {
    fn to_radians(&self) -> f64;
    fn to_degree(&self) -> f64;
}

impl DegRad for f64 {
    fn to_radians(&self) -> f64 {
        self * std::f64::consts::PI / 180.0
    }
    fn to_degree(&self) -> f64 {
        self * 180.0 / std::f64::consts::PI
    }
}

impl DegRad for Bearing {
    fn to_radians(&self) -> f64 {
        std::f64::consts::PI / 180.0 * (90.0 - self.angle)
    }

    fn to_degree(&self) -> f64 {
        0.0
    }
}

fn length(v: Vec<f64>) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_vector_folding() {
        let v: Vec<Location> = vec![
            Location { lat: 1.5, lon: 2.5 },
            Location { lat: 2.5, lon: 3.5 },
        ];
        assert_eq!(Location { lat: 2.0, lon: 3.0 }, fold_location_vec(v));
    }

    #[test]
    fn test_add_location() {
        let l1 = Location { lat: 5.0, lon: 6.0 };
        let l2 = Location { lat: 1.0, lon: 2.0 };
        assert_eq!(Location { lat: 6.0, lon: 8.0 }, l1 + l2);
    }

    #[test]
    fn test_location_interpolation() {
        let x = 1.5;
        let x1 = 1.0;
        let x2 = 2.0;

        let fx1 = Location {
            lat: 47.62,
            lon: -122.33,
        };
        let fx2 = Location {
            lat: 61.2,
            lon: -149.9,
        };
        assert_eq!(
            interpolate_location(x, x1, x2, fx1, fx2),
            Location {
                lat: 54.41,
                lon: -136.115
            }
        );
        let x = 1.9;

        let fx1 = Location {
            lat: 47.62,
            lon: -122.33,
        };
        let fx2 = Location {
            lat: 61.2,
            lon: -149.9,
        };
        assert_eq!(
            interpolate_location(x, x1, x2, fx1, fx2),
            Location {
                lat: 59.842,
                lon: -147.143
            }
        );
    }
}

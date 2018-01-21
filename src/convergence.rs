/// Type can check if iterative root-finding process has converged.
pub trait IsConverged {
    /// Indicate whether root-finding has converged.  The 'x_pre' and 'x_cur'
    /// are the previous and current iteration x values.  The 'f_cur' holds
    /// f(x_cur).
    fn is_converged(&self, x_pre: f64, x_cur: f64, f_cur: f64) -> bool;
}

/// SequenceDelta converges when the distance along the x-axis between
/// successive iterations becomes smaller than 'epsilon_abs'.
///
/// The test_pathology_microstep() shows a situation where this converges
/// prematurely.  Specifically, for a method like Newton-Raphson, a massive 
/// first derivative means taking only a small step along the x-axis even when 
/// far from the actual root.
pub struct SequenceDelta {
    epsilon_abs: f64,
}

impl SequenceDelta {
    pub fn new(epsilon_abs: f64) -> SequenceDelta {
        assert!(epsilon_abs > 0.0);
        assert!(epsilon_abs.is_finite());
        SequenceDelta { epsilon_abs }
    }
}

impl IsConverged for SequenceDelta {
    fn is_converged(&self, x_pre: f64, x_cur: f64, _f_cur: f64) -> bool {
        (x_pre - x_cur).abs() < self.epsilon_abs
    }
}

/// FnResidual converges when the residual is small: |f(x_cur)| < epsilon_abs.
///
/// Be aware that convergence can happen far from the actual root.  For example,
/// f(x)=-1e-7x+0.01 has the root at 100000, but with an epsilon_abs of 1e-3 we
/// would converge anywhere in the range [90000, 110000].
pub struct FnResidual {
    epsilon_abs: f64,
}

impl FnResidual {
    pub fn new(epsilon_abs: f64) -> FnResidual {
        assert!(epsilon_abs >= 0.0);
        assert!(epsilon_abs.is_finite());
        FnResidual { epsilon_abs }
    }
}

impl IsConverged for FnResidual {
    fn is_converged(&self, _x_pre: f64, _x_cur: f64, f_cur: f64) -> bool {
        f_cur.abs() < self.epsilon_abs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64;

    #[test]
    fn test_sequence_delta_convergence() {
        // too far apart
        let c = SequenceDelta::new(1e-9);
        let x_0 = 10.2;
        assert_eq!(false, c.is_converged(x_0, x_0 + 1e-8, 10.0));

        // just right
        assert_eq!(true, c.is_converged(x_0, x_0 + 5e-10, 10.0));
    }

    #[test]
    #[should_panic]
    fn test_sequence_delta_epsabs_zero() {
        let _ = SequenceDelta::new(0.0);
    }

    #[test]
    #[should_panic]
    fn test_sequence_delta_epsabs_negative() {
        let _ = SequenceDelta::new(-1.0);
    }

    #[test]
    #[should_panic]
    fn test_sequence_delta_epsabs_nan() {
        let _ = SequenceDelta::new(f64::NAN);
    }

    #[test]
    fn test_fn_residual_convergence() {
        let c = FnResidual::new(1e-3);
        assert_eq!(false, c.is_converged(0.0, 1e-10, 2e-3));
        assert_eq!(true, c.is_converged(0.0, 1e-10, 9e-4));
    }

    #[test]
    #[should_panic]
    fn test_fn_residual_epsabs_negative() {
        let _ = FnResidual::new(-1.0);
    }

    #[test]
    #[should_panic]
    fn test_fn_residual_epsabs_nan() {
        let _ = FnResidual::new(f64::NAN);
    }
}

//! Convergence criteria.
//!
//! This module defines the convergence criteria used to stop iterative root
//! finders.  Users can supply custom criteria by implementing the `IsConverged`
//! trait.
//!
//! This module also supplies two canned criteria:
//!
//! * DeltaX - stops when the steps along x-axis, |x_pre - x_cur|, gets small enough.
//! * FnResidual - stops when |f(x_cur)| gets small enough.
//!
//! It also provides a generic wrapper `DualCriteria` allowing two
//! IsConverged implementations to be combined.
//!
//! # Examples
//! ```
//! use rootfind::convergence::*;
//!
//! let c1 = DeltaX::new(1e-6);
//! let c2 = FnResidual::new(1e-9);
//! let finish = DualCriteria::new(&c1, &c2);
//!
//! // c1 satisfied but not c2
//! assert_eq!(finish.is_converged(0.1, 0.1+1e-7, 1e-8), false);
//!
//! // both required conditions satisfied
//! assert_eq!(finish.is_converged(0.1, 0.1+1e-7, 1e-12), true);
//! ```

/// Type can check if iterative root-finding process has converged.
pub trait IsConverged {
    /// Indicate whether root-finding has converged.
    ///
    /// The `x_pre` and `x_cur` are the previous and current iteration x values.
    /// The `f_cur` holds f(x_cur).
    fn is_converged(&self, x_pre: f64, x_cur: f64, f_cur: f64) -> bool;
}

/// DeltaX converges when the distance along the x-axis between successive
/// iterations becomes smaller than epsilon_abs.
///
/// The test_pathology_microstep() shows a situation where this converges
/// prematurely.  Specifically, for a method like Newton-Raphson, a massive
/// first derivative means taking only a small step along the x-axis even when
/// far from the actual root.
pub struct DeltaX {
    epsilon_abs: f64,
}

impl DeltaX {
    pub fn new(epsilon_abs: f64) -> DeltaX {
        assert!(epsilon_abs > 0.0);
        assert!(epsilon_abs.is_finite());
        DeltaX { epsilon_abs }
    }
}

impl IsConverged for DeltaX {
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

/// DualCriteria combines two IsConverged implementors.
///
/// Both must be true for convergence.
pub struct DualCriteria<'a, C1: 'a + IsConverged, C2: 'a + IsConverged> {
    c1: &'a C1,
    c2: &'a C2,
}

impl<'a, C1: 'a + IsConverged, C2: 'a + IsConverged> DualCriteria<'a, C1, C2> {
    pub fn new(c1: &'a C1, c2: &'a C2) -> DualCriteria<'a, C1, C2> {
        DualCriteria { c1, c2 }
    }
}

impl<'a, C1: IsConverged, C2: IsConverged> IsConverged for DualCriteria<'a, C1, C2> {
    fn is_converged(&self, x_pre: f64, x_cur: f64, f_cur: f64) -> bool {
        self.c1.is_converged(x_pre, x_cur, f_cur) && self.c2.is_converged(x_pre, x_cur, f_cur)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64;

    #[test]
    fn test_delta_x_convergence() {
        // too far apart
        let c = DeltaX::new(1e-9);
        let x_0 = 10.2;
        assert_eq!(false, c.is_converged(x_0, x_0 + 1e-8, 10.0));

        // just right
        assert_eq!(true, c.is_converged(x_0, x_0 + 5e-10, 10.0));
    }

    #[test]
    #[should_panic]
    fn test_delta_x_epsabs_zero() {
        let _ = DeltaX::new(0.0);
    }

    #[test]
    #[should_panic]
    fn test_delta_x_epsabs_negative() {
        let _ = DeltaX::new(-1.0);
    }

    #[test]
    #[should_panic]
    fn test_delta_x_epsabs_nan() {
        let _ = DeltaX::new(f64::NAN);
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

    #[test]
    fn test_dual_convergence() {
        let c1 = FnResidual::new(1e-4);
        let c2 = DeltaX::new(1e-9);
        let c = DualCriteria::new(&c1, &c2);

        // neither c1 nor c2
        let x_0 = -3.7;
        assert_eq!(false, c.is_converged(x_0, x_0 + 1.0, 113456.987));

        // c1 but not c2
        assert_eq!(false, c.is_converged(x_0, x_0 + 5e-10, 113456.987));

        // c2 but not c1
        assert_eq!(false, c.is_converged(x_0, x_0 + 1.0, 0.00008));

        // both c1 and c2
        assert_eq!(true, c.is_converged(0.0, 1e-10, 0.00008));
    }

}

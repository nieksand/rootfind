pub trait IsConverged {
    fn is_converged(&self, x_pre: f64, x_cur: f64, f_cur: f64) -> bool;
}

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
    fn test_sequence_delta_accuracy_zero() {
        let _ = SequenceDelta::new(0.0);
    }

    #[test]
    #[should_panic]
    fn test_sequence_delta_accuracy_negative() {
        let _ = SequenceDelta::new(-1.0);
    }

    #[test]
    #[should_panic]
    fn test_sequence_delta_accuracy_nan() {
        let _ = SequenceDelta::new(f64::NAN);
    }

}

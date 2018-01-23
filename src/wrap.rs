/// Trait evaluating: f(x) with x in R^1.
pub trait RealFnEval {
    fn eval_f(&self, x: f64) -> f64;
}

/// Trait evaluating the derivative: df(x) with x in R^1.
pub trait RealDfEval {
    fn eval_df(&self, x: f64) -> f64;
}

/// Trait evaluating the second derivative: d2f(x) with x in R^1.
pub trait RealD2fEval {
    fn eval_d2f(&self, x: f64) -> f64;
}

/// Wraps function to implement RealFnEval.
pub struct RealFn<'a, F>
where
    F: 'a + Fn(f64) -> f64,
{
    pub f: &'a F,
}

impl<'a, F> RealFn<'a, F>
where
    F: 'a + Fn(f64) -> f64,
{
    pub fn new(f: &'a F) -> RealFn<'a, F> {
        RealFn { f }
    }
}

impl<'a, F> RealFnEval for RealFn<'a, F>
where
    F: 'a + Fn(f64) -> f64,
{
    fn eval_f(&self, x: f64) -> f64 {
        (self.f)(x)
    }
}

/// Wraps functions to implement RealFnEval and RealDfEval.
pub struct RealFnAndFirst<'a, F1, F2>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
{
    pub f: &'a F1,
    pub df: &'a F2,
}

impl<'a, F1, F2> RealFnAndFirst<'a, F1, F2>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
{
    pub fn new(f: &'a F1, df: &'a F2) -> RealFnAndFirst<'a, F1, F2> {
        RealFnAndFirst { f, df }
    }
}

impl<'a, F1, F2> RealFnEval for RealFnAndFirst<'a, F1, F2>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
{
    fn eval_f(&self, x: f64) -> f64 {
        (self.f)(x)
    }
}

impl<'a, F1, F2> RealDfEval for RealFnAndFirst<'a, F1, F2>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
{
    fn eval_df(&self, x: f64) -> f64 {
        (self.df)(x)
    }
}

/// Wraps functions to implement RealFnEval, RealDfEval, and RealD2fEval.
pub struct RealFnAndFirstSecond<'a, F1, F2, F3>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
    F3: 'a + Fn(f64) -> f64,
{
    pub f: &'a F1,
    pub df: &'a F2,
    pub d2f: &'a F3,
}

impl<'a, F1, F2, F3> RealFnAndFirstSecond<'a, F1, F2, F3>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
    F3: 'a + Fn(f64) -> f64,
{
    pub fn new(f: &'a F1, df: &'a F2, d2f: &'a F3) -> RealFnAndFirstSecond<'a, F1, F2, F3> {
        RealFnAndFirstSecond { f, df, d2f }
    }
}

impl<'a, F1, F2, F3> RealFnEval for RealFnAndFirstSecond<'a, F1, F2, F3>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
    F3: 'a + Fn(f64) -> f64,
{
    fn eval_f(&self, x: f64) -> f64 {
        (self.f)(x)
    }
}

impl<'a, F1, F2, F3> RealDfEval for RealFnAndFirstSecond<'a, F1, F2, F3>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
    F3: 'a + Fn(f64) -> f64,
{
    fn eval_df(&self, x: f64) -> f64 {
        (self.df)(x)
    }
}

impl<'a, F1, F2, F3> RealD2fEval for RealFnAndFirstSecond<'a, F1, F2, F3>
where
    F1: 'a + Fn(f64) -> f64,
    F2: 'a + Fn(f64) -> f64,
    F3: 'a + Fn(f64) -> f64,
{
    fn eval_d2f(&self, x: f64) -> f64 {
        (self.d2f)(x)
    }
}

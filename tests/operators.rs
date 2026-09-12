//! The operators on the element types agree with the named methods, in every
//! ownership form.
//!
//! The borrowed forms delegate to the methods, so agreement there checks only
//! that each operator is wired to the right method. The owned forms reuse an
//! operand's map — `a + b` moves the smaller side's terms into the larger
//! side's map — and share no code with `SymFn::add`, so the sweep over pairs
//! of both sizes is what pins them.
//!
//! The coefficient types' operators are written over their `Ring` methods by
//! one macro, so each type's sweep checks that its impls call the right
//! method in every form.

use symfn::permutation::Perm;
use symfn::schubert::Schubert;
use symfn::{
    coproduct, partitions_of, Elementary, Forgotten, Homogeneous, Ht, Monomial, Partition,
    PowerSum, Rational, Ring, Schur, St, SymFn, SymTensor,
};
use symfn::{AFrac, Field, Frac, Guarded, GuardedRat, QAlgebra, QtPoly, Ratio};

/// Elements of one to three terms with coefficients of both signs, so a sum
/// has terms that cancel, terms that merge, and terms that are new, plus zero.
fn elements<C: Ring, B: SymFn<C>>() -> Vec<B> {
    let parts: Vec<Partition> = (0..=4).flat_map(partitions_of).collect();
    let mut out = Vec::new();
    for (i, p) in parts.iter().enumerate() {
        let c = [1, -2, 3][i % 3];
        out.push(B::monomial(p.clone(), C::from_i64(c)));
        if i + 2 < parts.len() {
            let mut f = B::monomial(p.clone(), C::from_i64(c));
            f.add_term(parts[i + 1].clone(), C::from_i64(-c));
            f.add_term(parts[i + 2].clone(), C::from_i64(2));
            out.push(f);
        }
    }
    out.push(B::zero());
    out
}

macro_rules! linear_operator_laws {
    ($module:ident, $ty:ident) => {
        mod $module {
            use super::*;

            #[test]
            fn every_ownership_form_of_plus_agrees_with_add() {
                for a in elements::<i64, $ty<i64>>() {
                    for b in elements::<i64, $ty<i64>>() {
                        let want = a.add(&b);
                        assert_eq!(&a + &b, want, "&a + &b at {a}, {b}");
                        assert_eq!(a.clone() + &b, want, "a + &b at {a}, {b}");
                        assert_eq!(&a + b.clone(), want, "&a + b at {a}, {b}");
                        assert_eq!(a.clone() + b.clone(), want, "a + b at {a}, {b}");
                        let mut c = a.clone();
                        c += &b;
                        assert_eq!(c, want, "c += &b at {a}, {b}");
                        let mut c = a.clone();
                        c += b.clone();
                        assert_eq!(c, want, "c += b at {a}, {b}");
                    }
                }
            }

            #[test]
            fn every_ownership_form_of_minus_agrees_with_sub() {
                for a in elements::<i64, $ty<i64>>() {
                    for b in elements::<i64, $ty<i64>>() {
                        let want = a.sub(&b);
                        assert_eq!(&a - &b, want, "&a - &b at {a}, {b}");
                        assert_eq!(a.clone() - &b, want, "a - &b at {a}, {b}");
                        assert_eq!(&a - b.clone(), want, "&a - b at {a}, {b}");
                        assert_eq!(a.clone() - b.clone(), want, "a - b at {a}, {b}");
                        let mut c = a.clone();
                        c -= &b;
                        assert_eq!(c, want, "c -= &b at {a}, {b}");
                        let mut c = a.clone();
                        c -= b.clone();
                        assert_eq!(c, want, "c -= b at {a}, {b}");
                    }
                    let zero = &a - &a;
                    assert!(zero.terms().is_empty(), "a - a stores a term at {a}");
                }
            }

            #[test]
            fn unary_minus_agrees_with_neg_and_is_an_involution() {
                for a in elements::<i64, $ty<i64>>() {
                    assert_eq!(-&a, a.neg(), "-&a at {a}");
                    assert_eq!(-a.clone(), a.neg(), "-a at {a}");
                    assert_eq!(-(-&a), a, "-(-a) at {a}");
                    assert!((&a + (-&a)).is_zero(), "a + (-a) at {a}");
                }
            }

            #[test]
            fn a_coefficient_on_the_right_agrees_with_scale() {
                for a in elements::<i64, $ty<i64>>() {
                    for c in [3i64, -1, 0] {
                        let want = a.scale(&c);
                        assert_eq!(&a * c, want, "&a * c at {a}, {c}");
                        assert_eq!(&a * &c, want, "&a * &c at {a}, {c}");
                        assert_eq!(a.clone() * c, want, "a * c at {a}, {c}");
                        assert_eq!(a.clone() * &c, want, "a * &c at {a}, {c}");
                        let mut d = a.clone();
                        d *= c;
                        assert_eq!(d, want, "d *= c at {a}, {c}");
                        let mut d = a.clone();
                        d *= &c;
                        assert_eq!(d, want, "d *= &c at {a}, {c}");
                    }
                    assert!((&a * 0).terms().is_empty(), "a * 0 stores a term at {a}");
                    assert_eq!(&a * -1, -&a, "a * -1 at {a}");
                }
            }
        }
    };
}

linear_operator_laws!(schur, Schur);
linear_operator_laws!(power_sum, PowerSum);
linear_operator_laws!(monomial, Monomial);
linear_operator_laws!(elementary, Elementary);
linear_operator_laws!(homogeneous, Homogeneous);
linear_operator_laws!(forgotten, Forgotten);
linear_operator_laws!(st, St);
linear_operator_laws!(ht, Ht);

macro_rules! product_operator_laws {
    ($module:ident, $ty:ident) => {
        mod $module {
            use super::*;

            #[test]
            fn every_ownership_form_of_times_agrees_with_mul() {
                for a in elements::<i64, $ty<i64>>() {
                    for b in elements::<i64, $ty<i64>>() {
                        let want = a.mul(&b);
                        assert_eq!(&a * &b, want, "&a * &b at {a}, {b}");
                        assert_eq!(a.clone() * &b, want, "a * &b at {a}, {b}");
                        assert_eq!(&a * b.clone(), want, "&a * b at {a}, {b}");
                        assert_eq!(a.clone() * b.clone(), want, "a * b at {a}, {b}");
                        let mut c = a.clone();
                        c *= &b;
                        assert_eq!(c, want, "c *= &b at {a}, {b}");
                        let mut c = a.clone();
                        c *= b.clone();
                        assert_eq!(c, want, "c *= b at {a}, {b}");
                    }
                }
            }
        }
    };
}

product_operator_laws!(schur_product, Schur);
product_operator_laws!(power_sum_product, PowerSum);
product_operator_laws!(monomial_product, Monomial);
product_operator_laws!(elementary_product, Elementary);
product_operator_laws!(homogeneous_product, Homogeneous);

/// The reduced-Kronecker product is the expensive one, so it is checked on a
/// pair rather than swept.
#[test]
fn the_character_basis_products_agree_with_mul() {
    let a: St<i64> = St::monomial(Partition::new([1]), 1);
    let b: St<i64> = St::monomial(Partition::new([2]), 1) - St::monomial(Partition::new([1, 1]), 2);
    assert_eq!(&a * &b, a.mul(&b), "st product");
    assert_eq!(a.clone() * b.clone(), a.mul(&b), "st product, owned");
    let a: Ht<i64> = Ht::monomial(Partition::new([1]), 1);
    let b: Ht<i64> = Ht::monomial(Partition::new([2]), 1) + Ht::monomial(Partition::new([1, 1]), 2);
    assert_eq!(&a * &b, a.mul(&b), "ht product");
    assert_eq!(a.clone() * b.clone(), a.mul(&b), "ht product, owned");
}

// `&a * &half` is the `Mul<&C>` impl under test, not a reference taken by habit.
#[allow(clippy::op_ref)]
#[test]
fn the_operators_hold_over_the_rationals() {
    let half = Rational::from_ratio(1, 2).unwrap();
    for a in elements::<Rational, Schur<Rational>>() {
        for b in elements::<Rational, Schur<Rational>>() {
            assert_eq!(&a + &b, a.add(&b), "&a + &b at {a}, {b}");
            assert_eq!(a.clone() - b.clone(), a.sub(&b), "a - b at {a}, {b}");
            assert_eq!(&a * &b, a.mul(&b), "&a * &b at {a}, {b}");
        }
        assert_eq!(&a * &half, a.scale(&half), "&a * 1/2 at {a}");
        assert_eq!((&a * &half) * Rational::from_i64(2), a, "(a/2)·2 at {a}");
    }
}

/// `2·s_3 − (s_3 + s_{21}) = s_3 − s_{21}`, the module doctest's value, as an
/// expression that mixes every operator and an integer literal.
#[test]
fn a_mixed_expression_reads_as_the_arithmetic_it_states() {
    let s =
        |parts: &[u32]| -> Schur<i64> { Schur::monomial(Partition::new(parts.iter().copied()), 1) };
    let got = &s(&[3]) * 2 - &s(&[2]) * &s(&[1]);
    let want = s(&[3]) - s(&[2, 1]);
    assert_eq!(got, want);
    assert_eq!(-got + &want, Schur::zero());
}

fn schubert_elements() -> Vec<Schubert<i64>> {
    let perms = [
        Perm::identity(),
        Perm::from_code(&[1]),
        Perm::from_code(&[0, 1]),
        Perm::from_code(&[2, 0, 1]),
        Perm::from_code(&[1, 1]),
    ];
    let mut out = Vec::new();
    for (i, w) in perms.iter().enumerate() {
        let c = [1, -2, 3][i % 3];
        out.push(Schubert::monomial(*w, c));
        if i + 1 < perms.len() {
            let mut f = Schubert::monomial(*w, c);
            f.add_term(perms[i + 1], &-c);
            out.push(f);
        }
    }
    out.push(Schubert::zero());
    out
}

// `a * 0` is the scale-to-zero path under test.
#[allow(clippy::erasing_op)]
#[test]
fn schubert_operators_agree_with_the_named_methods() {
    for a in schubert_elements() {
        for b in schubert_elements() {
            assert_eq!(&a + &b, a.add(&b), "&a + &b at {a:?}, {b:?}");
            assert_eq!(a.clone() + b.clone(), a.add(&b), "a + b at {a:?}, {b:?}");
            assert_eq!(&a - &b, a.sub(&b), "&a - &b at {a:?}, {b:?}");
            assert_eq!(&a - b.clone(), a.sub(&b), "&a - b at {a:?}, {b:?}");
            assert_eq!(&a * &b, a.mul(&b), "&a * &b at {a:?}, {b:?}");
            assert_eq!(a.clone() * &b, a.mul(&b), "a * &b at {a:?}, {b:?}");
            let mut c = a.clone();
            c += &b;
            c *= &b;
            assert_eq!(c, a.add(&b).mul(&b), "(c += b) *= b at {a:?}, {b:?}");
        }
        assert_eq!(-&a, a.scale(&-1), "-a at {a:?}");
        assert_eq!(&a * 3, a.scale(&3), "a * 3 at {a:?}");
        assert!((&a * 0).is_zero(), "a * 0 at {a:?}");
        assert!((&a - &a).terms().is_empty(), "a - a stores a term at {a:?}");
    }
}

// `t * 0` is the scale-to-zero path under test.
#[allow(clippy::erasing_op)]
#[test]
fn tensor_operators_agree_with_add_and_cancel_exactly() {
    let tensors: Vec<SymTensor<i64>> = elements::<i64, Schur<i64>>()
        .iter()
        .map(coproduct)
        .collect();
    for t in &tensors {
        for u in &tensors {
            assert_eq!(t + u, t.add(u), "&t + &u at {t:?}, {u:?}");
            assert_eq!(t.clone() + u.clone(), t.add(u), "t + u at {t:?}, {u:?}");
            assert_eq!((t - u) + u, *t, "(t - u) + u at {t:?}, {u:?}");
            let mut v = t.clone();
            v += u;
            v -= u;
            assert_eq!(v, *t, "(v += u) -= u at {t:?}, {u:?}");
        }
        assert_eq!(t * 2, t.add(t), "t * 2 at {t:?}");
        assert!((-t + t).terms().is_empty(), "-t + t stores a term at {t:?}");
        assert!((t.clone() * 0).is_zero(), "t * 0 at {t:?}");
    }
}

/// Every ownership form of `+`, `-`, `*`, unary `-` and the assigning forms
/// on a coefficient type, against the `Ring` method it calls, over every pair
/// of `$samples`.
macro_rules! ring_operator_laws {
    ($test:ident, $ty:ty, $samples:expr) => {
        #[test]
        #[allow(clippy::clone_on_copy)]
        fn $test() {
            let xs: Vec<$ty> = $samples;
            for a in &xs {
                for b in &xs {
                    let mut sum = a.clone();
                    Ring::add_assign(&mut sum, b);
                    let mut diff = a.clone();
                    Ring::sub_assign(&mut diff, b);
                    let prod = Ring::mul(a, b);
                    assert_eq!(a + b, sum, "&a + &b at {a:?}, {b:?}");
                    assert_eq!(a.clone() + b, sum, "a + &b at {a:?}, {b:?}");
                    assert_eq!(a + b.clone(), sum, "&a + b at {a:?}, {b:?}");
                    assert_eq!(a.clone() + b.clone(), sum, "a + b at {a:?}, {b:?}");
                    assert_eq!(a - b, diff, "&a - &b at {a:?}, {b:?}");
                    assert_eq!(a.clone() - b, diff, "a - &b at {a:?}, {b:?}");
                    assert_eq!(a - b.clone(), diff, "&a - b at {a:?}, {b:?}");
                    assert_eq!(a.clone() - b.clone(), diff, "a - b at {a:?}, {b:?}");
                    assert_eq!(a * b, prod, "&a * &b at {a:?}, {b:?}");
                    assert_eq!(a.clone() * b, prod, "a * &b at {a:?}, {b:?}");
                    assert_eq!(a * b.clone(), prod, "&a * b at {a:?}, {b:?}");
                    assert_eq!(a.clone() * b.clone(), prod, "a * b at {a:?}, {b:?}");
                    let mut c = a.clone();
                    c += b;
                    assert_eq!(c, sum, "c += &b at {a:?}, {b:?}");
                    let mut c = a.clone();
                    c += b.clone();
                    assert_eq!(c, sum, "c += b at {a:?}, {b:?}");
                    let mut c = a.clone();
                    c -= b;
                    assert_eq!(c, diff, "c -= &b at {a:?}, {b:?}");
                    let mut c = a.clone();
                    c -= b.clone();
                    assert_eq!(c, diff, "c -= b at {a:?}, {b:?}");
                    let mut c = a.clone();
                    c *= b;
                    assert_eq!(c, prod, "c *= &b at {a:?}, {b:?}");
                    let mut c = a.clone();
                    c *= b.clone();
                    assert_eq!(c, prod, "c *= b at {a:?}, {b:?}");
                }
                assert_eq!(-a, Ring::neg(a), "-&a at {a:?}");
                assert_eq!(-a.clone(), Ring::neg(a), "-a at {a:?}");
                assert!(Ring::is_zero(&(a - a)), "a - a at {a:?}");
            }
        }
    };
}

/// Zero, one, a variable, a scaled variable, and a three-term polynomial whose
/// terms cancel against the others in some sums.
fn qt_samples() -> Vec<QtPoly<i64>> {
    let mut f = QtPoly::term(0, 0, 1);
    f.add_term(1, 0, -1);
    f.add_term(1, 1, 3);
    vec![
        QtPoly::zero(),
        QtPoly::one(),
        QtPoly::q(),
        QtPoly::term(0, 1, -2),
        f,
    ]
}

ring_operator_laws!(
    qt_poly_operators_agree_with_the_ring_methods,
    QtPoly<i64>,
    qt_samples()
);
ring_operator_laws!(
    rational_operators_agree_with_the_ring_methods,
    Rational,
    vec![
        Rational::zero(),
        Rational::one(),
        Rational::new(-3, 2),
        Rational::new(5, 7)
    ]
);
ring_operator_laws!(
    guarded_operators_agree_with_the_ring_methods,
    Guarded,
    vec![Guarded(0), Guarded(3), Guarded(-5)]
);
ring_operator_laws!(
    guarded_rat_operators_agree_with_the_ring_methods,
    GuardedRat,
    vec![
        GuardedRat::zero(),
        GuardedRat::from_i64(-2),
        GuardedRat::from_i64(3).div_u128(4)
    ]
);
ring_operator_laws!(
    frac_operators_agree_with_the_ring_methods,
    Frac<i64>,
    qt_samples()
        .into_iter()
        .map(Frac::from_poly)
        .chain([Frac::inv_factor(1, 1), Frac::ratio(1, 0, 0, 1)])
        .collect()
);
ring_operator_laws!(
    afrac_operators_agree_with_the_ring_methods,
    AFrac<i64>,
    vec![
        AFrac::zero(),
        AFrac::one(),
        AFrac::from_coeffs(vec![1, -2]),
        AFrac::inv_linear(1, 1),
    ]
);
ring_operator_laws!(
    ratio_operators_agree_with_the_ring_methods,
    Ratio<i64>,
    qt_samples().into_iter().map(Ratio::from_poly).collect()
);

#[test]
fn rational_division_agrees_with_field_div_in_every_form() {
    let xs = [Rational::one(), Rational::new(-3, 2), Rational::new(5, 7)];
    for a in &xs {
        for b in &xs {
            let want = Field::div(a, b);
            assert_eq!(a / b, want, "&a / &b at {a:?}, {b:?}");
            assert_eq!(*a / b, want, "a / &b at {a:?}, {b:?}");
            assert_eq!(a / *b, want, "&a / b at {a:?}, {b:?}");
            assert_eq!(*a / *b, want, "a / b at {a:?}, {b:?}");
            let mut c = *a;
            c /= b;
            assert_eq!(c, want, "c /= &b at {a:?}, {b:?}");
            let mut c = *a;
            c /= *b;
            assert_eq!(c, want, "c /= b at {a:?}, {b:?}");
        }
    }
}

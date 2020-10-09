//
//  AppController.swift
//  NewInductanceCalc
//
//  Created by Peter Huber on 2020-10-08.
//

// Implementation of paper "Improvement in Calculation of the Self- and Mutual Inductance of Thin-Wall Solenoids and Disk Coils" (S. Babic & C. Akyel, IEEE). Note that these formulas calculate air-core inductances.

import Cocoa
import Accelerate

fileprivate let relError = 1.0E-8
fileprivate let absError = 1.0E-12

class AppController: NSObject {

    @IBAction func handleTest1(_ sender: Any) {
        
        let B = 0.1
        let Tk = T(beta: B)
        DLog("B: \(B) => T: \(Tk)")
        
        let a = 1.5
        let N = 100.0
        let innerR = 15.0 * meterPerInch
        let outerR = a * innerR
        
        let L = Ldc(N: N, innerR: innerR, outerR: outerR)
        let R = (innerR + outerR) / 2.0
        let test = L / (N * N * R)
        
        DLog("Alpha: \(a) => L: \(L) => L/(N*N*R): \(test)")
    }
    
    /// Return the self-inductance of a "thin-wall solenoid" (layer).
    /// - parameter N: number of turns
    /// - parameter R: average radius (m)
    /// - parameter ht: winding height (m)
    func Ltws(N:Double, R:Double, ht:Double) -> Double
    {
        // calculate "Beta" shape-factor
        let B = (ht / 2.0) / R
        
        return µ0 * π * N * N * R * R * T(beta: B) / ht
    }
    
    /// Return the value T(kI). Note that kI is a function of beta (only), so we use beta as the argument
    /// - parameter beta: Shape factor (thin-walled solenoids), equal to half the coil height over the mean radius of the coil
    func T(beta:Double) -> Double
    {
        let kI_squared = 1.0 / (1.0 + beta * beta)
        let kI = sqrt(kI_squared)
        let kI_cubed = pow(kI_squared, 1.5)
        let K = gsl_sf_ellint_Kcomp(kI, gsl_mode_t(GSL_PREC_DOUBLE))
        let E = gsl_sf_ellint_Ecomp(kI, gsl_mode_t(GSL_PREC_DOUBLE))
        
        return 4.0 / (3.0 * π * beta * kI_cubed) * ((2.0 * kI_squared - 1.0) * E + (1.0 - kI_squared) * K - kI_cubed)
    }
    
    /// Return the self-inductance of a "disc coil" (single disc)
    /// - parameter N: number of turns
    /// - parameter innerR: inner radius (m)
    /// - parameter outerR: outer radius (m)
    func Ldc(N:Double, innerR:Double, outerR:Double) -> Double
    {
        ZAssert(outerR >= innerR, message: "Outer radius must be greater than inner radius")
        
        let a = outerR / innerR
        
        DLog("S(alpha): \(S(alpha: a))")
        
        return µ0 * N * N * innerR * S(alpha: a) / (3.0 * (a - 1) * (a - 1))
    }

    /// Return the value S(kII). Note that kII is a function of alpha (only), so we take alpha as the argument
    /// - parameter alpha: Shape factor (disc coils), equal to the disc's outer radius over inner radius
    func S(alpha:Double) -> Double
    {
        let kII_squared = 4.0 * alpha / ((1.0 + alpha) * (1.0 + alpha))
        let kII_cubed = pow(kII_squared, 1.5)
        let kII = sqrt(kII_squared)
        let alpha_squared = alpha * alpha
        let alpha_cubed = alpha * alpha_squared
        let alpha_fourth = alpha_squared * alpha_squared
        let K = gsl_sf_ellint_Kcomp(kII, gsl_mode_t(GSL_PREC_DOUBLE))
        let E = gsl_sf_ellint_Ecomp(kII, gsl_mode_t(GSL_PREC_DOUBLE))
        let J1 = self.J1(alpha:alpha)
        let J2 = self.J2(alpha:alpha)
        
        let result = (alpha_cubed + 1.0) * (2.0 * G - 1.0) - π / 2.0 * log(2.0) - (alpha_cubed + 1) * K + (alpha + 1.0) / (2.0 * alpha * kII_squared) * (((alpha_fourth + 4.0 * alpha_cubed + 4.0 * alpha + 1.0) * kII_squared - 4.0 * alpha * (alpha_squared + 1.0)) * K + (-(alpha_fourth + 2.0 * alpha_cubed + 2.0 * alpha + 1.0) * kII_squared + 4.0 * alpha * (alpha_squared + 1.0)) * E) - alpha_cubed / 2.0 * J1 - J2
        
        return result
    }

    /// The J1 integral required to calculate the self-inductance of a disk-coil
    func J1(alpha:Double) -> Double
    {
        let quadrature = Quadrature(integrator: .qags(maxIntervals: 10), absoluteTolerance: absError, relativeTolerance: relError)
        
        let integrationResult = quadrature.integrate(over: 0.0...(π / 2.0)) { x in
            
            let numerator = sqrt(1.0 + alpha * alpha + 2.0 * alpha * cos(x)) + 1.0 + alpha * cos(x)
            let denom = sqrt(1.0 + alpha * alpha - 2.0 * alpha * sin(x)) + alpha * sin(x) - 1.0
            
            return log(numerator / denom)
        }
        
        switch integrationResult {
        
        // As of Xcode 11.4, this is the correct way to unwrap a tuple in a switch() pattern
        case .success(let resultTuple):
            let (result, _ /* absError */) = resultTuple
            return result
        
        case .failure(let error):
            ALog("Error calling integration routine. The error is: \(error)")
            return 0.0
        }
    }
    
    /// The J2 integral required to calculate the self-inductance of a disk-coil
    func J2(alpha:Double) -> Double
    {
        let quadrature = Quadrature(integrator: .qags(maxIntervals: 10), absoluteTolerance: absError, relativeTolerance: relError)
        
        let integrationResult = quadrature.integrate(over: 0.0...(π / 2.0)) { x in
            
            let innerTerm = sqrt(1.0 + alpha * alpha + 2.0 * alpha * cos(2.0 * x)) + alpha + cos(2.0 * x)
        
            return log(innerTerm)
        }
        
        switch integrationResult {
        
        // As of Xcode 11.4, this is the correct way to unwrap a tuple in a switch() pattern
        case .success(let resultTuple):
            let (result, _ /* absError */) = resultTuple
            return result
        
        case .failure(let error):
            ALog("Error calling integration routine. The error is: \(error)")
            return 0.0
        }
    }
}

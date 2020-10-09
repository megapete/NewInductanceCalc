//
//  AppController.swift
//  NewInductanceCalc
//
//  Created by Peter Huber on 2020-10-08.
//

import Cocoa
import Accelerate

fileprivate let relError = 1.0E-8
fileprivate let absError = 1.0E-12

class AppController: NSObject {

    @IBAction func handleTest1(_ sender: Any) {
        
        let B = 0.1
        let kI_squared = 1.0 / (1.0 + B * B)
        let Tk = T(kI_squared: kI_squared, beta: B)
        DLog("B: \(B) => T: \(Tk)")
    }
    
    func Ltws(N:Double, R:Double, ht:Double) -> Double
    {
        let B = (ht / 2.0) / R
        let kI_squared = 1.0 / (1.0 + B * B)
        
        return µ0 * π * N * N * R * R * T(kI_squared: kI_squared, beta: B) / ht
    }
    
    func T(kI_squared:Double, beta:Double) -> Double
    {
        let kI = sqrt(kI_squared)
        let kI_cubed = pow(kI_squared, 1.5)
        let K = gsl_sf_ellint_Kcomp(kI, gsl_mode_t(GSL_PREC_DOUBLE))
        let E = gsl_sf_ellint_Ecomp(kI, gsl_mode_t(GSL_PREC_DOUBLE))
        
        return 4.0 / (3.0 * π * beta * kI_cubed) * ((2.0 * kI_squared - 1.0) * E + (1.0 - kI_squared) * K - kI_cubed)
    }
    
    func J1(alpha:Double) -> Double
    {
        let quadrature = Quadrature(integrator: .qags(maxIntervals: 10), absoluteTolerance: absError, relativeTolerance: relError)
        
        let integrationResult = quadrature.integrate(over: 0.0...(π / 2.0)) { x in
            
            let numerator = sqrt(1.0 + alpha * alpha + 2.0 * alpha * cos(x)) + 1.0 + alpha * cos(x)
            let denom = sqrt(1.0 + alpha * alpha - 2.0 * alpha * sin(x) + alpha * sin(x) - 1)
            
            return log(numerator / denom)
        }
        
        switch integrationResult {
        
        // As of Xcode 11.4, this is the correct way to unwrap a tuple in a switch() pattern
        case .success(let resultTuple):
            let (result, _ /* absError */) = resultTuple
            return result * 2.0 / π
        
        case .failure(let error):
            ALog("Error calling integration routine. The error is: \(error)")
            return 0.0
        }
    }
    
}

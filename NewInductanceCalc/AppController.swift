//
//  AppController.swift
//  NewInductanceCalc
//
//  Created by Peter Huber on 2020-10-08.
//

import Cocoa

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
    
}

//
//  MatrixDecompositionTest.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2017/05/09.
//  Copyright © 2017年 Taketo Sano. All rights reserved.
//

import XCTest
import SwiftyMath
@testable import SwiftySolver

class MatrixEliminationTests: XCTestCase {
    
    typealias M = Matrix
    typealias M1 = Matrix1
    typealias M2 = Matrix2
    typealias M5<R: EuclideanRing> = Matrix<_5, _5, R>
    
    func testNormalize_Z() {
        let A = M1(-2)
        let B = M1(2)
        let E = MatrixEliminator.eliminate(target: A) 
        XCTAssertEqual(E.result, B)
    }
    
    func testNormalize_Q() {
        let A = M1(-3./1)
        let B = M1(1./1)
        let E = MatrixEliminator.eliminate(target: A)
        XCTAssertEqual(E.result, B)
    }

    func testFullRank() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        XCTAssertEqual(E.rank, 5)
    }
    
    func testRank4() {
        let A = M5(3, -5, -22, 20, 8, 6, -11, -50, 45, 18, -1, 2, 10, -9, -3, 3, -6, -30, 27, 10, -1, 2, 7, -6, -3)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        XCTAssertEqual(E.rank, 4)
    }
    
    func testFullRank_HNF() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .RowHermite)
        
        XCTAssertTrue(E.result.isIdentity)
    }
    
    func testRank4_HNF() {
        let A = M5(3, -5, -22, 20, 8, 6, -11, -50, 45, 18, -1, 2, 10, -9, -3, 3, -6, -30, 27, 10, -1, 2, 7, -6, -3)
        let E = MatrixEliminator.eliminate(target: A, form: .RowHermite)
        print(E.result.detailDescription)
    }
    
    
    func testLeftAndLeftInverse() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        XCTAssertEqual(E.left * E.leftInverse, M5.identity)
        XCTAssertEqual(E.leftInverse * E.left, M5.identity)
    }
    
    func testRightAndRightInverse() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        XCTAssertEqual(E.right * E.rightInverse, M5.identity)
        XCTAssertEqual(E.rightInverse * E.right, M5.identity)
    }
    
    func testLeftRestriction() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        let P = E.left
        XCTAssertEqual(E.left(restrictedToCols: 2 ..< 4), P.submatrix(colRange: 2 ..< 4))
        XCTAssertEqual(E.left(restrictedToRows: 2 ..< 4), P.submatrix(rowRange: 2 ..< 4))
    }
    
    func testLeftInverseRestriction() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        let P = E.leftInverse
        XCTAssertEqual(E.leftInverse(restrictedToCols: 2 ..< 4), P.submatrix(colRange: 2 ..< 4))
        XCTAssertEqual(E.leftInverse(restrictedToRows: 2 ..< 4), P.submatrix(rowRange: 2 ..< 4))
    }
    
    func testRightRestriction() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        let Q = E.right
        XCTAssertEqual(E.right(restrictedToCols: 2 ..< 4), Q.submatrix(colRange: 2 ..< 4))
        XCTAssertEqual(E.right(restrictedToRows: 2 ..< 4), Q.submatrix(rowRange: 2 ..< 4))
    }
    
    func testRightInverseRestriction() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        let Q = E.rightInverse
        XCTAssertEqual(E.rightInverse(restrictedToCols: 2 ..< 4), Q.submatrix(colRange: 2 ..< 4))
        XCTAssertEqual(E.rightInverse(restrictedToRows: 2 ..< 4), Q.submatrix(rowRange: 2 ..< 4))
    }
    
    func testPAQ() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        XCTAssertEqual(E.left * A * E.right, E.result)
    }
    
    func testZ55_regular() {
        let A = M5(2, -1, -2, -2, -3, 1, 2, -1, 1, -1, 2, -2, -4, -3, -6, 1, 7, 1, 5, 3, 1, -12, -6, -10, -11)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)

        XCTAssertEqual(E.result, M5.identity)
    }

    func testZ55_rank4() {
        let A = M5(3, -5, -22, 20, 8, 6, -11, -50, 45, 18, -1, 2, 10, -9, -3, 3, -6, -30, 27, 10, -1, 2, 7, -6, -3)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)

        XCTAssertEqual(E.result, M5(diagonal: [1,1,1,1,0]))
    }

    func testZ55_fullRankWithFactors() {
        let A = M5(-20, -7, -27, 2, 29, 17, 8, 14, -4, -10, 13, 8, 10, -4, -6, -9, -2, -14, 0, 16, 5, 0, 5, -1, -4)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)

        XCTAssertEqual(E.result, M5(diagonal: [1,1,1,2,60]))
    }

    func testZ55_rank3WithFactors() {
        let A = M5(4, 6, -18, -15, -46, -1, 0, 6, 4, 13, -13, -12, 36, 30, 97, -7, -6, 18, 15, 49, -6, -6, 18, 15, 48)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        XCTAssertEqual(E.result, M5(diagonal: [1,1,6]))
    }

    func testZ46_rank4WithFactors() {
        let A = M<_4, _6, 𝐙>(8, -6, 14, -10, -14, 6, 12, -8, 18, -18, -20, 8, -16, 7, -23, 22, 23, -7, 32, -17, 44, -49, -49, 17)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)

        XCTAssertEqual(E.result, M(size: (4,6), diagonal: [1,1,2,12]))
    }

    func testZ46_zero() {
        let A = M<_4, _6, 𝐙>.zero
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)
        
        XCTAssertEqual(E.result, M.zero)
    }

    func testQ55_regular() {
        let A = M5<𝐐>(-3./1, 0./1, 0./1, -9./2, 0./1, 10./3, 2./1, 0./1, -15./2, 6./1, -10./3, -2./1, 0./1, 15./2, -10./1, 0./1, 0./1, 3./4, -5./1, 0./1, 0./1, 0./1, 1./1, 0./1, 0./1)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)

        XCTAssertEqual(E.result, M.identity)
    }

    func testQ55_rank3() {
        let A = M5<𝐐>(1./1, 1./1, 0./1, 8./3, 10./3, -3./1, 0./1, 0./1, -3./1, -5./1, 2./1, 0./1, 10./3, 2./1, 16./3, 79./8, 0./1, 395./24, 79./8, 79./3, 7./2, 0./1, 35./6, 7./2, 28./3)
        let E = MatrixEliminator.eliminate(target: A, form: .Smith)

        XCTAssertEqual(E.result, M(diagonal: [1,1,1]))
    }
    
    func testQPolynomial() {
        typealias R = Polynomial<_x, 𝐐>
        
        let x = R.indeterminate
        let I = Matrix3<R>.identity
        let A = Matrix3([0, 2, 1,
                        -4, 6, 2,
                        4, -4, 0].map{ R($0)} )
        let P = x * I - A
        let e = MatrixEliminator.eliminate(target: P, form: .Smith)
        
        XCTAssertEqual(e.result, Matrix3<R>(diagonal: [R(1), x - R(2), (x - R(2)).pow(2)]))
    }
    
    public func testKernel() {
        let A = M2(1, 2, 1, 2)
        let E = MatrixEliminator.eliminate(target: A)
        let K = E.kernelMatrix
        
        XCTAssertTrue(K.size == (2, 1))
        XCTAssertTrue((A * K).isZero)

        let T = E.kernelTransitionMatrix
        XCTAssertEqual(T * K, DMatrix(size:(1, 1), grid: [1]))
    }
    
    func testKernel2() {
        let A = DMatrix(size: (6, 15), grid: [-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] )
        let E = MatrixEliminator.eliminate(target: A)
        let K = E.kernelMatrix
        
        XCTAssertTrue(K.size == (15, 10))
        XCTAssertEqual(A * K, DMatrix.zero(size: (6, 10)))
        
        let T = E.kernelTransitionMatrix
        XCTAssertTrue(T.size == (10, 15))
        XCTAssertEqual(T * K, DMatrix.identity(size: 10))
    }

    public func testImage() {
        let A = M2(2, 4, 2, 4)
        let E = MatrixEliminator.eliminate(target: A)
        let I = E.imageMatrix
        
        XCTAssertTrue(I.size == (2, 1))
        XCTAssertEqual(I.asArray, [2, 2])
        
        let T = E.imageTransitionMatrix
        XCTAssertTrue(T.size == (1, 2))
        XCTAssertEqual(T * I, DMatrix(size:(1, 1), grid: [2]))
    }
    
    public func testDet() {
        let A = Matrix4(3,-1,2,4,
                        2,1,1,3,
                        -2,0,3,-1,
                        0,-2,1,3)
        let E = MatrixEliminator.eliminate(target: A)
        XCTAssertEqual(E.determinant, 66)
    }
    
    public func testLinEq() {
        let A = M<_6, _4, 𝐙>(8, -6, 14, -10, -14, 6, 12, -8, 18, -18, -20, 8, -16, 7, -23, 22, 23, -7, 32, -17, 44, -49, -49, 17)
        let E = MatrixEliminator.eliminate(target: A)
        let y = A * Vector4(1,2,3,4)
        
        if let x = E.invert(y) {
            XCTAssertEqual(A * x, y)
        } else {
            XCTFail()
        }
        
        let y2 = ColVector<_6, 𝐙>(243996, -422477, 555238, -482263, 689731, 1363066)
        XCTAssertNil(E.invert(y2))
        
        let y3 = ColVector<_6, 𝐙>(520530, -901291, 1184519, -1028837, 1471438, 2907903)
        XCTAssertNil(E.invert(y3))
    }
}

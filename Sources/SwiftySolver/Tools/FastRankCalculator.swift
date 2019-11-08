//
//  RankCalculator.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/11/03.
//

import Foundation
import SwiftyMath

public final class FastRankCalculator {
    public typealias R = 𝐅₂ // TODO to general field
    typealias Row = RowAlignedMatrixData<R>.Row
    
    private let data: RowAlignedMatrixData<R>
    private var pivots: [MatrixComponent<R>]
    private var remainRows: [Int]
    private let slots: Int
    
    public static func calculateRank<n, m>(of A: Matrix<n, m, R>) -> Int {
        let pivots = MatrixPivotFinder.findPivots(of: A)
        return self.calculateRank(of: A, withPivots: pivots.pivots)
    }
    
    public static func calculateRank<n, m>(of A: Matrix<n, m, R>, withPivots pivots: [MatrixComponent<R>]) -> Int {
        let calc = FastRankCalculator(data: RowAlignedMatrixData(A), pivots: pivots)
        return calc.run()
    }
    
    private init(data: RowAlignedMatrixData<R>, pivots: [MatrixComponent<R>]) {
        self.data = data
        self.pivots = pivots
        self.remainRows = Set(0 ..< data.size.rows)
            .subtracting(pivots.map{ $0.0 })
            .subtracting(data.rows.enumerated().compactMap{ (i, l) in l.isEmpty ? i : nil })
            .sorted()
        self.slots = 8
    }

    private func run() -> Int {
        if remainRows.isEmpty {
            return pivots.count
        }

        var buffer = ContiguousArray(repeating: R.zero, count: data.size.cols * slots)
        buffer.withUnsafeMutableBufferPointer { p in
            self.run(buffer: p)
        }

        return pivots.count
    }
    
    private func run(buffer: UnsafeMutableBufferPointer<R>) {
        let N = remainRows.count
        var step = 1
        
        while step <= N {
            let newRows = Array(0 ..< slots).parallelCompactMap { i -> Row? in
                let row = self.randomRow(step, useOnlyUnits: step < N)
                
                self.eliminateRow(row, buffer: buffer, slot: i)
                
                return !row.isEmpty ? row : nil
            }
            
            var added = 0
            for row in newRows {
                eliminateRow(row, buffer: buffer, slot: 0)
                
                if row.isEmpty {
                    continue
                }
                
                let i = data.size.rows
                let c = row.headElement!
                
                data.append(row)
                pivots.append((i, c.col, c.value))
                added += 1
            }

            if added == 0 { // proceed to next step
                if step < N {
                    step = min(2 * step, N)
                } else {
                    break
                }
            }
        }

    }
    
    private func randomRow(_ k: Int, useOnlyUnits: Bool) -> Row {
        let row = Row()
        for _ in 0 ..< k {
            let randRange = useOnlyUnits ? 1 ..< 2 : 0 ..< 2
            let r = R(from: Int.random(in: randRange))
            if r.isZero { continue }
            
            let i = remainRows.randomElement()!
            RowAlignedMatrixData.addRow(data.row(i), into: row, multipliedBy: r)
        }
        
        return row
    }
    
    private func eliminateRow(_ row: Row, buffer: UnsafeMutableBufferPointer<R>, slot: Int) {
        let m = data.size.cols
        let s = slot * m
        
        for (j, a) in row {
            buffer[s + j] = a
        }
        
        for (i0, j0, u0) in pivots {
            let a = buffer[s + j0]
            if a.isZero {
                continue
            }
            
            let r = -a * u0.inverse!
            let pRow = data.row(i0)
            
            for (j, u) in pRow {
                buffer[s + j] = buffer[s + j] + r * u
            }
        }
        
        row.removeAll()
        
        for j in (0 ..< m).reversed() {
            let a = buffer[s + j]
            if a.isZero {
                continue
            }
            
            row.insertHead( (j, a) )
            buffer[s + j] = .zero
        }
    }
}

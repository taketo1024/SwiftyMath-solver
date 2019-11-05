// swift-tools-version:5.1
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SwiftySolver",
    products: [
        // Products define the executables and libraries produced by a package, and make them visible to other packages.
        .library(
            name: "SwiftySolver",
            targets: ["SwiftySolver"]),
    ],
    dependencies: [
        .package(url: "https://github.com/taketo1024/SwiftyMath.git", .branch("develop")),
    ],
    targets: [
        .target(
            name: "SwiftySolver",
            dependencies: ["SwiftyMath"]),
        .testTarget(
            name: "SwiftySolverTests",
            dependencies: ["SwiftySolver"]),
    ]
)

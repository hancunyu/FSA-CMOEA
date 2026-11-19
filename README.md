# FSA CMOEA

Source code and problem interfaces for FSA-CMOEA: Feasibility-Structure-Aware Evolutionary Optimization for Multi-Objective Optimal Power Flow.

## Contents

- algorithm/main contains FSA-CMOEA and its helper functions for IEEE30, IEEE57 and RWMOP40--46. The class retains its development name OPFPI_M3E.
- algorithm/IEEE118 contains the implementation used for the IEEE118 case.
- problems contains the IEEE30, IEEE57 and IEEE118 interfaces.
- run_fsa_case.m runs a single seeded example using the paper's default algorithm settings.

The archived algorithm and interface files are unchanged. source_sha256.txt records their checksums.

## Requirements

The revision experiments used MATLAB R2023b and MATPOWER 8.1.

Install PlatEMO from https://github.com/BIMK/PlatEMO and MATPOWER from https://github.com/MATPOWER/matpower separately. MATPOWER is required for the IEEE cases. RWMOP40--46 use the problem interfaces distributed with PlatEMO. dependency_sha256.txt identifies the PlatEMO core and RWMOP files used in the revision.

## Run

Start a fresh MATLAB session, change to this repository directory and set the two dependency paths for your computer. PLATEMO_ROOT points to the directory containing Algorithms and Problems.

~~~matlab
setenv('PLATEMO_ROOT', '/path/to/PlatEMO');
setenv('MATPOWER_ROOT', '/path/to/matpower8.1');
r = run_fsa_case("Case1", Seed=1);
save('case1_seed1.mat', 'r');
~~~

Supported case names are Case1--Case8, IEEE118 and RWMOP40--RWMOP46. Cases 1--4 use the IEEE 30-bus system, and Cases 5--8 use the IEEE 57-bus system. Each group contains two two-objective cases, one three-objective case and one four-objective case.

~~~matlab
r = run_fsa_case("Case5", Seed=1);
r = run_fsa_case("IEEE118", Seed=1);
r = run_fsa_case("RWMOP40", Seed=1);
~~~

Use a fresh MATLAB session when switching between the main and IEEE118 implementations because they contain helpers with shared names. Add this repository directory only, rather than all its subdirectories with genpath.

The wrapper uses N = 100 and MaxFE = 50000 and disables feasible seeding for IEEE30 and IEEE57. MaxFE can be changed through the corresponding name-value argument. Each result contains the final decisions, objectives, constraints and counted evaluations. The original implementations may write diagnostics under Data.

## Reproduction scope

These examples run individual seeded searches. The full historical run archive and metric reference sets are outside this source package. Some manuscript tables use an archived selection of 30 runs from 120 candidates. Exact reconstruction of those entries requires the selected-run identifiers, shared initial populations and reference sets.

Original source notices are retained. Refer to the upstream licenses for PlatEMO and MATPOWER.


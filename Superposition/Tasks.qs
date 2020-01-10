// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

namespace Quantum.Kata.Superposition {

    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;


    //////////////////////////////////////////////////////////////////
    // Welcome!
    //////////////////////////////////////////////////////////////////

    // "Superposition" quantum kata is a series of exercises designed
    // to get you familiar with programming in Q#.
    // It covers the following topics:
    //  - basic single-qubit and multi-qubit gates,
    //  - superposition,
    //  - flow control and recursion in Q#.

    // Each task is wrapped in one operation preceded by the description of the task.
    // Each task (except tasks in which you have to write a test) has a unit test associated with it,
    // which initially fails. Your goal is to fill in the blank (marked with // ... comment)
    // with some Q# code to make the failing test pass.

    // The tasks are given in approximate order of increasing difficulty; harder ones are marked with asterisks.

    // Task 1. Plus state
    // Input: a qubit in the |0⟩ state.
    // Goal: prepare a |+⟩ state on this qubit (|+⟩ = (|0⟩ + |1⟩) / sqrt(2)).
    operation PlusState (q : Qubit) : Unit {
        // Hadamard gate H will convert |0⟩ state to |+⟩ state.
        // Type the following: H(q);
        // Then rebuild the project and rerun the tests - T01_PlusState_Test should now pass!
        H(q);
    }


    // Task 2. Minus state
    // Input: a qubit in the |0⟩ state.
    // Goal: prepare a |-⟩ state on this qubit (|-⟩ = (|0⟩ - |1⟩) / sqrt(2)).
    operation MinusState (q : Qubit) : Unit {
        X(q);
        H(q);
    }


    // Task 3*. Unequal superposition
    // Inputs:
    //      1) a qubit in the |0⟩ state.
    //      2) angle alpha, in radians, represented as Double.
    // Goal: prepare a cos(alpha) * |0⟩ + sin(alpha) * |1⟩ state on this qubit.
    operation UnequalSuperposition (q : Qubit, alpha : Double) : Unit {
        // Hint: Experiment with rotation gates from the Microsoft.Quantum.Intrinsic namespace.
        // Note that all rotation operators rotate the state by _half_ of its angle argument.
        Ry(2.0 * alpha, q);
    }


    // Task 4. Superposition of all basis vectors on two qubits
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal:  create the following state on these qubits: (|00⟩ + |01⟩ + |10⟩ + |11⟩) / 2.
    operation AllBasisVectors_TwoQubits (qs : Qubit[]) : Unit {
        // The following lines enforce the constraints on the input that you are given.
        // You don't need to modify them. Feel free to remove them, this won't cause your code to fail.
        EqualityFactI(Length(qs), 2, "The array should have exactly 2 qubits.");
        H(qs[0]);
        H(qs[1]);
    }


    // Task 5. Superposition of basis vectors with phases
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal:  create the following state on these qubits: (|00⟩ + i*|01⟩ - |10⟩ - i*|11⟩) / 2.
    operation AllBasisVectorsWithPhases_TwoQubits (qs : Qubit[]) : Unit {
        // The following lines enforce the constraints on the input that you are given.
        // You don't need to modify them. Feel free to remove them, this won't cause your code to fail.
        EqualityFactI(Length(qs), 2, "The array should have exactly 2 qubits.");

        // Hint: Is this state separable?
        // Yes, |-> tensor |i>
        X(qs[0]);
        H(qs[0]);
        H(qs[1]);
        S(qs[1]);
    }


    // Task 6. Bell state
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal: create a Bell state |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2) on these qubits.
    operation BellState (qs : Qubit[]) : Unit {
        H(qs[0]);
        CNOT(qs[0], qs[1]);
    }


    // Task 7. All Bell states
    // Inputs:
    //      1) two qubits in |00⟩ state (stored in an array of length 2)
    //      2) an integer index
    // Goal: create one of the Bell states based on the value of index:
    //       0: |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2)
    //       1: |Φ⁻⟩ = (|00⟩ - |11⟩) / sqrt(2)
    //       2: |Ψ⁺⟩ = (|01⟩ + |10⟩) / sqrt(2)
    //       3: |Ψ⁻⟩ = (|01⟩ - |10⟩) / sqrt(2)
    operation AllBellStates (qs : Qubit[], index : Int) : Unit {
        BellState(qs);
        if (index % 2 == 1) {
            Z(qs[1]);
        }
        if (index / 2 == 1) {
            X(qs[1]);
        }
    }


    // Task 8. Greenberger–Horne–Zeilinger state
    // Input: N qubits in |0...0⟩ state.
    // Goal: create a GHZ state (|0...0⟩ + |1...1⟩) / sqrt(2) on these qubits.
    operation GHZ_State (qs : Qubit[]) : Unit {
        // Hint: N can be found as Length(qs).
        H(qs[0]);
        let N = Length(qs);
        for (i in 1..(N - 1)) {
            CNOT(qs[0], qs[i]);
        }
    }


    // Task 9. Superposition of all basis vectors
    // Input: N qubits in |0...0⟩ state.
    // Goal: create an equal superposition of all basis vectors from |0...0⟩ to |1...1⟩
    // (i.e. state (|0...0⟩ + ... + |1...1⟩) / sqrt(2^N) ).
    operation AllBasisVectorsSuperposition (qs : Qubit[]) : Unit {
        for (q in qs) {
            H(q);
        }
    }


    // Task 10. Superposition of all even or all odd numbers
    // Inputs:
    //      1) N qubits in |0...0⟩ state.
    //      2) A boolean isEven.
    // Goal: create a superposition of all even numbers on N qubits if isEven is true,
    //       or a superposition of all odd numbers on N qubits if isEven is false.
    // 
    // A basis state encodes an integer number using big-endian binary notation: 
    // state |01⟩ corresponds to the integer 1, and state |10⟩ - to the integer 2. 
    //
    // Example: for N = 2 and isEven = true the required state is (|00⟩ + |10⟩) / sqrt(2), 
    //      and for N = 2 and isEven = false - (|01⟩ + |11⟩) / sqrt(2).
    operation EvenOddNumbersSuperposition (qs : Qubit[], isEven : Bool) : Unit {
        AllBasisVectorsSuperposition(qs);
        H(qs[Length(qs) - 1]);
        if (not isEven) {
            X(qs[Length(qs) - 1]);
        }
    }


    // Task 11. |00⟩ + |01⟩ + |10⟩ state
    // Input: 2 qubits in |00⟩ state.
    // Goal: create the state (|00⟩ + |01⟩ + |10⟩) / sqrt(3) on these qubits.
    operation ThreeStates_TwoQubits (qs : Qubit[]) : Unit {
        let angle = ArcSin(Sqrt(2.0 / 3.0));
        Ry(2.0 * angle, qs[0]); // state is (|00> + 2 |10>) / Sqrt(3)
        Controlled H([qs[0]], qs[1]); // state is (|00> + |10> + |11>) / Sqrt(3)
        X(qs[0]); // state is (|10> + |00> + |01>) / Sqrt(3)
    }


    // Task 12*. Hardy State
    // Input: 2 qubits in |00⟩ state
    // Goal: create the state (3|00⟩ + |01⟩ + |10⟩ + |11⟩) / sqrt(12) on these qubits.
    operation Hardy_State (qs : Qubit[]) : Unit {
        // rotate by sqrt(2/12) to get probability of first qubit right
        let angle1 = ArcSin(Sqrt(2.0 / 12.0));
        Ry(2.0 * angle1, qs[0]);

        // hadamard second if first is |1>
        Controlled H([qs[0]], qs[1]);

        // rotate second by sqrt(1/10) if first is |0>
        let angle2 = ArcSin(Sqrt(1.0 / 10.0));
        X(qs[0]);
        Controlled Ry([qs[0]], (2.0 * angle2, qs[1]));
        X(qs[0]);
    }


    // Task 13. Superposition of |0...0⟩ and given bit string
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) bit string represented as Bool[]
    // Goal: create an equal superposition of |0...0⟩ and basis state given by the bit string.
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // You are guaranteed that the qubit array and the bit string have the same length.
    // You are guaranteed that the first bit of the bit string is true.
    // Example: for bit string = [true, false] the qubit state required is (|00⟩ + |10⟩) / sqrt(2).

    // Helper that Creates a superposition of given bitstring
    operation CreateBitString(qs: Qubit[], bits : Bool[]) : Unit is Adj + Ctl {
        EqualityFactI(Length(bits), Length(qs), "Arrays should have the same length");
        for (i in 0..(Length(bits) - 1)) {
            if (bits[i]) {
                X(qs[i]);
            }
        }
    }

    operation ZeroAndBitstringSuperposition (qs : Qubit[], bits : Bool[]) : Unit {
        // The following lines enforce the constraints on the input that you are given.
        // You don't need to modify them. Feel free to remove them, this won't cause your code to fail.
        EqualityFactI(Length(bits), Length(qs), "Arrays should have the same length");
        EqualityFactB(bits[0], true, "First bit of the input bit string should be set to true");

        using (ancilla = Qubit()) {
            H(ancilla);

            // In the |1> case of ancilla set qs to the desired bitstring
            Controlled CreateBitString([ancilla], (qs, bits));

            // Uncompute ancilla if resulting qs is input bits
            // NOTE: cannot use Hadamard because state is not a tensor product anymore (its entangled)
            (ControlledOnBitString(bits, X))(qs, ancilla);
        }
    }


    // Task 14. Superposition of two bit strings
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) two bit string represented as Bool[]s
    // Goal: create an equal superposition of two basis states given by the bit strings.
    //
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // Example: for bit strings [false, true, false] and [false, false, true]
    // the qubit state required is (|010⟩ + |001⟩) / sqrt(2).
    // You are guaranteed that the both bit strings have the same length as the qubit array,
    // and that the bit strings will differ in at least one bit.
    operation TwoBitstringSuperposition (qs : Qubit[], bits1 : Bool[], bits2 : Bool[]) : Unit is Ctl {
        using (anc = Qubit()) {
            H(anc);

            // In the |1> case of ancilla set to bits2
            Controlled CreateBitString([anc], (qs, bits2));

            // In the |0> case of ancilla set to bits1
            X(anc);
            Controlled CreateBitString([anc], (qs, bits1));
            X(anc);

            // Uncompute if resulting qs is input bits1
            (ControlledOnBitString(bits2, X))(qs, anc);
        }
    }


    // Task 15*. Superposition of four bit strings
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) four bit string represented as Bool[][] bits
    //         bits is an array of size 4 x N which describes the bit strings as follows:
    //         bits[i] describes the i-th bit string and has N elements.
    //         All four bit strings will be distinct.
    //
    // Goal: create an equal superposition of the four basis states given by the bit strings.
    //
    // Example: for N = 3 and bits = [[false, true, false], [true, false, false], [false, false, true], [true, true, false]]
    //          the state you need to prepare is (|010⟩ + |100⟩ + |001⟩ + |110⟩) / 2.
    operation FourBitstringSuperposition (qs : Qubit[], bits : Bool[][]) : Unit {
        // Hint: remember that you can allocate extra qubits.
        using (anc = Qubit()) {
            H(anc);

            // In the |1> case of ancilla set to equal superposition of bits2 and bits3
            Controlled TwoBitstringSuperposition([anc], (qs, bits[2], bits[3]));

            // In the |0> case of ancilla set to equal superposition of bits0 and bits1
            X(anc);
            Controlled TwoBitstringSuperposition([anc], (qs, bits[0], bits[1]));
            X(anc);

            // Uncompute if resulting qs is input bits2 or bits3
            (ControlledOnBitString(bits[2], X))(qs, anc);
            (ControlledOnBitString(bits[3], X))(qs, anc);
        }
    }


    // Task 16**. W state on 2ᵏ qubits
    // Input: N = 2ᵏ qubits in |0...0⟩ state.
    // Goal: create a W state (https://en.wikipedia.org/wiki/W_state) on these qubits.
    // W state is an equal superposition of all basis states on N qubits of Hamming weight 1.
    // Example: for N = 4, W state is (|1000⟩ + |0100⟩ + |0010⟩ + |0001⟩) / 2.
    operation WState_PowerOfTwo (qs : Qubit[]) : Unit {
        // O(log N)
        let N = Length(qs);
        if (N == 1) {
            // base case
            X(qs[0]);
        } else {
            let N2 = N / 2;

            // First compute the W state for the first half qubits
            WState_PowerOfTwo(qs[0..(N2 - 1)]);

            // Then use ancilla to conditionally move the |1> to the other half
            using (anc = Qubit()) {
                H(anc);

                // In |1> case of ancilla move all the ones to to the other half
                for (i in 0..(N2 - 1)) {
                    Controlled SWAP([anc], (qs[i], qs[i + N2]));
                }

                // Uncompute
                for (i in N2..(N - 1)) {
                    CNOT(qs[i], anc);
                }
            }
        }
    }


    // Task 17**. W state on arbitrary number of qubits
    // Input: N qubits in |0...0⟩ state (N is not necessarily a power of 2).
    // Goal: create a W state (https://en.wikipedia.org/wiki/W_state) on these qubits.
    // W state is an equal superposition of all basis states on N qubits of Hamming weight 1.
    // Example: for N = 3, W state is (|100⟩ + |010⟩ + |001⟩) / sqrt(3).
    operation WState_Arbitrary (qs : Qubit[]) : Unit is Adj + Ctl {
        // This one is a lamer O(N)
        let N = Length(qs);
        if (N == 1) {
            // Base case
            X(qs[0]);
        } else {
            // Rotate first qubit by angle so that |1> amplitude is correct
            let angle = ArcSin(Sqrt(1.0 / IntAsDouble(N)));
            Ry(2.0 * angle, qs[0]);

            // On case where first qubit is |0> recurse
            (ControlledOnInt(0, WState_Arbitrary))([qs[0]], qs[1..(N - 1)]);
        }
    }
}

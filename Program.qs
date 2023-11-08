namespace qrnd {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Bitwise;
    
    /// # Описание
    /// измерение значений (коллапсирование) кубитов в массиве (который рассматриваем как один регистр)
    /// и возврат числа (равного полученной двоичной последовательности)
    operation Measure(qubits: Qubit[]) : Int {
        let results = ForEach(M, qubits);
        let i = ResultArrayAsInt(results);
        return i;
    }

    /// # Описание
    /// вычисления знака переноса при арифметической операции сложения
    /// то есть трансформация вида |k>|0> -> |k>|sign(k+value)>
    operation Sign(target: Qubit[], sign: Qubit, value: Int) : Unit {
        let k = Length(target);
        let bools = IntAsBoolArray(2^k-value, k);
        use (qubits) = (Qubit[2]) {
            for idx in 0..k-1 {
                let carry = qubits[idx%2];
                let next = qubits[1-(idx%2)];
                // вычисляем следующее значение флага переноса разряда
                if(bools[idx]) {
                    // next = carry*target[idx]^carry^target[idx]
                    Controlled X([carry, target[idx]], next);
                    Controlled X([carry], next);
                    Controlled X([target[idx]], next);
                }
                else {
                    // next = carry*target[idx] = carry&target[idx]
                    Controlled X([carry, target[idx]], next);
                }
                Reset(carry);
            } 
            Controlled X([qubits[k%2]], sign);
            ResetAll(qubits);
        }
    }

    @EntryPoint()
    operation Main(n: Int, w: Int, tests: Int) : Unit {
        Message("Hello quantum world!");

        let N=2^n;
        let p=100*w/N;
        let q=100*(N-w)/N;
        
        Message($"n={n} w={w} N={N} p={p}% q={q}% tests={tests}");

        mutable counters = [0, size=2];
        use (x,y)=(Qubit[n],Qubit()){
            for _ in 1..tests {
                ApplyToEach(H,x);
                Sign(x,y,w);
                let res = Measure([y]);
                set counters w/= res <- counters[res]+1;
                ResetAll(x);
                Reset(y);
            }
        }

        let f = [w,N-w];

        mutable total = 0;
        for element in f {
            set total+=element;
        }
        for idx in 0..1 {
            let expect = tests*f[idx]/N;
            let fact = counters[idx];
            Message($"{idx}: ... expect={expect} ... fact={fact}");
        }
    }
}

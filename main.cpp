#include <algorithm>
#include <random>
#include <vector>
#include <ctime>
using namespace std;

std::random_device rd;  
std::mt19937 engine(rd());

#define START_TIME  start_time = clock();
#define FINISH_TIME printf("Running time: %.2f\n", (double)(clock() - start_time)/CLOCKS_PER_SEC);

void assert(bool flag, string message)
{
    if (!flag)
    {
        printf("%s\n", message.c_str());
        system("pause");
    }
}

const long long MAX_POWER = 62;
const long long MAX_MODULUS = LLONG_MAX; // = (1<<63) - 1

long long modAdd(unsigned long long a, unsigned long long b, unsigned long long q)
{    
    assert(q > 0 && q <= MAX_MODULUS, "modAdd: incorrect modulus size.");
    assert(a >= 0 && b >= 0 && a < q && b < q, "modAdd: incorrect parameter values.");
    unsigned long long res = a + b;
    if (res >= q) res -= q;
    return res;
}

long long modMult(unsigned long long a, unsigned long long b, unsigned long long q)
{
    assert(q > 0 && q <= MAX_MODULUS, "modMult: incorrect modulus size.");
    assert(a >= 0 && b >= 0 && a < q && b < q, "modMult: incorrect parameter values.");
    unsigned long long res = 0;
    while (a != 0) {
        if (a & 1) res = (res + b) % q;
        a >>= 1;
        b = (b << 1) % q;
    }
    return res;
}

long long modMultPow2(unsigned long long x, int pow, unsigned long long q)
{
    unsigned long long res;
    if (pow < 0)
    {
        pow = -pow;
        x = (x >> (pow - 1));
        int round = x % 2;
        res = ((x >> 1) + round) % q;
    }
    else
    {
        x = x % q;
        long long pow2 = ((long long) 1 << pow) % q;
        res = modMult(x, pow2, q);
    }
    return res;
}

typedef int plaintext; // 1 bit

class matrix
{
    friend class matrix;
public:
    matrix()
    {
    }

    matrix(int rows_in, int columns_in)
    {
        initialize(rows_in, columns_in);
    }

    void initialize(int rows_in, int columns_in)
    {
        rows = rows_in;
        columns = columns_in;
         
        m.resize(rows);
        for (int i = 0; i < rows; ++i)
            m[i].resize(columns, 0);
    }

    void randomize(long long l, long long r, long long q)
    {
        std::uniform_int_distribution<long long> dist(l, r);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < columns; ++j)
            {
                m[i][j] = dist(engine);
                if (m[i][j] < 0) m[i][j] += q;
            }
    }

    long long getCell(int row, int column) const
    {        
        assert(row > 0 && row <= rows && column > 0 && column <= columns, "getCell: wrong indices.");
        return m[row - 1][column - 1]; // note: rows and columns are numbered from 1
    }

    void setCell(int row, int column, long long v)
    {
        assert (row > 0 && row <= rows && column > 0 && column <= columns, "setCell: wrong indices.");
        assert(v >= 0, "setCell: positive value required.");
        m[row - 1][column - 1] = v;
    }

    void negate(long long q)
    {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < columns; ++j)
                if (m[i][j] != 0)
                {
                    m[i][j] = q - m[i][j];
                    assert(m[i][j] >= 0, "negate: values do not match.");
                }
    }

    void fillWithFractions(const long long si, int pow, long long q)
    {
        // pow = (l0 - l1)
        // take (si, si, si, ..., si)
        // multiply first element by 2^{pow}, second by 2^{pow + 1}, etc.
        // round each of the resulting elements
        assert(columns == 1, "fillWithFractions: must be a column vector.");
        for (int i = 0; i < rows; ++i)
            m[i][0] = modMultPow2(si, pow + i, q);
    }

    void addEncodedMessage(plaintext x, int n, int ell, long long q)
    {
        assert(columns == n && rows == n * ell, "addEncodedMessage: incorrect matrix dimenensions.");
        if (x == 0) return; // We're done.
        assert(x == 1, "addEncodedMessage: the plaintext must be 1 bit long.");
        assert(ell - 1 <= MAX_POWER, "addEncodedMessage: the value of ell too large.");
        for (int i = 0; i < ell; ++i)
            for (int j = 0; j < n; ++j)
            {
                int row = n * i + j;
                int column = j;
                long long r = ((long long) 1 << i);
                m[row][column] = modAdd(m[row][column], r, q);
            }
    }
    
    matrix mult(const matrix B, long long q) const
    {
        assert(columns == B.rows, "mult: columns != rows");
        matrix res(rows, B.columns);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < B.columns; ++j)
            {
                long long sum = 0;
                for (int k = 0; k < columns; ++k)
                    sum = modAdd(sum, modMult(m[i][k], B.m[k][j], q), q);
                res.m[i][j] = sum;
            }
        return res;
    }

    matrix add(const matrix B, long long q) const
    {
        assert(rows == B.rows && columns == B.columns, "add: matrix dimensions do not match.");
        matrix res(rows, columns);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < columns; ++j)
                res.m[i][j] = modAdd(m[i][j], B.m[i][j], q);
        return res;
    }

    matrix extractRowVector(int row) const
    {
        assert(row > 0 && row <= rows, "extractRowVector: wrong row number.");
        --row; // implementation indexed from 0, and the interface from 1
        matrix res(1, columns);
        for (int i = 0; i < columns; ++i)
            res.m[0][i] = m[row][i];
        return res;
    }

    matrix concatenate(const matrix B, bool horizontal) const
    {
        if (!horizontal)
        {
            assert(columns == B.columns, "concatenate: matrix dimensions do not match.");
            matrix res(rows + B.rows, columns);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < columns; ++j)
                    res.m[i][j] = m[i][j];
            for (int i = 0; i < B.rows; ++i)
                for (int j = 0; j < B.columns; ++j)
                    res.m[rows + i][j] = B.m[i][j];
            return res;
        }
        else // vertical
        {
            assert(rows == B.rows, "concatenate: matrix dimensions do not match.");
            matrix res(rows, columns + B.columns);
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < columns; ++j)
                    res.m[i][j] = m[i][j];
            for (int i = 0; i < B.rows; ++i)
                for (int j = 0; j < B.columns; ++j)
                    res.m[i][columns + j] = B.m[i][j];
            return res;
        }
    }

    matrix splitInHorizontalBitLayerConcatenation(int ell)
    {
        matrix res(rows, columns * ell);
        for (int layer = 0; layer < ell; ++layer)
        {
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < columns; ++j)
                    res.m[i][layer*columns + j] = (m[i][j] >> layer) % 2;
        }
        return res;
    }

private:
    int rows, columns;
    std::vector<std::vector<long long> > m;
};

typedef matrix secret_key;
typedef matrix row_vector;
typedef matrix scalar;
typedef matrix column_vector;

struct swhe_params
{
    int n, ell, beta;
    long long  q;
};

bool operator==(swhe_params pp0, swhe_params pp1)
{
    return (pp0.n == pp1.n
         && pp0.ell == pp1.ell
         && pp0.beta == pp1.beta
         && pp0.q == pp1.q);
}

struct ciphertext
{
    ciphertext() {}
    ciphertext(swhe_params pp, matrix m, int err)
        : pp(pp), m(m), error_bound(err) {}
    swhe_params pp;
    matrix      m;
    long long error_bound;
};

const double eps = 1e-10;

class SWHE
{
    friend class FHE;
public:
    SWHE(int n, int ell, int beta)
    {
        SetupKg(n, ell, beta);
    }

    SWHE(int n, int ell)
    {
        int beta = (int) ceil(sqrt(n) - eps);
        SetupKg(n, ell, beta);
    }

    ciphertext Enc(plaintext m)
    {
        matrix A_transpose(/*rows*/ pp.n * pp.ell, /*columns*/ pp.n - 1);
        column_vector e(/*rows*/ pp.n * pp.ell, /*columns*/ 1);

        A_transpose.randomize(0, pp.q - 1, pp.q);
        e.randomize(-pp.beta, pp.beta, pp.q);
        column_vector b = e.add(A_transpose.mult(s, pp.q), pp.q);
        matrix C = A_transpose.concatenate(b, /*horizontal*/ true);
        C.addEncodedMessage(m, pp.n, pp.ell, pp.q);        
        return ciphertext(pp, C, pp.beta);
    }

    plaintext Dec(ciphertext C)
    {
        row_vector c = C.m.extractRowVector(pp.n * pp.ell);
        scalar constCell(1, 1);
        constCell.setCell(1, 1, pp.q - 1);
        column_vector comb = s.concatenate(constCell, /*horizontal*/ false);
        comb.negate(pp.q);
        scalar res = c.mult(comb, pp.q);
        long long v = res.getCell(1, 1);
        if (pp.q / 4 < v && v < 3 * (pp.q / 4)) return 1;
        else return 0;
    }

private:
    void SetupKg(int n, int ell, int beta)
    {
        assert(ell <= MAX_POWER, "SWHE: ell value too large (must be less than 63)");
        pp.n = n;
        pp.ell = ell;
        pp.beta = beta;
        pp.q = ((long long) 1 << ell);
        s.initialize(/*rows*/ pp.n - 1, /*columns*/ 1);
        s.randomize(0, pp.q - 1, pp.q);
    }

    swhe_params pp;
    secret_key s;
};

ciphertext NOT(ciphertext C)
{
    ciphertext res = C;
    res.m.negate(res.pp.q);
    res.m.addEncodedMessage(/*x*/ 1, res.pp.n, res.pp.ell, res.pp.q);
    return res;
}

ciphertext XOR(ciphertext C1, ciphertext C2)
{
    assert(C1.pp == C2.pp, "XOR: the add ciphertexts have different public parameters.");
    // TODO: check that "C1.error_bound + C2.error_bound" is less than "C1.pp.q/4"
    matrix m = (C1.m).add(C2.m, C1.pp.q);
    return ciphertext(C1.pp, m, C1.error_bound + C2.error_bound);
}

ciphertext AND(ciphertext C1, ciphertext C2)
{
    assert(C1.pp == C2.pp, "AND: the add ciphertexts have different public parameters.");
    // TODO: check that "C1.error_bound * C2.error_bound" is less than "C1.pp.q/4"
    matrix a = C1.m.splitInHorizontalBitLayerConcatenation(C1.pp.ell);
    matrix b = a.mult(C2.m, C1.pp.q);
    return ciphertext(C1.pp, b, C1.error_bound * C2.error_bound);
}

struct evk_kswitch
{
    void init(swhe_params pp0_in, swhe_params pp1)
    {
        pp0 = pp0_in;
        C.resize(pp1.n);
    }
    swhe_params pp0;
    vector<matrix> C;
};

struct evk_hdec
{
    evk_hdec()
        : swhe1(5,5)
    {
    }
    void init(swhe_params pp0_in, swhe_params pp1_in, SWHE swhe1_in)
    {
        pp0 = pp0_in; pp1 = pp1_in; swhe1 = swhe1_in;
        e.resize(pp0.n);
        for (int n = 1; n < pp0.n; ++n)
        {
            e[n].resize(pp0.ell);
            for (int j = 0; j < pp0.ell; ++j)
                e[n][j].resize(pp0.q);
        }
    }
    SWHE swhe1;
    swhe_params pp1, pp0;
    vector<vector<vector<ciphertext> > > e;
};

struct evk
{
    // TODO: swhe1 is only passed for debugging, remove when done
    void init(swhe_params pp0, swhe_params pp1, SWHE swhe1)
    {
        eval_key_kswitch.init(pp0, pp1);
        eval_key_hdec.init(pp0, pp1, swhe1);
    }
    evk_kswitch eval_key_kswitch;
    evk_hdec    eval_key_hdec;
};

class FHE
{
public:
    FHE(int n, int l0, int l1)
        : swhe0(n, l0), swhe1(n, l1)
    {
        eval_key.init(swhe0.pp, swhe1.pp, swhe1);
        Kg_kswitch();
        Kg_hdec();
    }

    FHE(int n, int l0, int l1, int beta)
        : swhe0(n, l0, beta), swhe1(n, l1, beta)          
    {
        eval_key.init(swhe0.pp, swhe1.pp, swhe1);
        Kg_kswitch();
        Kg_hdec();
    }

    evk        get_evk()         { return eval_key; }
    ciphertext Enc(plaintext m)  { return swhe1.Enc(m); }
    plaintext  Dec(ciphertext C) { return swhe1.Dec(C); }

    // These two methods are used for testing KSwitch and HDec separately
    ciphertext __debug_Enc_with0(plaintext m)  { return swhe0.Enc(m); }
    plaintext  __debug_Dec_with0(ciphertext C) { return swhe0.Dec(C); }
private:
    void Kg_kswitch()
    {
        int n0 = swhe0.pp.n, n1 = swhe1.pp.n;
        int ell1 = swhe1.pp.ell;
        int beta0 = swhe0.pp.beta, ell0 = swhe0.pp.ell;
        long long q0 = swhe0.pp.q;

        for (int i = 1; i < n1; ++i)
        {
            matrix A_transpose(/*rows*/ ell1, /*columns*/ n0 - 1);
            column_vector e(/*rows*/ ell1, /*columns*/ 1);

            A_transpose.randomize(0, q0 - 1, q0);
            e.randomize(-beta0, beta0, q0);
            column_vector x(/*rows*/ ell1, /*columns*/ 1);
            x.fillWithFractions(swhe1.s.getCell(i, 1), ell0 - ell1, q0);
            column_vector b = x.add(e.add(A_transpose.mult(swhe0.s, q0), q0), q0);
            matrix C = A_transpose.concatenate(b, /*horizontal*/ true);            
            eval_key.eval_key_kswitch.C[i] = C;
        }
    }

    void Kg_hdec()
    {
        int n0 = swhe0.pp.n, ell0 = swhe0.pp.ell;
        long long q0 = swhe0.pp.q;

        for (int i = 1; i < n0; ++i)
            for (int j = 0; j < ell0; ++j)
            {
                long long x = modMult(swhe0.s.getCell(i, 1), ((long long) 1 << j), q0);
                for (long long k = 0; k < q0; ++k)
                    if (k == x)
                        eval_key.eval_key_hdec.e[i][j][k] = swhe1.Enc(1);
                    else 
                        eval_key.eval_key_hdec.e[i][j][k] = swhe1.Enc(0);
            }
    }

    evk  eval_key;
    SWHE swhe0, swhe1;
};

matrix KSwitch(evk_kswitch eval_key_kswitch, ciphertext C)
{
    int n0 = eval_key_kswitch.pp0.n, n1 = C.pp.n, ell1 = C.pp.ell;
    int ell0 = eval_key_kswitch.pp0.ell;
    long long q0 = eval_key_kswitch.pp0.q, q1 = C.pp.q;

    row_vector c = C.m.extractRowVector(n1 * ell1);
    long long r = modMultPow2(c.getCell(1, n1), ell0 - ell1, q0);
    row_vector res(/*rows*/ 1, /*columns*/ n0);
    res.setCell(1, n0, r);
    for (int i = 1; i < n1; ++i)
    {
        long long c_i = c.getCell(1, i);
        scalar tmp(1, 1);
        tmp.setCell(1, 1, c_i);
        row_vector a = tmp.splitInHorizontalBitLayerConcatenation(ell1);
        row_vector x = a.mult(eval_key_kswitch.C[i], q0);
        x.negate(q0);
        res = res.add(x, q0);
    }
    return res;
}

ciphertext HDec(evk_hdec eval_key_hdec, row_vector c)
{
    int n1 = eval_key_hdec.pp1.n, ell1 = eval_key_hdec.pp1.ell;
    long long q1 = eval_key_hdec.pp1.q, q0 = eval_key_hdec.pp0.q;
    int n0 = eval_key_hdec.pp0.n, ell0 = eval_key_hdec.pp0.ell;
    
    matrix empty_m(/*rows*/ n1*ell1, /*columns*/ n1);
    ciphertext empty_c(eval_key_hdec.pp1, empty_m, 0);
    vector<ciphertext> C(q0, empty_c);
    bool C_empty = true;

    /*
    printf("\nSummary:");
    for (int i = 1; i < n0; ++i)
        for (int j = 0; j < ell0; ++j)
        {
            printf("\n");
            for (int k = 0; k < q0; ++k)
                printf("%d", eval_key_hdec.swhe1.Dec(eval_key_hdec.e[i][j][k]));
        }
    */

    for (int i = 1; i < n0; ++i)
        for (int j = 0; j < ell0; ++j)
        {
            long long x = c.getCell(1, i);
            x = (x >> j) % 2;
            if (x == 1)
            {
                vector<ciphertext> C_star(q0, empty_c);
                if (C_empty)
                {
                    for (long long k = 0; k < q0; ++k)
                        C_star[k] = XOR(C_star[k], eval_key_hdec.e[i][j][k]);
                }
                else
                {
                    for (long long u = 0; u < q0; ++u)
                        for (long long v = 0; v < q0; ++v)
                        {
                            long long k = (u + v) % q0;
                            ciphertext C_prod = AND(C[u], eval_key_hdec.e[i][j][v]);
                            C_star[k] = XOR(C_star[k], C_prod);
                        }
                }
                C = C_star;
                C_empty = false;

                /*
                printf("\nIntermediate result:");
                for (int i = 0; i < q0; ++i)
                    printf("%d", eval_key_hdec.swhe1.Dec(C[i]));
                */
            }
        }
       
    /*
    printf("\nResult:");
    for (int i = 0; i < q0; ++i)
        printf("%d", eval_key_hdec.swhe1.Dec(C[i]));
    printf("\n");
    */

    ciphertext res(eval_key_hdec.pp1, empty_m, 0);
    long long c_n = c.getCell(1, n0);
    for (long long i = 0; i < q0; ++i)
    {
        long long x = (q0 + c_n - i) % q0;
        if (q0 /4 < x && x < (q0/4)*3)
            res = XOR(res, C[i]);
    }

    return res;
}

ciphertext Recrypt(evk eval_key, ciphertext C)
{
    matrix c = KSwitch(eval_key.eval_key_kswitch, C);
    return HDec(eval_key.eval_key_hdec, c);
}

int main()
{
    std::uniform_int_distribution<> dist(0, 1);

    clock_t start_time;
    {
        // FHE fhe(/*n*/ 2, /*ell0*/ 6, /*ell1*/ 20, /*beta*/ 1); // Passed 20 tests, each 9-37 sec (average ~18 sec)
        FHE fhe(/*n*/ 4, /*ell0*/ 6, /*ell1*/ 20, /*beta*/ 1); // Passed 3 tests, 520 - 1250 sec (average ~800 sec)
        evk eval_key = fhe.get_evk();
        for (int test = 0; test < 20; ++test)
        {
            int m0 = dist(engine), m1 = dist(engine);
            ciphertext C0 = fhe.Enc(m0), C1 = fhe.Enc(m1);

            START_TIME;
            ciphertext andC = AND(C0, C1);
            ciphertext nandC = NOT(andC);
            ciphertext rec_nandC = Recrypt(eval_key, nandC);
            FINISH_TIME;

            if (fhe.Dec(rec_nandC) == 1 - (m0 * m1))
                printf("Passed NAND time test for beta=1, n=2.\n");
            else
            {
                printf("Incorrect parameters!\n");
                break;
            }
        }
    }
    return 0;

    {
        FHE fhe(/*n*/ 10, /*ell0*/ 6, /*ell1*/ 2, /*beta*/ 0);
        evk eval_key = fhe.get_evk();
        int m0 = dist(engine), m1 = dist(engine);
        ciphertext C0 = fhe.Enc(m0), C1 = fhe.Enc(m1);

        START_TIME;
        ciphertext andC = AND(C0, C1);
        ciphertext nandC = NOT(andC);
        ciphertext rec_nandC = Recrypt(eval_key, nandC);
        FINISH_TIME;

        if (fhe.Dec(rec_nandC) == 1 - (m0 * m1))
            printf("Passed NAND time test for beta=0, n=10.\n");
    }

    printf("Testing the SWHE scheme: ");
    SWHE swhe(/*n*/ 10, /*ell*/ 8); // max ell is 62
    for (int i = 0; i < 10; ++i)
    {
        int m0 = dist(engine), m1 = dist(engine);
        ciphertext C0 = swhe.Enc(m0);
        ciphertext C1 = swhe.Enc(m1);
        ciphertext NOT_C0 = NOT(C0);
        ciphertext NOT_C1 = NOT(C1);
        ciphertext XOR_C0_C1 = XOR(C0, C1);
        ciphertext AND_C0_C1 = AND(C0, C1);
        bool ok = (swhe.Dec(C0) == m0
                && swhe.Dec(C1) == m1
                && swhe.Dec(NOT_C0) == (1 - m0)
                && swhe.Dec(NOT_C1) == (1 - m1)
                && swhe.Dec(XOR_C0_C1) == (m0 + m1) % 2
                && swhe.Dec(AND_C0_C1) == m0 * m1);
        if (ok) printf("+"); else printf("-");
    }
    printf("\n");

    // FHE fhe(/*n*/ 2, /*ell0*/ 4, /*ell1*/ 2, /*beta*/ 0);  // - Generally OK, sometimes fails [n=2 and c = (0, 8) would fail]
    FHE fhe(/*n*/ 3, /*ell0*/ 4, /*ell1*/ 2, /*beta*/ 0);  // - Generally OK, sometimes fails [n=3 and c = (0, 0, *) would fail]
    // FHE fhe(/*n*/ 4, /*ell0*/ 5, /*ell1*/ 2, /*beta*/ 0);  // - OK
    // FHE fhe(/*n*/ 10, /*ell0*/ 6, /*ell1*/ 2, /*beta*/ 0); // - OK, but slow
    // FHE fhe(/*n*/ 2, /*ell0*/ 9, /*ell1*/ 40, /*beta*/ 1); // - too slow
    evk eval_key = fhe.get_evk();

    int n0 = eval_key.eval_key_hdec.pp0.n, n1 = eval_key.eval_key_hdec.pp1.n;
    int ell0 = eval_key.eval_key_hdec.pp0.ell, ell1 = eval_key.eval_key_hdec.pp1.ell;
    bool ok = true;

    printf("Testing FHE KSwitch: ");
    for (int i = 0; i < 10; ++i)
    {
        int m = dist(engine);
        ciphertext C = fhe.Enc(m);
        plaintext m_chk = fhe.Dec(C);
        assert(m == m_chk, "Encryption is incorrect.");
        row_vector C_kswitch = KSwitch(eval_key.eval_key_kswitch, C);
        matrix C_kswitch_m = matrix(/*rows*/ n0 * ell0, /*columns*/ n0);
        for (int i = 1; i <= n0; ++i)
            C_kswitch_m.setCell(n0*ell0, i, C_kswitch.getCell(1, i));
        ciphertext C_kswitch_c(eval_key.eval_key_hdec.pp0, C_kswitch_m, 0); //TODO: set up the correct noise (instead of 0)
        ok = (fhe.__debug_Dec_with0(C_kswitch_c) == m);
        if (ok) printf("+"); else printf("-");
    }
    printf("\n");

    
    printf("Testing FHE HDec: ");
    for (int i = 0; i < 10; ++i)
    {
        int m = dist(engine);
        ciphertext C = fhe.__debug_Enc_with0(m);
        plaintext m_chk = fhe.__debug_Dec_with0(C);
        assert(m == m_chk, "Encryption is incorrect.");
        row_vector c = C.m.extractRowVector(n0 * ell0);
        ciphertext C_hdec = HDec(eval_key.eval_key_hdec, c);
        bool ok = (fhe.Dec(C_hdec) == m);
        if (ok) printf("+"); else printf("-");
    }
    printf("\n");

    printf("Testing the FHE scheme: ");
    for (int i = 0; i < 10; ++i)
    {
        int m0 = dist(engine), m1 = dist(engine);
        ciphertext C0 = fhe.Enc(m0);
        ciphertext C1 = fhe.Enc(m1); 
        ciphertext NOT_C0 = Recrypt(eval_key, NOT(C0));
        ciphertext NOT_C1 = Recrypt(eval_key, NOT(C1));
        ciphertext XOR_C0_C1 = Recrypt(eval_key, XOR(C0, C1));
        ciphertext AND_C0_C1 = Recrypt(eval_key, AND(C0, C1));
        bool ok = (fhe.Dec(C0) == m0
                && fhe.Dec(C1) == m1
                && fhe.Dec(NOT_C0) == (1 - m0)
                && fhe.Dec(NOT_C1) == (1 - m1)
                && fhe.Dec(XOR_C0_C1) == (m0 + m1) % 2
                && fhe.Dec(AND_C0_C1) == m0 * m1);
        if (ok) printf("+"); else printf("-");
    }
    printf("\n");

	return 0;
}

/*
beta = sqrt(n)

209268315701.47467
2 2 10 42 20.97152 4 min 18 MB
4369200843606.928
3 2 11 45 138.412032 73 min 147 MB
59373627899904.0
4 2 12 48 805.306368 17 h 864 MB

beta = 1

42707819530.9132
2 1 9 40 4.718592 1 min 8 MB
865013500350.4625
3 1 10 42 31.45728 15 min 59 MB
11958799564800.0
4 1 11 45 184.549376 4 h 349 MB

beta = 0

310864543.10519236
10 0 6 2 0.24576 1 min 1 MB
*/
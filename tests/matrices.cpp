// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/matrices.h>
#include <math_h/vectortransformations.h>
using namespace std;
using namespace MathTemplates;
const double epsilon = 0.0000000001;
TEST(Matrix, zero1)
{
    EXPECT_EQ(
        ZERO<1>(), line(0.)
    );
    EXPECT_EQ(
        ZERO<2>(), lines(
            desCartes(0., 0.),
            desCartes(0., 0.)
        )
    );
    EXPECT_EQ(
        ZERO<3>(), lines(
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.)
        )
    );
}
TEST(Matrix, zero2)
{
    EXPECT_EQ(
        ZERO<1>(), column(0.)
    );
    EXPECT_EQ(
        ZERO<2>(), columns(
            desCartes(0., 0.),
            desCartes(0., 0.)
        )
    );
    EXPECT_EQ(
        ZERO<3>(), columns(
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.),
            desCartes(0., 0., 0.)
        )
    );
}
TEST(Matrix, one1)
{
    EXPECT_EQ(
        ONE<1>(), line(1.)
    );
    EXPECT_EQ(
        ONE<2>(), lines(
            desCartes(1., 0.),
            desCartes(0., 1.)
        )
    );
    EXPECT_EQ(
        ONE<3>(), lines(
            desCartes(1., 0., 0.),
            desCartes(0., 1., 0.),
            desCartes(0., 0., 1.)
        )
    );
}
TEST(Matrix, one2)
{
    EXPECT_EQ(
        ONE<1>(), column(1.)
    );
    EXPECT_EQ(
        ONE<2>(), columns(
            desCartes(1., 0.),
            desCartes(0., 1.)
        )
    );
    EXPECT_EQ(
        ONE<3>(), columns(
            desCartes(1., 0., 0.),
            desCartes(0., 1., 0.),
            desCartes(0., 0., 1.)
        )
    );
}
TEST(Matrix, lines_columns)
{
    EXPECT_EQ(
        lines(
            desCartes(1, 2, 3),
            desCartes(4, 5, 6),
            desCartes(7, 8, 9)
        ),
        columns(
            desCartes(1, 4, 7),
            desCartes(2, 5, 8),
            desCartes(3, 6, 9)
        )
    );
    EXPECT_EQ(
        lines(
            desCartes(1, 2, 3),
            desCartes(4, 5, 6)
        ),
        columns(
            desCartes(1, 4),
            desCartes(2, 5),
            desCartes(3, 6)
        )
    );
}
TEST(Matrix, elements)
{
    const auto M=lines(
	desCartes(1,2,3,4),
	desCartes(5,6,7,8),
	desCartes(9,0,1,2)
    );
    const auto&e11=M.element<1,1>();
    const auto&e12=M.element<1,2>();
    const auto&e13=M.element<1,3>();
    const auto&e14=M.element<1,4>();
    const auto&e21=M.element<2,1>();
    const auto&e22=M.element<2,2>();
    const auto&e23=M.element<2,3>();
    const auto&e24=M.element<2,4>();
    const auto&e31=M.element<3,1>();
    const auto&e32=M.element<3,2>();
    const auto&e33=M.element<3,3>();
    const auto&e34=M.element<3,4>();
    EXPECT_EQ(e11,1);
    EXPECT_EQ(e12,2);
    EXPECT_EQ(e13,3);
    EXPECT_EQ(e14,4);
    EXPECT_EQ(e21,5);
    EXPECT_EQ(e22,6);
    EXPECT_EQ(e23,7);
    EXPECT_EQ(e24,8);
    EXPECT_EQ(e31,9);
    EXPECT_EQ(e32,0);
    EXPECT_EQ(e33,1);
    EXPECT_EQ(e34,2);
}
TEST(Matrix, mul1)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<1>(RG) * M(RG);
        EXPECT_TRUE(desCartes(0.).CloseTo(ZERO<1>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<1>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<1>() * (ONE<1>()*v), epsilon));
        const auto R1 = randomIsotropic<1>(RG).Rotations();
        const auto R2 = randomIsotropic<1>(RG).Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(Matrix, mul2)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<2>(RG) * M(RG);
        EXPECT_TRUE(zero().CloseTo(ZERO<2>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<2>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<2>() * (ONE<2>()*v), epsilon));
        const auto R1 = randomIsotropic<2>(RG).Rotations();
        const auto R2 = randomIsotropic<2>(RG).Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(Matrix, mul3)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<3>(RG) * M(RG);
        EXPECT_TRUE(Zero().CloseTo(ZERO<3>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<3>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<3>() * (ONE<3>()*v), epsilon));
        const auto R1 = randomIsotropic<3>(RG).Rotations();
        const auto R2 = randomIsotropic<3>(RG).Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(Matrix, mul32)
{
    RANDOM RG;
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto a = randomIsotropic<3>(RG) * M(RG);
        const auto b = randomIsotropic<3>(RG) * M(RG);
        const auto c = randomIsotropic<3>(RG) * M(RG);
        const auto d = randomIsotropic<3>(RG) * M(RG);
	EXPECT_EQ(
	    lines(a,b)*columns(c,d),
	    lines(desCartes(a*c,a*d),desCartes(b*c,b*d))
	);
    }
}
#ifdef ____optimized_version_of_matrices_h_____
TEST(Matrix, det1)
{
    const auto M=line(3.0);
    EXPECT_EQ(3,M.Determinant());
}
TEST(Matrix, det2)
{
    EXPECT_EQ(lines(
	desCartes(1.0,2.0),
	desCartes(3.0,4.0)
    ).Determinant(),-2);
    EXPECT_EQ(lines(
	desCartes(1.0,2.0),
	desCartes(2.0,4.0)
    ).Determinant(),0);
    EXPECT_EQ(lines(
	desCartes(1.0,0.0),
	desCartes(0.0,1.0)
    ).Determinant(),1);
}
TEST(Matrix, det3)
{
    EXPECT_EQ(lines(
	desCartes(1.,2.,3.),
	desCartes(3.,4.,5.),
	desCartes(4.,5.,6.)
    ).Determinant(),0);
    EXPECT_EQ(lines(
	desCartes(1,0,0),
	desCartes(0,1,0),
	desCartes(0,0,1)
    ).Determinant(),1);
    EXPECT_EQ(lines(
	desCartes(1.,0.,0.),
	desCartes(0.,1.,0.),
	desCartes(0.,0.,1.)
    ).Determinant(),1);
}
TEST(Matrix, minor3)
{
    auto M=lines(
	    desCartes(1.,2.,3.),
	    desCartes(3.,4.,5.),
	    desCartes(4.,5.,6.)
	).GetMinor<1,1>();
    EXPECT_EQ(
	M,
	lines(
	    desCartes(4.,5.),
	    desCartes(5.,6.)
	)
    );
        M=lines(
	    desCartes(1.,2.,3.),
	    desCartes(3.,4.,5.),
	    desCartes(4.,5.,6.)
	).GetMinor<1,2>();

    EXPECT_EQ(
	M,
	lines(
	    desCartes(3.,5.),
	    desCartes(4.,6.)
	)
    );
        M=lines(
	    desCartes(1.,2.,3.),
	    desCartes(3.,4.,5.),
	    desCartes(4.,5.,6.)
	).GetMinor<1,3>();

    EXPECT_EQ(
	M,
	lines(
	    desCartes(3.,4.),
	    desCartes(4.,5.)
	)
    );
        M=lines(
	    desCartes(1.,2.,3.),
	    desCartes(3.,4.,5.),
	    desCartes(4.,5.,6.)
	).GetMinor<2,1>();

    EXPECT_EQ(
	M,
	lines(
	    desCartes(2.,3.),
	    desCartes(5.,6.)
	)
    );
        M=lines(
	    desCartes(1.,2.,3.),
	    desCartes(3.,4.,5.),
	    desCartes(4.,5.,6.)
	).GetMinor<2,2>();

    EXPECT_EQ(
	M,
	lines(
	    desCartes(1.,3.),
	    desCartes(4.,6.)
	)
    );
        M=lines(
	    desCartes(1.,2.,3.),
	    desCartes(3.,4.,5.),
	    desCartes(4.,5.,6.)
	).GetMinor<2,3>();
    EXPECT_EQ(
	M,
	lines(
	    desCartes(1.,2.),
	    desCartes(4.,5.)
	)
    );
}
TEST(Matrix, InsertColumn)
{
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.,3.),
	    desCartes(3.,4.,5.),
	    desCartes(4.,5.,6.)
	),
	lines(
	    desCartes(1.,3.),
	    desCartes(3.,5.),
	    desCartes(4.,6.)
	).InsertColumn<2>(desCartes(2.,4.,5.))
    );
}
TEST(Matrix, InsertRow1)
{
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.)
	).InsertRow<1>(desCartes(0.,0.)),
	lines(
	    desCartes(0.,0.),
	    desCartes(1.,2.)
	)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.)
	).InsertRow<2>(desCartes(0.,0.)),
	lines(
	    desCartes(1.,2.),
	    desCartes(0.,0.)
	)
    );
}
TEST(Matrix, InsertRow2)
{
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.)
	).InsertRow<1>(desCartes(0.,0.)),
	lines(
	    desCartes(0.,0.),
	    desCartes(1.,2.),
	    desCartes(3.,4.)
	)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.)
	).InsertRow<2>(desCartes(0.,0.)),
	lines(
	    desCartes(1.,2.),
	    desCartes(0.,0.),
	    desCartes(3.,4.)
	)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.)
	).InsertRow<3>(desCartes(0.,0.)),
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(0.,0.)
	)
    );
}
TEST(Matrix, InsertRow3)
{
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(5.,6.)
	).InsertRow<1>(desCartes(0.,0.)),
	lines(
	    desCartes(0.,0.),
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(5.,6.)
	)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(5.,6.)
	).InsertRow<2>(desCartes(0.,0.)),
	lines(
	    desCartes(1.,2.),
	    desCartes(0.,0.),
	    desCartes(3.,4.),
	    desCartes(5.,6.)
	)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(5.,6.)
	).InsertRow<3>(desCartes(0.,0.)),
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(0.,0.),
	    desCartes(5.,6.)
	)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(5.,6.)
	).InsertRow<4>(desCartes(0.,0.)),
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.),
	    desCartes(5.,6.),
	    desCartes(0.,0.)
	)
    );
}
TEST(Matrix,Cramer1)
{
    EXPECT_EQ(line(2).Cramer(desCartes(6)),desCartes(3));
    EXPECT_EQ(line(1).Cramer(desCartes(1)),desCartes(1));
    EXPECT_EQ(line(1).Cramer(desCartes(0)),desCartes(0));
    EXPECT_ANY_THROW(line(0).Cramer(desCartes(1)));
}
TEST(Matrix,Cramer2)
{
    EXPECT_EQ(
	lines(
	    desCartes(1.,2.),
	    desCartes(3.,4.)
	).Cramer(desCartes(1.,1.)),
	desCartes(-1.,1.)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,0.),
	    desCartes(0.,1.)
	).Cramer(desCartes(1.,1.)),
	desCartes(1.,1.)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,0.),
	    desCartes(0.,1.)
	).Cramer(desCartes(2.,3.)),
	desCartes(2.,3.)
    );
    EXPECT_EQ(
	lines(
	    desCartes(0.,1.),
	    desCartes(1.,0.)
	).Cramer(desCartes(2.,3.)),
	desCartes(3.,2.)
    );
    EXPECT_ANY_THROW(
	lines(
	    desCartes(1.,1.),
	    desCartes(0.,0.)
	).Cramer(desCartes(1.,1.))
    );
}
TEST(Matrix,Cramer3)
{
    EXPECT_ANY_THROW(
	lines(
	    desCartes(1.,0.,1.),
	    desCartes(0.,1.,0.),
	    desCartes(0.,0.,0.)
	).Cramer(desCartes(1.,1.,1.))
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,0.,0.),
	    desCartes(0.,1.,0.),
	    desCartes(0.,0.,1.)
	).Cramer(desCartes(1.,1.,1.)),
	desCartes(1.,1.,1.)
    );
    EXPECT_EQ(
	lines(
	    desCartes(1.,0.,0.),
	    desCartes(0.,1.,0.),
	    desCartes(0.,0.,1.)
	).Cramer(desCartes(2.,3.,5.)),
	desCartes(2.,3.,5.)
    );
    EXPECT_EQ(
	lines(
	    desCartes(0.,0.,1.),
	    desCartes(0.,1.,0.),
	    desCartes(1.,0.,0.)
	).Cramer(desCartes(2.,3.,5.)),
	desCartes(5.,3.,2.)
    );
}
#endif

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
        ZERO<1>(), row(0.)
    );
    EXPECT_EQ(
        ZERO<2>(), rows(
            vec(0., 0.),
            vec(0., 0.)
        )
    );
    EXPECT_EQ(
        ZERO<3>(), rows(
            vec(0., 0., 0.),
            vec(0., 0., 0.),
            vec(0., 0., 0.)
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
            vec(0., 0.),
            vec(0., 0.)
        )
    );
    EXPECT_EQ(
        ZERO<3>(), columns(
            vec(0., 0., 0.),
            vec(0., 0., 0.),
            vec(0., 0., 0.)
        )
    );
}
TEST(Matrix, one1)
{
    EXPECT_EQ(
        ONE<1>(), row(1.)
    );
    EXPECT_EQ(
        ONE<2>(), rows(
            vec(1., 0.),
            vec(0., 1.)
        )
    );
    EXPECT_EQ(
        ONE<3>(), rows(
            vec(1., 0., 0.),
            vec(0., 1., 0.),
            vec(0., 0., 1.)
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
            vec(1., 0.),
            vec(0., 1.)
        )
    );
    EXPECT_EQ(
        ONE<3>(), columns(
            vec(1., 0., 0.),
            vec(0., 1., 0.),
            vec(0., 0., 1.)
        )
    );
}
TEST(Matrix, rows_columns)
{
    EXPECT_EQ(
        rows(
            vec(1, 2, 3),
            vec(4, 5, 6),
            vec(7, 8, 9)
        ),
        columns(
            vec(1, 4, 7),
            vec(2, 5, 8),
            vec(3, 6, 9)
        )
    );
    EXPECT_EQ(
        rows(
            vec(1, 2, 3),
            vec(4, 5, 6)
        ),
        columns(
            vec(1, 4),
            vec(2, 5),
            vec(3, 6)
        )
    );
}
TEST(Matrix, elements)
{
    const auto M=rows(
	vec(1,2,3,4),
	vec(5,6,7,8),
	vec(9,0,1,2)
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
    const auto&r1=M.row<1>();
    const auto&r2=M.row<2>();
    const auto&r3=M.row<3>();
    EXPECT_EQ(r1,vec(1,2,3,4));
    EXPECT_EQ(r2,vec(5,6,7,8));
    EXPECT_EQ(r3,vec(9,0,1,2));
}
TEST(Matrix, mul1)
{
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<1>() * M();
        EXPECT_TRUE(vec(0.).CloseTo(ZERO<1>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<1>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<1>() * (ONE<1>()*v), epsilon));
        const auto R1 = randomIsotropic<1>().Rotations();
        const auto R2 = randomIsotropic<1>().Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(Matrix, mul2)
{
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<2>() * M();
        EXPECT_TRUE(zero().CloseTo(ZERO<2>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<2>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<2>() * (ONE<2>()*v), epsilon));
        const auto R1 = randomIsotropic<2>().Rotations();
        const auto R2 = randomIsotropic<2>().Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(Matrix, mul3)
{
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto v = randomIsotropic<3>() * M();
        EXPECT_TRUE(Zero().CloseTo(ZERO<3>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<3>()*v, epsilon));
        EXPECT_TRUE(v.CloseTo(ONE<3>() * (ONE<3>()*v), epsilon));
        const auto R1 = randomIsotropic<3>().Rotations();
        const auto R2 = randomIsotropic<3>().Rotations();
        EXPECT_TRUE((R1 * (R2 * v)).CloseTo((R1 * R2)*v, epsilon));
    }
}
TEST(Matrix, mul32)
{
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto a = randomIsotropic<3>() * M();
        const auto b = randomIsotropic<3>() * M();
        const auto c = randomIsotropic<3>() * M();
        const auto d = randomIsotropic<3>() * M();
	EXPECT_EQ(
	    rows(a,b)*columns(c,d),
	    rows(vec(a*c,a*d),vec(b*c,b*d))
	);
    }
}
TEST(Matrix, AddColumns)
{
    EXPECT_TRUE(row(1).AddColumn(vec(4))==row(1,4));
    EXPECT_TRUE(row(1).AddColumns(columns(vec(4),vec(5)))==row(1,4,5));
    EXPECT_TRUE(row(1).AddColumns(columns(vec(4),vec(5),vec(6)))==row(1,4,5,6));
    EXPECT_TRUE(row(1,2).AddColumn(vec(4))==row(1,2,4));
    EXPECT_TRUE(row(1,2).AddColumns(columns(vec(4),vec(5)))==row(1,2,4,5));
    EXPECT_TRUE(row(1,2).AddColumns(columns(vec(4),vec(5),vec(6)))==row(1,2,4,5,6));

    EXPECT_TRUE(rows(vec(1),vec(2)).AddColumn(vec(4,5))
	    ==rows(vec(1,4),vec(2,5)));
    EXPECT_TRUE(rows(vec(1),vec(2)).AddColumns(columns(vec(4,5),vec(5,7)))
	    ==rows(vec(1,4,5),vec(2,5,7)));
    EXPECT_TRUE(rows(vec(1),vec(2)).AddColumns(columns(vec(4,5),vec(5,6),vec(6,7)))
	    ==rows(vec(1,4,5,6),vec(2,5,6,7)));
    EXPECT_TRUE(rows(vec(1,2),vec(2,3)).AddColumn(vec(4,5))
	    ==rows(vec(1,2,4),vec(2,3,5)));
    EXPECT_TRUE(rows(vec(1,2),vec(2,3)).AddColumns(columns(vec(4,5),vec(9,6)))
	    ==rows(vec(1,2,4,9),vec(2,3,5,6)));
    EXPECT_TRUE(rows(vec(1,2),vec(2,3)).AddColumns(columns(vec(4,5),vec(5,6),vec(6,7)))
	    ==rows(vec(1,2,4,5,6),vec(2,3,5,6,7)));

    EXPECT_TRUE(
	rows(vec(1,2),vec(3,4)).AddColumns(rows(vec(1,2),vec(3,4)))
	    ==rows(vec(1,2,1,2),vec(3,4,3,4))
    );
    EXPECT_TRUE(
	rows(vec(1,2)).AddColumns(row(3,4))
	    ==rows(vec(1,2,3,4))
    );
}
TEST(Matrix, AddRows)
{
    EXPECT_TRUE(rows(vec(1)).AddRow(vec(6))
	    ==rows(vec(1),vec(6)));
    EXPECT_TRUE(rows(vec(1)).AddRows(rows(vec(6),vec(7)))
	    ==rows(vec(1),vec(6),vec(7)));
    EXPECT_TRUE(rows(vec(1)).AddRows(rows(vec(6),vec(7),vec(8)))
	    ==rows(vec(1),vec(6),vec(7),vec(8)));

    EXPECT_TRUE(rows(vec(1,2)).AddRow(vec(6,8))
	    ==rows(vec(1,2),vec(6,8)));
    EXPECT_TRUE(rows(vec(1,2)).AddRows(rows(vec(6,8),vec(7,9)))
	    ==rows(vec(1,2),vec(6,8),vec(7,9)));
    EXPECT_TRUE(rows(vec(1,2)).AddRows(rows(vec(6,8),vec(7,9),vec(8,0)))
	    ==rows(vec(1,2),vec(6,8),vec(7,9),vec(8,0)));

    EXPECT_TRUE(rows(vec(1),vec(2)).AddRow(vec(6))
	    ==rows(vec(1),vec(2),vec(6)));
    EXPECT_TRUE(rows(vec(1),vec(2)).AddRows(rows(vec(6),vec(7)))
	    ==rows(vec(1),vec(2),vec(6),vec(7)));
    EXPECT_TRUE(rows(vec(1),vec(2)).AddRows(rows(vec(6),vec(7),vec(8)))
	    ==rows(vec(1),vec(2),vec(6),vec(7),vec(8)));

    EXPECT_TRUE(rows(vec(1,2),vec(3,4)).AddRow(vec(6,8))
	    ==rows(vec(1,2),vec(3,4),vec(6,8)));
    EXPECT_TRUE(rows(vec(1,2),vec(3,4)).AddRows(rows(vec(6,8),vec(7,9)))
	    ==rows(vec(1,2),vec(3,4),vec(6,8),vec(7,9)));
    EXPECT_TRUE(rows(vec(1,2),vec(3,4)).AddRows(rows(vec(6,8),vec(7,9),vec(8,0)))
	    ==rows(vec(1,2),vec(3,4),vec(6,8),vec(7,9),vec(8,0)));

    EXPECT_TRUE(
	rows(vec(1,2),vec(3,4)).AddRows(rows(vec(1,2),vec(3,4)))
	    ==rows(vec(1,2),vec(3,4),vec(1,2),vec(3,4))
    );
    EXPECT_TRUE(
	rows(vec(3,4)).AddRows(rows(vec(1,2),vec(3,4)))
	    ==rows(vec(3,4),vec(1,2),vec(3,4))
    );

}
#ifdef ____full_version_of_math_h_____
TEST(Matrix, det1)
{
    const auto M=row(3.0);
    EXPECT_EQ(3,M.Determinant());
}
TEST(Matrix, det2)
{
    EXPECT_EQ(rows(
	vec(1.0,2.0),
	vec(3.0,4.0)
    ).Determinant(),-2);
    EXPECT_EQ(rows(
	vec(1.0,2.0),
	vec(2.0,4.0)
    ).Determinant(),0);
    EXPECT_EQ(rows(
	vec(1.0,0.0),
	vec(0.0,1.0)
    ).Determinant(),1);
}
TEST(Matrix, det3)
{
    EXPECT_EQ(rows(
	vec(1.,2.,3.),
	vec(3.,4.,5.),
	vec(4.,5.,6.)
    ).Determinant(),0);
    EXPECT_EQ(rows(
	vec(1,0,0),
	vec(0,1,0),
	vec(0,0,1)
    ).Determinant(),1);
    EXPECT_EQ(rows(
	vec(1.,0.,0.),
	vec(0.,1.,0.),
	vec(0.,0.,1.)
    ).Determinant(),1);
}
TEST(Matrix, minor3)
{
    auto M=rows(
	    vec(1.,2.,3.),
	    vec(3.,4.,5.),
	    vec(4.,5.,6.)
	).GetMinor<1,1>();
    EXPECT_EQ(
	M,
	rows(
	    vec(4.,5.),
	    vec(5.,6.)
	)
    );
        M=rows(
	    vec(1.,2.,3.),
	    vec(3.,4.,5.),
	    vec(4.,5.,6.)
	).GetMinor<1,2>();

    EXPECT_EQ(
	M,
	rows(
	    vec(3.,5.),
	    vec(4.,6.)
	)
    );
        M=rows(
	    vec(1.,2.,3.),
	    vec(3.,4.,5.),
	    vec(4.,5.,6.)
	).GetMinor<1,3>();

    EXPECT_EQ(
	M,
	rows(
	    vec(3.,4.),
	    vec(4.,5.)
	)
    );
        M=rows(
	    vec(1.,2.,3.),
	    vec(3.,4.,5.),
	    vec(4.,5.,6.)
	).GetMinor<2,1>();

    EXPECT_EQ(
	M,
	rows(
	    vec(2.,3.),
	    vec(5.,6.)
	)
    );
        M=rows(
	    vec(1.,2.,3.),
	    vec(3.,4.,5.),
	    vec(4.,5.,6.)
	).GetMinor<2,2>();

    EXPECT_EQ(
	M,
	rows(
	    vec(1.,3.),
	    vec(4.,6.)
	)
    );
        M=rows(
	    vec(1.,2.,3.),
	    vec(3.,4.,5.),
	    vec(4.,5.,6.)
	).GetMinor<2,3>();
    EXPECT_EQ(
	M,
	rows(
	    vec(1.,2.),
	    vec(4.,5.)
	)
    );
}
TEST(Matrix, InsertColumns)
{
    EXPECT_EQ(
	rows(
	    vec(1.,2.,3.)
	),
	rows(
	    vec(1.,3.)
	).InsertColumn<2>(vec(2.))
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.,3.,4)
	),
	rows(
	    vec(1.,4.)
	).InsertColumns<2>(row(2.,3.))
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.,3.),
	    vec(3.,4.,5.),
	    vec(4.,5.,6.)
	),
	rows(
	    vec(1.,3.),
	    vec(3.,5.),
	    vec(4.,6.)
	).InsertColumn<2>(vec(2.,4.,5.))
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.,3.,4),
	    vec(4.,3.,2.,1)
	),
	rows(
	    vec(1.,4.),
	    vec(4.,1.)
	).InsertColumns<2>(rows(vec(2.,3.),vec(3.,2.)))
    );
}
TEST(Matrix, RemoveColumn)
{
    const auto M=rows(
	    vec(1.,2.,3.)
	).RemoveColumns<1,2>(),
	M2=rows(
	    vec(3.)
	);
    EXPECT_EQ(M,M2);
    const auto M3=rows(
	    vec(1.,2.,3.),
	    vec(4.,5.,6.)
	).RemoveColumns<1,2>(),
	M4=rows(
	    vec(3.),
	    vec(6.)
	);
    EXPECT_EQ(M3,M4);
}
TEST(Matrix, RemoveRow)
{
    const auto M=rows(
	    vec(1.,2.,3.),
	    vec(4.,5.,6.),
	    vec(7.,8.,9.),
	    vec(0.,1.,2.)
	).RemoveRows<2,2>(),
	M2=rows(
	    vec(1.,2.,3.),
	    vec(0.,1.,2.)
	);
    EXPECT_EQ(M,M2);
}
TEST(Matrix, InsertRows1)
{
    EXPECT_EQ(
	rows(
	    vec(1.,2.)
	).InsertRow<1>(vec(0.,0.)),
	rows(
	    vec(0.,0.),
	    vec(1.,2.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.)
	).InsertRows<1>(rows(vec(0.,0.),vec(3.,3.))),
	rows(
	    vec(0.,0.),
	    vec(3.,3.),
	    vec(1.,2.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.)
	).InsertRows<2>(rows(vec(0.,0.),vec(3.,3.))),
	rows(
	    vec(1.,2.),
	    vec(0.,0.),
	    vec(3.,3.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.)
	).InsertRow<2>(vec(0.,0.)),
	rows(
	    vec(1.,2.),
	    vec(0.,0.)
	)
    );
}
TEST(Matrix, InsertRows2)
{
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.)
	).InsertRow<1>(vec(0.,0.)),
	rows(
	    vec(0.,0.),
	    vec(1.,2.),
	    vec(3.,4.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.)
	).InsertRow<2>(vec(0.,0.)),
	rows(
	    vec(1.,2.),
	    vec(0.,0.),
	    vec(3.,4.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.)
	).InsertRow<3>(vec(0.,0.)),
	rows(
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(0.,0.)
	)
    );
}
TEST(Matrix, InsertRows3)
{
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(5.,6.)
	).InsertRow<1>(vec(0.,0.)),
	rows(
	    vec(0.,0.),
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(5.,6.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(5.,6.)
	).InsertRow<2>(vec(0.,0.)),
	rows(
	    vec(1.,2.),
	    vec(0.,0.),
	    vec(3.,4.),
	    vec(5.,6.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(5.,6.)
	).InsertRow<3>(vec(0.,0.)),
	rows(
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(0.,0.),
	    vec(5.,6.)
	)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(5.,6.)
	).InsertRow<4>(vec(0.,0.)),
	rows(
	    vec(1.,2.),
	    vec(3.,4.),
	    vec(5.,6.),
	    vec(0.,0.)
	)
    );
}
TEST(Matrix,Cramer1)
{
    EXPECT_EQ(row(2).Cramer(vec(6)),vec(3));
    EXPECT_EQ(row(1).Cramer(vec(1)),vec(1));
    EXPECT_EQ(row(1).Cramer(vec(0)),vec(0));
    EXPECT_ANY_THROW(row(0).Cramer(vec(1)));
}
TEST(Matrix,Cramer2)
{
    EXPECT_EQ(
	rows(
	    vec(1.,2.),
	    vec(3.,4.)
	).Cramer(vec(1.,1.)),
	vec(-1.,1.)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,0.),
	    vec(0.,1.)
	).Cramer(vec(1.,1.)),
	vec(1.,1.)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,0.),
	    vec(0.,1.)
	).Cramer(vec(2.,3.)),
	vec(2.,3.)
    );
    EXPECT_EQ(
	rows(
	    vec(0.,1.),
	    vec(1.,0.)
	).Cramer(vec(2.,3.)),
	vec(3.,2.)
    );
    EXPECT_ANY_THROW(
	rows(
	    vec(1.,1.),
	    vec(0.,0.)
	).Cramer(vec(1.,1.))
    );
}
TEST(Matrix,Cramer3)
{
    EXPECT_ANY_THROW(
	rows(
	    vec(1.,0.,1.),
	    vec(0.,1.,0.),
	    vec(0.,0.,0.)
	).Cramer(vec(1.,1.,1.))
    );
    EXPECT_EQ(
	rows(
	    vec(1.,0.,0.),
	    vec(0.,1.,0.),
	    vec(0.,0.,1.)
	).Cramer(vec(1.,1.,1.)),
	vec(1.,1.,1.)
    );
    EXPECT_EQ(
	rows(
	    vec(1.,0.,0.),
	    vec(0.,1.,0.),
	    vec(0.,0.,1.)
	).Cramer(vec(2.,3.,5.)),
	vec(2.,3.,5.)
    );
    EXPECT_EQ(
	rows(
	    vec(0.,0.,1.),
	    vec(0.,1.,0.),
	    vec(1.,0.,0.)
	).Cramer(vec(2.,3.,5.)),
	vec(5.,3.,2.)
    );
    EXPECT_ANY_THROW(
	rows(
	    vec(2.,0.,1.),
	    vec(0.,1.,0.),
	    vec(0.,0.,0.)
	).Cramer(vec(2.,3.,5.))
    );
}
TEST(Matrix, transponate1)
{
    EXPECT_EQ(
        row(1, 2,3).transponate(),
        column(1,2,3)
    );
}
TEST(Matrix, transponate2)
{
    EXPECT_EQ(
        rows(
            vec(1, 2),
            vec(3, 4),
            vec(5, 6)
        ).transponate(),
        columns(
            vec(1, 2),
            vec(3, 4),
            vec(5, 6)
        )
    );
    EXPECT_EQ(
        rows(
            vec(1, 2, 3),
            vec(4, 5, 6),
            vec(7, 8, 9)
        ).transponate(),
        rows(
            vec(1, 4, 7),
            vec(2, 5, 8),
            vec(3, 6, 9)
        )
    );
}
#endif
TEST(Matrix,diagonal)
{
    EXPECT_EQ(row(3.5).diagonal(),vec(3.5));
    EXPECT_EQ(
        rows(
            vec(1, 2, 3),
            vec(4, 5, 6),
            vec(7, 8, 9)
        ).diagonal(),
	vec(1, 5, 9)
    );
}

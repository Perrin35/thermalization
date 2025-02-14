OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.51337564) q[0];
sx q[0];
rz(1.4613419) q[0];
sx q[0];
rz(10.339633) q[0];
rz(-0.86923161) q[1];
sx q[1];
rz(-2.4183122) q[1];
sx q[1];
rz(1.7887315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76465121) q[0];
sx q[0];
rz(-0.94263715) q[0];
sx q[0];
rz(-2.2215722) q[0];
rz(-pi) q[1];
rz(-0.74341821) q[2];
sx q[2];
rz(-1.8525436) q[2];
sx q[2];
rz(-2.8670058) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1654403) q[1];
sx q[1];
rz(-0.70164499) q[1];
sx q[1];
rz(1.8681504) q[1];
rz(0.21463359) q[3];
sx q[3];
rz(-2.1691536) q[3];
sx q[3];
rz(2.2325688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4065518) q[2];
sx q[2];
rz(-2.328379) q[2];
sx q[2];
rz(1.0068033) q[2];
rz(-1.462228) q[3];
sx q[3];
rz(-1.4454602) q[3];
sx q[3];
rz(-1.538895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8524858) q[0];
sx q[0];
rz(-0.91133457) q[0];
sx q[0];
rz(-0.12595969) q[0];
rz(-1.1360629) q[1];
sx q[1];
rz(-1.2964396) q[1];
sx q[1];
rz(-2.0819285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.179716) q[0];
sx q[0];
rz(-1.5700514) q[0];
sx q[0];
rz(3.1396833) q[0];
x q[1];
rz(2.4074209) q[2];
sx q[2];
rz(-2.3245272) q[2];
sx q[2];
rz(1.2881607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48908614) q[1];
sx q[1];
rz(-2.2678657) q[1];
sx q[1];
rz(0.98784312) q[1];
x q[2];
rz(1.6206762) q[3];
sx q[3];
rz(-1.2119519) q[3];
sx q[3];
rz(1.9045496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58553592) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(0.6297099) q[2];
rz(-1.2398047) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(-2.0201717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12104163) q[0];
sx q[0];
rz(-0.23153767) q[0];
sx q[0];
rz(-1.9112021) q[0];
rz(2.4379099) q[1];
sx q[1];
rz(-0.92862248) q[1];
sx q[1];
rz(1.5603265) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4245371) q[0];
sx q[0];
rz(-1.5858264) q[0];
sx q[0];
rz(1.1771955) q[0];
rz(-0.73810021) q[2];
sx q[2];
rz(-1.7154652) q[2];
sx q[2];
rz(1.6044782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81654233) q[1];
sx q[1];
rz(-2.5245164) q[1];
sx q[1];
rz(-2.3338334) q[1];
rz(-pi) q[2];
rz(-1.1915474) q[3];
sx q[3];
rz(-0.82189188) q[3];
sx q[3];
rz(2.3683907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8065717) q[2];
sx q[2];
rz(-1.0806012) q[2];
sx q[2];
rz(-0.054277167) q[2];
rz(1.0862167) q[3];
sx q[3];
rz(-0.79038668) q[3];
sx q[3];
rz(-0.49200341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739968) q[0];
sx q[0];
rz(-0.01570276) q[0];
sx q[0];
rz(0.98454654) q[0];
rz(1.7355851) q[1];
sx q[1];
rz(-1.5407298) q[1];
sx q[1];
rz(2.9071009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9391044) q[0];
sx q[0];
rz(-1.8456158) q[0];
sx q[0];
rz(2.2780134) q[0];
rz(-pi) q[1];
rz(1.8580397) q[2];
sx q[2];
rz(-1.8063288) q[2];
sx q[2];
rz(-1.4969431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.381466) q[1];
sx q[1];
rz(-2.1153808) q[1];
sx q[1];
rz(-1.8167956) q[1];
rz(-2.1470039) q[3];
sx q[3];
rz(-1.5458917) q[3];
sx q[3];
rz(-1.2664317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24388127) q[2];
sx q[2];
rz(-2.477024) q[2];
sx q[2];
rz(2.393874) q[2];
rz(-1.0866577) q[3];
sx q[3];
rz(-0.99435157) q[3];
sx q[3];
rz(0.35191107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20027593) q[0];
sx q[0];
rz(-2.2714697) q[0];
sx q[0];
rz(-1.5927607) q[0];
rz(-3.0039655) q[1];
sx q[1];
rz(-2.177114) q[1];
sx q[1];
rz(2.4639938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0726822) q[0];
sx q[0];
rz(-0.015791647) q[0];
sx q[0];
rz(-2.6591405) q[0];
x q[1];
rz(2.4816072) q[2];
sx q[2];
rz(-2.8980245) q[2];
sx q[2];
rz(-3.0483831) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.69285821) q[1];
sx q[1];
rz(-3.0046607) q[1];
sx q[1];
rz(-1.6137692) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26664385) q[3];
sx q[3];
rz(-0.75481774) q[3];
sx q[3];
rz(2.9637869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64348334) q[2];
sx q[2];
rz(-0.52916873) q[2];
sx q[2];
rz(0.37437487) q[2];
rz(-0.74786413) q[3];
sx q[3];
rz(-1.4948083) q[3];
sx q[3];
rz(-1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0102274) q[0];
sx q[0];
rz(-1.3656728) q[0];
sx q[0];
rz(0.22450547) q[0];
rz(-0.35047105) q[1];
sx q[1];
rz(-1.1843362) q[1];
sx q[1];
rz(1.9559466) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910584) q[0];
sx q[0];
rz(-1.6259999) q[0];
sx q[0];
rz(1.4226705) q[0];
rz(0.19059981) q[2];
sx q[2];
rz(-2.6707044) q[2];
sx q[2];
rz(2.6550271) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.396186) q[1];
sx q[1];
rz(-2.2377659) q[1];
sx q[1];
rz(2.1275669) q[1];
rz(-pi) q[2];
rz(0.66817254) q[3];
sx q[3];
rz(-2.1134317) q[3];
sx q[3];
rz(0.99110161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4993569) q[2];
sx q[2];
rz(-0.88358742) q[2];
sx q[2];
rz(-2.9767766) q[2];
rz(-0.76588255) q[3];
sx q[3];
rz(-0.82847786) q[3];
sx q[3];
rz(-0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0484475) q[0];
sx q[0];
rz(-0.026938139) q[0];
sx q[0];
rz(0.41937605) q[0];
rz(-1.5834437) q[1];
sx q[1];
rz(-2.7132468) q[1];
sx q[1];
rz(1.3126119) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93363491) q[0];
sx q[0];
rz(-2.5267753) q[0];
sx q[0];
rz(1.2864248) q[0];
rz(2.3840226) q[2];
sx q[2];
rz(-2.295864) q[2];
sx q[2];
rz(0.49882327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53322843) q[1];
sx q[1];
rz(-1.3238412) q[1];
sx q[1];
rz(-0.19503959) q[1];
x q[2];
rz(-0.26793865) q[3];
sx q[3];
rz(-1.4456141) q[3];
sx q[3];
rz(-2.1942558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7340362) q[2];
sx q[2];
rz(-2.4702256) q[2];
sx q[2];
rz(2.7835795) q[2];
rz(0.19896209) q[3];
sx q[3];
rz(-2.3301061) q[3];
sx q[3];
rz(0.097804047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1426549) q[0];
sx q[0];
rz(-2.1079347) q[0];
sx q[0];
rz(2.348483) q[0];
rz(2.2456887) q[1];
sx q[1];
rz(-1.221012) q[1];
sx q[1];
rz(-2.6633247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5332141) q[0];
sx q[0];
rz(-1.694931) q[0];
sx q[0];
rz(2.690655) q[0];
rz(-pi) q[1];
rz(1.2094648) q[2];
sx q[2];
rz(-1.981297) q[2];
sx q[2];
rz(-2.362006) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7859808) q[1];
sx q[1];
rz(-2.0352239) q[1];
sx q[1];
rz(2.8465438) q[1];
rz(-0.10981126) q[3];
sx q[3];
rz(-1.0581019) q[3];
sx q[3];
rz(1.5231665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.521495) q[2];
sx q[2];
rz(-1.5767117) q[2];
sx q[2];
rz(-3.133797) q[2];
rz(1.9589641) q[3];
sx q[3];
rz(-1.7595485) q[3];
sx q[3];
rz(1.8462605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1496534) q[0];
sx q[0];
rz(-2.6872334) q[0];
sx q[0];
rz(0.40865189) q[0];
rz(-2.0463792) q[1];
sx q[1];
rz(-1.6588255) q[1];
sx q[1];
rz(0.85831395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82018554) q[0];
sx q[0];
rz(-1.5044065) q[0];
sx q[0];
rz(1.4231349) q[0];
x q[1];
rz(-1.8377531) q[2];
sx q[2];
rz(-2.6745195) q[2];
sx q[2];
rz(-0.89887039) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.146119) q[1];
sx q[1];
rz(-1.766537) q[1];
sx q[1];
rz(-0.11222307) q[1];
rz(-0.32118843) q[3];
sx q[3];
rz(-0.62231718) q[3];
sx q[3];
rz(2.8701607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3506763) q[2];
sx q[2];
rz(-0.75296593) q[2];
sx q[2];
rz(-0.6616627) q[2];
rz(0.800313) q[3];
sx q[3];
rz(-1.1626264) q[3];
sx q[3];
rz(0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8372339) q[0];
sx q[0];
rz(-2.9999314) q[0];
sx q[0];
rz(-0.71453553) q[0];
rz(0.47572687) q[1];
sx q[1];
rz(-0.58554244) q[1];
sx q[1];
rz(0.48019662) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6960951) q[0];
sx q[0];
rz(-1.0871743) q[0];
sx q[0];
rz(-2.4726443) q[0];
rz(1.1473898) q[2];
sx q[2];
rz(-0.69726476) q[2];
sx q[2];
rz(0.84691511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.811753) q[1];
sx q[1];
rz(-1.8564714) q[1];
sx q[1];
rz(2.4655254) q[1];
rz(-pi) q[2];
rz(-2.8586724) q[3];
sx q[3];
rz(-1.1503387) q[3];
sx q[3];
rz(0.69129163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6314038) q[2];
sx q[2];
rz(-1.1978585) q[2];
sx q[2];
rz(0.27604827) q[2];
rz(-2.7035642) q[3];
sx q[3];
rz(-1.5166401) q[3];
sx q[3];
rz(0.62694222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20090564) q[0];
sx q[0];
rz(-1.154366) q[0];
sx q[0];
rz(0.05703297) q[0];
rz(3.1055462) q[1];
sx q[1];
rz(-1.6135975) q[1];
sx q[1];
rz(1.5560908) q[1];
rz(2.0644319) q[2];
sx q[2];
rz(-1.3282599) q[2];
sx q[2];
rz(-1.9400772) q[2];
rz(-0.48595002) q[3];
sx q[3];
rz(-0.28578352) q[3];
sx q[3];
rz(-1.8430229) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

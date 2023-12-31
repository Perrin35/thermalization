OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(2.8154362) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.125995) q[0];
sx q[0];
rz(-1.3552357) q[0];
sx q[0];
rz(2.9247345) q[0];
rz(-pi) q[1];
rz(1.3165663) q[2];
sx q[2];
rz(-0.98346522) q[2];
sx q[2];
rz(1.6531144) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5427248) q[1];
sx q[1];
rz(-1.3398783) q[1];
sx q[1];
rz(-2.8938107) q[1];
rz(-0.78674973) q[3];
sx q[3];
rz(-1.5610715) q[3];
sx q[3];
rz(1.7973289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(0.88511434) q[2];
rz(0.72201133) q[3];
sx q[3];
rz(-1.4530028) q[3];
sx q[3];
rz(3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-2.0425178) q[0];
rz(0.66501578) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(0.87759334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646334) q[0];
sx q[0];
rz(-1.3898802) q[0];
sx q[0];
rz(-0.69676708) q[0];
rz(-pi) q[1];
rz(-1.7052824) q[2];
sx q[2];
rz(-1.9324537) q[2];
sx q[2];
rz(-2.550617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6192012) q[1];
sx q[1];
rz(-3.0295105) q[1];
sx q[1];
rz(-1.2680608) q[1];
x q[2];
rz(-1.1869933) q[3];
sx q[3];
rz(-1.882453) q[3];
sx q[3];
rz(0.66031885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39891222) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(1.1304643) q[2];
rz(-1.8418664) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14625064) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(-2.8515942) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(3.0677632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89798112) q[0];
sx q[0];
rz(-1.6323166) q[0];
sx q[0];
rz(-3.0520526) q[0];
rz(0.29088144) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(1.3128624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5362629) q[1];
sx q[1];
rz(-1.3932091) q[1];
sx q[1];
rz(3.0847343) q[1];
rz(1.6894433) q[3];
sx q[3];
rz(-1.9199748) q[3];
sx q[3];
rz(-0.94482869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6510216) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(-0.64669615) q[2];
rz(-2.0329287) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9639503) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(1.1886764) q[0];
rz(2.1229318) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(1.4368988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3387412) q[0];
sx q[0];
rz(-0.84206284) q[0];
sx q[0];
rz(1.6943504) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3827299) q[2];
sx q[2];
rz(-0.43735158) q[2];
sx q[2];
rz(-2.8749089) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.97836271) q[1];
sx q[1];
rz(-1.1095188) q[1];
sx q[1];
rz(-0.39050885) q[1];
rz(-pi) q[2];
rz(3.0152263) q[3];
sx q[3];
rz(-1.4411981) q[3];
sx q[3];
rz(1.6882997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.556095) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(2.0193224) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68200237) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(-0.37297747) q[0];
rz(0.22398082) q[1];
sx q[1];
rz(-1.9517027) q[1];
sx q[1];
rz(-1.3164828) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3217357) q[0];
sx q[0];
rz(-2.5279191) q[0];
sx q[0];
rz(-0.86921285) q[0];
rz(-pi) q[1];
rz(0.26222783) q[2];
sx q[2];
rz(-1.8143468) q[2];
sx q[2];
rz(1.3893931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97536918) q[1];
sx q[1];
rz(-1.6227286) q[1];
sx q[1];
rz(-0.012253461) q[1];
rz(-2.5943807) q[3];
sx q[3];
rz(-1.1037165) q[3];
sx q[3];
rz(-2.3576749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(-0.80580795) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(-1.0114975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7508115) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(-0.090963013) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(-1.3202753) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0108171) q[0];
sx q[0];
rz(-2.5512619) q[0];
sx q[0];
rz(-0.70461313) q[0];
x q[1];
rz(0.57029057) q[2];
sx q[2];
rz(-1.7208793) q[2];
sx q[2];
rz(2.1652086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4016061) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(-2.2120038) q[1];
x q[2];
rz(-2.8517013) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(2.1260726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0992574) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(-0.48505923) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(1.6962359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27959529) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(2.5860508) q[0];
rz(0.034596054) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(-1.3909891) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334675) q[0];
sx q[0];
rz(-1.9017178) q[0];
sx q[0];
rz(2.5740741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7752635) q[2];
sx q[2];
rz(-1.804367) q[2];
sx q[2];
rz(-2.4457268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3607189) q[1];
sx q[1];
rz(-2.2370173) q[1];
sx q[1];
rz(-0.16155508) q[1];
rz(-pi) q[2];
rz(-0.075918003) q[3];
sx q[3];
rz(-0.77709353) q[3];
sx q[3];
rz(-2.0457552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81327072) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(-0.39548809) q[2];
rz(-1.8528806) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(-2.18816) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(1.4377726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9490818) q[0];
sx q[0];
rz(-1.3882982) q[0];
sx q[0];
rz(-2.0927246) q[0];
x q[1];
rz(-0.17022325) q[2];
sx q[2];
rz(-1.7049689) q[2];
sx q[2];
rz(-1.2302878) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21282141) q[1];
sx q[1];
rz(-1.0911687) q[1];
sx q[1];
rz(-0.70478435) q[1];
rz(-2.0182761) q[3];
sx q[3];
rz(-1.5527225) q[3];
sx q[3];
rz(-1.7759089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64547223) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(2.9628741) q[2];
rz(0.86137613) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9361967) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(2.0478915) q[0];
rz(2.4049092) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(1.0587943) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98204389) q[0];
sx q[0];
rz(-1.0315572) q[0];
sx q[0];
rz(-0.3190785) q[0];
rz(-2.043622) q[2];
sx q[2];
rz(-1.9377922) q[2];
sx q[2];
rz(2.8851913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4150548) q[1];
sx q[1];
rz(-1.448436) q[1];
sx q[1];
rz(1.1916257) q[1];
rz(2.6685647) q[3];
sx q[3];
rz(-1.932945) q[3];
sx q[3];
rz(-0.87272296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(-2.7311834) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8695628) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(0.29125443) q[0];
rz(-2.5323396) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3647389) q[0];
sx q[0];
rz(-0.59989444) q[0];
sx q[0];
rz(2.5584496) q[0];
rz(-1.8068238) q[2];
sx q[2];
rz(-1.126295) q[2];
sx q[2];
rz(-2.1399463) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0309279) q[1];
sx q[1];
rz(-1.0193766) q[1];
sx q[1];
rz(1.5196147) q[1];
rz(-1.7581975) q[3];
sx q[3];
rz(-0.43863505) q[3];
sx q[3];
rz(-0.24500971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6802406) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(2.0001901) q[2];
rz(-1.5348148) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-1.776406) q[2];
sx q[2];
rz(-1.793135) q[2];
sx q[2];
rz(-1.8829913) q[2];
rz(-1.7678123) q[3];
sx q[3];
rz(-1.9867298) q[3];
sx q[3];
rz(0.17476535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

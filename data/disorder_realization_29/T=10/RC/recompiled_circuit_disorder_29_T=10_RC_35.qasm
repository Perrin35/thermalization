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
rz(-0.32615647) q[0];
rz(-1.1905319) q[1];
sx q[1];
rz(-1.3500554) q[1];
sx q[1];
rz(-1.5989369) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.125995) q[0];
sx q[0];
rz(-1.7863569) q[0];
sx q[0];
rz(-2.9247345) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7803696) q[2];
sx q[2];
rz(-0.63399705) q[2];
sx q[2];
rz(1.0499954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76347683) q[1];
sx q[1];
rz(-0.33707481) q[1];
sx q[1];
rz(-2.377305) q[1];
x q[2];
rz(0.013734038) q[3];
sx q[3];
rz(-0.78679689) q[3];
sx q[3];
rz(-0.2168344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(2.2564783) q[2];
rz(-0.72201133) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(1.0990748) q[0];
rz(2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(-0.87759334) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7769593) q[0];
sx q[0];
rz(-1.7517125) q[0];
sx q[0];
rz(2.4448256) q[0];
rz(-pi) q[1];
rz(-1.7052824) q[2];
sx q[2];
rz(-1.209139) q[2];
sx q[2];
rz(2.550617) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3493537) q[1];
sx q[1];
rz(-1.5374447) q[1];
sx q[1];
rz(1.4637714) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33438501) q[3];
sx q[3];
rz(-1.2063724) q[3];
sx q[3];
rz(2.3543166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7426804) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(-2.0111283) q[2];
rz(-1.8418664) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.995342) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(2.8515942) q[0];
rz(0.6668123) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(3.0677632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683257) q[0];
sx q[0];
rz(-0.10859117) q[0];
sx q[0];
rz(0.6032087) q[0];
x q[1];
rz(-1.2286085) q[2];
sx q[2];
rz(-1.295919) q[2];
sx q[2];
rz(-0.16155044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9170345) q[1];
sx q[1];
rz(-0.18637603) q[1];
sx q[1];
rz(1.264155) q[1];
rz(-1.6894433) q[3];
sx q[3];
rz(-1.9199748) q[3];
sx q[3];
rz(0.94482869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4905711) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(-1.1086639) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(-1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9639503) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(1.1886764) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(1.7046938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98722092) q[0];
sx q[0];
rz(-2.4043596) q[0];
sx q[0];
rz(-0.1371951) q[0];
x q[1];
rz(2.8145153) q[2];
sx q[2];
rz(-1.8665258) q[2];
sx q[2];
rz(1.1277652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9096133) q[1];
sx q[1];
rz(-0.59514272) q[1];
sx q[1];
rz(0.91722782) q[1];
x q[2];
rz(-3.0152263) q[3];
sx q[3];
rz(-1.7003945) q[3];
sx q[3];
rz(-1.453293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58549762) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-1.1222703) q[2];
rz(-2.1155817) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595903) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(0.37297747) q[0];
rz(2.9176118) q[1];
sx q[1];
rz(-1.9517027) q[1];
sx q[1];
rz(1.3164828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1236498) q[0];
sx q[0];
rz(-2.0262449) q[0];
sx q[0];
rz(-0.42670593) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3189795) q[2];
sx q[2];
rz(-1.8251112) q[2];
sx q[2];
rz(-3.0248259) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59606325) q[1];
sx q[1];
rz(-1.5585594) q[1];
sx q[1];
rz(1.6227325) q[1];
rz(-2.5943807) q[3];
sx q[3];
rz(-1.1037165) q[3];
sx q[3];
rz(-2.3576749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.10107723) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(-2.3357847) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(1.0114975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(0.090963013) q[0];
rz(0.85995752) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.8213173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108171) q[0];
sx q[0];
rz(-2.5512619) q[0];
sx q[0];
rz(-2.4369795) q[0];
rz(2.5713021) q[2];
sx q[2];
rz(-1.4207134) q[2];
sx q[2];
rz(2.1652086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7399866) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(2.2120038) q[1];
x q[2];
rz(-0.91464197) q[3];
sx q[3];
rz(-1.802889) q[3];
sx q[3];
rz(-2.7616012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0992574) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(-2.5406204) q[2];
rz(2.6565334) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(-1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27959529) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(-2.5860508) q[0];
rz(-3.1069966) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(1.7506036) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3334675) q[0];
sx q[0];
rz(-1.2398749) q[0];
sx q[0];
rz(2.5740741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7752635) q[2];
sx q[2];
rz(-1.804367) q[2];
sx q[2];
rz(0.69586588) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2512867) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(-2.2433953) q[1];
rz(-0.77565149) q[3];
sx q[3];
rz(-1.6240048) q[3];
sx q[3];
rz(-2.720811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.81327072) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(-2.7461046) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46273461) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(-1.6495552) q[0];
rz(0.95343268) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(1.7038201) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4576365) q[0];
sx q[0];
rz(-0.55011612) q[0];
sx q[0];
rz(-1.9253299) q[0];
x q[1];
rz(2.9713694) q[2];
sx q[2];
rz(-1.4366237) q[2];
sx q[2];
rz(1.2302878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21282141) q[1];
sx q[1];
rz(-2.050424) q[1];
sx q[1];
rz(0.70478435) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1215454) q[3];
sx q[3];
rz(-1.123395) q[3];
sx q[3];
rz(-2.945154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(-2.9628741) q[2];
rz(2.2802165) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(-0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20539595) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(-2.0478915) q[0];
rz(0.73668346) q[1];
sx q[1];
rz(-1.8700347) q[1];
sx q[1];
rz(1.0587943) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7208784) q[0];
sx q[0];
rz(-1.2982561) q[0];
sx q[0];
rz(1.0084624) q[0];
x q[1];
rz(-1.0979707) q[2];
sx q[2];
rz(-1.9377922) q[2];
sx q[2];
rz(-2.8851913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79566369) q[1];
sx q[1];
rz(-1.1946031) q[1];
sx q[1];
rz(-3.0099807) q[1];
rz(0.47302795) q[3];
sx q[3];
rz(-1.932945) q[3];
sx q[3];
rz(0.87272296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8686691) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(2.7311834) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8695628) q[0];
sx q[0];
rz(-0.27150387) q[0];
sx q[0];
rz(-0.29125443) q[0];
rz(-2.5323396) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(-1.7094918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0390022) q[0];
sx q[0];
rz(-2.0615091) q[0];
sx q[0];
rz(-1.9309994) q[0];
x q[1];
rz(0.45639313) q[2];
sx q[2];
rz(-2.6420339) q[2];
sx q[2];
rz(1.6299786) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4869625) q[1];
sx q[1];
rz(-1.527206) q[1];
sx q[1];
rz(-0.55200465) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1389144) q[3];
sx q[3];
rz(-1.6500041) q[3];
sx q[3];
rz(1.6457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6802406) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(2.0001901) q[2];
rz(-1.6067778) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(-2.6721568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5466945) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-1.3651866) q[2];
sx q[2];
rz(-1.3484577) q[2];
sx q[2];
rz(1.2586013) q[2];
rz(-2.7244444) q[3];
sx q[3];
rz(-2.683831) q[3];
sx q[3];
rz(2.8575069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

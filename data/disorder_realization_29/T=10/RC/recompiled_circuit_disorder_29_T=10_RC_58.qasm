OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0260789) q[0];
sx q[0];
rz(-1.6576515) q[0];
sx q[0];
rz(-2.8154362) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.125995) q[0];
sx q[0];
rz(-1.7863569) q[0];
sx q[0];
rz(-0.21685812) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8250263) q[2];
sx q[2];
rz(-2.1581274) q[2];
sx q[2];
rz(1.6531144) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1118288) q[1];
sx q[1];
rz(-1.811869) q[1];
sx q[1];
rz(1.3328711) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5845675) q[3];
sx q[3];
rz(-2.3574986) q[3];
sx q[3];
rz(0.23628326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2797543) q[2];
sx q[2];
rz(-0.97936169) q[2];
sx q[2];
rz(-2.2564783) q[2];
rz(-2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9279813) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-2.0425178) q[0];
rz(2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(2.2639993) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35560247) q[0];
sx q[0];
rz(-2.2539833) q[0];
sx q[0];
rz(-1.3366633) q[0];
rz(-pi) q[1];
rz(-1.7052824) q[2];
sx q[2];
rz(-1.209139) q[2];
sx q[2];
rz(-0.59097564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21786015) q[1];
sx q[1];
rz(-1.6777615) q[1];
sx q[1];
rz(3.1080493) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8072076) q[3];
sx q[3];
rz(-1.9352203) q[3];
sx q[3];
rz(0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7426804) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(2.0111283) q[2];
rz(-1.2997262) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.6931504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(0.28999844) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(-3.0677632) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4742972) q[0];
sx q[0];
rz(-1.6601666) q[0];
sx q[0];
rz(1.6325634) q[0];
rz(0.29088144) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(1.3128624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0970043) q[1];
sx q[1];
rz(-1.6267596) q[1];
sx q[1];
rz(-1.3929277) q[1];
rz(-pi) q[2];
rz(0.35145268) q[3];
sx q[3];
rz(-1.4593399) q[3];
sx q[3];
rz(0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6510216) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(-2.0329287) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17764238) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(-1.9529163) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(-1.7046938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.456159) q[0];
sx q[0];
rz(-1.478727) q[0];
sx q[0];
rz(2.409056) q[0];
rz(-2.8145153) q[2];
sx q[2];
rz(-1.2750669) q[2];
sx q[2];
rz(1.1277652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1632299) q[1];
sx q[1];
rz(-1.1095188) q[1];
sx q[1];
rz(-2.7510838) q[1];
rz(-pi) q[2];
rz(0.12636633) q[3];
sx q[3];
rz(-1.7003945) q[3];
sx q[3];
rz(1.6882997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58549762) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68200237) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(-0.37297747) q[0];
rz(2.9176118) q[1];
sx q[1];
rz(-1.9517027) q[1];
sx q[1];
rz(-1.8251098) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017942863) q[0];
sx q[0];
rz(-1.1153478) q[0];
sx q[0];
rz(2.7148867) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8793648) q[2];
sx q[2];
rz(-1.3272459) q[2];
sx q[2];
rz(1.3893931) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59606325) q[1];
sx q[1];
rz(-1.5830333) q[1];
sx q[1];
rz(1.6227325) q[1];
x q[2];
rz(-0.54721197) q[3];
sx q[3];
rz(-1.1037165) q[3];
sx q[3];
rz(2.3576749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(0.80580795) q[2];
rz(-0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7508115) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(-3.0506296) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(1.3202753) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66619191) q[0];
sx q[0];
rz(-1.1328508) q[0];
sx q[0];
rz(-1.9802718) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57029057) q[2];
sx q[2];
rz(-1.4207134) q[2];
sx q[2];
rz(-2.1652086) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4016061) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(0.92958881) q[1];
rz(-pi) q[2];
rz(2.8517013) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(1.0155201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0423353) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(2.5406204) q[2];
rz(-2.6565334) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(-1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8619974) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(-2.5860508) q[0];
rz(3.1069966) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(-1.7506036) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23404183) q[0];
sx q[0];
rz(-2.4939135) q[0];
sx q[0];
rz(2.5729022) q[0];
x q[1];
rz(0.36632914) q[2];
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
x q[0];
rz(-0.78087378) q[1];
sx q[1];
rz(-2.2370173) q[1];
sx q[1];
rz(-2.9800376) q[1];
rz(-pi) q[2];
rz(0.075918003) q[3];
sx q[3];
rz(-2.3644991) q[3];
sx q[3];
rz(1.0958375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3283219) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(1.288712) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(-2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46273461) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(-0.95343268) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(1.7038201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8673082) q[0];
sx q[0];
rz(-2.0831997) q[0];
sx q[0];
rz(0.209765) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67280025) q[2];
sx q[2];
rz(-2.9252508) q[2];
sx q[2];
rz(2.1397482) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86086035) q[1];
sx q[1];
rz(-2.3128465) q[1];
sx q[1];
rz(0.6764722) q[1];
x q[2];
rz(1.5290456) q[3];
sx q[3];
rz(-2.6937727) q[3];
sx q[3];
rz(-2.8988422) q[3];
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
rz(-0.86137613) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(0.3716968) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9361967) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(-2.0478915) q[0];
rz(-0.73668346) q[1];
sx q[1];
rz(-1.8700347) q[1];
sx q[1];
rz(-1.0587943) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5878384) q[0];
sx q[0];
rz(-0.61843473) q[0];
sx q[0];
rz(-1.0879602) q[0];
x q[1];
rz(-2.043622) q[2];
sx q[2];
rz(-1.2038004) q[2];
sx q[2];
rz(-2.8851913) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4150548) q[1];
sx q[1];
rz(-1.6931567) q[1];
sx q[1];
rz(1.1916257) q[1];
x q[2];
rz(-2.4478854) q[3];
sx q[3];
rz(-0.5872763) q[3];
sx q[3];
rz(-1.8380084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(0.41040928) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(-0.60925305) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(-1.4321009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2925443) q[0];
sx q[0];
rz(-1.2546854) q[0];
sx q[0];
rz(2.6228117) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6851995) q[2];
sx q[2];
rz(-2.6420339) q[2];
sx q[2];
rz(1.511614) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6546302) q[1];
sx q[1];
rz(-1.6143867) q[1];
sx q[1];
rz(0.55200465) q[1];
rz(3.0544153) q[3];
sx q[3];
rz(-2.0012337) q[3];
sx q[3];
rz(-0.038539683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6802406) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(1.1414026) q[2];
rz(-1.5348148) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(-0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5466945) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(-0.36874157) q[1];
sx q[1];
rz(-1.2422961) q[1];
sx q[1];
rz(3.0098343) q[1];
rz(0.73451191) q[2];
sx q[2];
rz(-2.8399158) q[2];
sx q[2];
rz(0.5010571) q[2];
rz(1.7678123) q[3];
sx q[3];
rz(-1.1548629) q[3];
sx q[3];
rz(-2.9668273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

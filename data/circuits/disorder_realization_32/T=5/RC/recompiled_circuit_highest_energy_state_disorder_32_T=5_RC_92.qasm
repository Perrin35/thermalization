OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.028539874) q[0];
sx q[0];
rz(3.9588504) q[0];
sx q[0];
rz(9.1398653) q[0];
rz(-2.3218396) q[1];
sx q[1];
rz(-2.2567891) q[1];
sx q[1];
rz(1.1362145) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4856271) q[0];
sx q[0];
rz(-0.95447094) q[0];
sx q[0];
rz(-3.0462395) q[0];
rz(0.66464632) q[2];
sx q[2];
rz(-0.79865361) q[2];
sx q[2];
rz(-0.89671521) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2023102) q[1];
sx q[1];
rz(-1.9946672) q[1];
sx q[1];
rz(1.2594628) q[1];
rz(-pi) q[2];
rz(-0.44228999) q[3];
sx q[3];
rz(-1.7768404) q[3];
sx q[3];
rz(1.1161436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39584407) q[2];
sx q[2];
rz(-2.6555847) q[2];
sx q[2];
rz(-2.8493472) q[2];
rz(2.5203943) q[3];
sx q[3];
rz(-2.0665593) q[3];
sx q[3];
rz(2.9228389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.578823) q[0];
sx q[0];
rz(-1.1172453) q[0];
sx q[0];
rz(-0.086409464) q[0];
rz(2.7482765) q[1];
sx q[1];
rz(-1.6875279) q[1];
sx q[1];
rz(-1.6578065) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42272767) q[0];
sx q[0];
rz(-1.8047835) q[0];
sx q[0];
rz(-2.2586593) q[0];
rz(1.9772524) q[2];
sx q[2];
rz(-0.33824846) q[2];
sx q[2];
rz(1.2814613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.377903) q[1];
sx q[1];
rz(-2.2202146) q[1];
sx q[1];
rz(-0.19291921) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8350164) q[3];
sx q[3];
rz(-1.0162899) q[3];
sx q[3];
rz(-2.4305024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1958486) q[2];
sx q[2];
rz(-1.5986779) q[2];
sx q[2];
rz(-3.1405385) q[2];
rz(-0.70683181) q[3];
sx q[3];
rz(-2.2471434) q[3];
sx q[3];
rz(-2.8640174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23866776) q[0];
sx q[0];
rz(-0.9280197) q[0];
sx q[0];
rz(2.6631885) q[0];
rz(-2.2967285) q[1];
sx q[1];
rz(-1.076315) q[1];
sx q[1];
rz(-2.1873651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060392875) q[0];
sx q[0];
rz(-2.0234479) q[0];
sx q[0];
rz(-1.1552108) q[0];
rz(-pi) q[1];
rz(0.83691544) q[2];
sx q[2];
rz(-2.0544312) q[2];
sx q[2];
rz(-1.4326295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3500757) q[1];
sx q[1];
rz(-1.9261987) q[1];
sx q[1];
rz(-2.3018964) q[1];
x q[2];
rz(-1.30458) q[3];
sx q[3];
rz(-2.2053218) q[3];
sx q[3];
rz(-2.3887544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2039631) q[2];
sx q[2];
rz(-0.340168) q[2];
sx q[2];
rz(0.071268737) q[2];
rz(-1.010262) q[3];
sx q[3];
rz(-1.9805757) q[3];
sx q[3];
rz(2.8205813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45692745) q[0];
sx q[0];
rz(-1.2926956) q[0];
sx q[0];
rz(0.36351031) q[0];
rz(-2.2495296) q[1];
sx q[1];
rz(-1.0328247) q[1];
sx q[1];
rz(1.2537289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3433468) q[0];
sx q[0];
rz(-1.1805184) q[0];
sx q[0];
rz(-2.0623931) q[0];
x q[1];
rz(0.10102306) q[2];
sx q[2];
rz(-1.2083078) q[2];
sx q[2];
rz(-2.6355138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2334785) q[1];
sx q[1];
rz(-1.5292948) q[1];
sx q[1];
rz(2.8913778) q[1];
rz(-pi) q[2];
rz(0.39616743) q[3];
sx q[3];
rz(-1.2855144) q[3];
sx q[3];
rz(-2.0821636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38972798) q[2];
sx q[2];
rz(-0.77406293) q[2];
sx q[2];
rz(1.5285726) q[2];
rz(-2.6914237) q[3];
sx q[3];
rz(-1.7894141) q[3];
sx q[3];
rz(1.2990797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12440974) q[0];
sx q[0];
rz(-2.9876509) q[0];
sx q[0];
rz(-2.2907139) q[0];
rz(1.4022013) q[1];
sx q[1];
rz(-1.5504928) q[1];
sx q[1];
rz(0.15241399) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35203347) q[0];
sx q[0];
rz(-1.5769582) q[0];
sx q[0];
rz(1.8101109) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1609116) q[2];
sx q[2];
rz(-0.64233795) q[2];
sx q[2];
rz(2.040524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4153395) q[1];
sx q[1];
rz(-1.2881973) q[1];
sx q[1];
rz(-2.9690109) q[1];
x q[2];
rz(-2.3984473) q[3];
sx q[3];
rz(-0.51880793) q[3];
sx q[3];
rz(1.102953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5695345) q[2];
sx q[2];
rz(-0.92430884) q[2];
sx q[2];
rz(-2.7657236) q[2];
rz(2.87319) q[3];
sx q[3];
rz(-1.0207876) q[3];
sx q[3];
rz(1.0387897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1035476) q[0];
sx q[0];
rz(-0.91829848) q[0];
sx q[0];
rz(3.015836) q[0];
rz(0.045299564) q[1];
sx q[1];
rz(-0.89729512) q[1];
sx q[1];
rz(0.50517857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4046458) q[0];
sx q[0];
rz(-2.8789799) q[0];
sx q[0];
rz(-0.24554952) q[0];
x q[1];
rz(1.366682) q[2];
sx q[2];
rz(-2.6417126) q[2];
sx q[2];
rz(1.4171039) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1152079) q[1];
sx q[1];
rz(-1.7167101) q[1];
sx q[1];
rz(-0.77512957) q[1];
rz(-pi) q[2];
rz(2.807928) q[3];
sx q[3];
rz(-1.1265605) q[3];
sx q[3];
rz(0.21008106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81856412) q[2];
sx q[2];
rz(-0.33107859) q[2];
sx q[2];
rz(2.8884812) q[2];
rz(-0.74462914) q[3];
sx q[3];
rz(-1.587845) q[3];
sx q[3];
rz(-2.871992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.700106) q[0];
sx q[0];
rz(-1.2865257) q[0];
sx q[0];
rz(0.043524608) q[0];
rz(1.173191) q[1];
sx q[1];
rz(-1.2930608) q[1];
sx q[1];
rz(-0.80165577) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35079682) q[0];
sx q[0];
rz(-2.9614095) q[0];
sx q[0];
rz(-0.58191802) q[0];
x q[1];
rz(-2.413676) q[2];
sx q[2];
rz(-2.2322828) q[2];
sx q[2];
rz(-0.9582516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1420586) q[1];
sx q[1];
rz(-2.2323066) q[1];
sx q[1];
rz(2.9444047) q[1];
x q[2];
rz(0.67569895) q[3];
sx q[3];
rz(-2.7609918) q[3];
sx q[3];
rz(2.8322647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7511071) q[2];
sx q[2];
rz(-2.0483053) q[2];
sx q[2];
rz(0.92457479) q[2];
rz(-3.0604002) q[3];
sx q[3];
rz(-1.3982747) q[3];
sx q[3];
rz(-0.72807062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28295949) q[0];
sx q[0];
rz(-1.2456243) q[0];
sx q[0];
rz(1.0754732) q[0];
rz(-1.7643499) q[1];
sx q[1];
rz(-1.7694764) q[1];
sx q[1];
rz(-2.4864054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0690883) q[0];
sx q[0];
rz(-0.069511183) q[0];
sx q[0];
rz(-2.4609835) q[0];
x q[1];
rz(-0.78777625) q[2];
sx q[2];
rz(-1.3609741) q[2];
sx q[2];
rz(0.94888806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1671518) q[1];
sx q[1];
rz(-0.87572423) q[1];
sx q[1];
rz(-1.4822757) q[1];
rz(-pi) q[2];
rz(-1.9094134) q[3];
sx q[3];
rz(-1.5985672) q[3];
sx q[3];
rz(-2.9737792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7738771) q[2];
sx q[2];
rz(-2.750062) q[2];
sx q[2];
rz(1.9850622) q[2];
rz(1.1441506) q[3];
sx q[3];
rz(-2.5113228) q[3];
sx q[3];
rz(-1.8088079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82041204) q[0];
sx q[0];
rz(-0.92083609) q[0];
sx q[0];
rz(-1.17571) q[0];
rz(1.8020449) q[1];
sx q[1];
rz(-2.1422155) q[1];
sx q[1];
rz(2.7340926) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20243199) q[0];
sx q[0];
rz(-0.5078041) q[0];
sx q[0];
rz(2.5277471) q[0];
x q[1];
rz(-1.3027329) q[2];
sx q[2];
rz(-1.1020793) q[2];
sx q[2];
rz(0.38662903) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5989227) q[1];
sx q[1];
rz(-0.58784584) q[1];
sx q[1];
rz(1.6076615) q[1];
rz(-pi) q[2];
rz(0.27176933) q[3];
sx q[3];
rz(-1.6394233) q[3];
sx q[3];
rz(-0.29952279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0743865) q[2];
sx q[2];
rz(-1.1585453) q[2];
sx q[2];
rz(2.1161533) q[2];
rz(-2.6440559) q[3];
sx q[3];
rz(-2.189744) q[3];
sx q[3];
rz(-0.40646762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24791726) q[0];
sx q[0];
rz(-1.9837288) q[0];
sx q[0];
rz(-1.3854618) q[0];
rz(-2.0939743) q[1];
sx q[1];
rz(-1.7966929) q[1];
sx q[1];
rz(-0.97283831) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5248643) q[0];
sx q[0];
rz(-1.4848733) q[0];
sx q[0];
rz(1.026154) q[0];
x q[1];
rz(-0.49373105) q[2];
sx q[2];
rz(-1.0876473) q[2];
sx q[2];
rz(-1.8918623) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9094338) q[1];
sx q[1];
rz(-1.9240018) q[1];
sx q[1];
rz(-1.722419) q[1];
rz(-pi) q[2];
rz(-2.0165958) q[3];
sx q[3];
rz(-0.74927038) q[3];
sx q[3];
rz(0.42271915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4086548) q[2];
sx q[2];
rz(-2.813377) q[2];
sx q[2];
rz(2.0396063) q[2];
rz(1.7105626) q[3];
sx q[3];
rz(-2.3165063) q[3];
sx q[3];
rz(1.903532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86193209) q[0];
sx q[0];
rz(-1.4623549) q[0];
sx q[0];
rz(-2.097492) q[0];
rz(-0.73061371) q[1];
sx q[1];
rz(-0.63332557) q[1];
sx q[1];
rz(0.81061737) q[1];
rz(-1.8865449) q[2];
sx q[2];
rz(-0.98636711) q[2];
sx q[2];
rz(-1.2489088) q[2];
rz(1.295119) q[3];
sx q[3];
rz(-1.9119921) q[3];
sx q[3];
rz(0.32031624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

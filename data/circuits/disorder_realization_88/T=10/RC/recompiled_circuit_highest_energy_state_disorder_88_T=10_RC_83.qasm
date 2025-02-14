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
rz(-0.1470546) q[0];
sx q[0];
rz(-1.4015863) q[0];
sx q[0];
rz(0.92010486) q[0];
rz(3.0609581) q[1];
sx q[1];
rz(-0.57890761) q[1];
sx q[1];
rz(2.1644367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.282519) q[0];
sx q[0];
rz(-2.7096097) q[0];
sx q[0];
rz(0.52406128) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0394568) q[2];
sx q[2];
rz(-0.13291026) q[2];
sx q[2];
rz(-1.6249715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5861533) q[1];
sx q[1];
rz(-1.4358913) q[1];
sx q[1];
rz(0.6026661) q[1];
rz(-pi) q[2];
rz(0.36698384) q[3];
sx q[3];
rz(-2.4189848) q[3];
sx q[3];
rz(2.5740576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5620293) q[2];
sx q[2];
rz(-1.2144438) q[2];
sx q[2];
rz(-2.5207632) q[2];
rz(2.6260455) q[3];
sx q[3];
rz(-1.7190869) q[3];
sx q[3];
rz(-0.8425042) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2466549) q[0];
sx q[0];
rz(-1.6064914) q[0];
sx q[0];
rz(-2.5625693) q[0];
rz(2.31965) q[1];
sx q[1];
rz(-1.6252981) q[1];
sx q[1];
rz(-0.45375219) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762208) q[0];
sx q[0];
rz(-0.95673086) q[0];
sx q[0];
rz(-1.2943511) q[0];
rz(1.4582107) q[2];
sx q[2];
rz(-2.0346918) q[2];
sx q[2];
rz(-0.91700208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.380341) q[1];
sx q[1];
rz(-2.5793992) q[1];
sx q[1];
rz(0.23558771) q[1];
rz(-pi) q[2];
rz(-2.6964397) q[3];
sx q[3];
rz(-1.0694519) q[3];
sx q[3];
rz(-1.8568764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4072676) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(0.56488758) q[2];
rz(2.7867553) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.7179276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722209) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(0.55150223) q[0];
rz(-0.36477271) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(0.50484467) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0009815) q[0];
sx q[0];
rz(-1.7303082) q[0];
sx q[0];
rz(0.41054748) q[0];
rz(-pi) q[1];
rz(-2.7288923) q[2];
sx q[2];
rz(-1.6566212) q[2];
sx q[2];
rz(-0.43666652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7271115) q[1];
sx q[1];
rz(-1.0950146) q[1];
sx q[1];
rz(-2.4610795) q[1];
rz(2.4931413) q[3];
sx q[3];
rz(-0.57944991) q[3];
sx q[3];
rz(-3.0045829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(3.0943387) q[2];
rz(-0.83682483) q[3];
sx q[3];
rz(-2.4182726) q[3];
sx q[3];
rz(1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.717201) q[0];
sx q[0];
rz(-1.0294788) q[0];
sx q[0];
rz(1.8151872) q[0];
rz(-1.4855509) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(0.63492376) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6328652) q[0];
sx q[0];
rz(-1.5888402) q[0];
sx q[0];
rz(3.1403072) q[0];
rz(-pi) q[1];
rz(-2.1702386) q[2];
sx q[2];
rz(-0.58341372) q[2];
sx q[2];
rz(2.3177878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30764515) q[1];
sx q[1];
rz(-1.9059423) q[1];
sx q[1];
rz(3.0651211) q[1];
x q[2];
rz(-1.0403955) q[3];
sx q[3];
rz(-1.0113071) q[3];
sx q[3];
rz(-2.6992269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9012458) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(1.3725613) q[2];
rz(1.9339804) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(-0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8127301) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(-0.10511705) q[0];
rz(2.8016727) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-1.0629268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2824616) q[0];
sx q[0];
rz(-0.71159186) q[0];
sx q[0];
rz(3.1200073) q[0];
rz(-pi) q[1];
rz(1.8196443) q[2];
sx q[2];
rz(-0.82468678) q[2];
sx q[2];
rz(-2.7052372) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.059710709) q[1];
sx q[1];
rz(-3.0813327) q[1];
sx q[1];
rz(2.9032272) q[1];
rz(-pi) q[2];
rz(-0.0023554597) q[3];
sx q[3];
rz(-1.071953) q[3];
sx q[3];
rz(2.3209907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0792599) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(-1.8049392) q[2];
rz(0.054232728) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17276758) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(1.2139976) q[0];
rz(-2.0955775) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(-2.2878343) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19148286) q[0];
sx q[0];
rz(-1.5111418) q[0];
sx q[0];
rz(1.4364987) q[0];
x q[1];
rz(-0.60889027) q[2];
sx q[2];
rz(-0.85431803) q[2];
sx q[2];
rz(2.4520055) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5750372) q[1];
sx q[1];
rz(-2.2872529) q[1];
sx q[1];
rz(-0.8030007) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7925909) q[3];
sx q[3];
rz(-1.7141388) q[3];
sx q[3];
rz(0.81763148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0115016) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(0.74756527) q[2];
rz(-0.93303624) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(0.30195495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1099243) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(-0.64144301) q[0];
rz(-2.172442) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(-2.0279121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3308709) q[0];
sx q[0];
rz(-1.4735876) q[0];
sx q[0];
rz(-1.9770798) q[0];
rz(-pi) q[1];
x q[1];
rz(0.895787) q[2];
sx q[2];
rz(-0.75715827) q[2];
sx q[2];
rz(-0.96162187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42150527) q[1];
sx q[1];
rz(-1.0522138) q[1];
sx q[1];
rz(-0.20629136) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6051172) q[3];
sx q[3];
rz(-2.0908818) q[3];
sx q[3];
rz(1.2421364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.574719) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(2.3835772) q[2];
rz(0.30592439) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(0.44736403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4161943) q[0];
sx q[0];
rz(-2.5328126) q[0];
sx q[0];
rz(-2.388227) q[0];
rz(-2.2196409) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(-3.0070378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0908703) q[0];
sx q[0];
rz(-0.36733741) q[0];
sx q[0];
rz(1.3274756) q[0];
rz(-pi) q[1];
rz(1.8250663) q[2];
sx q[2];
rz(-0.94506028) q[2];
sx q[2];
rz(2.3601687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2462897) q[1];
sx q[1];
rz(-0.091769204) q[1];
sx q[1];
rz(0.20905881) q[1];
x q[2];
rz(-2.5882072) q[3];
sx q[3];
rz(-2.3883005) q[3];
sx q[3];
rz(2.2224421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3677463) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(1.194225) q[2];
rz(-1.9959244) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(1.0660508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0399748) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(-2.7591144) q[0];
rz(2.3756012) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.8006178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1857698) q[0];
sx q[0];
rz(-2.9583724) q[0];
sx q[0];
rz(2.630156) q[0];
rz(-pi) q[1];
rz(2.3843896) q[2];
sx q[2];
rz(-1.5067889) q[2];
sx q[2];
rz(-0.369095) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7241269) q[1];
sx q[1];
rz(-1.7977906) q[1];
sx q[1];
rz(-2.3844152) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2148592) q[3];
sx q[3];
rz(-1.4934469) q[3];
sx q[3];
rz(2.429395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0535023) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(-3.0255393) q[2];
rz(-1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.457823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5905404) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(2.9569448) q[0];
rz(-2.2231936) q[1];
sx q[1];
rz(-1.5946439) q[1];
sx q[1];
rz(1.7041697) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0032063403) q[0];
sx q[0];
rz(-1.0995563) q[0];
sx q[0];
rz(-0.22881656) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.080731656) q[2];
sx q[2];
rz(-2.417832) q[2];
sx q[2];
rz(-2.2249976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9574163) q[1];
sx q[1];
rz(-0.97089689) q[1];
sx q[1];
rz(-0.23317144) q[1];
rz(-pi) q[2];
rz(-2.734951) q[3];
sx q[3];
rz(-2.4235117) q[3];
sx q[3];
rz(-2.6585916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.664428) q[2];
sx q[2];
rz(-1.2457341) q[2];
sx q[2];
rz(0.63875443) q[2];
rz(2.3264558) q[3];
sx q[3];
rz(-0.4709979) q[3];
sx q[3];
rz(2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(-1.9736171) q[2];
sx q[2];
rz(-0.30277534) q[2];
sx q[2];
rz(1.4390611) q[2];
rz(0.15624795) q[3];
sx q[3];
rz(-2.8946946) q[3];
sx q[3];
rz(1.3794086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

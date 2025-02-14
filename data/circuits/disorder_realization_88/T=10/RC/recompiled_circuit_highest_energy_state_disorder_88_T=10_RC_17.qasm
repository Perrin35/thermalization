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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8492964) q[0];
sx q[0];
rz(-1.1998645) q[0];
sx q[0];
rz(1.7975259) q[0];
rz(-1.5844272) q[2];
sx q[2];
rz(-1.4385828) q[2];
sx q[2];
rz(1.4135828) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0339573) q[1];
sx q[1];
rz(-0.9743685) q[1];
sx q[1];
rz(-1.4075085) q[1];
rz(-1.8771873) q[3];
sx q[3];
rz(-2.236111) q[3];
sx q[3];
rz(-0.09395919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(2.5207632) q[2];
rz(2.6260455) q[3];
sx q[3];
rz(-1.7190869) q[3];
sx q[3];
rz(2.2990885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(0.57902336) q[0];
rz(0.82194263) q[1];
sx q[1];
rz(-1.6252981) q[1];
sx q[1];
rz(-2.6878405) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4849783) q[0];
sx q[0];
rz(-1.3458283) q[0];
sx q[0];
rz(-2.5091835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4582107) q[2];
sx q[2];
rz(-2.0346918) q[2];
sx q[2];
rz(-2.2245906) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4848043) q[1];
sx q[1];
rz(-2.1156807) q[1];
sx q[1];
rz(-1.4247895) q[1];
rz(-pi) q[2];
rz(2.1165127) q[3];
sx q[3];
rz(-1.9580055) q[3];
sx q[3];
rz(-0.51147616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4072676) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(2.5767051) q[2];
rz(2.7867553) q[3];
sx q[3];
rz(-2.1653039) q[3];
sx q[3];
rz(1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9722209) q[0];
sx q[0];
rz(-0.2114978) q[0];
sx q[0];
rz(-2.5900904) q[0];
rz(-2.7768199) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(-0.50484467) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0802949) q[0];
sx q[0];
rz(-0.43879959) q[0];
sx q[0];
rz(0.3831692) q[0];
rz(-pi) q[1];
rz(-0.41270035) q[2];
sx q[2];
rz(-1.6566212) q[2];
sx q[2];
rz(-2.7049261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7271115) q[1];
sx q[1];
rz(-2.046578) q[1];
sx q[1];
rz(2.4610795) q[1];
rz(-pi) q[2];
rz(-2.4931413) q[3];
sx q[3];
rz(-2.5621427) q[3];
sx q[3];
rz(-3.0045829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7079033) q[2];
sx q[2];
rz(-1.7283551) q[2];
sx q[2];
rz(3.0943387) q[2];
rz(-0.83682483) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(2.1164472) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.717201) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(1.3264054) q[0];
rz(1.6560417) q[1];
sx q[1];
rz(-2.0573719) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43760532) q[0];
sx q[0];
rz(-0.018089596) q[0];
sx q[0];
rz(1.4996858) q[0];
rz(0.97135404) q[2];
sx q[2];
rz(-0.58341372) q[2];
sx q[2];
rz(-0.82380481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.078770854) q[1];
sx q[1];
rz(-2.7981542) q[1];
sx q[1];
rz(-1.7867286) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67948273) q[3];
sx q[3];
rz(-2.3906997) q[3];
sx q[3];
rz(0.39284947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2403468) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(1.7690313) q[2];
rz(1.2076123) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(-2.8670368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8127301) q[0];
sx q[0];
rz(-2.6565318) q[0];
sx q[0];
rz(-3.0364756) q[0];
rz(0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-2.0786659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.88763) q[0];
sx q[0];
rz(-2.2821873) q[0];
sx q[0];
rz(-1.5894058) q[0];
rz(-pi) q[1];
rz(-0.26039298) q[2];
sx q[2];
rz(-0.77885926) q[2];
sx q[2];
rz(0.078140251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0818819) q[1];
sx q[1];
rz(-0.060259911) q[1];
sx q[1];
rz(-0.23836542) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5751198) q[3];
sx q[3];
rz(-0.49884847) q[3];
sx q[3];
rz(-0.81567848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0792599) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(1.8049392) q[2];
rz(0.054232728) q[3];
sx q[3];
rz(-1.1905866) q[3];
sx q[3];
rz(-2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17276758) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(-1.927595) q[0];
rz(-1.0460151) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(2.2878343) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9501098) q[0];
sx q[0];
rz(-1.6304509) q[0];
sx q[0];
rz(-1.7050939) q[0];
x q[1];
rz(-2.1519203) q[2];
sx q[2];
rz(-0.90384353) q[2];
sx q[2];
rz(-0.12573952) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5484838) q[1];
sx q[1];
rz(-0.9973155) q[1];
sx q[1];
rz(-2.4683263) q[1];
rz(-2.7925909) q[3];
sx q[3];
rz(-1.4274538) q[3];
sx q[3];
rz(2.3239612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13009109) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(2.3940274) q[2];
rz(0.93303624) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1099243) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(-2.5001496) q[0];
rz(-0.96915069) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(-1.1136805) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9819558) q[0];
sx q[0];
rz(-2.7244719) q[0];
sx q[0];
rz(-1.3288767) q[0];
rz(0.53345726) q[2];
sx q[2];
rz(-1.0048303) q[2];
sx q[2];
rz(-0.12803687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2526452) q[1];
sx q[1];
rz(-1.749649) q[1];
sx q[1];
rz(-2.0986544) q[1];
rz(-pi) q[2];
rz(1.6051172) q[3];
sx q[3];
rz(-1.0507108) q[3];
sx q[3];
rz(-1.8994562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5668737) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(0.75801545) q[2];
rz(-2.8356683) q[3];
sx q[3];
rz(-2.7829792) q[3];
sx q[3];
rz(2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7253983) q[0];
sx q[0];
rz(-2.5328126) q[0];
sx q[0];
rz(0.7533657) q[0];
rz(-0.92195177) q[1];
sx q[1];
rz(-1.4267068) q[1];
sx q[1];
rz(-3.0070378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20924231) q[0];
sx q[0];
rz(-1.9268231) q[0];
sx q[0];
rz(3.0491475) q[0];
rz(-2.806611) q[2];
sx q[2];
rz(-0.6689531) q[2];
sx q[2];
rz(-1.9425336) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8953029) q[1];
sx q[1];
rz(-3.0498235) q[1];
sx q[1];
rz(-2.9325338) q[1];
x q[2];
rz(0.55338545) q[3];
sx q[3];
rz(-0.7532922) q[3];
sx q[3];
rz(-2.2224421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3677463) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(-1.194225) q[2];
rz(1.9959244) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(-1.0660508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0399748) q[0];
sx q[0];
rz(-2.6823253) q[0];
sx q[0];
rz(-0.38247821) q[0];
rz(0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.3409748) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670973) q[0];
sx q[0];
rz(-1.7303559) q[0];
sx q[0];
rz(-1.480353) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6587615) q[2];
sx q[2];
rz(-0.81552699) q[2];
sx q[2];
rz(-1.2620827) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.080440259) q[1];
sx q[1];
rz(-2.3576479) q[1];
sx q[1];
rz(-0.32439167) q[1];
rz(-pi) q[2];
rz(-1.789647) q[3];
sx q[3];
rz(-2.7776981) q[3];
sx q[3];
rz(2.4879251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0880903) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(-0.11605334) q[2];
rz(1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(-2.9569448) q[0];
rz(-2.2231936) q[1];
sx q[1];
rz(-1.5946439) q[1];
sx q[1];
rz(1.7041697) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4686615) q[0];
sx q[0];
rz(-1.7743006) q[0];
sx q[0];
rz(-1.0887906) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.499648) q[2];
sx q[2];
rz(-0.84991036) q[2];
sx q[2];
rz(2.3325553) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.25333693) q[1];
sx q[1];
rz(-1.3789021) q[1];
sx q[1];
rz(-2.1835421) q[1];
x q[2];
rz(-1.9035133) q[3];
sx q[3];
rz(-0.92192382) q[3];
sx q[3];
rz(2.1391265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47716466) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(0.81513682) q[3];
sx q[3];
rz(-0.4709979) q[3];
sx q[3];
rz(0.42069978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4169793) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(-0.78085113) q[1];
sx q[1];
rz(-1.4518705) q[1];
sx q[1];
rz(-0.89444583) q[1];
rz(-1.1679756) q[2];
sx q[2];
rz(-2.8388173) q[2];
sx q[2];
rz(-1.7025316) q[2];
rz(2.9853447) q[3];
sx q[3];
rz(-0.24689804) q[3];
sx q[3];
rz(-1.762184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

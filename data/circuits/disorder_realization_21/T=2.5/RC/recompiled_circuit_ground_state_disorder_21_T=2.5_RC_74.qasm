OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2709687) q[0];
sx q[0];
rz(-0.55611098) q[0];
sx q[0];
rz(2.1882353) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(-2.1058197) q[1];
sx q[1];
rz(0.9486202) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23466388) q[0];
sx q[0];
rz(-1.0415823) q[0];
sx q[0];
rz(2.0509023) q[0];
rz(-pi) q[1];
rz(-1.5679915) q[2];
sx q[2];
rz(-1.8082779) q[2];
sx q[2];
rz(1.592092) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.030101209) q[1];
sx q[1];
rz(-1.9989357) q[1];
sx q[1];
rz(2.4011022) q[1];
rz(2.7953732) q[3];
sx q[3];
rz(-1.0448714) q[3];
sx q[3];
rz(0.093737515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3775776) q[2];
sx q[2];
rz(-0.074187584) q[2];
sx q[2];
rz(2.3589676) q[2];
rz(-0.11893663) q[3];
sx q[3];
rz(-1.0340034) q[3];
sx q[3];
rz(-2.9940166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.9453732) q[0];
sx q[0];
rz(-1.1213028) q[0];
sx q[0];
rz(-0.25783208) q[0];
rz(0.091015426) q[1];
sx q[1];
rz(-1.0992522) q[1];
sx q[1];
rz(1.6450504) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5071867) q[0];
sx q[0];
rz(-2.1226465) q[0];
sx q[0];
rz(0.080291434) q[0];
x q[1];
rz(1.9159258) q[2];
sx q[2];
rz(-1.218443) q[2];
sx q[2];
rz(2.8690668) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1196668) q[1];
sx q[1];
rz(-0.41321427) q[1];
sx q[1];
rz(1.6456804) q[1];
rz(-0.5037276) q[3];
sx q[3];
rz(-0.73299512) q[3];
sx q[3];
rz(-0.021857787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1220793) q[2];
sx q[2];
rz(-1.8668819) q[2];
sx q[2];
rz(-0.40461928) q[2];
rz(0.98958611) q[3];
sx q[3];
rz(-1.3000969) q[3];
sx q[3];
rz(-2.5185481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1125546) q[0];
sx q[0];
rz(-2.9974388) q[0];
sx q[0];
rz(-2.3032904) q[0];
rz(2.2932032) q[1];
sx q[1];
rz(-1.219607) q[1];
sx q[1];
rz(0.99260509) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780571) q[0];
sx q[0];
rz(-2.40487) q[0];
sx q[0];
rz(-0.20918734) q[0];
rz(-1.6279334) q[2];
sx q[2];
rz(-1.9046202) q[2];
sx q[2];
rz(-1.305507) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.61286825) q[1];
sx q[1];
rz(-2.6515685) q[1];
sx q[1];
rz(-2.571066) q[1];
rz(0.35723585) q[3];
sx q[3];
rz(-0.30870507) q[3];
sx q[3];
rz(-1.5451408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8695716) q[2];
sx q[2];
rz(-1.1676936) q[2];
sx q[2];
rz(-2.0167548) q[2];
rz(2.520842) q[3];
sx q[3];
rz(-2.1706332) q[3];
sx q[3];
rz(3.0264405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.89897412) q[0];
sx q[0];
rz(-1.1599351) q[0];
sx q[0];
rz(0.17380357) q[0];
rz(-2.4558892) q[1];
sx q[1];
rz(-1.4861636) q[1];
sx q[1];
rz(-2.3033843) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4259151) q[0];
sx q[0];
rz(-2.0260149) q[0];
sx q[0];
rz(1.0550189) q[0];
x q[1];
rz(-2.0897267) q[2];
sx q[2];
rz(-1.0609259) q[2];
sx q[2];
rz(3.1367658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2357852) q[1];
sx q[1];
rz(-2.2407994) q[1];
sx q[1];
rz(1.3895821) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6082786) q[3];
sx q[3];
rz(-1.1082543) q[3];
sx q[3];
rz(1.9328062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59914261) q[2];
sx q[2];
rz(-1.6833545) q[2];
sx q[2];
rz(0.35169265) q[2];
rz(1.1951949) q[3];
sx q[3];
rz(-1.9545133) q[3];
sx q[3];
rz(0.024141969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46566063) q[0];
sx q[0];
rz(-0.52023482) q[0];
sx q[0];
rz(-2.7886673) q[0];
rz(2.832761) q[1];
sx q[1];
rz(-2.1052723) q[1];
sx q[1];
rz(-0.028506361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4052947) q[0];
sx q[0];
rz(-2.5319244) q[0];
sx q[0];
rz(-1.2627312) q[0];
rz(-pi) q[1];
rz(-0.38178954) q[2];
sx q[2];
rz(-1.2068527) q[2];
sx q[2];
rz(2.292143) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18896477) q[1];
sx q[1];
rz(-2.3213534) q[1];
sx q[1];
rz(-2.8431975) q[1];
x q[2];
rz(2.0621544) q[3];
sx q[3];
rz(-2.6565928) q[3];
sx q[3];
rz(-0.62937832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9286524) q[2];
sx q[2];
rz(-3.0781367) q[2];
sx q[2];
rz(0.97079903) q[2];
rz(2.094723) q[3];
sx q[3];
rz(-1.4528843) q[3];
sx q[3];
rz(2.651732) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30188072) q[0];
sx q[0];
rz(-1.3588926) q[0];
sx q[0];
rz(-3.0506328) q[0];
rz(1.2184527) q[1];
sx q[1];
rz(-0.33088845) q[1];
sx q[1];
rz(-0.63922304) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7122593) q[0];
sx q[0];
rz(-2.7875226) q[0];
sx q[0];
rz(-2.7430915) q[0];
x q[1];
rz(-3.0193954) q[2];
sx q[2];
rz(-1.7079884) q[2];
sx q[2];
rz(2.4047763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4358232) q[1];
sx q[1];
rz(-0.87638777) q[1];
sx q[1];
rz(2.4468876) q[1];
rz(-pi) q[2];
rz(-1.5766673) q[3];
sx q[3];
rz(-0.3265872) q[3];
sx q[3];
rz(1.1297376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0794534) q[2];
sx q[2];
rz(-0.97525758) q[2];
sx q[2];
rz(0.47387588) q[2];
rz(-1.8114932) q[3];
sx q[3];
rz(-0.88089839) q[3];
sx q[3];
rz(2.7661095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57629958) q[0];
sx q[0];
rz(-1.0813035) q[0];
sx q[0];
rz(2.9190049) q[0];
rz(-1.0635771) q[1];
sx q[1];
rz(-1.56366) q[1];
sx q[1];
rz(1.9482013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759068) q[0];
sx q[0];
rz(-0.75729174) q[0];
sx q[0];
rz(-2.971847) q[0];
rz(-pi) q[1];
rz(1.7980099) q[2];
sx q[2];
rz(-1.6889204) q[2];
sx q[2];
rz(2.3051777) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4970253) q[1];
sx q[1];
rz(-0.46176061) q[1];
sx q[1];
rz(-0.30739947) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3037549) q[3];
sx q[3];
rz(-1.7283261) q[3];
sx q[3];
rz(-0.13539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0003164) q[2];
sx q[2];
rz(-2.2681984) q[2];
sx q[2];
rz(0.6401965) q[2];
rz(-0.021942465) q[3];
sx q[3];
rz(-1.4098189) q[3];
sx q[3];
rz(2.791361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6496277) q[0];
sx q[0];
rz(-1.3977298) q[0];
sx q[0];
rz(0.52038991) q[0];
rz(2.261816) q[1];
sx q[1];
rz(-1.9089411) q[1];
sx q[1];
rz(-0.89281503) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8917903) q[0];
sx q[0];
rz(-0.098488657) q[0];
sx q[0];
rz(-2.6633224) q[0];
rz(-0.73142902) q[2];
sx q[2];
rz(-0.9584223) q[2];
sx q[2];
rz(2.2867672) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5984447) q[1];
sx q[1];
rz(-0.81564995) q[1];
sx q[1];
rz(1.1786908) q[1];
rz(-1.9939984) q[3];
sx q[3];
rz(-1.0931921) q[3];
sx q[3];
rz(1.3582317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2616547) q[2];
sx q[2];
rz(-1.1038021) q[2];
sx q[2];
rz(-0.086816303) q[2];
rz(2.1964729) q[3];
sx q[3];
rz(-0.34543959) q[3];
sx q[3];
rz(2.0115578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9574808) q[0];
sx q[0];
rz(-3.0511973) q[0];
sx q[0];
rz(3.0775253) q[0];
rz(-1.13569) q[1];
sx q[1];
rz(-1.4402729) q[1];
sx q[1];
rz(-2.6293829) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9001246) q[0];
sx q[0];
rz(-1.1607803) q[0];
sx q[0];
rz(-1.2154237) q[0];
rz(1.9718599) q[2];
sx q[2];
rz(-2.3980015) q[2];
sx q[2];
rz(2.0498073) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54359084) q[1];
sx q[1];
rz(-2.8683122) q[1];
sx q[1];
rz(-0.32229801) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5867611) q[3];
sx q[3];
rz(-1.9796238) q[3];
sx q[3];
rz(-2.7328934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3247437) q[2];
sx q[2];
rz(-2.4090359) q[2];
sx q[2];
rz(2.0762439) q[2];
rz(1.6522853) q[3];
sx q[3];
rz(-0.95732147) q[3];
sx q[3];
rz(0.092441946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5548993) q[0];
sx q[0];
rz(-2.7296992) q[0];
sx q[0];
rz(0.33083415) q[0];
rz(-1.2384442) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(2.8668561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3150644) q[0];
sx q[0];
rz(-1.2312268) q[0];
sx q[0];
rz(2.0281606) q[0];
x q[1];
rz(0.11604734) q[2];
sx q[2];
rz(-0.2751285) q[2];
sx q[2];
rz(-0.76062084) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7995791) q[1];
sx q[1];
rz(-1.7106317) q[1];
sx q[1];
rz(-2.0796156) q[1];
rz(-pi) q[2];
rz(-1.4624743) q[3];
sx q[3];
rz(-0.59039718) q[3];
sx q[3];
rz(0.27449671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9773418) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(-1.3573307) q[2];
rz(0.50061289) q[3];
sx q[3];
rz(-1.5784135) q[3];
sx q[3];
rz(1.64465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5464583) q[0];
sx q[0];
rz(-1.5820137) q[0];
sx q[0];
rz(-0.59376846) q[0];
rz(3.0733227) q[1];
sx q[1];
rz(-1.2208114) q[1];
sx q[1];
rz(0.023718871) q[1];
rz(2.7082607) q[2];
sx q[2];
rz(-2.7637252) q[2];
sx q[2];
rz(-2.6890474) q[2];
rz(1.3832292) q[3];
sx q[3];
rz(-1.3959342) q[3];
sx q[3];
rz(0.29831553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

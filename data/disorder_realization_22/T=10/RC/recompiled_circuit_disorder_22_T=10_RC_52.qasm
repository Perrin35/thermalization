OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(4.099457) q[0];
sx q[0];
rz(9.2803331) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(0.97775835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5241961) q[0];
sx q[0];
rz(-2.8688736) q[0];
sx q[0];
rz(-0.32193907) q[0];
rz(0.30774967) q[2];
sx q[2];
rz(-2.0648742) q[2];
sx q[2];
rz(2.5623164) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6482918) q[1];
sx q[1];
rz(-0.19953218) q[1];
sx q[1];
rz(-2.3471911) q[1];
rz(-pi) q[2];
rz(1.3917189) q[3];
sx q[3];
rz(-2.350051) q[3];
sx q[3];
rz(-0.66034987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4686761) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(0.19876924) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4345877) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(0.36112753) q[0];
rz(-1.7065642) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(-0.8180058) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0571787) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(-1.2731228) q[0];
rz(3.0188574) q[2];
sx q[2];
rz(-1.2018179) q[2];
sx q[2];
rz(-2.4462647) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9243014) q[1];
sx q[1];
rz(-1.1326619) q[1];
sx q[1];
rz(-2.2469254) q[1];
x q[2];
rz(0.057095842) q[3];
sx q[3];
rz(-1.9340056) q[3];
sx q[3];
rz(2.925194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8890185) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(1.4206295) q[2];
rz(-1.8255) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(2.7533598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(0.87093583) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(-0.2972163) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57097829) q[0];
sx q[0];
rz(-1.1665205) q[0];
sx q[0];
rz(1.5017609) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88163968) q[2];
sx q[2];
rz(-2.214553) q[2];
sx q[2];
rz(1.8635441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9416562) q[1];
sx q[1];
rz(-1.7574649) q[1];
sx q[1];
rz(-3.093064) q[1];
rz(-pi) q[2];
rz(2.2753784) q[3];
sx q[3];
rz(-0.88170393) q[3];
sx q[3];
rz(-2.6967274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7685984) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-0.42759839) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(-0.52743131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4363842) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(-0.77392459) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(-2.4598222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1958774) q[0];
sx q[0];
rz(-1.3099075) q[0];
sx q[0];
rz(2.9257141) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21708023) q[2];
sx q[2];
rz(-0.97765572) q[2];
sx q[2];
rz(2.8341688) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7766429) q[1];
sx q[1];
rz(-0.72307359) q[1];
sx q[1];
rz(-0.88100453) q[1];
x q[2];
rz(1.9346312) q[3];
sx q[3];
rz(-1.2483276) q[3];
sx q[3];
rz(-2.913167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(-2.0751674) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(-0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(0.56030309) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(1.6220185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9902089) q[0];
sx q[0];
rz(-2.4575893) q[0];
sx q[0];
rz(-2.3955976) q[0];
rz(2.5298169) q[2];
sx q[2];
rz(-2.0423186) q[2];
sx q[2];
rz(-0.39786354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2677742) q[1];
sx q[1];
rz(-2.3582637) q[1];
sx q[1];
rz(-0.28247139) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7468466) q[3];
sx q[3];
rz(-2.5213443) q[3];
sx q[3];
rz(2.8125417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4524298) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(-2.9186644) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(-3.070014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(1.0700595) q[0];
rz(-1.3609715) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(1.2840575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43858466) q[0];
sx q[0];
rz(-2.2692338) q[0];
sx q[0];
rz(0.75011487) q[0];
rz(-pi) q[1];
rz(1.1805004) q[2];
sx q[2];
rz(-1.830606) q[2];
sx q[2];
rz(0.63894546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8821841) q[1];
sx q[1];
rz(-0.55667294) q[1];
sx q[1];
rz(-2.5747719) q[1];
rz(-1.0809903) q[3];
sx q[3];
rz(-2.5488857) q[3];
sx q[3];
rz(2.5268775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47508919) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(-0.63759032) q[2];
rz(0.83667886) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56629431) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(-0.56754011) q[0];
rz(0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(-2.2033851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0465614) q[0];
sx q[0];
rz(-1.9782269) q[0];
sx q[0];
rz(-2.808232) q[0];
rz(-pi) q[1];
rz(1.0023408) q[2];
sx q[2];
rz(-2.217514) q[2];
sx q[2];
rz(-1.9288174) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.36233703) q[1];
sx q[1];
rz(-1.2149097) q[1];
sx q[1];
rz(-0.87902714) q[1];
rz(-1.8754962) q[3];
sx q[3];
rz(-0.88677553) q[3];
sx q[3];
rz(1.5414433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2618835) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(-1.3860469) q[2];
rz(-1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.26738527) q[0];
sx q[0];
rz(-2.8081942) q[0];
sx q[0];
rz(1.7077131) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(-2.3103255) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48152637) q[0];
sx q[0];
rz(-0.43018451) q[0];
sx q[0];
rz(-0.35920401) q[0];
rz(-pi) q[1];
rz(2.2057461) q[2];
sx q[2];
rz(-1.6747464) q[2];
sx q[2];
rz(-0.39633358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2483406) q[1];
sx q[1];
rz(-2.2447531) q[1];
sx q[1];
rz(2.2158951) q[1];
rz(-2.3723699) q[3];
sx q[3];
rz(-2.8260494) q[3];
sx q[3];
rz(-0.10570082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5635809) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-2.1772299) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7678541) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(-1.7550229) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(2.0358553) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10406636) q[0];
sx q[0];
rz(-2.6797047) q[0];
sx q[0];
rz(2.877263) q[0];
rz(-1.4377016) q[2];
sx q[2];
rz(-2.4705952) q[2];
sx q[2];
rz(-2.4457096) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3016262) q[1];
sx q[1];
rz(-1.7558388) q[1];
sx q[1];
rz(1.806083) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35131485) q[3];
sx q[3];
rz(-1.2993386) q[3];
sx q[3];
rz(-2.9815477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5332807) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(-2.9917955) q[2];
rz(1.7685361) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6479284) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(0.1272442) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(0.1677992) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509335) q[0];
sx q[0];
rz(-1.5899842) q[0];
sx q[0];
rz(0.030254342) q[0];
rz(-pi) q[1];
rz(2.5334353) q[2];
sx q[2];
rz(-1.6981352) q[2];
sx q[2];
rz(-1.5978158) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.64393109) q[1];
sx q[1];
rz(-1.5776458) q[1];
sx q[1];
rz(2.5531406) q[1];
rz(-pi) q[2];
rz(1.4526669) q[3];
sx q[3];
rz(-1.4483293) q[3];
sx q[3];
rz(-0.99692217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(0.090027697) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5933843) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(-2.2446752) q[2];
sx q[2];
rz(-0.88340906) q[2];
sx q[2];
rz(-0.39426343) q[2];
rz(-0.17100632) q[3];
sx q[3];
rz(-2.1033559) q[3];
sx q[3];
rz(-0.81567473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

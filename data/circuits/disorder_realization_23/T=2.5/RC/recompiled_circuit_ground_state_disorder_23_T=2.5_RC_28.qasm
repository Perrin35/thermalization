OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.58148122) q[0];
sx q[0];
rz(-1.886263) q[0];
sx q[0];
rz(-0.49129593) q[0];
rz(2.8331941) q[1];
sx q[1];
rz(-1.3003132) q[1];
sx q[1];
rz(-3.0829859) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97293905) q[0];
sx q[0];
rz(-2.0117451) q[0];
sx q[0];
rz(-0.95885528) q[0];
x q[1];
rz(-0.7841447) q[2];
sx q[2];
rz(-1.9040739) q[2];
sx q[2];
rz(3.0212457) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72926988) q[1];
sx q[1];
rz(-0.34360057) q[1];
sx q[1];
rz(-2.6471828) q[1];
rz(-pi) q[2];
x q[2];
rz(0.079777579) q[3];
sx q[3];
rz(-1.8892131) q[3];
sx q[3];
rz(1.4246045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0356174) q[2];
sx q[2];
rz(-0.79088894) q[2];
sx q[2];
rz(-2.7604575) q[2];
rz(-3.1086339) q[3];
sx q[3];
rz(-1.4573174) q[3];
sx q[3];
rz(-1.0491252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9829262) q[0];
sx q[0];
rz(-0.52678147) q[0];
sx q[0];
rz(-0.43989936) q[0];
rz(3.0797709) q[1];
sx q[1];
rz(-1.6114176) q[1];
sx q[1];
rz(2.8672112) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44654057) q[0];
sx q[0];
rz(-2.1529004) q[0];
sx q[0];
rz(2.5869408) q[0];
rz(-pi) q[1];
rz(2.0437061) q[2];
sx q[2];
rz(-2.6691872) q[2];
sx q[2];
rz(-2.1771569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3267131) q[1];
sx q[1];
rz(-0.45456262) q[1];
sx q[1];
rz(1.3978895) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7663244) q[3];
sx q[3];
rz(-0.57884848) q[3];
sx q[3];
rz(1.8414094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2939833) q[2];
sx q[2];
rz(-3.1276438) q[2];
sx q[2];
rz(3.0713165) q[2];
rz(-0.62486068) q[3];
sx q[3];
rz(-1.8127541) q[3];
sx q[3];
rz(0.074987324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845508) q[0];
sx q[0];
rz(-1.5032737) q[0];
sx q[0];
rz(-2.7497838) q[0];
rz(-0.99569544) q[1];
sx q[1];
rz(-0.83637801) q[1];
sx q[1];
rz(0.044274274) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05737722) q[0];
sx q[0];
rz(-0.65859933) q[0];
sx q[0];
rz(1.0530217) q[0];
x q[1];
rz(-2.2703553) q[2];
sx q[2];
rz(-2.2164564) q[2];
sx q[2];
rz(1.7019367) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4489608) q[1];
sx q[1];
rz(-1.7794183) q[1];
sx q[1];
rz(-1.3814998) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6674394) q[3];
sx q[3];
rz(-0.97580384) q[3];
sx q[3];
rz(2.4293824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6280881) q[2];
sx q[2];
rz(-0.91614437) q[2];
sx q[2];
rz(1.8872895) q[2];
rz(-0.11145505) q[3];
sx q[3];
rz(-0.97193757) q[3];
sx q[3];
rz(-2.286262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8348389) q[0];
sx q[0];
rz(-0.67973891) q[0];
sx q[0];
rz(0.87598959) q[0];
rz(0.39729473) q[1];
sx q[1];
rz(-1.5700424) q[1];
sx q[1];
rz(1.3321336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5529163) q[0];
sx q[0];
rz(-2.0400926) q[0];
sx q[0];
rz(1.3568715) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1213673) q[2];
sx q[2];
rz(-0.75006333) q[2];
sx q[2];
rz(-1.6353861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5343473) q[1];
sx q[1];
rz(-1.4922472) q[1];
sx q[1];
rz(-0.15395366) q[1];
rz(-pi) q[2];
rz(-1.9537266) q[3];
sx q[3];
rz(-0.57102244) q[3];
sx q[3];
rz(1.8117222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9925655) q[2];
sx q[2];
rz(-1.488648) q[2];
sx q[2];
rz(1.3943025) q[2];
rz(-0.99327883) q[3];
sx q[3];
rz(-1.9781878) q[3];
sx q[3];
rz(1.2629898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8472327) q[0];
sx q[0];
rz(-2.9672186) q[0];
sx q[0];
rz(-2.2827523) q[0];
rz(2.7186588) q[1];
sx q[1];
rz(-2.6507288) q[1];
sx q[1];
rz(1.4609969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.925526) q[0];
sx q[0];
rz(-0.51617814) q[0];
sx q[0];
rz(0.25052089) q[0];
x q[1];
rz(-1.5654353) q[2];
sx q[2];
rz(-0.44525422) q[2];
sx q[2];
rz(-1.7536947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2475869) q[1];
sx q[1];
rz(-1.3787706) q[1];
sx q[1];
rz(-0.2094261) q[1];
x q[2];
rz(1.8247174) q[3];
sx q[3];
rz(-0.51957209) q[3];
sx q[3];
rz(0.40031734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4451663) q[2];
sx q[2];
rz(-2.7327765) q[2];
sx q[2];
rz(1.4985296) q[2];
rz(1.1284575) q[3];
sx q[3];
rz(-2.0339637) q[3];
sx q[3];
rz(-2.4988153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6978825) q[0];
sx q[0];
rz(-1.8118129) q[0];
sx q[0];
rz(-0.71769303) q[0];
rz(-0.36092654) q[1];
sx q[1];
rz(-0.91054994) q[1];
sx q[1];
rz(-2.8275729) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4249202) q[0];
sx q[0];
rz(-2.0659528) q[0];
sx q[0];
rz(1.6177482) q[0];
rz(0.038617066) q[2];
sx q[2];
rz(-0.84682276) q[2];
sx q[2];
rz(-2.4702009) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2870658) q[1];
sx q[1];
rz(-2.5662133) q[1];
sx q[1];
rz(2.1772312) q[1];
x q[2];
rz(1.7526423) q[3];
sx q[3];
rz(-2.3324128) q[3];
sx q[3];
rz(1.2004747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6641984) q[2];
sx q[2];
rz(-1.3496642) q[2];
sx q[2];
rz(-2.9273709) q[2];
rz(2.2231806) q[3];
sx q[3];
rz(-0.94614202) q[3];
sx q[3];
rz(-1.7414179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.65900954) q[0];
sx q[0];
rz(-0.56203401) q[0];
sx q[0];
rz(2.5157978) q[0];
rz(-1.91045) q[1];
sx q[1];
rz(-1.8742671) q[1];
sx q[1];
rz(-0.81987461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5441262) q[0];
sx q[0];
rz(-2.1334927) q[0];
sx q[0];
rz(-1.94348) q[0];
rz(2.6447222) q[2];
sx q[2];
rz(-0.77503237) q[2];
sx q[2];
rz(-0.97737038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7149799) q[1];
sx q[1];
rz(-2.8016114) q[1];
sx q[1];
rz(1.3424385) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9514337) q[3];
sx q[3];
rz(-3.0453322) q[3];
sx q[3];
rz(-0.13842336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95255533) q[2];
sx q[2];
rz(-1.2574715) q[2];
sx q[2];
rz(-0.38743585) q[2];
rz(2.6881325) q[3];
sx q[3];
rz(-2.267434) q[3];
sx q[3];
rz(2.9816154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.3979935) q[0];
sx q[0];
rz(-1.1389808) q[0];
sx q[0];
rz(2.0177662) q[0];
rz(-2.3881312) q[1];
sx q[1];
rz(-1.1898142) q[1];
sx q[1];
rz(-1.168728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982464) q[0];
sx q[0];
rz(-0.69730824) q[0];
sx q[0];
rz(1.4278317) q[0];
rz(-pi) q[1];
rz(3.1356461) q[2];
sx q[2];
rz(-2.0281254) q[2];
sx q[2];
rz(1.7944687) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67186165) q[1];
sx q[1];
rz(-0.97588343) q[1];
sx q[1];
rz(0.39284006) q[1];
x q[2];
rz(-2.2271676) q[3];
sx q[3];
rz(-1.9468938) q[3];
sx q[3];
rz(2.8825783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1184825) q[2];
sx q[2];
rz(-2.1239943) q[2];
sx q[2];
rz(2.6712096) q[2];
rz(1.6810301) q[3];
sx q[3];
rz(-1.8813671) q[3];
sx q[3];
rz(-2.1605261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39101446) q[0];
sx q[0];
rz(-1.3704726) q[0];
sx q[0];
rz(1.6455261) q[0];
rz(-1.1356614) q[1];
sx q[1];
rz(-1.2874425) q[1];
sx q[1];
rz(0.48887238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58502349) q[0];
sx q[0];
rz(-0.14358601) q[0];
sx q[0];
rz(-1.344259) q[0];
x q[1];
rz(2.082294) q[2];
sx q[2];
rz(-1.4867721) q[2];
sx q[2];
rz(2.3531599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1294293) q[1];
sx q[1];
rz(-1.5485073) q[1];
sx q[1];
rz(3.0921357) q[1];
x q[2];
rz(-2.5219157) q[3];
sx q[3];
rz(-1.6256728) q[3];
sx q[3];
rz(0.60407274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44498542) q[2];
sx q[2];
rz(-2.5312436) q[2];
sx q[2];
rz(-2.0390017) q[2];
rz(-1.6164814) q[3];
sx q[3];
rz(-0.81086719) q[3];
sx q[3];
rz(1.6703687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1333604) q[0];
sx q[0];
rz(-0.52953774) q[0];
sx q[0];
rz(1.8089005) q[0];
rz(-2.7545199) q[1];
sx q[1];
rz(-1.8862855) q[1];
sx q[1];
rz(2.0852087) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1176193) q[0];
sx q[0];
rz(-1.7101544) q[0];
sx q[0];
rz(-1.3246714) q[0];
rz(-pi) q[1];
rz(-2.7235673) q[2];
sx q[2];
rz(-1.4565832) q[2];
sx q[2];
rz(-0.25164652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0478967) q[1];
sx q[1];
rz(-2.2655089) q[1];
sx q[1];
rz(-1.0662778) q[1];
rz(-pi) q[2];
rz(-0.67513533) q[3];
sx q[3];
rz(-2.0656485) q[3];
sx q[3];
rz(2.2371045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.010783823) q[2];
sx q[2];
rz(-0.84785145) q[2];
sx q[2];
rz(1.7424142) q[2];
rz(2.0684659) q[3];
sx q[3];
rz(-1.9173887) q[3];
sx q[3];
rz(0.46260241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.241093) q[0];
sx q[0];
rz(-1.4971965) q[0];
sx q[0];
rz(1.5370488) q[0];
rz(0.72147876) q[1];
sx q[1];
rz(-0.57054467) q[1];
sx q[1];
rz(1.9346938) q[1];
rz(-1.9611364) q[2];
sx q[2];
rz(-1.9705638) q[2];
sx q[2];
rz(0.86314291) q[2];
rz(0.7714959) q[3];
sx q[3];
rz(-1.4998927) q[3];
sx q[3];
rz(-3.0580228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

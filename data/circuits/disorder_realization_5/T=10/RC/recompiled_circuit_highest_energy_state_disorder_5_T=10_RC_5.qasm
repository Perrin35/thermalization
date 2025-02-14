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
rz(1.16601) q[0];
sx q[0];
rz(2.2628885) q[0];
sx q[0];
rz(9.2287697) q[0];
rz(1.6485543) q[1];
sx q[1];
rz(-2.2033446) q[1];
sx q[1];
rz(-0.87747639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24841079) q[0];
sx q[0];
rz(-2.0201004) q[0];
sx q[0];
rz(0.23907678) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67812596) q[2];
sx q[2];
rz(-2.6161751) q[2];
sx q[2];
rz(-3.0209288) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0044549) q[1];
sx q[1];
rz(-0.52912213) q[1];
sx q[1];
rz(2.1251388) q[1];
rz(-pi) q[2];
rz(2.0095413) q[3];
sx q[3];
rz(-2.6225704) q[3];
sx q[3];
rz(-1.2281017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.33285546) q[2];
sx q[2];
rz(-2.4337807) q[2];
sx q[2];
rz(0.037192496) q[2];
rz(-0.91341364) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(-1.6325525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(-0.2221701) q[0];
rz(2.9479058) q[1];
sx q[1];
rz(-1.8733459) q[1];
sx q[1];
rz(-0.2598612) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3701137) q[0];
sx q[0];
rz(-2.6278017) q[0];
sx q[0];
rz(2.0822078) q[0];
rz(-1.031135) q[2];
sx q[2];
rz(-2.9339613) q[2];
sx q[2];
rz(0.59484824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7227962) q[1];
sx q[1];
rz(-1.7823311) q[1];
sx q[1];
rz(1.3822894) q[1];
rz(-pi) q[2];
rz(-1.0845029) q[3];
sx q[3];
rz(-0.87275617) q[3];
sx q[3];
rz(1.8956888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3746609) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(-2.984821) q[2];
rz(-0.80592704) q[3];
sx q[3];
rz(-2.0511878) q[3];
sx q[3];
rz(-2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2369775) q[0];
sx q[0];
rz(-3.0293063) q[0];
sx q[0];
rz(-1.9554546) q[0];
rz(3.0778432) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(-0.41248163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97226364) q[0];
sx q[0];
rz(-2.0893851) q[0];
sx q[0];
rz(2.5725468) q[0];
rz(2.9846014) q[2];
sx q[2];
rz(-0.89688939) q[2];
sx q[2];
rz(1.6403331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3220249) q[1];
sx q[1];
rz(-1.7734999) q[1];
sx q[1];
rz(2.8611819) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2354911) q[3];
sx q[3];
rz(-2.0482488) q[3];
sx q[3];
rz(-0.64732823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5731262) q[2];
sx q[2];
rz(-0.44034475) q[2];
sx q[2];
rz(-1.1661412) q[2];
rz(2.3447573) q[3];
sx q[3];
rz(-1.7723869) q[3];
sx q[3];
rz(-0.78127512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.9154938) q[0];
sx q[0];
rz(-1.508536) q[0];
sx q[0];
rz(1.728212) q[0];
rz(-2.6426897) q[1];
sx q[1];
rz(-1.3830802) q[1];
sx q[1];
rz(-1.2938719) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32484937) q[0];
sx q[0];
rz(-1.6559282) q[0];
sx q[0];
rz(3.1146634) q[0];
rz(3.0656112) q[2];
sx q[2];
rz(-1.2442028) q[2];
sx q[2];
rz(-0.031215515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.368694) q[1];
sx q[1];
rz(-0.062090896) q[1];
sx q[1];
rz(-3.0903878) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0741229) q[3];
sx q[3];
rz(-1.6327792) q[3];
sx q[3];
rz(0.27674473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6125907) q[2];
sx q[2];
rz(-0.68098536) q[2];
sx q[2];
rz(-2.9717818) q[2];
rz(0.42810193) q[3];
sx q[3];
rz(-1.4295108) q[3];
sx q[3];
rz(2.9812109) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47996461) q[0];
sx q[0];
rz(-0.39230883) q[0];
sx q[0];
rz(-2.3062134) q[0];
rz(2.5021878) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(-1.0909874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24966403) q[0];
sx q[0];
rz(-0.87452379) q[0];
sx q[0];
rz(1.554717) q[0];
rz(-2.8418737) q[2];
sx q[2];
rz(-1.4106361) q[2];
sx q[2];
rz(1.0417787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9613232) q[1];
sx q[1];
rz(-1.1774447) q[1];
sx q[1];
rz(-2.250227) q[1];
rz(-pi) q[2];
rz(-2.8150878) q[3];
sx q[3];
rz(-1.0708263) q[3];
sx q[3];
rz(0.96159354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7731446) q[2];
sx q[2];
rz(-1.1870563) q[2];
sx q[2];
rz(3.1375569) q[2];
rz(-1.076738) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(-0.18690404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3183597) q[0];
sx q[0];
rz(-1.7802745) q[0];
sx q[0];
rz(0.63364345) q[0];
rz(-0.40114316) q[1];
sx q[1];
rz(-0.70924962) q[1];
sx q[1];
rz(-2.4493682) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.322583) q[0];
sx q[0];
rz(-0.87912175) q[0];
sx q[0];
rz(-0.54499596) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31816407) q[2];
sx q[2];
rz(-0.79472322) q[2];
sx q[2];
rz(2.367521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3347098) q[1];
sx q[1];
rz(-0.63155424) q[1];
sx q[1];
rz(1.1401661) q[1];
rz(0.18544894) q[3];
sx q[3];
rz(-1.7191965) q[3];
sx q[3];
rz(-0.68391358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9358518) q[2];
sx q[2];
rz(-2.3571099) q[2];
sx q[2];
rz(1.240823) q[2];
rz(1.1848909) q[3];
sx q[3];
rz(-1.3013867) q[3];
sx q[3];
rz(1.544781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.009509) q[0];
sx q[0];
rz(-1.3030095) q[0];
sx q[0];
rz(-1.3707772) q[0];
rz(2.7235183) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(2.6549227) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91557568) q[0];
sx q[0];
rz(-2.5737983) q[0];
sx q[0];
rz(1.3894807) q[0];
rz(-2.1213221) q[2];
sx q[2];
rz(-1.2472868) q[2];
sx q[2];
rz(-0.24187096) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8312953) q[1];
sx q[1];
rz(-2.8716209) q[1];
sx q[1];
rz(-0.84862535) q[1];
rz(-1.6177931) q[3];
sx q[3];
rz(-1.1011657) q[3];
sx q[3];
rz(-2.3835973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.439094) q[2];
sx q[2];
rz(-2.0926026) q[2];
sx q[2];
rz(0.19719633) q[2];
rz(-2.6321453) q[3];
sx q[3];
rz(-0.26064894) q[3];
sx q[3];
rz(-2.313224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8615123) q[0];
sx q[0];
rz(-1.0741638) q[0];
sx q[0];
rz(1.7751088) q[0];
rz(-0.81661433) q[1];
sx q[1];
rz(-2.0520703) q[1];
sx q[1];
rz(2.8588967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0007083) q[0];
sx q[0];
rz(-1.3907897) q[0];
sx q[0];
rz(-0.96830826) q[0];
rz(-pi) q[1];
rz(2.8300072) q[2];
sx q[2];
rz(-1.4805613) q[2];
sx q[2];
rz(0.923522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1360838) q[1];
sx q[1];
rz(-0.89673954) q[1];
sx q[1];
rz(-0.23917178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6201309) q[3];
sx q[3];
rz(-1.2288501) q[3];
sx q[3];
rz(-1.4777044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.035159811) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(2.1584568) q[2];
rz(1.8631009) q[3];
sx q[3];
rz(-0.92517868) q[3];
sx q[3];
rz(0.92888752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66016692) q[0];
sx q[0];
rz(-1.7919414) q[0];
sx q[0];
rz(2.8283258) q[0];
rz(-2.616864) q[1];
sx q[1];
rz(-0.26291651) q[1];
sx q[1];
rz(-1.252334) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9577873) q[0];
sx q[0];
rz(-2.5982214) q[0];
sx q[0];
rz(2.9068391) q[0];
rz(-pi) q[1];
rz(0.60151327) q[2];
sx q[2];
rz(-2.2819715) q[2];
sx q[2];
rz(1.2486697) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58432594) q[1];
sx q[1];
rz(-2.1030104) q[1];
sx q[1];
rz(2.5126712) q[1];
rz(-0.081324025) q[3];
sx q[3];
rz(-0.81800753) q[3];
sx q[3];
rz(-2.6214056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9423882) q[2];
sx q[2];
rz(-2.51666) q[2];
sx q[2];
rz(1.3814629) q[2];
rz(2.840461) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(-2.7605831) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23406601) q[0];
sx q[0];
rz(-1.0028239) q[0];
sx q[0];
rz(-1.3475077) q[0];
rz(-0.8849591) q[1];
sx q[1];
rz(-0.62565175) q[1];
sx q[1];
rz(1.7431097) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8690344) q[0];
sx q[0];
rz(-1.5814651) q[0];
sx q[0];
rz(1.19633) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86549564) q[2];
sx q[2];
rz(-1.4303606) q[2];
sx q[2];
rz(0.27117929) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8191479) q[1];
sx q[1];
rz(-1.5335173) q[1];
sx q[1];
rz(-3.0611866) q[1];
x q[2];
rz(1.9246401) q[3];
sx q[3];
rz(-0.5914591) q[3];
sx q[3];
rz(2.9194057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.053190319) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(-2.1957446) q[2];
rz(0.69019067) q[3];
sx q[3];
rz(-1.918957) q[3];
sx q[3];
rz(-1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90203862) q[0];
sx q[0];
rz(-1.5382465) q[0];
sx q[0];
rz(2.4334346) q[0];
rz(-3.0802849) q[1];
sx q[1];
rz(-0.61534449) q[1];
sx q[1];
rz(-1.4475488) q[1];
rz(-2.8715677) q[2];
sx q[2];
rz(-0.30120987) q[2];
sx q[2];
rz(0.78200151) q[2];
rz(-1.3907018) q[3];
sx q[3];
rz(-1.4994499) q[3];
sx q[3];
rz(-2.3105619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

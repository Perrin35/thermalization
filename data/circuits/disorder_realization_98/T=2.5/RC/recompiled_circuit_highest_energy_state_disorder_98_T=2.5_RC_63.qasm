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
rz(-0.12914817) q[0];
sx q[0];
rz(-1.669786) q[0];
sx q[0];
rz(0.9653402) q[0];
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(1.7641915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5063874) q[0];
sx q[0];
rz(-1.1381519) q[0];
sx q[0];
rz(-0.80308993) q[0];
rz(-0.032261511) q[2];
sx q[2];
rz(-2.2164377) q[2];
sx q[2];
rz(1.0021068) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4272312) q[1];
sx q[1];
rz(-2.9032482) q[1];
sx q[1];
rz(2.4803216) q[1];
rz(-0.018947424) q[3];
sx q[3];
rz(-0.80120443) q[3];
sx q[3];
rz(0.5294746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2237902) q[2];
sx q[2];
rz(-1.4897572) q[2];
sx q[2];
rz(-1.6415143) q[2];
rz(2.3598059) q[3];
sx q[3];
rz(-2.2512071) q[3];
sx q[3];
rz(0.75209832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.016945275) q[0];
sx q[0];
rz(-0.77777672) q[0];
sx q[0];
rz(0.97050226) q[0];
rz(1.3264725) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(2.2296947) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2925948) q[0];
sx q[0];
rz(-0.98750293) q[0];
sx q[0];
rz(-1.3720157) q[0];
x q[1];
rz(0.85477306) q[2];
sx q[2];
rz(-2.8840384) q[2];
sx q[2];
rz(1.4285029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.32102958) q[1];
sx q[1];
rz(-2.2707175) q[1];
sx q[1];
rz(-0.053348347) q[1];
rz(-1.6025869) q[3];
sx q[3];
rz(-0.45511757) q[3];
sx q[3];
rz(-1.8546752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9087002) q[2];
sx q[2];
rz(-1.1129271) q[2];
sx q[2];
rz(-0.98928893) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-1.0003072) q[3];
sx q[3];
rz(-2.4205128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.99783889) q[0];
sx q[0];
rz(-1.9793352) q[0];
sx q[0];
rz(-1.0379399) q[0];
rz(-2.2833917) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(0.16955489) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35992453) q[0];
sx q[0];
rz(-1.5667652) q[0];
sx q[0];
rz(1.4272593) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2409147) q[2];
sx q[2];
rz(-0.49374106) q[2];
sx q[2];
rz(-2.7639219) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.819014) q[1];
sx q[1];
rz(-0.59016229) q[1];
sx q[1];
rz(2.1165089) q[1];
rz(1.9603332) q[3];
sx q[3];
rz(-2.7523605) q[3];
sx q[3];
rz(1.3078944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1797336) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(2.7355984) q[2];
rz(0.77754846) q[3];
sx q[3];
rz(-2.8113139) q[3];
sx q[3];
rz(0.8146666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6430214) q[0];
sx q[0];
rz(-1.2914456) q[0];
sx q[0];
rz(0.92798573) q[0];
rz(-0.058874933) q[1];
sx q[1];
rz(-1.6935655) q[1];
sx q[1];
rz(-1.7108797) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74364036) q[0];
sx q[0];
rz(-0.40983202) q[0];
sx q[0];
rz(1.3751956) q[0];
rz(-pi) q[1];
rz(2.7785382) q[2];
sx q[2];
rz(-1.7573234) q[2];
sx q[2];
rz(-2.8721953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2637564) q[1];
sx q[1];
rz(-2.0621513) q[1];
sx q[1];
rz(-2.8884535) q[1];
rz(1.9742613) q[3];
sx q[3];
rz(-2.472252) q[3];
sx q[3];
rz(2.739213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6835988) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(1.9963473) q[2];
rz(-0.84732071) q[3];
sx q[3];
rz(-1.3342131) q[3];
sx q[3];
rz(1.5020717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1640846) q[0];
sx q[0];
rz(-1.1825528) q[0];
sx q[0];
rz(0.384828) q[0];
rz(0.79967868) q[1];
sx q[1];
rz(-0.5189907) q[1];
sx q[1];
rz(1.5481366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.528397) q[0];
sx q[0];
rz(-1.3269182) q[0];
sx q[0];
rz(-2.8923678) q[0];
rz(-pi) q[1];
rz(-2.9351531) q[2];
sx q[2];
rz(-1.9597133) q[2];
sx q[2];
rz(-0.85064519) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0210798) q[1];
sx q[1];
rz(-1.9200293) q[1];
sx q[1];
rz(-1.2362739) q[1];
rz(0.37171419) q[3];
sx q[3];
rz(-1.0519042) q[3];
sx q[3];
rz(-2.1154161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7657713) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(3.0613464) q[2];
rz(-2.5806228) q[3];
sx q[3];
rz(-1.6601945) q[3];
sx q[3];
rz(-2.5188353) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(1.2605865) q[0];
rz(0.27443019) q[1];
sx q[1];
rz(-1.471289) q[1];
sx q[1];
rz(0.71896499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2043122) q[0];
sx q[0];
rz(-2.215909) q[0];
sx q[0];
rz(0.44207032) q[0];
rz(-1.5814797) q[2];
sx q[2];
rz(-0.74655338) q[2];
sx q[2];
rz(3.106046) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.58108854) q[1];
sx q[1];
rz(-1.6916654) q[1];
sx q[1];
rz(-0.64847364) q[1];
x q[2];
rz(1.1996693) q[3];
sx q[3];
rz(-1.7748482) q[3];
sx q[3];
rz(0.95641092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6545973) q[2];
sx q[2];
rz(-1.8751112) q[2];
sx q[2];
rz(-0.60301644) q[2];
rz(1.6675789) q[3];
sx q[3];
rz(-1.0242198) q[3];
sx q[3];
rz(2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5176158) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(-2.9543167) q[0];
rz(-0.97224832) q[1];
sx q[1];
rz(-0.15521237) q[1];
sx q[1];
rz(-2.7211199) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0213326) q[0];
sx q[0];
rz(-1.1545404) q[0];
sx q[0];
rz(0.72159213) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37897972) q[2];
sx q[2];
rz(-1.0595269) q[2];
sx q[2];
rz(0.015470964) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5125888) q[1];
sx q[1];
rz(-1.2930451) q[1];
sx q[1];
rz(1.1729097) q[1];
rz(-pi) q[2];
rz(-0.99655788) q[3];
sx q[3];
rz(-1.7645932) q[3];
sx q[3];
rz(-2.7791948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(-0.25071684) q[2];
rz(1.8935253) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(-0.59984508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2774169) q[0];
sx q[0];
rz(-0.59159788) q[0];
sx q[0];
rz(0.94171062) q[0];
rz(0.038453728) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(-3.0252735) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.250688) q[0];
sx q[0];
rz(-0.90849344) q[0];
sx q[0];
rz(-0.57470365) q[0];
x q[1];
rz(2.2519926) q[2];
sx q[2];
rz(-1.7187013) q[2];
sx q[2];
rz(-3.0513632) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6623936) q[1];
sx q[1];
rz(-0.31649193) q[1];
sx q[1];
rz(2.5280158) q[1];
x q[2];
rz(-2.8781901) q[3];
sx q[3];
rz(-2.7399094) q[3];
sx q[3];
rz(2.4954737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2945127) q[2];
sx q[2];
rz(-0.81081644) q[2];
sx q[2];
rz(1.026356) q[2];
rz(-1.2096842) q[3];
sx q[3];
rz(-2.1116202) q[3];
sx q[3];
rz(-2.1799555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.472979) q[0];
sx q[0];
rz(-2.0675779) q[0];
sx q[0];
rz(-0.37856722) q[0];
rz(1.1096795) q[1];
sx q[1];
rz(-0.30866426) q[1];
sx q[1];
rz(-3.1139156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6768764) q[0];
sx q[0];
rz(-0.72273556) q[0];
sx q[0];
rz(0.56761543) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9566628) q[2];
sx q[2];
rz(-1.7835296) q[2];
sx q[2];
rz(2.8367868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3301446) q[1];
sx q[1];
rz(-1.3864349) q[1];
sx q[1];
rz(-1.9346647) q[1];
x q[2];
rz(2.1975193) q[3];
sx q[3];
rz(-2.6289796) q[3];
sx q[3];
rz(-0.98820283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3045584) q[2];
sx q[2];
rz(-0.95418945) q[2];
sx q[2];
rz(0.41998106) q[2];
rz(0.84152451) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(-0.37874547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3751752) q[0];
sx q[0];
rz(-3.0247122) q[0];
sx q[0];
rz(0.28512678) q[0];
rz(0.20052234) q[1];
sx q[1];
rz(-1.9384117) q[1];
sx q[1];
rz(2.3060422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052976089) q[0];
sx q[0];
rz(-1.1167913) q[0];
sx q[0];
rz(0.80843057) q[0];
rz(-1.5773238) q[2];
sx q[2];
rz(-2.2178429) q[2];
sx q[2];
rz(1.7056291) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0855779) q[1];
sx q[1];
rz(-1.9025355) q[1];
sx q[1];
rz(0.47420926) q[1];
rz(-0.16160139) q[3];
sx q[3];
rz(-2.9029663) q[3];
sx q[3];
rz(-2.0493226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42651132) q[2];
sx q[2];
rz(-2.0501523) q[2];
sx q[2];
rz(-0.69768989) q[2];
rz(-2.6155124) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(-1.6349767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1220916) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(-1.1493692) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(-1.1463852) q[2];
sx q[2];
rz(-2.41832) q[2];
sx q[2];
rz(0.24857749) q[2];
rz(-2.5935843) q[3];
sx q[3];
rz(-1.8262134) q[3];
sx q[3];
rz(-1.2977737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.030169686) q[0];
sx q[0];
rz(4.5017894) q[0];
sx q[0];
rz(5.9202249) q[0];
rz(1.9650004) q[1];
sx q[1];
rz(-1.3624374) q[1];
sx q[1];
rz(-1.6972313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5267619) q[0];
sx q[0];
rz(-0.20069622) q[0];
sx q[0];
rz(2.654128) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28791223) q[2];
sx q[2];
rz(-1.4511564) q[2];
sx q[2];
rz(-1.5900915) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3807371) q[1];
sx q[1];
rz(-2.0049913) q[1];
sx q[1];
rz(-3.0369706) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6300171) q[3];
sx q[3];
rz(-2.539336) q[3];
sx q[3];
rz(-1.8730375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.968367) q[2];
sx q[2];
rz(-2.3543365) q[2];
sx q[2];
rz(-3.0464029) q[2];
rz(-1.3683052) q[3];
sx q[3];
rz(-0.94822001) q[3];
sx q[3];
rz(-0.29909721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66265166) q[0];
sx q[0];
rz(-2.3421685) q[0];
sx q[0];
rz(0.69315243) q[0];
rz(2.8043546) q[1];
sx q[1];
rz(-1.8353029) q[1];
sx q[1];
rz(2.3487924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1900729) q[0];
sx q[0];
rz(-1.8197104) q[0];
sx q[0];
rz(-2.3547257) q[0];
x q[1];
rz(-2.9117965) q[2];
sx q[2];
rz(-1.0655128) q[2];
sx q[2];
rz(1.8577223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0568367) q[1];
sx q[1];
rz(-2.3178219) q[1];
sx q[1];
rz(0.39577248) q[1];
x q[2];
rz(0.81962193) q[3];
sx q[3];
rz(-1.2030798) q[3];
sx q[3];
rz(1.468827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1482131) q[2];
sx q[2];
rz(-0.17017636) q[2];
sx q[2];
rz(1.6896089) q[2];
rz(2.4827237) q[3];
sx q[3];
rz(-1.6784607) q[3];
sx q[3];
rz(-2.7448867) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1228929) q[0];
sx q[0];
rz(-2.1227699) q[0];
sx q[0];
rz(-0.10957154) q[0];
rz(0.84596363) q[1];
sx q[1];
rz(-0.94596091) q[1];
sx q[1];
rz(2.3931961) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7331241) q[0];
sx q[0];
rz(-1.4811902) q[0];
sx q[0];
rz(2.1277079) q[0];
rz(-pi) q[1];
rz(0.72785925) q[2];
sx q[2];
rz(-2.4185491) q[2];
sx q[2];
rz(-0.41280189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8019164) q[1];
sx q[1];
rz(-1.5640364) q[1];
sx q[1];
rz(-1.6338145) q[1];
rz(-pi) q[2];
rz(1.2975537) q[3];
sx q[3];
rz(-0.58876172) q[3];
sx q[3];
rz(-1.2514284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0261859) q[2];
sx q[2];
rz(-1.4226961) q[2];
sx q[2];
rz(0.7126979) q[2];
rz(-1.6648434) q[3];
sx q[3];
rz(-1.8524757) q[3];
sx q[3];
rz(0.084499806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0067714) q[0];
sx q[0];
rz(-2.0016142) q[0];
sx q[0];
rz(-0.086932927) q[0];
rz(2.8839819) q[1];
sx q[1];
rz(-1.0898217) q[1];
sx q[1];
rz(-1.2923406) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9168388) q[0];
sx q[0];
rz(-1.3319022) q[0];
sx q[0];
rz(-0.80015667) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7316225) q[2];
sx q[2];
rz(-1.8984744) q[2];
sx q[2];
rz(-0.92152061) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72491652) q[1];
sx q[1];
rz(-1.9179763) q[1];
sx q[1];
rz(2.2031839) q[1];
rz(-2.7601808) q[3];
sx q[3];
rz(-2.7578045) q[3];
sx q[3];
rz(0.82510766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8452235) q[2];
sx q[2];
rz(-1.951428) q[2];
sx q[2];
rz(2.2947218) q[2];
rz(-1.3891247) q[3];
sx q[3];
rz(-1.4877157) q[3];
sx q[3];
rz(2.5510767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1883989) q[0];
sx q[0];
rz(-2.0162835) q[0];
sx q[0];
rz(0.97287792) q[0];
rz(-0.94213048) q[1];
sx q[1];
rz(-2.3140488) q[1];
sx q[1];
rz(-3.1408302) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1551092) q[0];
sx q[0];
rz(-2.0707978) q[0];
sx q[0];
rz(-0.9802823) q[0];
rz(2.7050651) q[2];
sx q[2];
rz(-0.68219705) q[2];
sx q[2];
rz(1.4774708) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9474802) q[1];
sx q[1];
rz(-1.0092217) q[1];
sx q[1];
rz(-0.33163867) q[1];
x q[2];
rz(-0.59511472) q[3];
sx q[3];
rz(-1.6356118) q[3];
sx q[3];
rz(-2.7832372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7611277) q[2];
sx q[2];
rz(-1.4030115) q[2];
sx q[2];
rz(0.38499704) q[2];
rz(-2.4132686) q[3];
sx q[3];
rz(-0.60898048) q[3];
sx q[3];
rz(2.716632) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4094792) q[0];
sx q[0];
rz(-2.3347169) q[0];
sx q[0];
rz(-2.683486) q[0];
rz(-1.7181646) q[1];
sx q[1];
rz(-0.79045311) q[1];
sx q[1];
rz(-0.093553392) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5377184) q[0];
sx q[0];
rz(-0.97386003) q[0];
sx q[0];
rz(-0.56054795) q[0];
rz(-pi) q[1];
rz(0.92624591) q[2];
sx q[2];
rz(-2.8684596) q[2];
sx q[2];
rz(2.5018843) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1653891) q[1];
sx q[1];
rz(-1.7862583) q[1];
sx q[1];
rz(2.5991751) q[1];
rz(0.5597975) q[3];
sx q[3];
rz(-1.5526532) q[3];
sx q[3];
rz(1.1813576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2576922) q[2];
sx q[2];
rz(-0.42282405) q[2];
sx q[2];
rz(-0.60919961) q[2];
rz(2.5404663) q[3];
sx q[3];
rz(-1.6060035) q[3];
sx q[3];
rz(0.48457423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6236098) q[0];
sx q[0];
rz(-1.4731151) q[0];
sx q[0];
rz(-1.3129039) q[0];
rz(3.0598705) q[1];
sx q[1];
rz(-1.3879958) q[1];
sx q[1];
rz(1.0645701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9814822) q[0];
sx q[0];
rz(-0.89799268) q[0];
sx q[0];
rz(0.84168169) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5470869) q[2];
sx q[2];
rz(-0.61653149) q[2];
sx q[2];
rz(1.7414684) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.895911) q[1];
sx q[1];
rz(-2.5616922) q[1];
sx q[1];
rz(1.7008476) q[1];
rz(-pi) q[2];
rz(-1.0113297) q[3];
sx q[3];
rz(-2.4593393) q[3];
sx q[3];
rz(-0.38328612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99819034) q[2];
sx q[2];
rz(-2.7394131) q[2];
sx q[2];
rz(1.9833938) q[2];
rz(-2.4573333) q[3];
sx q[3];
rz(-1.8161769) q[3];
sx q[3];
rz(-1.6392596) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96712464) q[0];
sx q[0];
rz(-3.0376349) q[0];
sx q[0];
rz(0.52484584) q[0];
rz(-1.9606494) q[1];
sx q[1];
rz(-0.8371822) q[1];
sx q[1];
rz(-2.4139074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1215661) q[0];
sx q[0];
rz(-1.354573) q[0];
sx q[0];
rz(-0.89631594) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0725756) q[2];
sx q[2];
rz(-0.91337088) q[2];
sx q[2];
rz(1.1595961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39036738) q[1];
sx q[1];
rz(-1.8341244) q[1];
sx q[1];
rz(-0.45231426) q[1];
x q[2];
rz(-0.10730174) q[3];
sx q[3];
rz(-2.2224217) q[3];
sx q[3];
rz(2.1880423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22964792) q[2];
sx q[2];
rz(-1.0934528) q[2];
sx q[2];
rz(1.3463119) q[2];
rz(0.66156578) q[3];
sx q[3];
rz(-1.9965568) q[3];
sx q[3];
rz(-3.055011) q[3];
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
rz(-0.10053703) q[0];
sx q[0];
rz(-2.4424398) q[0];
sx q[0];
rz(0.93851411) q[0];
rz(1.8427294) q[1];
sx q[1];
rz(-2.2524736) q[1];
sx q[1];
rz(-0.13359698) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4105281) q[0];
sx q[0];
rz(-0.34240568) q[0];
sx q[0];
rz(-2.8591741) q[0];
rz(1.8739545) q[2];
sx q[2];
rz(-2.1253307) q[2];
sx q[2];
rz(1.6078439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6386968) q[1];
sx q[1];
rz(-1.9303444) q[1];
sx q[1];
rz(-1.4877968) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4997098) q[3];
sx q[3];
rz(-1.2716023) q[3];
sx q[3];
rz(3.0282216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1453229) q[2];
sx q[2];
rz(-1.0066373) q[2];
sx q[2];
rz(-2.2553196) q[2];
rz(-0.26850548) q[3];
sx q[3];
rz(-2.0566745) q[3];
sx q[3];
rz(1.6296384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71235424) q[0];
sx q[0];
rz(-0.3952643) q[0];
sx q[0];
rz(-0.23671737) q[0];
rz(-2.9780544) q[1];
sx q[1];
rz(-1.6103368) q[1];
sx q[1];
rz(2.9530361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8343413) q[0];
sx q[0];
rz(-2.3545958) q[0];
sx q[0];
rz(2.045162) q[0];
rz(-1.8256906) q[2];
sx q[2];
rz(-2.0370649) q[2];
sx q[2];
rz(2.3026274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.982354) q[1];
sx q[1];
rz(-0.51430632) q[1];
sx q[1];
rz(1.4532754) q[1];
rz(-pi) q[2];
rz(1.7018407) q[3];
sx q[3];
rz(-2.0072972) q[3];
sx q[3];
rz(-0.97697645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4208372) q[2];
sx q[2];
rz(-1.5034224) q[2];
sx q[2];
rz(-2.3756964) q[2];
rz(0.01865538) q[3];
sx q[3];
rz(-1.3466287) q[3];
sx q[3];
rz(0.53001058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50273773) q[0];
sx q[0];
rz(-1.1549594) q[0];
sx q[0];
rz(1.4422944) q[0];
rz(-2.8522708) q[1];
sx q[1];
rz(-2.0850291) q[1];
sx q[1];
rz(-2.6883968) q[1];
rz(2.917541) q[2];
sx q[2];
rz(-0.55006156) q[2];
sx q[2];
rz(2.8203996) q[2];
rz(-1.7073333) q[3];
sx q[3];
rz(-2.668181) q[3];
sx q[3];
rz(-2.6287813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

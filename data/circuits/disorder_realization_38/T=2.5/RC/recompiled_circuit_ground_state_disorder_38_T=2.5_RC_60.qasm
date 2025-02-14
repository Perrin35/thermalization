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
rz(-1.7813959) q[0];
sx q[0];
rz(2.7786322) q[0];
rz(1.9650004) q[1];
sx q[1];
rz(-1.3624374) q[1];
sx q[1];
rz(1.4443614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5267619) q[0];
sx q[0];
rz(-2.9408964) q[0];
sx q[0];
rz(0.48746462) q[0];
rz(-pi) q[1];
rz(1.4460725) q[2];
sx q[2];
rz(-1.2850003) q[2];
sx q[2];
rz(-0.054626183) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23420424) q[1];
sx q[1];
rz(-1.475913) q[1];
sx q[1];
rz(1.1345052) q[1];
rz(2.1722342) q[3];
sx q[3];
rz(-1.5372608) q[3];
sx q[3];
rz(2.7905317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1732257) q[2];
sx q[2];
rz(-2.3543365) q[2];
sx q[2];
rz(3.0464029) q[2];
rz(1.7732874) q[3];
sx q[3];
rz(-0.94822001) q[3];
sx q[3];
rz(2.8424954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66265166) q[0];
sx q[0];
rz(-0.79942411) q[0];
sx q[0];
rz(0.69315243) q[0];
rz(2.8043546) q[1];
sx q[1];
rz(-1.8353029) q[1];
sx q[1];
rz(-0.79280028) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13986282) q[0];
sx q[0];
rz(-2.324449) q[0];
sx q[0];
rz(0.344622) q[0];
rz(-pi) q[1];
rz(1.9614001) q[2];
sx q[2];
rz(-2.5906569) q[2];
sx q[2];
rz(1.40755) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.08475595) q[1];
sx q[1];
rz(-0.82377071) q[1];
sx q[1];
rz(0.39577248) q[1];
rz(-pi) q[2];
rz(-2.0846689) q[3];
sx q[3];
rz(-0.82026359) q[3];
sx q[3];
rz(-0.46945527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1482131) q[2];
sx q[2];
rz(-2.9714163) q[2];
sx q[2];
rz(-1.4519838) q[2];
rz(-0.65886894) q[3];
sx q[3];
rz(-1.4631319) q[3];
sx q[3];
rz(-0.39670593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1228929) q[0];
sx q[0];
rz(-2.1227699) q[0];
sx q[0];
rz(3.0320211) q[0];
rz(-0.84596363) q[1];
sx q[1];
rz(-2.1956317) q[1];
sx q[1];
rz(-0.74839655) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8365252) q[0];
sx q[0];
rz(-2.5782707) q[0];
sx q[0];
rz(1.402424) q[0];
x q[1];
rz(2.559021) q[2];
sx q[2];
rz(-1.114985) q[2];
sx q[2];
rz(-1.3945902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9108991) q[1];
sx q[1];
rz(-1.5077796) q[1];
sx q[1];
rz(3.1348193) q[1];
rz(-pi) q[2];
rz(-0.17828618) q[3];
sx q[3];
rz(-1.0066089) q[3];
sx q[3];
rz(0.92638864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1154068) q[2];
sx q[2];
rz(-1.4226961) q[2];
sx q[2];
rz(2.4288948) q[2];
rz(1.6648434) q[3];
sx q[3];
rz(-1.289117) q[3];
sx q[3];
rz(0.084499806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0067714) q[0];
sx q[0];
rz(-2.0016142) q[0];
sx q[0];
rz(-3.0546597) q[0];
rz(2.8839819) q[1];
sx q[1];
rz(-1.0898217) q[1];
sx q[1];
rz(-1.2923406) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9168388) q[0];
sx q[0];
rz(-1.8096905) q[0];
sx q[0];
rz(-0.80015667) q[0];
x q[1];
rz(2.7013515) q[2];
sx q[2];
rz(-2.7778565) q[2];
sx q[2];
rz(-0.4546356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0513689) q[1];
sx q[1];
rz(-2.1601292) q[1];
sx q[1];
rz(2.7199171) q[1];
x q[2];
rz(-2.7830026) q[3];
sx q[3];
rz(-1.7106283) q[3];
sx q[3];
rz(-0.38966013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8452235) q[2];
sx q[2];
rz(-1.1901647) q[2];
sx q[2];
rz(2.2947218) q[2];
rz(1.752468) q[3];
sx q[3];
rz(-1.4877157) q[3];
sx q[3];
rz(2.5510767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95319372) q[0];
sx q[0];
rz(-2.0162835) q[0];
sx q[0];
rz(2.1687147) q[0];
rz(2.1994622) q[1];
sx q[1];
rz(-0.82754389) q[1];
sx q[1];
rz(-0.00076248893) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98648345) q[0];
sx q[0];
rz(-1.0707948) q[0];
sx q[0];
rz(0.9802823) q[0];
rz(-2.5070287) q[2];
sx q[2];
rz(-1.3009614) q[2];
sx q[2];
rz(2.8874804) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9474802) q[1];
sx q[1];
rz(-2.132371) q[1];
sx q[1];
rz(0.33163867) q[1];
x q[2];
rz(-1.4925748) q[3];
sx q[3];
rz(-0.97710246) q[3];
sx q[3];
rz(-1.8853312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7611277) q[2];
sx q[2];
rz(-1.4030115) q[2];
sx q[2];
rz(0.38499704) q[2];
rz(-0.72832406) q[3];
sx q[3];
rz(-0.60898048) q[3];
sx q[3];
rz(0.42496067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4094792) q[0];
sx q[0];
rz(-0.80687579) q[0];
sx q[0];
rz(2.683486) q[0];
rz(1.4234281) q[1];
sx q[1];
rz(-2.3511395) q[1];
sx q[1];
rz(-3.0480393) q[1];
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
rz(-1.3504998) q[2];
sx q[2];
rz(-1.7335906) q[2];
sx q[2];
rz(-0.3045813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8878758) q[1];
sx q[1];
rz(-2.5619526) q[1];
sx q[1];
rz(2.7405906) q[1];
rz(-pi) q[2];
rz(-1.549386) q[3];
sx q[3];
rz(-1.011102) q[3];
sx q[3];
rz(-2.7407854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2576922) q[2];
sx q[2];
rz(-2.7187686) q[2];
sx q[2];
rz(0.60919961) q[2];
rz(-2.5404663) q[3];
sx q[3];
rz(-1.6060035) q[3];
sx q[3];
rz(2.6570184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5179829) q[0];
sx q[0];
rz(-1.4731151) q[0];
sx q[0];
rz(1.8286888) q[0];
rz(-3.0598705) q[1];
sx q[1];
rz(-1.3879958) q[1];
sx q[1];
rz(-1.0645701) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9186606) q[0];
sx q[0];
rz(-2.1188508) q[0];
sx q[0];
rz(-2.3231044) q[0];
x q[1];
rz(-1.5470869) q[2];
sx q[2];
rz(-2.5250612) q[2];
sx q[2];
rz(1.4001242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.895911) q[1];
sx q[1];
rz(-0.57990042) q[1];
sx q[1];
rz(1.4407451) q[1];
rz(2.1737868) q[3];
sx q[3];
rz(-1.2295614) q[3];
sx q[3];
rz(-1.6399217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99819034) q[2];
sx q[2];
rz(-0.4021796) q[2];
sx q[2];
rz(1.1581988) q[2];
rz(-0.68425933) q[3];
sx q[3];
rz(-1.3254157) q[3];
sx q[3];
rz(1.5023331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96712464) q[0];
sx q[0];
rz(-0.10395771) q[0];
sx q[0];
rz(-2.6167468) q[0];
rz(-1.1809433) q[1];
sx q[1];
rz(-0.8371822) q[1];
sx q[1];
rz(2.4139074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.380881) q[0];
sx q[0];
rz(-2.2267939) q[0];
sx q[0];
rz(2.8674346) q[0];
rz(2.0725756) q[2];
sx q[2];
rz(-0.91337088) q[2];
sx q[2];
rz(-1.9819966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7512253) q[1];
sx q[1];
rz(-1.3074682) q[1];
sx q[1];
rz(-2.6892784) q[1];
x q[2];
rz(-1.7102881) q[3];
sx q[3];
rz(-2.4824641) q[3];
sx q[3];
rz(1.1293251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22964792) q[2];
sx q[2];
rz(-1.0934528) q[2];
sx q[2];
rz(1.7952807) q[2];
rz(-0.66156578) q[3];
sx q[3];
rz(-1.1450359) q[3];
sx q[3];
rz(0.0865817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10053703) q[0];
sx q[0];
rz(-2.4424398) q[0];
sx q[0];
rz(-0.93851411) q[0];
rz(-1.2988633) q[1];
sx q[1];
rz(-0.88911903) q[1];
sx q[1];
rz(0.13359698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4105281) q[0];
sx q[0];
rz(-2.799187) q[0];
sx q[0];
rz(2.8591741) q[0];
x q[1];
rz(-2.6924437) q[2];
sx q[2];
rz(-2.5173016) q[2];
sx q[2];
rz(-2.1438832) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6386968) q[1];
sx q[1];
rz(-1.9303444) q[1];
sx q[1];
rz(1.4877968) q[1];
rz(-pi) q[2];
rz(-1.9383921) q[3];
sx q[3];
rz(-0.9617241) q[3];
sx q[3];
rz(-1.6742791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1453229) q[2];
sx q[2];
rz(-1.0066373) q[2];
sx q[2];
rz(-0.88627306) q[2];
rz(0.26850548) q[3];
sx q[3];
rz(-2.0566745) q[3];
sx q[3];
rz(-1.6296384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71235424) q[0];
sx q[0];
rz(-0.3952643) q[0];
sx q[0];
rz(0.23671737) q[0];
rz(0.16353823) q[1];
sx q[1];
rz(-1.5312559) q[1];
sx q[1];
rz(0.18855655) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.057373) q[0];
sx q[0];
rz(-1.9002267) q[0];
sx q[0];
rz(-0.84217846) q[0];
rz(1.3159021) q[2];
sx q[2];
rz(-2.0370649) q[2];
sx q[2];
rz(2.3026274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30912599) q[1];
sx q[1];
rz(-1.6285076) q[1];
sx q[1];
rz(1.0594491) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27310009) q[3];
sx q[3];
rz(-2.6870603) q[3];
sx q[3];
rz(2.4668193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7207555) q[2];
sx q[2];
rz(-1.5034224) q[2];
sx q[2];
rz(-0.76589626) q[2];
rz(0.01865538) q[3];
sx q[3];
rz(-1.3466287) q[3];
sx q[3];
rz(-2.6115821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50273773) q[0];
sx q[0];
rz(-1.9866332) q[0];
sx q[0];
rz(-1.6992983) q[0];
rz(2.8522708) q[1];
sx q[1];
rz(-1.0565636) q[1];
sx q[1];
rz(0.4531959) q[1];
rz(-2.917541) q[2];
sx q[2];
rz(-2.5915311) q[2];
sx q[2];
rz(-0.32119309) q[2];
rz(-1.1011685) q[3];
sx q[3];
rz(-1.508699) q[3];
sx q[3];
rz(-1.1796622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

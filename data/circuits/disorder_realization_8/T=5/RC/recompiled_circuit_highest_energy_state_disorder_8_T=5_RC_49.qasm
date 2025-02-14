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
rz(0.60628211) q[0];
sx q[0];
rz(4.192306) q[0];
sx q[0];
rz(11.769073) q[0];
rz(0.94766131) q[1];
sx q[1];
rz(-1.892765) q[1];
sx q[1];
rz(-2.7864784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7240802) q[0];
sx q[0];
rz(-1.8204429) q[0];
sx q[0];
rz(1.5019519) q[0];
rz(3.1217977) q[2];
sx q[2];
rz(-1.0934085) q[2];
sx q[2];
rz(2.2098178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4754063) q[1];
sx q[1];
rz(-0.99129035) q[1];
sx q[1];
rz(1.866706) q[1];
rz(-pi) q[2];
rz(-0.45052455) q[3];
sx q[3];
rz(-1.1199711) q[3];
sx q[3];
rz(2.6032347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.015410272) q[2];
sx q[2];
rz(-2.5842032) q[2];
sx q[2];
rz(1.1861447) q[2];
rz(-1.0960389) q[3];
sx q[3];
rz(-2.3604184) q[3];
sx q[3];
rz(-1.621915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.010051589) q[0];
sx q[0];
rz(-3.0942656) q[0];
sx q[0];
rz(1.9897687) q[0];
rz(-0.25543073) q[1];
sx q[1];
rz(-1.7841745) q[1];
sx q[1];
rz(2.7327322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9329703) q[0];
sx q[0];
rz(-1.044036) q[0];
sx q[0];
rz(2.7072858) q[0];
rz(-pi) q[1];
rz(2.2636072) q[2];
sx q[2];
rz(-1.4804941) q[2];
sx q[2];
rz(-2.754359) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0989895) q[1];
sx q[1];
rz(-0.87638646) q[1];
sx q[1];
rz(2.8907772) q[1];
x q[2];
rz(-2.3395038) q[3];
sx q[3];
rz(-0.89381605) q[3];
sx q[3];
rz(-1.9888376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46201181) q[2];
sx q[2];
rz(-2.2363594) q[2];
sx q[2];
rz(-0.26642624) q[2];
rz(-0.5599432) q[3];
sx q[3];
rz(-2.4661049) q[3];
sx q[3];
rz(-2.717836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68488455) q[0];
sx q[0];
rz(-0.86181462) q[0];
sx q[0];
rz(-1.7488712) q[0];
rz(-2.4640153) q[1];
sx q[1];
rz(-1.1303439) q[1];
sx q[1];
rz(2.5052501) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22998097) q[0];
sx q[0];
rz(-1.9634569) q[0];
sx q[0];
rz(2.5225181) q[0];
x q[1];
rz(-0.88865033) q[2];
sx q[2];
rz(-1.3633308) q[2];
sx q[2];
rz(-0.51817521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96752466) q[1];
sx q[1];
rz(-1.8764157) q[1];
sx q[1];
rz(-2.4449062) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40319188) q[3];
sx q[3];
rz(-0.82357126) q[3];
sx q[3];
rz(0.43063569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5431518) q[2];
sx q[2];
rz(-1.8676912) q[2];
sx q[2];
rz(-2.4178823) q[2];
rz(1.9321457) q[3];
sx q[3];
rz(-2.3085322) q[3];
sx q[3];
rz(2.3256653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.74056149) q[0];
sx q[0];
rz(-0.36442345) q[0];
sx q[0];
rz(-2.5357699) q[0];
rz(1.9762899) q[1];
sx q[1];
rz(-1.106768) q[1];
sx q[1];
rz(-0.63444) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5215817) q[0];
sx q[0];
rz(-1.6130035) q[0];
sx q[0];
rz(-2.6642534) q[0];
x q[1];
rz(1.4342821) q[2];
sx q[2];
rz(-2.1079328) q[2];
sx q[2];
rz(0.31202635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8603348) q[1];
sx q[1];
rz(-1.0737281) q[1];
sx q[1];
rz(-0.71038891) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9021413) q[3];
sx q[3];
rz(-0.2718018) q[3];
sx q[3];
rz(-0.1258752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67924649) q[2];
sx q[2];
rz(-2.1470368) q[2];
sx q[2];
rz(-1.7328978) q[2];
rz(0.6399706) q[3];
sx q[3];
rz(-2.7100115) q[3];
sx q[3];
rz(-1.4225167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6022279) q[0];
sx q[0];
rz(-0.27502763) q[0];
sx q[0];
rz(-1.6590903) q[0];
rz(-1.1461343) q[1];
sx q[1];
rz(-0.71123743) q[1];
sx q[1];
rz(1.3292638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70721165) q[0];
sx q[0];
rz(-2.2298621) q[0];
sx q[0];
rz(-2.4102609) q[0];
rz(-0.48945697) q[2];
sx q[2];
rz(-2.202919) q[2];
sx q[2];
rz(-2.7165627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.043221323) q[1];
sx q[1];
rz(-0.30128208) q[1];
sx q[1];
rz(-0.92219086) q[1];
rz(-0.063956694) q[3];
sx q[3];
rz(-1.8956479) q[3];
sx q[3];
rz(2.9717546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9596404) q[2];
sx q[2];
rz(-2.1722309) q[2];
sx q[2];
rz(-2.1592965) q[2];
rz(2.9082409) q[3];
sx q[3];
rz(-1.8459277) q[3];
sx q[3];
rz(2.6006234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0503466) q[0];
sx q[0];
rz(-0.27882689) q[0];
sx q[0];
rz(-2.4009551) q[0];
rz(-2.8431173) q[1];
sx q[1];
rz(-1.9918171) q[1];
sx q[1];
rz(-2.1508335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6388004) q[0];
sx q[0];
rz(-0.37257344) q[0];
sx q[0];
rz(0.034046455) q[0];
rz(-pi) q[1];
rz(-1.9022983) q[2];
sx q[2];
rz(-2.0306892) q[2];
sx q[2];
rz(-2.1761328) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1220529) q[1];
sx q[1];
rz(-0.96457997) q[1];
sx q[1];
rz(-2.6531924) q[1];
rz(0.43277432) q[3];
sx q[3];
rz(-0.7996586) q[3];
sx q[3];
rz(1.8679269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6159181) q[2];
sx q[2];
rz(-2.2423223) q[2];
sx q[2];
rz(-0.57596469) q[2];
rz(1.9383664) q[3];
sx q[3];
rz(-1.7267745) q[3];
sx q[3];
rz(-1.5244213) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7791157) q[0];
sx q[0];
rz(-0.285382) q[0];
sx q[0];
rz(-1.9051636) q[0];
rz(1.3800887) q[1];
sx q[1];
rz(-0.80685902) q[1];
sx q[1];
rz(2.7469452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7056859) q[0];
sx q[0];
rz(-1.5360695) q[0];
sx q[0];
rz(1.3484363) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2038266) q[2];
sx q[2];
rz(-2.870976) q[2];
sx q[2];
rz(-1.1039162) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4952364) q[1];
sx q[1];
rz(-0.64351057) q[1];
sx q[1];
rz(1.8724426) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5731845) q[3];
sx q[3];
rz(-0.34332192) q[3];
sx q[3];
rz(-1.1487701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75957647) q[2];
sx q[2];
rz(-2.3754109) q[2];
sx q[2];
rz(-1.7570599) q[2];
rz(-1.1197439) q[3];
sx q[3];
rz(-1.9789968) q[3];
sx q[3];
rz(-0.82974452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2914332) q[0];
sx q[0];
rz(-2.9740574) q[0];
sx q[0];
rz(2.7450388) q[0];
rz(1.5605759) q[1];
sx q[1];
rz(-2.719559) q[1];
sx q[1];
rz(2.5975554) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2185036) q[0];
sx q[0];
rz(-0.83280665) q[0];
sx q[0];
rz(1.6138747) q[0];
x q[1];
rz(-1.1625353) q[2];
sx q[2];
rz(-0.28851337) q[2];
sx q[2];
rz(0.66513326) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1431191) q[1];
sx q[1];
rz(-2.0859189) q[1];
sx q[1];
rz(-1.4004986) q[1];
x q[2];
rz(2.0062449) q[3];
sx q[3];
rz(-1.4822591) q[3];
sx q[3];
rz(-3.1380668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.80514002) q[2];
sx q[2];
rz(-0.84638941) q[2];
sx q[2];
rz(-1.7479755) q[2];
rz(-2.9584068) q[3];
sx q[3];
rz(-1.8905247) q[3];
sx q[3];
rz(-2.2209871) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2038323) q[0];
sx q[0];
rz(-0.60619175) q[0];
sx q[0];
rz(0.071618557) q[0];
rz(2.3764745) q[1];
sx q[1];
rz(-1.7444997) q[1];
sx q[1];
rz(-0.58806932) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70635722) q[0];
sx q[0];
rz(-1.0093413) q[0];
sx q[0];
rz(-2.7538408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6438573) q[2];
sx q[2];
rz(-0.68314766) q[2];
sx q[2];
rz(2.0367931) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13473171) q[1];
sx q[1];
rz(-2.2344338) q[1];
sx q[1];
rz(-1.4506542) q[1];
x q[2];
rz(1.16251) q[3];
sx q[3];
rz(-0.65925778) q[3];
sx q[3];
rz(1.4548649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4582943) q[2];
sx q[2];
rz(-2.8647162) q[2];
sx q[2];
rz(-1.5923306) q[2];
rz(-1.8706627) q[3];
sx q[3];
rz(-1.1829605) q[3];
sx q[3];
rz(2.8196238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1147605) q[0];
sx q[0];
rz(-2.8883567) q[0];
sx q[0];
rz(-2.7993171) q[0];
rz(0.5510785) q[1];
sx q[1];
rz(-1.2197878) q[1];
sx q[1];
rz(-0.50318998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7321271) q[0];
sx q[0];
rz(-1.9163462) q[0];
sx q[0];
rz(-2.5045128) q[0];
rz(-pi) q[1];
rz(-1.5120498) q[2];
sx q[2];
rz(-2.4565832) q[2];
sx q[2];
rz(1.0420805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.79784517) q[1];
sx q[1];
rz(-0.78727951) q[1];
sx q[1];
rz(2.3549453) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14519174) q[3];
sx q[3];
rz(-0.54779966) q[3];
sx q[3];
rz(2.5992268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2912579) q[2];
sx q[2];
rz(-1.0185654) q[2];
sx q[2];
rz(-0.84602082) q[2];
rz(-0.033898354) q[3];
sx q[3];
rz(-0.97550052) q[3];
sx q[3];
rz(0.94023824) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7593183) q[0];
sx q[0];
rz(-1.2713534) q[0];
sx q[0];
rz(-0.8937339) q[0];
rz(-1.8779124) q[1];
sx q[1];
rz(-2.3448941) q[1];
sx q[1];
rz(5/(14*pi)) q[1];
rz(-0.47847099) q[2];
sx q[2];
rz(-0.84979117) q[2];
sx q[2];
rz(-0.69862943) q[2];
rz(-1.229677) q[3];
sx q[3];
rz(-2.1003953) q[3];
sx q[3];
rz(0.25672124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

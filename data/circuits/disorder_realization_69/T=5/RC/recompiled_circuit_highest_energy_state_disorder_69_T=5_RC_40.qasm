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
rz(0.34613553) q[0];
sx q[0];
rz(-1.5032285) q[0];
sx q[0];
rz(0.7315973) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(-0.43391689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2303378) q[0];
sx q[0];
rz(-0.83090913) q[0];
sx q[0];
rz(-3.0940542) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97831867) q[2];
sx q[2];
rz(-1.0772395) q[2];
sx q[2];
rz(1.6101642) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.190769) q[1];
sx q[1];
rz(-2.1956823) q[1];
sx q[1];
rz(2.6527219) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1103575) q[3];
sx q[3];
rz(-1.3797576) q[3];
sx q[3];
rz(-1.4174549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4291541) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(-2.1377371) q[2];
rz(2.6058274) q[3];
sx q[3];
rz(-2.2771213) q[3];
sx q[3];
rz(-3.1404449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67343229) q[0];
sx q[0];
rz(-0.68244451) q[0];
sx q[0];
rz(0.72714192) q[0];
rz(3.0739821) q[1];
sx q[1];
rz(-1.3550242) q[1];
sx q[1];
rz(2.0179857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1557187) q[0];
sx q[0];
rz(-0.56502461) q[0];
sx q[0];
rz(2.5955517) q[0];
rz(-1.0826473) q[2];
sx q[2];
rz(-0.13859525) q[2];
sx q[2];
rz(0.32735936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20642463) q[1];
sx q[1];
rz(-1.4713396) q[1];
sx q[1];
rz(1.8518492) q[1];
rz(2.0455849) q[3];
sx q[3];
rz(-1.1752306) q[3];
sx q[3];
rz(0.46985782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8565389) q[2];
sx q[2];
rz(-1.5736138) q[2];
sx q[2];
rz(0.48173586) q[2];
rz(0.5528062) q[3];
sx q[3];
rz(-1.0779287) q[3];
sx q[3];
rz(3.0520181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8168617) q[0];
sx q[0];
rz(-2.6664) q[0];
sx q[0];
rz(0.09356308) q[0];
rz(-1.3803253) q[1];
sx q[1];
rz(-0.9535791) q[1];
sx q[1];
rz(0.63724744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475315) q[0];
sx q[0];
rz(-1.1406745) q[0];
sx q[0];
rz(-0.57292666) q[0];
rz(-0.31934505) q[2];
sx q[2];
rz(-2.508906) q[2];
sx q[2];
rz(-1.3038532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90594572) q[1];
sx q[1];
rz(-1.5448346) q[1];
sx q[1];
rz(1.9683377) q[1];
rz(-pi) q[2];
rz(-1.6138541) q[3];
sx q[3];
rz(-1.8612487) q[3];
sx q[3];
rz(1.2424948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20643413) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(-0.68378249) q[2];
rz(2.911496) q[3];
sx q[3];
rz(-1.7347387) q[3];
sx q[3];
rz(-0.42207119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0919331) q[0];
sx q[0];
rz(-2.6347418) q[0];
sx q[0];
rz(-1.7012713) q[0];
rz(-2.6662042) q[1];
sx q[1];
rz(-0.63440591) q[1];
sx q[1];
rz(-0.58116523) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0345349) q[0];
sx q[0];
rz(-2.2904582) q[0];
sx q[0];
rz(-0.11373539) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0175242) q[2];
sx q[2];
rz(-0.5047732) q[2];
sx q[2];
rz(3.1044132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6789279) q[1];
sx q[1];
rz(-1.2691255) q[1];
sx q[1];
rz(-2.660391) q[1];
x q[2];
rz(-1.7293594) q[3];
sx q[3];
rz(-2.057029) q[3];
sx q[3];
rz(-0.18058807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0873969) q[2];
sx q[2];
rz(-2.9554458) q[2];
sx q[2];
rz(-0.52097121) q[2];
rz(-1.4969131) q[3];
sx q[3];
rz(-1.4119586) q[3];
sx q[3];
rz(1.514667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5533376) q[0];
sx q[0];
rz(-0.27565685) q[0];
sx q[0];
rz(1.1302554) q[0];
rz(-2.0626119) q[1];
sx q[1];
rz(-1.1993473) q[1];
sx q[1];
rz(2.030453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2030226) q[0];
sx q[0];
rz(-0.49020728) q[0];
sx q[0];
rz(1.2128279) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52972858) q[2];
sx q[2];
rz(-1.5324943) q[2];
sx q[2];
rz(1.2438347) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20566733) q[1];
sx q[1];
rz(-1.7280518) q[1];
sx q[1];
rz(-0.92604154) q[1];
rz(-2.2258198) q[3];
sx q[3];
rz(-2.0213599) q[3];
sx q[3];
rz(0.612606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7378716) q[2];
sx q[2];
rz(-2.626494) q[2];
sx q[2];
rz(2.1007288) q[2];
rz(-2.3216085) q[3];
sx q[3];
rz(-1.9085725) q[3];
sx q[3];
rz(-2.5177054) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89989221) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(-0.65648055) q[0];
rz(-1.3735324) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(-1.1318644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5053018) q[0];
sx q[0];
rz(-1.3735008) q[0];
sx q[0];
rz(2.4407766) q[0];
rz(-pi) q[1];
rz(0.64986367) q[2];
sx q[2];
rz(-0.78956008) q[2];
sx q[2];
rz(-2.2672578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.060185) q[1];
sx q[1];
rz(-2.2814732) q[1];
sx q[1];
rz(-0.036545444) q[1];
x q[2];
rz(1.5664212) q[3];
sx q[3];
rz(-1.4203912) q[3];
sx q[3];
rz(-0.37585092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.446283) q[2];
sx q[2];
rz(-2.535847) q[2];
sx q[2];
rz(-2.821935) q[2];
rz(0.22932912) q[3];
sx q[3];
rz(-2.1181483) q[3];
sx q[3];
rz(-2.2314609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0730154) q[0];
sx q[0];
rz(-1.9972766) q[0];
sx q[0];
rz(-1.2314433) q[0];
rz(3.0629509) q[1];
sx q[1];
rz(-1.4631203) q[1];
sx q[1];
rz(0.15422779) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.94119) q[0];
sx q[0];
rz(-1.909167) q[0];
sx q[0];
rz(1.75941) q[0];
x q[1];
rz(1.9811822) q[2];
sx q[2];
rz(-1.850381) q[2];
sx q[2];
rz(1.2321763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7423305) q[1];
sx q[1];
rz(-0.62472099) q[1];
sx q[1];
rz(2.4769251) q[1];
x q[2];
rz(1.7835791) q[3];
sx q[3];
rz(-1.9733505) q[3];
sx q[3];
rz(-2.8496131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8057033) q[2];
sx q[2];
rz(-1.615639) q[2];
sx q[2];
rz(1.4914782) q[2];
rz(-1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(-1.570805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0989646) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(-1.6866823) q[0];
rz(0.97575724) q[1];
sx q[1];
rz(-2.0819596) q[1];
sx q[1];
rz(1.3386493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0473235) q[0];
sx q[0];
rz(-1.1828831) q[0];
sx q[0];
rz(2.3045425) q[0];
rz(-pi) q[1];
rz(-0.33390121) q[2];
sx q[2];
rz(-2.0508162) q[2];
sx q[2];
rz(2.1154108) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6852831) q[1];
sx q[1];
rz(-1.3585618) q[1];
sx q[1];
rz(-2.2077399) q[1];
x q[2];
rz(1.9429132) q[3];
sx q[3];
rz(-0.94148472) q[3];
sx q[3];
rz(-3.1222285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51155382) q[2];
sx q[2];
rz(-2.126882) q[2];
sx q[2];
rz(0.83941984) q[2];
rz(1.2342341) q[3];
sx q[3];
rz(-1.0890361) q[3];
sx q[3];
rz(0.031410005) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3473174) q[0];
sx q[0];
rz(-1.3823771) q[0];
sx q[0];
rz(-0.41912249) q[0];
rz(1.6436815) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(-2.7083414) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.142665) q[0];
sx q[0];
rz(-1.6892528) q[0];
sx q[0];
rz(1.1112369) q[0];
rz(-1.5574607) q[2];
sx q[2];
rz(-0.91382256) q[2];
sx q[2];
rz(-1.8767534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9135838) q[1];
sx q[1];
rz(-2.5435547) q[1];
sx q[1];
rz(-2.0889455) q[1];
x q[2];
rz(2.8178704) q[3];
sx q[3];
rz(-0.60327134) q[3];
sx q[3];
rz(-2.3761185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35932943) q[2];
sx q[2];
rz(-1.9908345) q[2];
sx q[2];
rz(-2.9456054) q[2];
rz(-1.1228784) q[3];
sx q[3];
rz(-0.71176088) q[3];
sx q[3];
rz(-1.4518729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48542431) q[0];
sx q[0];
rz(-0.66101414) q[0];
sx q[0];
rz(-2.0458903) q[0];
rz(-0.51086673) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(-2.9313472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1258233) q[0];
sx q[0];
rz(-2.5298241) q[0];
sx q[0];
rz(1.2248898) q[0];
rz(-pi) q[1];
rz(1.7201682) q[2];
sx q[2];
rz(-2.4813642) q[2];
sx q[2];
rz(2.3200043) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89344674) q[1];
sx q[1];
rz(-0.67091828) q[1];
sx q[1];
rz(0.72709658) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9821591) q[3];
sx q[3];
rz(-1.0108915) q[3];
sx q[3];
rz(-0.066889569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0290587) q[2];
sx q[2];
rz(-1.5052648) q[2];
sx q[2];
rz(-2.0564334) q[2];
rz(-0.30760136) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(-2.4881081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6899684) q[0];
sx q[0];
rz(-1.987048) q[0];
sx q[0];
rz(-1.2183627) q[0];
rz(-0.65144173) q[1];
sx q[1];
rz(-0.94964288) q[1];
sx q[1];
rz(-2.3655187) q[1];
rz(-1.5522926) q[2];
sx q[2];
rz(-0.69102364) q[2];
sx q[2];
rz(-0.24204248) q[2];
rz(-0.92523706) q[3];
sx q[3];
rz(-0.74151562) q[3];
sx q[3];
rz(0.61552463) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

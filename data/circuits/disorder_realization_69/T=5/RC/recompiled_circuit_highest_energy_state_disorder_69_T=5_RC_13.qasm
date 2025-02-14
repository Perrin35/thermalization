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
rz(2.7076758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.514115) q[0];
sx q[0];
rz(-1.6058996) q[0];
sx q[0];
rz(-0.83034626) q[0];
rz(-2.5662759) q[2];
sx q[2];
rz(-1.0566718) q[2];
sx q[2];
rz(-2.7935087) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.190769) q[1];
sx q[1];
rz(-0.94591037) q[1];
sx q[1];
rz(-0.48887078) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1103575) q[3];
sx q[3];
rz(-1.7618351) q[3];
sx q[3];
rz(-1.7241378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71243858) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(2.1377371) q[2];
rz(-2.6058274) q[3];
sx q[3];
rz(-2.2771213) q[3];
sx q[3];
rz(-0.0011477688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4681604) q[0];
sx q[0];
rz(-0.68244451) q[0];
sx q[0];
rz(-0.72714192) q[0];
rz(-3.0739821) q[1];
sx q[1];
rz(-1.3550242) q[1];
sx q[1];
rz(-2.0179857) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7794117) q[0];
sx q[0];
rz(-1.0955278) q[0];
sx q[0];
rz(-1.2527466) q[0];
rz(-pi) q[1];
rz(3.0762663) q[2];
sx q[2];
rz(-1.6931173) q[2];
sx q[2];
rz(-2.9768012) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0331081) q[1];
sx q[1];
rz(-0.29769167) q[1];
sx q[1];
rz(1.9161403) q[1];
rz(-pi) q[2];
rz(2.7026342) q[3];
sx q[3];
rz(-1.135313) q[3];
sx q[3];
rz(1.8451231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2850538) q[2];
sx q[2];
rz(-1.5736138) q[2];
sx q[2];
rz(2.6598568) q[2];
rz(0.5528062) q[3];
sx q[3];
rz(-2.063664) q[3];
sx q[3];
rz(0.089574561) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8168617) q[0];
sx q[0];
rz(-2.6664) q[0];
sx q[0];
rz(-3.0480296) q[0];
rz(-1.7612673) q[1];
sx q[1];
rz(-0.9535791) q[1];
sx q[1];
rz(2.5043452) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475315) q[0];
sx q[0];
rz(-1.1406745) q[0];
sx q[0];
rz(0.57292666) q[0];
x q[1];
rz(-0.60814823) q[2];
sx q[2];
rz(-1.7575193) q[2];
sx q[2];
rz(0.0063467912) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5384851) q[1];
sx q[1];
rz(-2.7432495) q[1];
sx q[1];
rz(1.5038236) q[1];
rz(-pi) q[2];
rz(-2.8508857) q[3];
sx q[3];
rz(-1.5295431) q[3];
sx q[3];
rz(-2.8009529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.20643413) q[2];
sx q[2];
rz(-2.5863402) q[2];
sx q[2];
rz(-0.68378249) q[2];
rz(-2.911496) q[3];
sx q[3];
rz(-1.7347387) q[3];
sx q[3];
rz(0.42207119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049659599) q[0];
sx q[0];
rz(-2.6347418) q[0];
sx q[0];
rz(-1.7012713) q[0];
rz(-2.6662042) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(-2.5604274) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0770532) q[0];
sx q[0];
rz(-2.4145899) q[0];
sx q[0];
rz(-1.6995656) q[0];
x q[1];
rz(2.0330795) q[2];
sx q[2];
rz(-1.7812742) q[2];
sx q[2];
rz(2.0049948) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46266471) q[1];
sx q[1];
rz(-1.8724672) q[1];
sx q[1];
rz(2.660391) q[1];
rz(-pi) q[2];
rz(1.4122333) q[3];
sx q[3];
rz(-1.0845636) q[3];
sx q[3];
rz(-2.9610046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0541957) q[2];
sx q[2];
rz(-2.9554458) q[2];
sx q[2];
rz(-2.6206214) q[2];
rz(-1.6446796) q[3];
sx q[3];
rz(-1.729634) q[3];
sx q[3];
rz(-1.6269256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5882551) q[0];
sx q[0];
rz(-2.8659358) q[0];
sx q[0];
rz(-2.0113373) q[0];
rz(2.0626119) q[1];
sx q[1];
rz(-1.9422453) q[1];
sx q[1];
rz(2.030453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1905907) q[0];
sx q[0];
rz(-1.4050806) q[0];
sx q[0];
rz(-2.034305) q[0];
x q[1];
rz(3.0658998) q[2];
sx q[2];
rz(-0.53097979) q[2];
sx q[2];
rz(2.7493283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1597978) q[1];
sx q[1];
rz(-0.66098958) q[1];
sx q[1];
rz(1.3128408) q[1];
x q[2];
rz(-0.89966117) q[3];
sx q[3];
rz(-2.3658345) q[3];
sx q[3];
rz(0.44246261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40372103) q[2];
sx q[2];
rz(-2.626494) q[2];
sx q[2];
rz(-1.0408638) q[2];
rz(-0.8199842) q[3];
sx q[3];
rz(-1.2330202) q[3];
sx q[3];
rz(-2.5177054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(0.89989221) q[0];
sx q[0];
rz(-2.5427759) q[0];
sx q[0];
rz(-0.65648055) q[0];
rz(-1.7680602) q[1];
sx q[1];
rz(-2.3479925) q[1];
sx q[1];
rz(-1.1318644) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098399408) q[0];
sx q[0];
rz(-0.88623673) q[0];
sx q[0];
rz(-1.3149904) q[0];
x q[1];
rz(-2.4651338) q[2];
sx q[2];
rz(-2.0148811) q[2];
sx q[2];
rz(-1.9537587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1161728) q[1];
sx q[1];
rz(-2.4301404) q[1];
sx q[1];
rz(1.6132212) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9911861) q[3];
sx q[3];
rz(-1.5664706) q[3];
sx q[3];
rz(-1.195601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.446283) q[2];
sx q[2];
rz(-2.535847) q[2];
sx q[2];
rz(-2.821935) q[2];
rz(2.9122635) q[3];
sx q[3];
rz(-1.0234443) q[3];
sx q[3];
rz(0.91013175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0685773) q[0];
sx q[0];
rz(-1.9972766) q[0];
sx q[0];
rz(-1.9101494) q[0];
rz(0.07864174) q[1];
sx q[1];
rz(-1.6784724) q[1];
sx q[1];
rz(-2.9873649) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.94119) q[0];
sx q[0];
rz(-1.2324257) q[0];
sx q[0];
rz(1.75941) q[0];
x q[1];
rz(1.1604105) q[2];
sx q[2];
rz(-1.2912116) q[2];
sx q[2];
rz(-1.9094163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5104766) q[1];
sx q[1];
rz(-2.0492024) q[1];
sx q[1];
rz(-1.9892742) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4108415) q[3];
sx q[3];
rz(-1.7663398) q[3];
sx q[3];
rz(-1.1943749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33588931) q[2];
sx q[2];
rz(-1.615639) q[2];
sx q[2];
rz(1.4914782) q[2];
rz(1.4132181) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(-1.570805) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0989646) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(-1.6866823) q[0];
rz(-0.97575724) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(1.3386493) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8737301) q[0];
sx q[0];
rz(-2.3288245) q[0];
sx q[0];
rz(-2.1186747) q[0];
rz(2.1326124) q[2];
sx q[2];
rz(-0.57719165) q[2];
sx q[2];
rz(1.6704337) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2690791) q[1];
sx q[1];
rz(-0.95035205) q[1];
sx q[1];
rz(-0.2618813) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9429132) q[3];
sx q[3];
rz(-0.94148472) q[3];
sx q[3];
rz(-0.019364186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51155382) q[2];
sx q[2];
rz(-1.0147107) q[2];
sx q[2];
rz(-0.83941984) q[2];
rz(1.9073585) q[3];
sx q[3];
rz(-2.0525565) q[3];
sx q[3];
rz(0.031410005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79427528) q[0];
sx q[0];
rz(-1.7592156) q[0];
sx q[0];
rz(0.41912249) q[0];
rz(1.4979111) q[1];
sx q[1];
rz(-1.3587911) q[1];
sx q[1];
rz(0.43325123) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99892761) q[0];
sx q[0];
rz(-1.6892528) q[0];
sx q[0];
rz(-1.1112369) q[0];
rz(1.584132) q[2];
sx q[2];
rz(-0.91382256) q[2];
sx q[2];
rz(-1.8767534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8318787) q[1];
sx q[1];
rz(-2.0818748) q[1];
sx q[1];
rz(-2.8161777) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57859666) q[3];
sx q[3];
rz(-1.3893327) q[3];
sx q[3];
rz(-0.53574783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35932943) q[2];
sx q[2];
rz(-1.1507582) q[2];
sx q[2];
rz(-0.19598728) q[2];
rz(-1.1228784) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(-1.6897197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48542431) q[0];
sx q[0];
rz(-0.66101414) q[0];
sx q[0];
rz(2.0458903) q[0];
rz(2.6307259) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(-2.9313472) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71111403) q[0];
sx q[0];
rz(-1.0000044) q[0];
sx q[0];
rz(-2.9080703) q[0];
rz(-pi) q[1];
rz(0.91598454) q[2];
sx q[2];
rz(-1.4794) q[2];
sx q[2];
rz(0.86751988) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.044503511) q[1];
sx q[1];
rz(-1.0877481) q[1];
sx q[1];
rz(-2.0562857) q[1];
rz(1.1594335) q[3];
sx q[3];
rz(-1.0108915) q[3];
sx q[3];
rz(-0.066889569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0290587) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(1.0851592) q[2];
rz(-0.30760136) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(-2.4881081) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6899684) q[0];
sx q[0];
rz(-1.1545447) q[0];
sx q[0];
rz(1.92323) q[0];
rz(2.4901509) q[1];
sx q[1];
rz(-0.94964288) q[1];
sx q[1];
rz(-2.3655187) q[1];
rz(-0.87985676) q[2];
sx q[2];
rz(-1.5590038) q[2];
sx q[2];
rz(-1.8270983) q[2];
rz(2.6379587) q[3];
sx q[3];
rz(-1.0009652) q[3];
sx q[3];
rz(2.9611369) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

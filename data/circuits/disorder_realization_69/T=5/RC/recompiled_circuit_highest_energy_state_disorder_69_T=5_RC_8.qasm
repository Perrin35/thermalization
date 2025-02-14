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
rz(4.7799568) q[0];
sx q[0];
rz(10.156375) q[0];
rz(0.44959623) q[1];
sx q[1];
rz(-2.9540588) q[1];
sx q[1];
rz(-2.7076758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62747763) q[0];
sx q[0];
rz(-1.5356931) q[0];
sx q[0];
rz(-0.83034626) q[0];
rz(-2.5662759) q[2];
sx q[2];
rz(-1.0566718) q[2];
sx q[2];
rz(0.34808394) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.190769) q[1];
sx q[1];
rz(-0.94591037) q[1];
sx q[1];
rz(-2.6527219) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7619261) q[3];
sx q[3];
rz(-1.5401296) q[3];
sx q[3];
rz(0.14740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71243858) q[2];
sx q[2];
rz(-1.466208) q[2];
sx q[2];
rz(-1.0038556) q[2];
rz(-0.53576523) q[3];
sx q[3];
rz(-2.2771213) q[3];
sx q[3];
rz(-3.1404449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4681604) q[0];
sx q[0];
rz(-2.4591481) q[0];
sx q[0];
rz(2.4144507) q[0];
rz(0.06761059) q[1];
sx q[1];
rz(-1.3550242) q[1];
sx q[1];
rz(1.1236069) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0824995) q[0];
sx q[0];
rz(-1.8525665) q[0];
sx q[0];
rz(0.49651329) q[0];
rz(-pi) q[1];
rz(1.4482165) q[2];
sx q[2];
rz(-1.6356339) q[2];
sx q[2];
rz(1.4139869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1084845) q[1];
sx q[1];
rz(-2.843901) q[1];
sx q[1];
rz(1.9161403) q[1];
x q[2];
rz(-2.3109834) q[3];
sx q[3];
rz(-0.6081444) q[3];
sx q[3];
rz(0.45765314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8565389) q[2];
sx q[2];
rz(-1.5736138) q[2];
sx q[2];
rz(0.48173586) q[2];
rz(2.5887865) q[3];
sx q[3];
rz(-1.0779287) q[3];
sx q[3];
rz(-3.0520181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.1880136) q[1];
sx q[1];
rz(0.63724744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15055949) q[0];
sx q[0];
rz(-0.70165082) q[0];
sx q[0];
rz(-0.70233624) q[0];
rz(-pi) q[1];
x q[1];
rz(1.797051) q[2];
sx q[2];
rz(-2.166894) q[2];
sx q[2];
rz(-1.6929733) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67575021) q[1];
sx q[1];
rz(-1.1733965) q[1];
sx q[1];
rz(-3.1134362) q[1];
x q[2];
rz(2.9985688) q[3];
sx q[3];
rz(-0.29353729) q[3];
sx q[3];
rz(-2.0484201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.20643413) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(-2.4578102) q[2];
rz(2.911496) q[3];
sx q[3];
rz(-1.7347387) q[3];
sx q[3];
rz(-0.42207119) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0919331) q[0];
sx q[0];
rz(-2.6347418) q[0];
sx q[0];
rz(1.7012713) q[0];
rz(-2.6662042) q[1];
sx q[1];
rz(-0.63440591) q[1];
sx q[1];
rz(2.5604274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5388881) q[0];
sx q[0];
rz(-1.4853444) q[0];
sx q[0];
rz(-0.84792015) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9072806) q[2];
sx q[2];
rz(-2.0221124) q[2];
sx q[2];
rz(-2.6036604) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6789279) q[1];
sx q[1];
rz(-1.8724672) q[1];
sx q[1];
rz(0.48120166) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.851296) q[3];
sx q[3];
rz(-2.6321342) q[3];
sx q[3];
rz(-0.51028937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0873969) q[2];
sx q[2];
rz(-0.18614686) q[2];
sx q[2];
rz(2.6206214) q[2];
rz(-1.4969131) q[3];
sx q[3];
rz(-1.729634) q[3];
sx q[3];
rz(-1.514667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5533376) q[0];
sx q[0];
rz(-2.8659358) q[0];
sx q[0];
rz(-1.1302554) q[0];
rz(-1.0789808) q[1];
sx q[1];
rz(-1.9422453) q[1];
sx q[1];
rz(2.030453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6040627) q[0];
sx q[0];
rz(-1.1141234) q[0];
sx q[0];
rz(2.9567493) q[0];
rz(1.5264185) q[2];
sx q[2];
rz(-1.0414972) q[2];
sx q[2];
rz(0.30454301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1597978) q[1];
sx q[1];
rz(-0.66098958) q[1];
sx q[1];
rz(-1.8287519) q[1];
rz(2.2258198) q[3];
sx q[3];
rz(-1.1202328) q[3];
sx q[3];
rz(-2.5289867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7378716) q[2];
sx q[2];
rz(-0.51509866) q[2];
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
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89989221) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(0.65648055) q[0];
rz(-1.7680602) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(-2.0097282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8478125) q[0];
sx q[0];
rz(-0.72350693) q[0];
sx q[0];
rz(-0.30059867) q[0];
rz(-pi) q[1];
rz(-0.64986367) q[2];
sx q[2];
rz(-2.3520326) q[2];
sx q[2];
rz(-2.2672578) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0254198) q[1];
sx q[1];
rz(-2.4301404) q[1];
sx q[1];
rz(1.5283714) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9911861) q[3];
sx q[3];
rz(-1.5664706) q[3];
sx q[3];
rz(1.195601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.446283) q[2];
sx q[2];
rz(-0.60574564) q[2];
sx q[2];
rz(-0.31965762) q[2];
rz(-2.9122635) q[3];
sx q[3];
rz(-2.1181483) q[3];
sx q[3];
rz(0.91013175) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0685773) q[0];
sx q[0];
rz(-1.1443161) q[0];
sx q[0];
rz(-1.2314433) q[0];
rz(-0.07864174) q[1];
sx q[1];
rz(-1.6784724) q[1];
sx q[1];
rz(2.9873649) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.94119) q[0];
sx q[0];
rz(-1.909167) q[0];
sx q[0];
rz(1.3821826) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8381589) q[2];
sx q[2];
rz(-1.1772441) q[2];
sx q[2];
rz(2.9224666) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7423305) q[1];
sx q[1];
rz(-0.62472099) q[1];
sx q[1];
rz(-2.4769251) q[1];
x q[2];
rz(-0.46040543) q[3];
sx q[3];
rz(-0.45259991) q[3];
sx q[3];
rz(2.3456338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33588931) q[2];
sx q[2];
rz(-1.615639) q[2];
sx q[2];
rz(-1.4914782) q[2];
rz(-1.4132181) q[3];
sx q[3];
rz(-1.3946984) q[3];
sx q[3];
rz(-1.570805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042628057) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(-1.4549103) q[0];
rz(-0.97575724) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(1.3386493) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2678626) q[0];
sx q[0];
rz(-2.3288245) q[0];
sx q[0];
rz(-1.0229179) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8076914) q[2];
sx q[2];
rz(-1.0907764) q[2];
sx q[2];
rz(2.1154108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83709675) q[1];
sx q[1];
rz(-0.66667914) q[1];
sx q[1];
rz(1.9183938) q[1];
rz(2.4782031) q[3];
sx q[3];
rz(-1.8691392) q[3];
sx q[3];
rz(-1.364352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51155382) q[2];
sx q[2];
rz(-2.126882) q[2];
sx q[2];
rz(0.83941984) q[2];
rz(1.9073585) q[3];
sx q[3];
rz(-1.0890361) q[3];
sx q[3];
rz(3.1101826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3473174) q[0];
sx q[0];
rz(-1.7592156) q[0];
sx q[0];
rz(-2.7224702) q[0];
rz(1.4979111) q[1];
sx q[1];
rz(-1.3587911) q[1];
sx q[1];
rz(-2.7083414) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99892761) q[0];
sx q[0];
rz(-1.4523399) q[0];
sx q[0];
rz(-2.0303557) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4845759) q[2];
sx q[2];
rz(-1.5602367) q[2];
sx q[2];
rz(-2.8274909) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2280088) q[1];
sx q[1];
rz(-2.5435547) q[1];
sx q[1];
rz(1.0526471) q[1];
rz(-1.3550536) q[3];
sx q[3];
rz(-2.1387055) q[3];
sx q[3];
rz(1.9892094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7822632) q[2];
sx q[2];
rz(-1.9908345) q[2];
sx q[2];
rz(2.9456054) q[2];
rz(-2.0187142) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(-1.4518729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48542431) q[0];
sx q[0];
rz(-2.4805785) q[0];
sx q[0];
rz(-1.0957023) q[0];
rz(2.6307259) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(0.21024545) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4097262) q[0];
sx q[0];
rz(-1.3748226) q[0];
sx q[0];
rz(-2.1541697) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2256081) q[2];
sx q[2];
rz(-1.6621926) q[2];
sx q[2];
rz(-2.2740728) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0970891) q[1];
sx q[1];
rz(-2.0538446) q[1];
sx q[1];
rz(-2.0562857) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1594335) q[3];
sx q[3];
rz(-2.1307011) q[3];
sx q[3];
rz(-3.0747031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0290587) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(1.0851592) q[2];
rz(0.30760136) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(2.4881081) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6899684) q[0];
sx q[0];
rz(-1.1545447) q[0];
sx q[0];
rz(1.92323) q[0];
rz(-0.65144173) q[1];
sx q[1];
rz(-0.94964288) q[1];
sx q[1];
rz(-2.3655187) q[1];
rz(1.5522926) q[2];
sx q[2];
rz(-2.450569) q[2];
sx q[2];
rz(2.8995502) q[2];
rz(2.2163556) q[3];
sx q[3];
rz(-0.74151562) q[3];
sx q[3];
rz(0.61552463) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

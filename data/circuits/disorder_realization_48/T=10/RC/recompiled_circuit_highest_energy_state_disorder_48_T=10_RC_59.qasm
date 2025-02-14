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
rz(-1.7464632) q[0];
sx q[0];
rz(-0.75463086) q[0];
sx q[0];
rz(-2.057743) q[0];
rz(2.2639182) q[1];
sx q[1];
rz(4.35507) q[1];
sx q[1];
rz(9.9014643) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86171976) q[0];
sx q[0];
rz(-1.828816) q[0];
sx q[0];
rz(1.1610384) q[0];
rz(0.4557336) q[2];
sx q[2];
rz(-1.991716) q[2];
sx q[2];
rz(-1.2088261) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.013416524) q[1];
sx q[1];
rz(-2.8833296) q[1];
sx q[1];
rz(-1.5034666) q[1];
x q[2];
rz(-0.4313978) q[3];
sx q[3];
rz(-1.7713303) q[3];
sx q[3];
rz(-1.3667184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1276663) q[2];
sx q[2];
rz(-1.650859) q[2];
sx q[2];
rz(-0.15065436) q[2];
rz(-1.539218) q[3];
sx q[3];
rz(-0.37709245) q[3];
sx q[3];
rz(-0.25130513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056331228) q[0];
sx q[0];
rz(-1.769861) q[0];
sx q[0];
rz(0.37013176) q[0];
rz(-2.6689957) q[1];
sx q[1];
rz(-0.31817803) q[1];
sx q[1];
rz(-0.95742375) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8214398) q[0];
sx q[0];
rz(-1.6509045) q[0];
sx q[0];
rz(-1.0715241) q[0];
rz(-pi) q[1];
rz(1.1195807) q[2];
sx q[2];
rz(-1.4133845) q[2];
sx q[2];
rz(-1.5296641) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4366001) q[1];
sx q[1];
rz(-0.7448405) q[1];
sx q[1];
rz(3.0158494) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9022835) q[3];
sx q[3];
rz(-0.9210862) q[3];
sx q[3];
rz(-3.0420764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0051673278) q[2];
sx q[2];
rz(-0.74446669) q[2];
sx q[2];
rz(-1.0235419) q[2];
rz(1.7510022) q[3];
sx q[3];
rz(-1.7460881) q[3];
sx q[3];
rz(-1.8243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.35109529) q[0];
sx q[0];
rz(-2.1320765) q[0];
sx q[0];
rz(0.56280953) q[0];
rz(-1.5142745) q[1];
sx q[1];
rz(-1.7104251) q[1];
sx q[1];
rz(-0.079158457) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25185967) q[0];
sx q[0];
rz(-1.2556228) q[0];
sx q[0];
rz(2.6227975) q[0];
rz(-1.6116033) q[2];
sx q[2];
rz(-2.0611603) q[2];
sx q[2];
rz(0.23863579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97387633) q[1];
sx q[1];
rz(-1.1493681) q[1];
sx q[1];
rz(-0.80444471) q[1];
rz(-pi) q[2];
rz(0.7003673) q[3];
sx q[3];
rz(-2.1383965) q[3];
sx q[3];
rz(1.8710473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2341653) q[2];
sx q[2];
rz(-1.558446) q[2];
sx q[2];
rz(-1.3321715) q[2];
rz(2.7857156) q[3];
sx q[3];
rz(-1.1938813) q[3];
sx q[3];
rz(-0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6818162) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(2.0303149) q[0];
rz(0.92012826) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(0.2535325) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0363844) q[0];
sx q[0];
rz(-1.0691133) q[0];
sx q[0];
rz(1.9355477) q[0];
rz(-pi) q[1];
rz(1.4561256) q[2];
sx q[2];
rz(-1.7301699) q[2];
sx q[2];
rz(1.3068975) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56139466) q[1];
sx q[1];
rz(-0.78402482) q[1];
sx q[1];
rz(0.26303388) q[1];
x q[2];
rz(-2.0959804) q[3];
sx q[3];
rz(-0.56852698) q[3];
sx q[3];
rz(2.9884058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8186875) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(-1.6820071) q[2];
rz(0.5736351) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(0.044376317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6673073) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(1.7942418) q[0];
rz(-0.59066331) q[1];
sx q[1];
rz(-1.2202411) q[1];
sx q[1];
rz(-1.4264872) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9085981) q[0];
sx q[0];
rz(-2.5979266) q[0];
sx q[0];
rz(-2.3999016) q[0];
rz(-2.6408004) q[2];
sx q[2];
rz(-2.4363849) q[2];
sx q[2];
rz(1.0358182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8914681) q[1];
sx q[1];
rz(-1.824675) q[1];
sx q[1];
rz(-1.5230383) q[1];
x q[2];
rz(-2.9699202) q[3];
sx q[3];
rz(-2.7816157) q[3];
sx q[3];
rz(-1.5300919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23919375) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(1.1725461) q[2];
rz(-0.086056195) q[3];
sx q[3];
rz(-2.1890169) q[3];
sx q[3];
rz(-2.9854767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758783) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-2.2789047) q[0];
rz(-0.24738303) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(-2.6695796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.056045) q[0];
sx q[0];
rz(-1.9084832) q[0];
sx q[0];
rz(2.4187536) q[0];
x q[1];
rz(-2.9905867) q[2];
sx q[2];
rz(-1.0698294) q[2];
sx q[2];
rz(-2.6789064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1122163) q[1];
sx q[1];
rz(-0.35307717) q[1];
sx q[1];
rz(-0.53332163) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2788762) q[3];
sx q[3];
rz(-2.4987767) q[3];
sx q[3];
rz(-1.1228648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(1.7426573) q[2];
rz(1.6861964) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(2.5808891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48651925) q[0];
sx q[0];
rz(-1.9955248) q[0];
sx q[0];
rz(2.6229677) q[0];
rz(-0.71867603) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(1.8124883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78120056) q[0];
sx q[0];
rz(-2.5787163) q[0];
sx q[0];
rz(0.84873523) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7460004) q[2];
sx q[2];
rz(-2.8602798) q[2];
sx q[2];
rz(-0.11872053) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80793984) q[1];
sx q[1];
rz(-0.83631714) q[1];
sx q[1];
rz(0.044313852) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1217478) q[3];
sx q[3];
rz(-2.1571419) q[3];
sx q[3];
rz(1.3278409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79080498) q[2];
sx q[2];
rz(-2.4420276) q[2];
sx q[2];
rz(-0.038185509) q[2];
rz(-0.83032483) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(-3.0939046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6787978) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(-1.3702673) q[0];
rz(0.38617745) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(0.6699627) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7350175) q[0];
sx q[0];
rz(-1.1018714) q[0];
sx q[0];
rz(1.0107299) q[0];
x q[1];
rz(-1.6514844) q[2];
sx q[2];
rz(-0.85570691) q[2];
sx q[2];
rz(2.2632368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.252287) q[1];
sx q[1];
rz(-0.60599594) q[1];
sx q[1];
rz(1.3713981) q[1];
rz(-pi) q[2];
rz(-3.1042323) q[3];
sx q[3];
rz(-1.3589763) q[3];
sx q[3];
rz(2.7461014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2204444) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(-2.0849126) q[2];
rz(1.4165261) q[3];
sx q[3];
rz(-1.8959911) q[3];
sx q[3];
rz(1.8511124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0678299) q[0];
sx q[0];
rz(-1.0404328) q[0];
sx q[0];
rz(1.006806) q[0];
rz(-2.6489068) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(2.3115092) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5093436) q[0];
sx q[0];
rz(-1.2918988) q[0];
sx q[0];
rz(1.6053902) q[0];
rz(1.3557498) q[2];
sx q[2];
rz(-2.1587662) q[2];
sx q[2];
rz(2.7933592) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3398847) q[1];
sx q[1];
rz(-2.9045068) q[1];
sx q[1];
rz(2.4380142) q[1];
rz(-pi) q[2];
rz(-1.5437627) q[3];
sx q[3];
rz(-2.6984331) q[3];
sx q[3];
rz(-1.5704607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3229708) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(1.0969561) q[2];
rz(-1.8638301) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(2.4975615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85035664) q[0];
sx q[0];
rz(-1.1024029) q[0];
sx q[0];
rz(0.37503234) q[0];
rz(-0.11861435) q[1];
sx q[1];
rz(-1.7053441) q[1];
sx q[1];
rz(-2.2172701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345006) q[0];
sx q[0];
rz(-2.1372927) q[0];
sx q[0];
rz(1.2248216) q[0];
x q[1];
rz(-0.67258851) q[2];
sx q[2];
rz(-1.9832356) q[2];
sx q[2];
rz(-0.46063603) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1076775) q[1];
sx q[1];
rz(-2.7795791) q[1];
sx q[1];
rz(3.0605143) q[1];
rz(-0.44845805) q[3];
sx q[3];
rz(-1.8844386) q[3];
sx q[3];
rz(-2.3109316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2529926) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(-2.2484153) q[2];
rz(-2.0231953) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(2.4848188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3746344) q[0];
sx q[0];
rz(-1.081291) q[0];
sx q[0];
rz(-1.3670856) q[0];
rz(0.18192667) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(0.091808783) q[2];
sx q[2];
rz(-1.8845857) q[2];
sx q[2];
rz(1.2086445) q[2];
rz(-2.4891067) q[3];
sx q[3];
rz(-2.5204044) q[3];
sx q[3];
rz(1.6922097) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

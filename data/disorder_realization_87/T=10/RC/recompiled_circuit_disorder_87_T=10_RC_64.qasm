OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(-0.72609225) q[0];
sx q[0];
rz(-0.2015764) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0694323) q[0];
sx q[0];
rz(-1.9919112) q[0];
sx q[0];
rz(0.62299563) q[0];
x q[1];
rz(-0.61203981) q[2];
sx q[2];
rz(-1.1571552) q[2];
sx q[2];
rz(-0.94758247) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1569251) q[1];
sx q[1];
rz(-2.942454) q[1];
sx q[1];
rz(-1.9763293) q[1];
rz(-0.097221656) q[3];
sx q[3];
rz(-1.3029649) q[3];
sx q[3];
rz(-2.3107446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(-0.086159555) q[2];
rz(2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(-1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.8151059) q[0];
rz(-1.2558698) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(-0.27145162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6730246) q[0];
sx q[0];
rz(-0.66267555) q[0];
sx q[0];
rz(-2.2492692) q[0];
rz(-pi) q[1];
rz(2.62918) q[2];
sx q[2];
rz(-0.57704848) q[2];
sx q[2];
rz(1.7713828) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27122341) q[1];
sx q[1];
rz(-1.1961812) q[1];
sx q[1];
rz(-0.80925525) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7852887) q[3];
sx q[3];
rz(-2.2998527) q[3];
sx q[3];
rz(-2.9475398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0945956) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(-2.5644152) q[2];
rz(-2.2180637) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(-1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24401027) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(-0.14257167) q[0];
rz(1.3525195) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(0.20908633) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3554879) q[0];
sx q[0];
rz(-2.3007563) q[0];
sx q[0];
rz(-0.79746042) q[0];
rz(-pi) q[1];
x q[1];
rz(2.694391) q[2];
sx q[2];
rz(-2.3502091) q[2];
sx q[2];
rz(-0.42391047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0002034) q[1];
sx q[1];
rz(-0.63265002) q[1];
sx q[1];
rz(-2.3246517) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7861869) q[3];
sx q[3];
rz(-1.8672018) q[3];
sx q[3];
rz(-0.6682369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7769988) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(-2.7374632) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(-2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(2.1787815) q[0];
rz(-2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(3.1406291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92361802) q[0];
sx q[0];
rz(-2.1305363) q[0];
sx q[0];
rz(-2.2773507) q[0];
rz(-pi) q[1];
rz(0.12400603) q[2];
sx q[2];
rz(-1.5702855) q[2];
sx q[2];
rz(-2.309547) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2767267) q[1];
sx q[1];
rz(-0.58502561) q[1];
sx q[1];
rz(1.3328972) q[1];
rz(-pi) q[2];
rz(-2.5721781) q[3];
sx q[3];
rz(-1.4953519) q[3];
sx q[3];
rz(2.8358592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0756388) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(-3.0299419) q[2];
rz(2.3305437) q[3];
sx q[3];
rz(-0.45447293) q[3];
sx q[3];
rz(-0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1902996) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(2.7850889) q[0];
rz(-2.6351392) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(-0.26062632) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168419) q[0];
sx q[0];
rz(-1.6126313) q[0];
sx q[0];
rz(-0.14194685) q[0];
rz(-0.79029681) q[2];
sx q[2];
rz(-2.7919263) q[2];
sx q[2];
rz(-0.5860354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.534429) q[1];
sx q[1];
rz(-0.57463127) q[1];
sx q[1];
rz(-3.1092005) q[1];
rz(-pi) q[2];
rz(-0.77633206) q[3];
sx q[3];
rz(-0.20271248) q[3];
sx q[3];
rz(-0.67938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(2.9122706) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(1.0361766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.649958) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4947858) q[0];
sx q[0];
rz(-1.8456869) q[0];
sx q[0];
rz(1.6199714) q[0];
x q[1];
rz(-2.6238407) q[2];
sx q[2];
rz(-1.8038097) q[2];
sx q[2];
rz(-2.0875967) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1097088) q[1];
sx q[1];
rz(-1.110853) q[1];
sx q[1];
rz(-1.0201395) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3544975) q[3];
sx q[3];
rz(-1.6445451) q[3];
sx q[3];
rz(-2.4019965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.002939) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(-0.027739851) q[2];
rz(0.49267832) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7572927) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(3.0199155) q[0];
rz(-1.1514459) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-0.40245232) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0284751) q[0];
sx q[0];
rz(-2.3987781) q[0];
sx q[0];
rz(1.9755367) q[0];
rz(-pi) q[1];
rz(-1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(0.41551057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23558815) q[1];
sx q[1];
rz(-1.8530122) q[1];
sx q[1];
rz(1.7547592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4943069) q[3];
sx q[3];
rz(-0.34919958) q[3];
sx q[3];
rz(-2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(2.2793615) q[2];
rz(-2.6640653) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7858793) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(0.73295897) q[0];
rz(3.0015302) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(-1.0345116) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0089793423) q[0];
sx q[0];
rz(-1.7404557) q[0];
sx q[0];
rz(-1.3404113) q[0];
rz(-pi) q[1];
rz(2.5402252) q[2];
sx q[2];
rz(-0.12430087) q[2];
sx q[2];
rz(2.3483495) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4701177) q[1];
sx q[1];
rz(-0.56709328) q[1];
sx q[1];
rz(-1.5017897) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41820742) q[3];
sx q[3];
rz(-1.4718664) q[3];
sx q[3];
rz(-1.6855155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2404279) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(0.77587664) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4416606) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(3.124776) q[0];
rz(0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-0.7787849) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.461207) q[0];
sx q[0];
rz(-1.2347504) q[0];
sx q[0];
rz(1.8340322) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0428033) q[2];
sx q[2];
rz(-1.325377) q[2];
sx q[2];
rz(2.5185042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3719337) q[1];
sx q[1];
rz(-1.9366591) q[1];
sx q[1];
rz(2.8918173) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3498464) q[3];
sx q[3];
rz(-1.6337992) q[3];
sx q[3];
rz(-0.27730478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(0.71643913) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.7313892) q[3];
sx q[3];
rz(-1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.107782) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-2.9901436) q[0];
rz(-2.9653213) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(-2.418628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.003119) q[0];
sx q[0];
rz(-1.6771294) q[0];
sx q[0];
rz(-1.4195725) q[0];
rz(-1.1638219) q[2];
sx q[2];
rz(-2.7343035) q[2];
sx q[2];
rz(-2.2182857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1180229) q[1];
sx q[1];
rz(-1.5155063) q[1];
sx q[1];
rz(1.0667332) q[1];
rz(-pi) q[2];
rz(2.74182) q[3];
sx q[3];
rz(-2.3120566) q[3];
sx q[3];
rz(-1.5012036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4621949) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(-2.6160713) q[2];
rz(-2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(1.0992959) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(2.2407871) q[2];
sx q[2];
rz(-1.1849891) q[2];
sx q[2];
rz(1.6580788) q[2];
rz(2.8156149) q[3];
sx q[3];
rz(-1.3578284) q[3];
sx q[3];
rz(1.1708543) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

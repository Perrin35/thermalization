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
rz(1.243408) q[0];
sx q[0];
rz(-1.5771447) q[0];
sx q[0];
rz(0.91140437) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(-2.8652006) q[1];
sx q[1];
rz(-0.82958329) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88032915) q[0];
sx q[0];
rz(-1.9682878) q[0];
sx q[0];
rz(-0.3573768) q[0];
rz(2.8943823) q[2];
sx q[2];
rz(-1.6610378) q[2];
sx q[2];
rz(1.6177141) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.098886996) q[1];
sx q[1];
rz(-1.8708085) q[1];
sx q[1];
rz(0.76619958) q[1];
x q[2];
rz(1.7237648) q[3];
sx q[3];
rz(-1.139183) q[3];
sx q[3];
rz(-0.6237517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36898819) q[2];
sx q[2];
rz(-2.172894) q[2];
sx q[2];
rz(-1.9770835) q[2];
rz(-2.4146967) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(-0.070778457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-2.4788719) q[0];
sx q[0];
rz(2.8501999) q[0];
rz(0.60028752) q[1];
sx q[1];
rz(-2.184506) q[1];
sx q[1];
rz(-2.9765863) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3032119) q[0];
sx q[0];
rz(-1.8329805) q[0];
sx q[0];
rz(-3.068046) q[0];
rz(-pi) q[1];
rz(-1.1218614) q[2];
sx q[2];
rz(-1.672272) q[2];
sx q[2];
rz(-1.6368653) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7126727) q[1];
sx q[1];
rz(-2.4578252) q[1];
sx q[1];
rz(3.0316975) q[1];
rz(-pi) q[2];
rz(-0.24834541) q[3];
sx q[3];
rz(-1.6237061) q[3];
sx q[3];
rz(-1.2179045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0642285) q[2];
sx q[2];
rz(-2.1891429) q[2];
sx q[2];
rz(-0.23951086) q[2];
rz(-0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(-3.0467564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8851682) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(-1.5688131) q[0];
rz(-0.39372152) q[1];
sx q[1];
rz(-2.5865793) q[1];
sx q[1];
rz(-0.47600019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6248019) q[0];
sx q[0];
rz(-1.778852) q[0];
sx q[0];
rz(-2.1387908) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3838686) q[2];
sx q[2];
rz(-1.6003967) q[2];
sx q[2];
rz(-1.5160949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.16533599) q[1];
sx q[1];
rz(-1.9812728) q[1];
sx q[1];
rz(-0.93595502) q[1];
x q[2];
rz(-2.1952403) q[3];
sx q[3];
rz(-2.7331684) q[3];
sx q[3];
rz(-0.32865903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8126882) q[2];
sx q[2];
rz(-2.4220971) q[2];
sx q[2];
rz(-2.1997814) q[2];
rz(-1.4224667) q[3];
sx q[3];
rz(-2.7673281) q[3];
sx q[3];
rz(-0.67371887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.6085251) q[0];
rz(-2.7948921) q[1];
sx q[1];
rz(-2.0997212) q[1];
sx q[1];
rz(2.9516721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6081734) q[0];
sx q[0];
rz(-2.3629239) q[0];
sx q[0];
rz(-2.7649568) q[0];
rz(-pi) q[1];
rz(1.0713351) q[2];
sx q[2];
rz(-0.53218666) q[2];
sx q[2];
rz(-2.4376873) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.080481361) q[1];
sx q[1];
rz(-2.5330392) q[1];
sx q[1];
rz(-1.4707277) q[1];
x q[2];
rz(1.1433825) q[3];
sx q[3];
rz(-1.9460925) q[3];
sx q[3];
rz(2.0209598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3612264) q[2];
sx q[2];
rz(-2.2915338) q[2];
sx q[2];
rz(2.7596607) q[2];
rz(-3.0319038) q[3];
sx q[3];
rz(-1.0332801) q[3];
sx q[3];
rz(3.0221353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(3.0383976) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(0.84684816) q[0];
rz(-0.13941828) q[1];
sx q[1];
rz(-2.3250695) q[1];
sx q[1];
rz(0.74904186) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7199166) q[0];
sx q[0];
rz(-2.4329209) q[0];
sx q[0];
rz(0.22765521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7125569) q[2];
sx q[2];
rz(-2.0516011) q[2];
sx q[2];
rz(-2.9713809) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3440123) q[1];
sx q[1];
rz(-1.3770119) q[1];
sx q[1];
rz(-1.4414201) q[1];
rz(-pi) q[2];
rz(0.11744515) q[3];
sx q[3];
rz(-1.0721888) q[3];
sx q[3];
rz(1.487446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2060812) q[2];
sx q[2];
rz(-1.8847909) q[2];
sx q[2];
rz(2.0743267) q[2];
rz(-0.71077985) q[3];
sx q[3];
rz(-1.4830517) q[3];
sx q[3];
rz(1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34273219) q[0];
sx q[0];
rz(-1.6628168) q[0];
sx q[0];
rz(-0.96858281) q[0];
rz(0.69106483) q[1];
sx q[1];
rz(-1.7993118) q[1];
sx q[1];
rz(1.3074646) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.867315) q[0];
sx q[0];
rz(-2.4558892) q[0];
sx q[0];
rz(-0.4088485) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5619823) q[2];
sx q[2];
rz(-2.7977849) q[2];
sx q[2];
rz(1.0506309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1561574) q[1];
sx q[1];
rz(-1.623454) q[1];
sx q[1];
rz(1.9582204) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9581466) q[3];
sx q[3];
rz(-1.8357092) q[3];
sx q[3];
rz(-1.9062756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4751733) q[2];
sx q[2];
rz(-0.23491493) q[2];
sx q[2];
rz(1.8443745) q[2];
rz(1.3951067) q[3];
sx q[3];
rz(-1.8078943) q[3];
sx q[3];
rz(1.6102128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7995826) q[0];
sx q[0];
rz(-2.3260249) q[0];
sx q[0];
rz(-0.94014257) q[0];
rz(-0.6483342) q[1];
sx q[1];
rz(-1.9377361) q[1];
sx q[1];
rz(-2.8727093) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9342182) q[0];
sx q[0];
rz(-1.3257799) q[0];
sx q[0];
rz(0.55565683) q[0];
rz(-pi) q[1];
rz(-2.4388695) q[2];
sx q[2];
rz(-1.3303192) q[2];
sx q[2];
rz(0.81288494) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9536679) q[1];
sx q[1];
rz(-1.8005193) q[1];
sx q[1];
rz(0.98120316) q[1];
x q[2];
rz(2.712065) q[3];
sx q[3];
rz(-1.8801573) q[3];
sx q[3];
rz(1.7842899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1342643) q[2];
sx q[2];
rz(-1.8625883) q[2];
sx q[2];
rz(0.68592611) q[2];
rz(2.4053597) q[3];
sx q[3];
rz(-1.0400306) q[3];
sx q[3];
rz(-0.22549103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4418884) q[0];
sx q[0];
rz(-1.9206973) q[0];
sx q[0];
rz(0.7487444) q[0];
rz(-0.092983149) q[1];
sx q[1];
rz(-0.75983202) q[1];
sx q[1];
rz(-1.071113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37866805) q[0];
sx q[0];
rz(-1.0707741) q[0];
sx q[0];
rz(2.0471694) q[0];
x q[1];
rz(-2.0290211) q[2];
sx q[2];
rz(-2.1734218) q[2];
sx q[2];
rz(-2.9650709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0784356) q[1];
sx q[1];
rz(-2.401135) q[1];
sx q[1];
rz(2.034212) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0893965) q[3];
sx q[3];
rz(-1.9813271) q[3];
sx q[3];
rz(-3.0018788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6894655) q[2];
sx q[2];
rz(-0.63673821) q[2];
sx q[2];
rz(-0.20991906) q[2];
rz(-3.012015) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(-1.5773704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.121948) q[0];
sx q[0];
rz(-2.8804998) q[0];
sx q[0];
rz(3.1280532) q[0];
rz(-2.6059222) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(-3.1279235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5807728) q[0];
sx q[0];
rz(-0.81308621) q[0];
sx q[0];
rz(-1.2570501) q[0];
rz(-pi) q[1];
rz(2.7510468) q[2];
sx q[2];
rz(-1.1154113) q[2];
sx q[2];
rz(0.18965882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7175258) q[1];
sx q[1];
rz(-2.2560511) q[1];
sx q[1];
rz(-1.1478893) q[1];
rz(-pi) q[2];
rz(1.1982273) q[3];
sx q[3];
rz(-2.2856224) q[3];
sx q[3];
rz(-1.4095588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51836625) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(-1.7784485) q[2];
rz(-2.3220883) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(-0.68377408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118947) q[0];
sx q[0];
rz(-2.2830257) q[0];
sx q[0];
rz(1.6534506) q[0];
rz(0.34291521) q[1];
sx q[1];
rz(-1.5577134) q[1];
sx q[1];
rz(-0.75103474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7900811) q[0];
sx q[0];
rz(-1.292391) q[0];
sx q[0];
rz(-0.06601575) q[0];
x q[1];
rz(-1.1840759) q[2];
sx q[2];
rz(-0.49418517) q[2];
sx q[2];
rz(-1.3665716) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6760867) q[1];
sx q[1];
rz(-1.4526396) q[1];
sx q[1];
rz(-2.5398269) q[1];
rz(-0.033174935) q[3];
sx q[3];
rz(-1.4806595) q[3];
sx q[3];
rz(-1.4226923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.670383) q[2];
sx q[2];
rz(-2.1977916) q[2];
sx q[2];
rz(2.75441) q[2];
rz(-2.6662628) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(-2.3502684) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.249007) q[0];
sx q[0];
rz(-2.1795166) q[0];
sx q[0];
rz(-2.4695061) q[0];
rz(0.36920209) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(1.3849707) q[2];
sx q[2];
rz(-2.2792049) q[2];
sx q[2];
rz(-0.21273108) q[2];
rz(1.9561097) q[3];
sx q[3];
rz(-0.29373071) q[3];
sx q[3];
rz(3.1391524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(-2.1699177) q[0];
sx q[0];
rz(1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(-0.08637698) q[1];
sx q[1];
rz(-3.1187305) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1296596) q[0];
sx q[0];
rz(-1.4225905) q[0];
sx q[0];
rz(-0.11797842) q[0];
rz(-2.4542698) q[2];
sx q[2];
rz(-2.5097373) q[2];
sx q[2];
rz(-0.99422115) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7658246) q[1];
sx q[1];
rz(-1.4210912) q[1];
sx q[1];
rz(1.9181812) q[1];
rz(-0.25306563) q[3];
sx q[3];
rz(-0.96732891) q[3];
sx q[3];
rz(-0.30482182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-3.0549333) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-3.0564953) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3413977) q[0];
sx q[0];
rz(-1.5872103) q[0];
sx q[0];
rz(-2.5192758) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30479635) q[2];
sx q[2];
rz(-2.1509503) q[2];
sx q[2];
rz(-3.0457029) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18499204) q[1];
sx q[1];
rz(-1.4586095) q[1];
sx q[1];
rz(-0.086602028) q[1];
x q[2];
rz(-0.99382932) q[3];
sx q[3];
rz(-2.5964063) q[3];
sx q[3];
rz(-0.59515566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8791085) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(1.4651728) q[2];
rz(-1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6572606) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(-2.8088645) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9860486) q[0];
sx q[0];
rz(-1.4882003) q[0];
sx q[0];
rz(0.47604576) q[0];
rz(-pi) q[1];
rz(2.2134174) q[2];
sx q[2];
rz(-0.60545063) q[2];
sx q[2];
rz(1.6124992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82921925) q[1];
sx q[1];
rz(-0.81554283) q[1];
sx q[1];
rz(-1.6014598) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82377394) q[3];
sx q[3];
rz(-0.71411055) q[3];
sx q[3];
rz(-2.2658474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-0.34003568) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(-1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760647) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(2.4705825) q[0];
rz(-1.7975851) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(-2.9342594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57732338) q[0];
sx q[0];
rz(-3.0992357) q[0];
sx q[0];
rz(1.3576515) q[0];
x q[1];
rz(-3.1111654) q[2];
sx q[2];
rz(-1.7446339) q[2];
sx q[2];
rz(-2.0174842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85992766) q[1];
sx q[1];
rz(-1.7237701) q[1];
sx q[1];
rz(-1.8106736) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51165032) q[3];
sx q[3];
rz(-1.5582065) q[3];
sx q[3];
rz(-2.1714641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8308782) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(2.3332398) q[2];
rz(2.3338142) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5508674) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(2.2461058) q[0];
rz(-0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(0.53363824) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14347178) q[0];
sx q[0];
rz(-0.40092418) q[0];
sx q[0];
rz(-0.801416) q[0];
x q[1];
rz(0.23405481) q[2];
sx q[2];
rz(-0.78163994) q[2];
sx q[2];
rz(2.3603338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4662019) q[1];
sx q[1];
rz(-1.4247822) q[1];
sx q[1];
rz(-0.22202613) q[1];
x q[2];
rz(0.20547262) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(-1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9020033) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(-2.5732102) q[2];
rz(1.9249632) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(-1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5414808) q[0];
sx q[0];
rz(-1.913537) q[0];
sx q[0];
rz(-0.39370763) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7447694) q[2];
sx q[2];
rz(-1.4962215) q[2];
sx q[2];
rz(-0.65149388) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5914454) q[1];
sx q[1];
rz(-2.1131556) q[1];
sx q[1];
rz(-1.28246) q[1];
rz(-pi) q[2];
x q[2];
rz(0.011990487) q[3];
sx q[3];
rz(-1.1569835) q[3];
sx q[3];
rz(3.0179265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(2.693434) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3649243) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(1.8572042) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(0.67214322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1561688) q[0];
sx q[0];
rz(-0.16616136) q[0];
sx q[0];
rz(-1.6155433) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1291814) q[2];
sx q[2];
rz(-1.9600944) q[2];
sx q[2];
rz(-2.0814975) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8183648) q[1];
sx q[1];
rz(-1.0867026) q[1];
sx q[1];
rz(-1.1546043) q[1];
x q[2];
rz(2.4403205) q[3];
sx q[3];
rz(-1.6810732) q[3];
sx q[3];
rz(0.14227223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(0.59857541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(2.9833802) q[0];
rz(1.0868866) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(2.3349082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5935532) q[0];
sx q[0];
rz(-2.002945) q[0];
sx q[0];
rz(3.0120864) q[0];
rz(-pi) q[1];
rz(0.20571795) q[2];
sx q[2];
rz(-1.4652025) q[2];
sx q[2];
rz(0.92668698) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9629509) q[1];
sx q[1];
rz(-2.2370403) q[1];
sx q[1];
rz(-1.9402177) q[1];
rz(-pi) q[2];
rz(-0.032012352) q[3];
sx q[3];
rz(-1.8714928) q[3];
sx q[3];
rz(0.34919448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5640101) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-0.99964833) q[2];
rz(0.10351652) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(0.34067571) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7332471) q[0];
sx q[0];
rz(-2.354963) q[0];
sx q[0];
rz(1.1128845) q[0];
rz(-pi) q[1];
rz(1.2286439) q[2];
sx q[2];
rz(-0.59235448) q[2];
sx q[2];
rz(0.60281384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.454168) q[1];
sx q[1];
rz(-1.3892281) q[1];
sx q[1];
rz(-2.1329692) q[1];
x q[2];
rz(-0.19823234) q[3];
sx q[3];
rz(-1.4054338) q[3];
sx q[3];
rz(0.73393047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-0.49794751) q[2];
rz(-2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(2.6224459) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-2.8630032) q[0];
rz(0.57922286) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(-3.0648807) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0522703) q[0];
sx q[0];
rz(-2.8025586) q[0];
sx q[0];
rz(-1.9498755) q[0];
rz(-pi) q[1];
rz(1.8329343) q[2];
sx q[2];
rz(-0.58912504) q[2];
sx q[2];
rz(0.72074705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0753239) q[1];
sx q[1];
rz(-1.9362209) q[1];
sx q[1];
rz(1.5931904) q[1];
x q[2];
rz(-2.6807908) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(-0.27289665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3672093) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(2.6860766) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(3.07807) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1800304) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(2.5554399) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(-1.099844) q[2];
sx q[2];
rz(-1.8269314) q[2];
sx q[2];
rz(-1.6614428) q[2];
rz(-1.9361817) q[3];
sx q[3];
rz(-1.972651) q[3];
sx q[3];
rz(-2.7046711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
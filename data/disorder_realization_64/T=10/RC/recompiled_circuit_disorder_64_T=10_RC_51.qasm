OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7862608) q[0];
sx q[0];
rz(-0.064602764) q[0];
sx q[0];
rz(0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(-1.3316766) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0639122) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(2.1190686) q[0];
rz(-1.9900436) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-0.56564769) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0845619) q[1];
sx q[1];
rz(-1.7006526) q[1];
sx q[1];
rz(-0.76743857) q[1];
rz(-pi) q[2];
rz(-2.0256151) q[3];
sx q[3];
rz(-2.4664719) q[3];
sx q[3];
rz(-1.7636553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33064476) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(-0.60418207) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169096) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(-2.0781793) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(-3.1399472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7749274) q[0];
sx q[0];
rz(-1.5257376) q[0];
sx q[0];
rz(-0.010618322) q[0];
x q[1];
rz(1.3130982) q[2];
sx q[2];
rz(-0.94444599) q[2];
sx q[2];
rz(1.4305654) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.590608) q[1];
sx q[1];
rz(-1.814517) q[1];
sx q[1];
rz(0.79382146) q[1];
x q[2];
rz(-1.315829) q[3];
sx q[3];
rz(-1.4062738) q[3];
sx q[3];
rz(-2.4428575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.985618) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(2.4334811) q[2];
rz(1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(0.8272585) q[0];
rz(0.0050841252) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(-1.089383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3141146) q[0];
sx q[0];
rz(-0.70686045) q[0];
sx q[0];
rz(0.88721888) q[0];
x q[1];
rz(-1.9956279) q[2];
sx q[2];
rz(-1.0152738) q[2];
sx q[2];
rz(-0.59031634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.032420302) q[1];
sx q[1];
rz(-1.4854327) q[1];
sx q[1];
rz(1.9568155) q[1];
x q[2];
rz(0.79077625) q[3];
sx q[3];
rz(-1.9867992) q[3];
sx q[3];
rz(-2.5115867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0777145) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(-2.2568978) q[2];
rz(-1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(-1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(-2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(-0.35983905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7774178) q[0];
sx q[0];
rz(-0.67876498) q[0];
sx q[0];
rz(-1.2749519) q[0];
x q[1];
rz(2.4904576) q[2];
sx q[2];
rz(-1.6614117) q[2];
sx q[2];
rz(0.091094253) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29405669) q[1];
sx q[1];
rz(-1.9901853) q[1];
sx q[1];
rz(1.1555175) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8538586) q[3];
sx q[3];
rz(-0.44970185) q[3];
sx q[3];
rz(-0.85292294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(-0.7652258) q[2];
rz(-0.75677538) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9969479) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.7255406) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(1.0353154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1128164) q[0];
sx q[0];
rz(-1.1338286) q[0];
sx q[0];
rz(0.90941888) q[0];
rz(-pi) q[1];
rz(1.6100699) q[2];
sx q[2];
rz(-2.0619259) q[2];
sx q[2];
rz(2.479535) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.1086515) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(-0.095468949) q[1];
rz(-2.4953793) q[3];
sx q[3];
rz(-0.34020243) q[3];
sx q[3];
rz(-2.7681805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(-1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-2.9343658) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80400318) q[0];
sx q[0];
rz(-0.9540671) q[0];
sx q[0];
rz(-1.3998652) q[0];
x q[1];
rz(0.60116641) q[2];
sx q[2];
rz(-1.1784369) q[2];
sx q[2];
rz(-2.5208134) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0761557) q[1];
sx q[1];
rz(-1.3216615) q[1];
sx q[1];
rz(-2.8912828) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3239922) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(-0.51018754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(0.72171372) q[2];
rz(-1.2747814) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-0.26003626) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24467829) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(0.73202837) q[0];
rz(-0.0094982068) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(-0.2917372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4371944) q[0];
sx q[0];
rz(-1.6219553) q[0];
sx q[0];
rz(2.722446) q[0];
rz(-pi) q[1];
rz(1.8740011) q[2];
sx q[2];
rz(-1.469194) q[2];
sx q[2];
rz(0.50819699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8768423) q[1];
sx q[1];
rz(-1.2683588) q[1];
sx q[1];
rz(1.7251863) q[1];
rz(-2.6351356) q[3];
sx q[3];
rz(-0.582687) q[3];
sx q[3];
rz(-1.2014233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-0.83731246) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(2.4777381) q[0];
rz(-3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(0.95867872) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3631358) q[0];
sx q[0];
rz(-1.0582557) q[0];
sx q[0];
rz(0.97495671) q[0];
rz(2.8616222) q[2];
sx q[2];
rz(-2.0915871) q[2];
sx q[2];
rz(2.1998646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.465185) q[1];
sx q[1];
rz(-2.0945815) q[1];
sx q[1];
rz(0.1233867) q[1];
x q[2];
rz(-2.557425) q[3];
sx q[3];
rz(-1.468588) q[3];
sx q[3];
rz(2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0083996) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(-0.91910249) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(2.6760496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532928) q[0];
sx q[0];
rz(-2.4294937) q[0];
sx q[0];
rz(-0.13029356) q[0];
rz(0.12840694) q[2];
sx q[2];
rz(-1.1616716) q[2];
sx q[2];
rz(0.94305925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8852639) q[1];
sx q[1];
rz(-2.0205106) q[1];
sx q[1];
rz(1.1880258) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2113308) q[3];
sx q[3];
rz(-1.3350447) q[3];
sx q[3];
rz(-0.93709968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(-2.9368029) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(-2.8695316) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64514226) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-2.0196594) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(-2.5591992) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56349194) q[0];
sx q[0];
rz(-2.2390215) q[0];
sx q[0];
rz(-2.8816954) q[0];
x q[1];
rz(2.9899644) q[2];
sx q[2];
rz(-1.0851589) q[2];
sx q[2];
rz(0.033657311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.13233391) q[1];
sx q[1];
rz(-1.3887822) q[1];
sx q[1];
rz(0.11972129) q[1];
rz(-pi) q[2];
rz(0.21721812) q[3];
sx q[3];
rz(-1.5981042) q[3];
sx q[3];
rz(2.581493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(-2.3790512) q[2];
rz(1.4108346) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.0653771) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-1.8680686) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
rz(2.0740261) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

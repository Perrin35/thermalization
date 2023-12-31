OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(-1.6834747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9561477) q[0];
sx q[0];
rz(-1.3584134) q[0];
sx q[0];
rz(2.3953715) q[0];
rz(-pi) q[1];
rz(-0.027878472) q[2];
sx q[2];
rz(-0.61402245) q[2];
sx q[2];
rz(-0.30464722) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7588501) q[1];
sx q[1];
rz(-1.9470125) q[1];
sx q[1];
rz(1.43169) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17957844) q[3];
sx q[3];
rz(-0.60185963) q[3];
sx q[3];
rz(1.7367401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(2.9620985) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(-0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935788) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5886473) q[0];
sx q[0];
rz(-1.5536904) q[0];
sx q[0];
rz(0.018789142) q[0];
rz(1.0500533) q[2];
sx q[2];
rz(-0.69486952) q[2];
sx q[2];
rz(-2.3827202) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3719912) q[1];
sx q[1];
rz(-0.76906119) q[1];
sx q[1];
rz(3.0197057) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96315893) q[3];
sx q[3];
rz(-0.31664407) q[3];
sx q[3];
rz(2.7505927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-0.31323355) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(-2.5405163) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-1.0916969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9515581) q[0];
sx q[0];
rz(-1.0150195) q[0];
sx q[0];
rz(-0.65727289) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7269084) q[2];
sx q[2];
rz(-2.1846002) q[2];
sx q[2];
rz(2.2874122) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97869067) q[1];
sx q[1];
rz(-2.2802417) q[1];
sx q[1];
rz(2.6497926) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46488071) q[3];
sx q[3];
rz(-2.998623) q[3];
sx q[3];
rz(-2.3425992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(2.8312347) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6757641) q[0];
sx q[0];
rz(-2.7390263) q[0];
sx q[0];
rz(-0.34253828) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15375806) q[2];
sx q[2];
rz(-0.69090828) q[2];
sx q[2];
rz(1.549987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6314108) q[1];
sx q[1];
rz(-1.5522172) q[1];
sx q[1];
rz(1.2127962) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6944828) q[3];
sx q[3];
rz(-1.5354772) q[3];
sx q[3];
rz(-0.24231054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(3.0854026) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.5475387) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-0.87019428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4207626) q[0];
sx q[0];
rz(-1.2569068) q[0];
sx q[0];
rz(-1.785196) q[0];
rz(-1.7587897) q[2];
sx q[2];
rz(-2.2586939) q[2];
sx q[2];
rz(2.5685513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2493077) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(-2.409163) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16964511) q[3];
sx q[3];
rz(-2.1262453) q[3];
sx q[3];
rz(2.5559705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(2.4678521) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.649883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0817889) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(-2.9102737) q[0];
rz(-pi) q[1];
rz(2.4620373) q[2];
sx q[2];
rz(-1.9050042) q[2];
sx q[2];
rz(0.92781767) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7548435) q[1];
sx q[1];
rz(-2.1347087) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(-pi) q[2];
rz(-0.10156472) q[3];
sx q[3];
rz(-1.2079117) q[3];
sx q[3];
rz(-2.9543119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(0.89404026) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82350746) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(3.1076028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720374) q[0];
sx q[0];
rz(-0.79622686) q[0];
sx q[0];
rz(1.7982593) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(0.63461441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6778292) q[1];
sx q[1];
rz(-2.4619108) q[1];
sx q[1];
rz(-1.6512647) q[1];
rz(-pi) q[2];
rz(2.5266685) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(-1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(-2.774033) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(1.4454909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5515585) q[0];
sx q[0];
rz(-0.43338767) q[0];
sx q[0];
rz(-1.5831468) q[0];
rz(-0.58480279) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(0.98758299) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3280914) q[1];
sx q[1];
rz(-0.5760759) q[1];
sx q[1];
rz(-2.8495795) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3950855) q[3];
sx q[3];
rz(-2.3254447) q[3];
sx q[3];
rz(-1.5619123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(1.5173222) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-2.8318185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71000242) q[0];
sx q[0];
rz(-0.54034034) q[0];
sx q[0];
rz(2.8938328) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1321399) q[2];
sx q[2];
rz(-1.3938892) q[2];
sx q[2];
rz(-1.1284019) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.579638) q[1];
sx q[1];
rz(-1.6288174) q[1];
sx q[1];
rz(-1.3647563) q[1];
rz(-2.7731032) q[3];
sx q[3];
rz(-0.62928761) q[3];
sx q[3];
rz(-0.027241782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(-1.2072198) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(-0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.4996128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9046281) q[0];
sx q[0];
rz(-1.6722073) q[0];
sx q[0];
rz(1.8629835) q[0];
x q[1];
rz(0.76403107) q[2];
sx q[2];
rz(-1.521763) q[2];
sx q[2];
rz(-1.3438091) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91883509) q[1];
sx q[1];
rz(-1.7822052) q[1];
sx q[1];
rz(-2.0445776) q[1];
rz(1.7435944) q[3];
sx q[3];
rz(-1.5683335) q[3];
sx q[3];
rz(-0.60009225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88400921) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(2.3161841) q[2];
sx q[2];
rz(-2.5054629) q[2];
sx q[2];
rz(2.7873743) q[2];
rz(0.55340135) q[3];
sx q[3];
rz(-2.3407866) q[3];
sx q[3];
rz(1.8368807) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

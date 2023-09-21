OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(6.0072748) q[0];
sx q[0];
rz(10.732565) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6558134) q[0];
sx q[0];
rz(-1.2278779) q[0];
sx q[0];
rz(1.9302084) q[0];
x q[1];
rz(-0.76531305) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(3.090976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.03304122) q[1];
sx q[1];
rz(-1.3290977) q[1];
sx q[1];
rz(2.7944195) q[1];
rz(1.0899815) q[3];
sx q[3];
rz(-2.7365723) q[3];
sx q[3];
rz(2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(1.4663565) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(-2.9247608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8502055) q[0];
sx q[0];
rz(-0.7190401) q[0];
sx q[0];
rz(1.1262116) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47302834) q[2];
sx q[2];
rz(-1.066726) q[2];
sx q[2];
rz(0.31993983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.195897) q[1];
sx q[1];
rz(-2.2356114) q[1];
sx q[1];
rz(-2.871454) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87644491) q[3];
sx q[3];
rz(-1.0507686) q[3];
sx q[3];
rz(0.84393535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(-1.2878093) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-2.2089829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9661449) q[0];
sx q[0];
rz(-1.517059) q[0];
sx q[0];
rz(1.2530112) q[0];
rz(2.9595397) q[2];
sx q[2];
rz(-1.7943873) q[2];
sx q[2];
rz(1.455866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0886791) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(-1.710379) q[1];
rz(2.2266085) q[3];
sx q[3];
rz(-0.6647771) q[3];
sx q[3];
rz(1.8613601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.4979699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.362975) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(1.0233364) q[0];
rz(-pi) q[1];
rz(1.285032) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(2.6464268) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5768891) q[1];
sx q[1];
rz(-2.7895045) q[1];
sx q[1];
rz(-2.4114154) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1180531) q[3];
sx q[3];
rz(-1.1467883) q[3];
sx q[3];
rz(0.31183576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.74636373) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(-0.70181075) q[2];
rz(-0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(-1.6197846) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002776) q[0];
sx q[0];
rz(-1.321723) q[0];
sx q[0];
rz(1.8427909) q[0];
rz(-1.5830718) q[2];
sx q[2];
rz(-2.2222812) q[2];
sx q[2];
rz(2.2001681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8480307) q[1];
sx q[1];
rz(-1.8082431) q[1];
sx q[1];
rz(0.32229396) q[1];
x q[2];
rz(-2.7311677) q[3];
sx q[3];
rz(-0.99567185) q[3];
sx q[3];
rz(-0.23469532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(-2.0416416) q[2];
rz(2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(-0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0681756) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(-1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-3.0117603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79486217) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(-2.5559588) q[0];
rz(-pi) q[1];
rz(0.99545698) q[2];
sx q[2];
rz(-1.2577004) q[2];
sx q[2];
rz(0.42195937) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1631158) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(0.58660581) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31452175) q[3];
sx q[3];
rz(-0.57146996) q[3];
sx q[3];
rz(2.275327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(1.9355109) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(2.4553305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9335564) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(2.3963388) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58418602) q[2];
sx q[2];
rz(-2.2224732) q[2];
sx q[2];
rz(0.78648957) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2346748) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(-0.075637416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2495743) q[3];
sx q[3];
rz(-1.1676844) q[3];
sx q[3];
rz(1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4454322) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(1.5015645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38599309) q[0];
sx q[0];
rz(-1.0781243) q[0];
sx q[0];
rz(-3.1194869) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6325475) q[2];
sx q[2];
rz(-1.5973063) q[2];
sx q[2];
rz(-2.3960631) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.20128076) q[1];
sx q[1];
rz(-2.8021325) q[1];
sx q[1];
rz(-2.8708354) q[1];
x q[2];
rz(-2.5128965) q[3];
sx q[3];
rz(-1.0843715) q[3];
sx q[3];
rz(1.6823671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35187307) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.8224576) q[2];
rz(-1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(-2.7583292) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(0.35167545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4708913) q[0];
sx q[0];
rz(-1.1180709) q[0];
sx q[0];
rz(-0.41505138) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69182379) q[2];
sx q[2];
rz(-1.3335388) q[2];
sx q[2];
rz(1.1292063) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.44605276) q[1];
sx q[1];
rz(-0.76247588) q[1];
sx q[1];
rz(1.6474849) q[1];
x q[2];
rz(-1.6230691) q[3];
sx q[3];
rz(-1.6802603) q[3];
sx q[3];
rz(2.0930406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7982771) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(-1.8593672) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(0.19432755) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949085) q[0];
sx q[0];
rz(-2.2241728) q[0];
sx q[0];
rz(2.9772467) q[0];
rz(-2.2718272) q[2];
sx q[2];
rz(-1.1681721) q[2];
sx q[2];
rz(-2.082777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1397814) q[1];
sx q[1];
rz(-1.0011295) q[1];
sx q[1];
rz(-0.12057481) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7887573) q[3];
sx q[3];
rz(-0.86943227) q[3];
sx q[3];
rz(2.4234114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.4355961) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(2.1736017) q[2];
sx q[2];
rz(-1.597076) q[2];
sx q[2];
rz(1.1521641) q[2];
rz(1.773949) q[3];
sx q[3];
rz(-0.33645867) q[3];
sx q[3];
rz(-2.4154739) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

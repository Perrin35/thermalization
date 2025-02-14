OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.59392053) q[0];
sx q[0];
rz(-2.5424818) q[0];
sx q[0];
rz(0.17679086) q[0];
rz(-1.085936) q[1];
sx q[1];
rz(-0.50958264) q[1];
sx q[1];
rz(-3.0546319) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6330948) q[0];
sx q[0];
rz(-1.0035536) q[0];
sx q[0];
rz(-1.1192516) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7638788) q[2];
sx q[2];
rz(-1.0079441) q[2];
sx q[2];
rz(-0.70729296) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1198169) q[1];
sx q[1];
rz(-1.4939185) q[1];
sx q[1];
rz(-1.6989811) q[1];
rz(-pi) q[2];
rz(-1.2779034) q[3];
sx q[3];
rz(-2.170142) q[3];
sx q[3];
rz(-1.0224316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95734346) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(-0.59172612) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(2.2653968) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0210719) q[0];
sx q[0];
rz(-0.35728917) q[0];
sx q[0];
rz(1.0907809) q[0];
rz(-0.72129321) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(0.24478197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2782841) q[0];
sx q[0];
rz(-2.0119036) q[0];
sx q[0];
rz(-0.37190227) q[0];
rz(-pi) q[1];
rz(0.4366283) q[2];
sx q[2];
rz(-2.2712913) q[2];
sx q[2];
rz(-1.8281405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3192057) q[1];
sx q[1];
rz(-0.97619769) q[1];
sx q[1];
rz(0.90984224) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7177204) q[3];
sx q[3];
rz(-1.2770481) q[3];
sx q[3];
rz(2.5054641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81405866) q[2];
sx q[2];
rz(-0.56698292) q[2];
sx q[2];
rz(2.2954693) q[2];
rz(-2.7766679) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(2.164771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038079809) q[0];
sx q[0];
rz(-0.80779034) q[0];
sx q[0];
rz(0.14347759) q[0];
rz(1.6733276) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-3.0751244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9212355) q[0];
sx q[0];
rz(-2.4548303) q[0];
sx q[0];
rz(-0.14089091) q[0];
rz(-pi) q[1];
rz(-2.8861553) q[2];
sx q[2];
rz(-2.0991) q[2];
sx q[2];
rz(-2.262399) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4926949) q[1];
sx q[1];
rz(-1.6872511) q[1];
sx q[1];
rz(-1.4728738) q[1];
rz(-pi) q[2];
rz(-0.18174882) q[3];
sx q[3];
rz(-2.3929993) q[3];
sx q[3];
rz(2.4570217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0299783) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(0.31271333) q[2];
rz(2.8800268) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(2.3436782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12544352) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(-3.0545767) q[0];
rz(1.7866987) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(3.0932025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3589879) q[0];
sx q[0];
rz(-1.6986956) q[0];
sx q[0];
rz(-2.3652698) q[0];
rz(-pi) q[1];
rz(0.87746967) q[2];
sx q[2];
rz(-0.23261586) q[2];
sx q[2];
rz(-2.2777429) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34710566) q[1];
sx q[1];
rz(-2.2655267) q[1];
sx q[1];
rz(1.6415651) q[1];
x q[2];
rz(-0.88416962) q[3];
sx q[3];
rz(-0.83809747) q[3];
sx q[3];
rz(0.17894408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4309569) q[2];
sx q[2];
rz(-0.29971665) q[2];
sx q[2];
rz(-0.44130138) q[2];
rz(2.273061) q[3];
sx q[3];
rz(-1.7732311) q[3];
sx q[3];
rz(-1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271673) q[0];
sx q[0];
rz(-1.704957) q[0];
sx q[0];
rz(2.4145678) q[0];
rz(-1.4432888) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(-1.7194933) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073339065) q[0];
sx q[0];
rz(-1.4311106) q[0];
sx q[0];
rz(2.5951067) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.431972) q[2];
sx q[2];
rz(-2.8292252) q[2];
sx q[2];
rz(-1.4210977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.369097) q[1];
sx q[1];
rz(-0.14892347) q[1];
sx q[1];
rz(-0.61222525) q[1];
rz(-2.5636219) q[3];
sx q[3];
rz(-0.9909361) q[3];
sx q[3];
rz(-2.487929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(1.5809853) q[2];
rz(1.3382781) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1230028) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(-0.71989584) q[0];
rz(-0.24049354) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(2.8954411) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3703281) q[0];
sx q[0];
rz(-1.851558) q[0];
sx q[0];
rz(-1.595003) q[0];
rz(-pi) q[1];
rz(1.9754378) q[2];
sx q[2];
rz(-0.85524594) q[2];
sx q[2];
rz(3.0457979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0726021) q[1];
sx q[1];
rz(-1.700987) q[1];
sx q[1];
rz(2.4058002) q[1];
x q[2];
rz(0.4972552) q[3];
sx q[3];
rz(-2.3905633) q[3];
sx q[3];
rz(-2.8145144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4973732) q[2];
sx q[2];
rz(-2.5770598) q[2];
sx q[2];
rz(-0.79968828) q[2];
rz(2.6175446) q[3];
sx q[3];
rz(-2.7553813) q[3];
sx q[3];
rz(-3.1246429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2972357) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(0.921184) q[0];
rz(-1.720403) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(-2.1465819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1608022) q[0];
sx q[0];
rz(-1.1931927) q[0];
sx q[0];
rz(0.23997216) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0757881) q[2];
sx q[2];
rz(-1.1514246) q[2];
sx q[2];
rz(-0.21682993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.095671818) q[1];
sx q[1];
rz(-2.0883882) q[1];
sx q[1];
rz(-1.0817097) q[1];
x q[2];
rz(2.7217676) q[3];
sx q[3];
rz(-1.1683162) q[3];
sx q[3];
rz(-2.4703006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93817389) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(0.43756884) q[2];
rz(2.3214052) q[3];
sx q[3];
rz(-1.5929675) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95847982) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(1.0108277) q[0];
rz(-0.084835947) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(2.5929677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6609782) q[0];
sx q[0];
rz(-1.3950431) q[0];
sx q[0];
rz(3.1094167) q[0];
x q[1];
rz(-0.43253501) q[2];
sx q[2];
rz(-1.2706332) q[2];
sx q[2];
rz(-1.6520239) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8884436) q[1];
sx q[1];
rz(-2.4427919) q[1];
sx q[1];
rz(-0.22380016) q[1];
rz(-2.2574378) q[3];
sx q[3];
rz(-0.13152371) q[3];
sx q[3];
rz(0.92708528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7944472) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(2.6084206) q[2];
rz(-2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086640373) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(2.3593498) q[0];
rz(0.070146322) q[1];
sx q[1];
rz(-0.47824305) q[1];
sx q[1];
rz(-0.22629647) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5180888) q[0];
sx q[0];
rz(-1.6395578) q[0];
sx q[0];
rz(-3.070773) q[0];
x q[1];
rz(-1.3391206) q[2];
sx q[2];
rz(-0.69139987) q[2];
sx q[2];
rz(-1.1065799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.037776064) q[1];
sx q[1];
rz(-1.0352219) q[1];
sx q[1];
rz(1.6984254) q[1];
rz(-pi) q[2];
rz(-1.6582011) q[3];
sx q[3];
rz(-2.8897396) q[3];
sx q[3];
rz(2.9853068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0349064) q[2];
sx q[2];
rz(-0.28768134) q[2];
sx q[2];
rz(1.7314343) q[2];
rz(-2.8790224) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866289) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(2.9578399) q[0];
rz(-2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(-0.62350887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79270169) q[0];
sx q[0];
rz(-2.0306132) q[0];
sx q[0];
rz(-1.2280653) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5207768) q[2];
sx q[2];
rz(-1.0631655) q[2];
sx q[2];
rz(1.0342645) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.19621721) q[1];
sx q[1];
rz(-1.5935285) q[1];
sx q[1];
rz(-3.0356221) q[1];
rz(-1.9454089) q[3];
sx q[3];
rz(-1.5419898) q[3];
sx q[3];
rz(-2.0458986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3925675) q[2];
sx q[2];
rz(-0.6124658) q[2];
sx q[2];
rz(0.037671063) q[2];
rz(2.7231349) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(-0.17879626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2035718) q[0];
sx q[0];
rz(-0.9104712) q[0];
sx q[0];
rz(-0.18785432) q[0];
rz(0.45166311) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(-0.53244932) q[2];
sx q[2];
rz(-1.3146853) q[2];
sx q[2];
rz(2.5934861) q[2];
rz(-2.5463811) q[3];
sx q[3];
rz(-0.32929904) q[3];
sx q[3];
rz(2.3198447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5476721) q[0];
sx q[0];
rz(-0.5991109) q[0];
sx q[0];
rz(-0.17679086) q[0];
rz(2.0556567) q[1];
sx q[1];
rz(-2.63201) q[1];
sx q[1];
rz(-0.086960763) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3172042) q[0];
sx q[0];
rz(-1.9476711) q[0];
sx q[0];
rz(0.61610846) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8463507) q[2];
sx q[2];
rz(-0.59165162) q[2];
sx q[2];
rz(-1.0585143) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0217758) q[1];
sx q[1];
rz(-1.6476742) q[1];
sx q[1];
rz(1.4426115) q[1];
rz(-pi) q[2];
rz(1.8636892) q[3];
sx q[3];
rz(-0.97145069) q[3];
sx q[3];
rz(1.0224316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1842492) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(-2.5498665) q[2];
rz(-2.7033778) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(0.87619585) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1205207) q[0];
sx q[0];
rz(-0.35728917) q[0];
sx q[0];
rz(-2.0508118) q[0];
rz(2.4202994) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(0.24478197) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1274968) q[0];
sx q[0];
rz(-1.9056221) q[0];
sx q[0];
rz(2.0398519) q[0];
rz(-0.82142395) q[2];
sx q[2];
rz(-1.9000179) q[2];
sx q[2];
rz(3.1067348) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87585014) q[1];
sx q[1];
rz(-2.2835554) q[1];
sx q[1];
rz(-2.4044988) q[1];
rz(-0.29676389) q[3];
sx q[3];
rz(-1.4302084) q[3];
sx q[3];
rz(2.1641017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.81405866) q[2];
sx q[2];
rz(-0.56698292) q[2];
sx q[2];
rz(0.84612334) q[2];
rz(-0.36492473) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(0.97682166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038079809) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(0.14347759) q[0];
rz(1.4682651) q[1];
sx q[1];
rz(-1.9915308) q[1];
sx q[1];
rz(0.066468261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398427) q[0];
sx q[0];
rz(-2.2494613) q[0];
sx q[0];
rz(1.6854273) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8861553) q[2];
sx q[2];
rz(-1.0424926) q[2];
sx q[2];
rz(2.262399) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4926949) q[1];
sx q[1];
rz(-1.4543415) q[1];
sx q[1];
rz(1.4728738) q[1];
x q[2];
rz(0.74031728) q[3];
sx q[3];
rz(-1.6941287) q[3];
sx q[3];
rz(-0.75238673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0299783) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(0.31271333) q[2];
rz(0.2615658) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(-0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12544352) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(3.0545767) q[0];
rz(-1.7866987) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(0.048390128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780086) q[0];
sx q[0];
rz(-2.3391294) q[0];
sx q[0];
rz(-1.7491231) q[0];
x q[1];
rz(1.3905647) q[2];
sx q[2];
rz(-1.7186621) q[2];
sx q[2];
rz(0.027005349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.794487) q[1];
sx q[1];
rz(-2.2655267) q[1];
sx q[1];
rz(-1.6415651) q[1];
x q[2];
rz(0.6137621) q[3];
sx q[3];
rz(-2.1832972) q[3];
sx q[3];
rz(-0.70741725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4309569) q[2];
sx q[2];
rz(-2.841876) q[2];
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
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61442536) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(-0.72702485) q[0];
rz(1.4432888) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(1.7194933) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073339065) q[0];
sx q[0];
rz(-1.710482) q[0];
sx q[0];
rz(-0.54648593) q[0];
x q[1];
rz(-0.70962064) q[2];
sx q[2];
rz(-2.8292252) q[2];
sx q[2];
rz(-1.720495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15502587) q[1];
sx q[1];
rz(-1.4490713) q[1];
sx q[1];
rz(-1.4847859) q[1];
rz(-0.57797076) q[3];
sx q[3];
rz(-0.9909361) q[3];
sx q[3];
rz(-0.65366369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(-1.5606073) q[2];
rz(1.8033146) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(0.54792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018589858) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(0.71989584) q[0];
rz(-2.9010991) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(-2.8954411) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6841078) q[0];
sx q[0];
rz(-0.2817758) q[0];
sx q[0];
rz(3.0578567) q[0];
rz(-pi) q[1];
rz(1.9754378) q[2];
sx q[2];
rz(-0.85524594) q[2];
sx q[2];
rz(3.0457979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5009821) q[1];
sx q[1];
rz(-2.3965008) q[1];
sx q[1];
rz(0.19265811) q[1];
rz(-1.9897377) q[3];
sx q[3];
rz(-2.2139858) q[3];
sx q[3];
rz(-2.1758428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4973732) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(-0.79968828) q[2];
rz(-0.52404809) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2972357) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(2.2204087) q[0];
rz(-1.4211897) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(-0.99501077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9807905) q[0];
sx q[0];
rz(-1.9483999) q[0];
sx q[0];
rz(2.9016205) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.424355) q[2];
sx q[2];
rz(-2.7173923) q[2];
sx q[2];
rz(2.7643124) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0459208) q[1];
sx q[1];
rz(-2.0883882) q[1];
sx q[1];
rz(1.0817097) q[1];
rz(-pi) q[2];
rz(2.7217676) q[3];
sx q[3];
rz(-1.9732765) q[3];
sx q[3];
rz(2.4703006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2034188) q[2];
sx q[2];
rz(-2.9439681) q[2];
sx q[2];
rz(2.7040238) q[2];
rz(-2.3214052) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95847982) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(-2.130765) q[0];
rz(3.0567567) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(2.5929677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6609782) q[0];
sx q[0];
rz(-1.7465495) q[0];
sx q[0];
rz(0.032175933) q[0];
rz(-pi) q[1];
rz(-0.43253501) q[2];
sx q[2];
rz(-1.2706332) q[2];
sx q[2];
rz(-1.6520239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1774166) q[1];
sx q[1];
rz(-2.2488222) q[1];
sx q[1];
rz(1.3864338) q[1];
rz(-pi) q[2];
rz(0.88415481) q[3];
sx q[3];
rz(-0.13152371) q[3];
sx q[3];
rz(-2.2145074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7944472) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(0.53317201) q[2];
rz(-2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086640373) q[0];
sx q[0];
rz(-2.7821879) q[0];
sx q[0];
rz(-2.3593498) q[0];
rz(-3.0714463) q[1];
sx q[1];
rz(-0.47824305) q[1];
sx q[1];
rz(2.9152962) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057581456) q[0];
sx q[0];
rz(-1.5001443) q[0];
sx q[0];
rz(1.5018626) q[0];
rz(-pi) q[1];
rz(0.89266291) q[2];
sx q[2];
rz(-1.423866) q[2];
sx q[2];
rz(-0.28444296) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4676241) q[1];
sx q[1];
rz(-1.6804763) q[1];
sx q[1];
rz(-2.6024271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3198648) q[3];
sx q[3];
rz(-1.5925515) q[3];
sx q[3];
rz(-1.8117429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(-2.8790224) q[3];
sx q[3];
rz(-1.5841443) q[3];
sx q[3];
rz(3.0098651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054963741) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(2.9578399) q[0];
rz(0.19206583) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(0.62350887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.469891) q[0];
sx q[0];
rz(-2.5755223) q[0];
sx q[0];
rz(-2.5453687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97098668) q[2];
sx q[2];
rz(-2.1040593) q[2];
sx q[2];
rz(-2.2704934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5564881) q[1];
sx q[1];
rz(-0.10837238) q[1];
sx q[1];
rz(-0.21173112) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1961837) q[3];
sx q[3];
rz(-1.5419898) q[3];
sx q[3];
rz(-1.0956941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7490251) q[2];
sx q[2];
rz(-0.6124658) q[2];
sx q[2];
rz(-0.037671063) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(-2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9380209) q[0];
sx q[0];
rz(-0.9104712) q[0];
sx q[0];
rz(-0.18785432) q[0];
rz(-0.45166311) q[1];
sx q[1];
rz(-1.3126806) q[1];
sx q[1];
rz(-1.5246593) q[1];
rz(-0.47623604) q[2];
sx q[2];
rz(-0.58544896) q[2];
sx q[2];
rz(0.61665012) q[2];
rz(0.27576294) q[3];
sx q[3];
rz(-1.388474) q[3];
sx q[3];
rz(0.17920517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

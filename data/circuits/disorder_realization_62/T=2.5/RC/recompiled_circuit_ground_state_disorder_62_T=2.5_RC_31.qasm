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
rz(2.9648018) q[0];
rz(2.0556567) q[1];
sx q[1];
rz(3.6511753) q[1];
sx q[1];
rz(9.3378172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8243884) q[0];
sx q[0];
rz(-1.9476711) q[0];
sx q[0];
rz(-2.5254842) q[0];
rz(-pi) q[1];
rz(-2.5702417) q[2];
sx q[2];
rz(-1.4077912) q[2];
sx q[2];
rz(-0.75955078) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0217758) q[1];
sx q[1];
rz(-1.6476742) q[1];
sx q[1];
rz(-1.6989811) q[1];
x q[2];
rz(-0.39985184) q[3];
sx q[3];
rz(-2.4824871) q[3];
sx q[3];
rz(2.6100998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.95734346) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(2.5498665) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(2.2653968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1205207) q[0];
sx q[0];
rz(-0.35728917) q[0];
sx q[0];
rz(-1.0907809) q[0];
rz(2.4202994) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(0.24478197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1274968) q[0];
sx q[0];
rz(-1.2359706) q[0];
sx q[0];
rz(-1.1017407) q[0];
rz(-pi) q[1];
rz(-2.0356947) q[2];
sx q[2];
rz(-2.3361578) q[2];
sx q[2];
rz(1.9400846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.822387) q[1];
sx q[1];
rz(-0.97619769) q[1];
sx q[1];
rz(0.90984224) q[1];
x q[2];
rz(-1.7177204) q[3];
sx q[3];
rz(-1.2770481) q[3];
sx q[3];
rz(-0.6361286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81405866) q[2];
sx q[2];
rz(-0.56698292) q[2];
sx q[2];
rz(2.2954693) q[2];
rz(0.36492473) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(-0.97682166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1035128) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(-1.4682651) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-3.0751244) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9003948) q[0];
sx q[0];
rz(-1.6599492) q[0];
sx q[0];
rz(-2.4597079) q[0];
rz(-2.1135751) q[2];
sx q[2];
rz(-1.7907871) q[2];
sx q[2];
rz(2.5808711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.089515226) q[1];
sx q[1];
rz(-1.4735392) q[1];
sx q[1];
rz(-3.0245824) q[1];
x q[2];
rz(0.74031728) q[3];
sx q[3];
rz(-1.4474639) q[3];
sx q[3];
rz(0.75238673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0299783) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(-0.31271333) q[2];
rz(-0.2615658) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(-2.3436782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12544352) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(0.087015986) q[0];
rz(-1.7866987) q[1];
sx q[1];
rz(-1.5785297) q[1];
sx q[1];
rz(-0.048390128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78260471) q[0];
sx q[0];
rz(-1.6986956) q[0];
sx q[0];
rz(0.77632287) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87746967) q[2];
sx q[2];
rz(-2.9089768) q[2];
sx q[2];
rz(-2.2777429) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.4573867) q[1];
sx q[1];
rz(-0.69772875) q[1];
sx q[1];
rz(3.0569424) q[1];
rz(-0.88416962) q[3];
sx q[3];
rz(-0.83809747) q[3];
sx q[3];
rz(-2.9626486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71063572) q[2];
sx q[2];
rz(-0.29971665) q[2];
sx q[2];
rz(-2.7002913) q[2];
rz(2.273061) q[3];
sx q[3];
rz(-1.3683616) q[3];
sx q[3];
rz(1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61442536) q[0];
sx q[0];
rz(-1.704957) q[0];
sx q[0];
rz(-0.72702485) q[0];
rz(-1.4432888) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(-1.7194933) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2724394) q[0];
sx q[0];
rz(-0.56229385) q[0];
sx q[0];
rz(2.8773688) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3634113) q[2];
sx q[2];
rz(-1.8060914) q[2];
sx q[2];
rz(-2.4547142) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4053012) q[1];
sx q[1];
rz(-1.6561688) q[1];
sx q[1];
rz(-0.12217223) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2660013) q[3];
sx q[3];
rz(-2.3471222) q[3];
sx q[3];
rz(-0.21847413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21165851) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(-1.5809853) q[2];
rz(1.3382781) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(-0.54792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(3.1230028) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(0.71989584) q[0];
rz(-2.9010991) q[1];
sx q[1];
rz(-1.980282) q[1];
sx q[1];
rz(-0.24615157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7712646) q[0];
sx q[0];
rz(-1.2900347) q[0];
sx q[0];
rz(1.5465897) q[0];
x q[1];
rz(1.9754378) q[2];
sx q[2];
rz(-2.2863467) q[2];
sx q[2];
rz(-3.0457979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5009821) q[1];
sx q[1];
rz(-0.74509186) q[1];
sx q[1];
rz(0.19265811) q[1];
x q[2];
rz(-2.4544958) q[3];
sx q[3];
rz(-1.9023484) q[3];
sx q[3];
rz(2.2755663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64421946) q[2];
sx q[2];
rz(-2.5770598) q[2];
sx q[2];
rz(0.79968828) q[2];
rz(2.6175446) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(-0.016949765) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2972357) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(0.921184) q[0];
rz(-1.720403) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(-0.99501077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1608022) q[0];
sx q[0];
rz(-1.1931927) q[0];
sx q[0];
rz(2.9016205) q[0];
x q[1];
rz(-0.065804577) q[2];
sx q[2];
rz(-1.9901681) q[2];
sx q[2];
rz(-0.21682993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.095671818) q[1];
sx q[1];
rz(-2.0883882) q[1];
sx q[1];
rz(1.0817097) q[1];
x q[2];
rz(-1.1345484) q[3];
sx q[3];
rz(-1.9552257) q[3];
sx q[3];
rz(-2.4151797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2034188) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(-0.43756884) q[2];
rz(-2.3214052) q[3];
sx q[3];
rz(-1.5929675) q[3];
sx q[3];
rz(-3.0062413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95847982) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(-1.0108277) q[0];
rz(0.084835947) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(0.54862499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6626604) q[0];
sx q[0];
rz(-2.9629483) q[0];
sx q[0];
rz(1.3915855) q[0];
x q[1];
rz(0.63603129) q[2];
sx q[2];
rz(-2.6205728) q[2];
sx q[2];
rz(-0.65083671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2531491) q[1];
sx q[1];
rz(-0.69880077) q[1];
sx q[1];
rz(2.9177925) q[1];
rz(-pi) q[2];
rz(2.2574378) q[3];
sx q[3];
rz(-0.13152371) q[3];
sx q[3];
rz(-0.92708528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34714547) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(-2.6084206) q[2];
rz(2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(-2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0549523) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(-2.3593498) q[0];
rz(-3.0714463) q[1];
sx q[1];
rz(-2.6633496) q[1];
sx q[1];
rz(-2.9152962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5180888) q[0];
sx q[0];
rz(-1.5020348) q[0];
sx q[0];
rz(3.070773) q[0];
rz(-0.18780577) q[2];
sx q[2];
rz(-0.90131288) q[2];
sx q[2];
rz(-1.4037496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1038166) q[1];
sx q[1];
rz(-2.1063707) q[1];
sx q[1];
rz(1.4431672) q[1];
x q[2];
rz(0.022458301) q[3];
sx q[3];
rz(-1.3199255) q[3];
sx q[3];
rz(-2.89507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0349064) q[2];
sx q[2];
rz(-0.28768134) q[2];
sx q[2];
rz(-1.4101583) q[2];
rz(2.8790224) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(-0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866289) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(-2.9578399) q[0];
rz(2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(0.62350887) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62105084) q[0];
sx q[0];
rz(-1.8767002) q[0];
sx q[0];
rz(2.6575179) q[0];
x q[1];
rz(-2.170606) q[2];
sx q[2];
rz(-1.0375334) q[2];
sx q[2];
rz(2.2704934) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9453754) q[1];
sx q[1];
rz(-1.5480642) q[1];
sx q[1];
rz(-3.0356221) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6493813) q[3];
sx q[3];
rz(-0.37566638) q[3];
sx q[3];
rz(-0.54822719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(0.037671063) q[2];
rz(2.7231349) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.2035718) q[0];
sx q[0];
rz(-2.2311214) q[0];
sx q[0];
rz(2.9537383) q[0];
rz(0.45166311) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(0.47623604) q[2];
sx q[2];
rz(-2.5561437) q[2];
sx q[2];
rz(-2.5249425) q[2];
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

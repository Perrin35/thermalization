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
rz(3.6511753) q[1];
sx q[1];
rz(9.3378172) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6330948) q[0];
sx q[0];
rz(-1.0035536) q[0];
sx q[0];
rz(-1.1192516) q[0];
rz(-pi) q[1];
rz(-2.8463507) q[2];
sx q[2];
rz(-2.549941) q[2];
sx q[2];
rz(-1.0585143) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.086395212) q[1];
sx q[1];
rz(-2.9922303) q[1];
sx q[1];
rz(-1.0284852) q[1];
rz(-2.7417408) q[3];
sx q[3];
rz(-0.65910554) q[3];
sx q[3];
rz(2.6100998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1842492) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(-0.59172612) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(-0.87619585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0210719) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(2.0508118) q[0];
rz(-0.72129321) q[1];
sx q[1];
rz(-1.6585766) q[1];
sx q[1];
rz(-0.24478197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8633085) q[0];
sx q[0];
rz(-2.0119036) q[0];
sx q[0];
rz(2.7696904) q[0];
x q[1];
rz(-0.4366283) q[2];
sx q[2];
rz(-2.2712913) q[2];
sx q[2];
rz(1.8281405) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3192057) q[1];
sx q[1];
rz(-2.165395) q[1];
sx q[1];
rz(0.90984224) q[1];
rz(1.7177204) q[3];
sx q[3];
rz(-1.8645446) q[3];
sx q[3];
rz(-0.6361286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.327534) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(0.84612334) q[2];
rz(-2.7766679) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(2.164771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038079809) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(1.4682651) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-0.066468261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9003948) q[0];
sx q[0];
rz(-1.4816435) q[0];
sx q[0];
rz(0.68188473) q[0];
rz(1.9793545) q[2];
sx q[2];
rz(-2.5600932) q[2];
sx q[2];
rz(1.7844019) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79083645) q[1];
sx q[1];
rz(-2.9895824) q[1];
sx q[1];
rz(-0.69614567) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18174882) q[3];
sx q[3];
rz(-2.3929993) q[3];
sx q[3];
rz(-2.4570217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1116144) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(-2.8288793) q[2];
rz(0.2615658) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0161491) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(0.087015986) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(3.0932025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66358405) q[0];
sx q[0];
rz(-2.3391294) q[0];
sx q[0];
rz(-1.7491231) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87746967) q[2];
sx q[2];
rz(-2.9089768) q[2];
sx q[2];
rz(0.86384976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.684206) q[1];
sx q[1];
rz(-0.69772875) q[1];
sx q[1];
rz(-0.084650234) q[1];
rz(0.86080024) q[3];
sx q[3];
rz(-1.080092) q[3];
sx q[3];
rz(-1.2482289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71063572) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(-2.7002913) q[2];
rz(-2.273061) q[3];
sx q[3];
rz(-1.3683616) q[3];
sx q[3];
rz(-1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.8836421) q[1];
sx q[1];
rz(-1.4220994) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0682536) q[0];
sx q[0];
rz(-1.710482) q[0];
sx q[0];
rz(2.5951067) q[0];
x q[1];
rz(2.431972) q[2];
sx q[2];
rz(-0.31236744) q[2];
sx q[2];
rz(1.720495) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9865668) q[1];
sx q[1];
rz(-1.6925214) q[1];
sx q[1];
rz(1.6568068) q[1];
rz(-2.234455) q[3];
sx q[3];
rz(-2.0454413) q[3];
sx q[3];
rz(-1.8812219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(1.5606073) q[2];
rz(1.8033146) q[3];
sx q[3];
rz(-0.18733297) q[3];
sx q[3];
rz(-0.54792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018589858) q[0];
sx q[0];
rz(-0.81273166) q[0];
sx q[0];
rz(-2.4216968) q[0];
rz(-2.9010991) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(0.24615157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3703281) q[0];
sx q[0];
rz(-1.2900347) q[0];
sx q[0];
rz(1.595003) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9754378) q[2];
sx q[2];
rz(-2.2863467) q[2];
sx q[2];
rz(0.095794769) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0689905) q[1];
sx q[1];
rz(-1.4406057) q[1];
sx q[1];
rz(2.4058002) q[1];
x q[2];
rz(-0.68709685) q[3];
sx q[3];
rz(-1.2392443) q[3];
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
rz(-2.3419044) q[2];
rz(-0.52404809) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(3.1246429) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84435695) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(2.2204087) q[0];
rz(1.4211897) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(-0.99501077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667082) q[0];
sx q[0];
rz(-0.44429438) q[0];
sx q[0];
rz(-2.1106476) q[0];
x q[1];
rz(-1.7172377) q[2];
sx q[2];
rz(-2.7173923) q[2];
sx q[2];
rz(-2.7643124) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.095671818) q[1];
sx q[1];
rz(-2.0883882) q[1];
sx q[1];
rz(2.059883) q[1];
rz(-2.7217676) q[3];
sx q[3];
rz(-1.9732765) q[3];
sx q[3];
rz(0.67129204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93817389) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(0.43756884) q[2];
rz(2.3214052) q[3];
sx q[3];
rz(-1.5929675) q[3];
sx q[3];
rz(3.0062413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.95847982) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(-2.130765) q[0];
rz(3.0567567) q[1];
sx q[1];
rz(-1.9761706) q[1];
sx q[1];
rz(-0.54862499) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6609782) q[0];
sx q[0];
rz(-1.3950431) q[0];
sx q[0];
rz(-0.032175933) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2422416) q[2];
sx q[2];
rz(-1.9827843) q[2];
sx q[2];
rz(0.054445353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2531491) q[1];
sx q[1];
rz(-0.69880077) q[1];
sx q[1];
rz(0.22380016) q[1];
x q[2];
rz(3.0579257) q[3];
sx q[3];
rz(-1.4691969) q[3];
sx q[3];
rz(-1.6179832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.34714547) q[2];
sx q[2];
rz(-2.562279) q[2];
sx q[2];
rz(-0.53317201) q[2];
rz(2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(0.50895154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086640373) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(-0.78224283) q[0];
rz(-3.0714463) q[1];
sx q[1];
rz(-2.6633496) q[1];
sx q[1];
rz(0.22629647) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6235038) q[0];
sx q[0];
rz(-1.5020348) q[0];
sx q[0];
rz(-3.070773) q[0];
rz(-pi) q[1];
rz(1.3391206) q[2];
sx q[2];
rz(-2.4501928) q[2];
sx q[2];
rz(-1.1065799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8574684) q[1];
sx q[1];
rz(-0.54912607) q[1];
sx q[1];
rz(2.930307) q[1];
rz(1.6582011) q[3];
sx q[3];
rz(-2.8897396) q[3];
sx q[3];
rz(-2.9853068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1066863) q[2];
sx q[2];
rz(-0.28768134) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(2.8790224) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(3.0098651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866289) q[0];
sx q[0];
rz(-0.2121191) q[0];
sx q[0];
rz(0.18375272) q[0];
rz(0.19206583) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(0.62350887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348891) q[0];
sx q[0];
rz(-2.0306132) q[0];
sx q[0];
rz(1.9135273) q[0];
rz(-2.170606) q[2];
sx q[2];
rz(-1.0375334) q[2];
sx q[2];
rz(-0.87109921) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9453754) q[1];
sx q[1];
rz(-1.5480642) q[1];
sx q[1];
rz(-0.10597056) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6493813) q[3];
sx q[3];
rz(-2.7659263) q[3];
sx q[3];
rz(2.5933655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(-3.1039216) q[2];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.6091433) q[2];
sx q[2];
rz(-1.8269074) q[2];
sx q[2];
rz(-0.54810654) q[2];
rz(1.3814817) q[3];
sx q[3];
rz(-1.2997205) q[3];
sx q[3];
rz(-1.442853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

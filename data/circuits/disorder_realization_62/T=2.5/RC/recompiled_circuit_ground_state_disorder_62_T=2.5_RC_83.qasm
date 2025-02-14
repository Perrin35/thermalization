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
rz(3.7407036) q[0];
sx q[0];
rz(9.6015688) q[0];
rz(-1.085936) q[1];
sx q[1];
rz(-0.50958264) q[1];
sx q[1];
rz(0.086960763) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22573839) q[0];
sx q[0];
rz(-2.432352) q[0];
sx q[0];
rz(-2.5410557) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3777138) q[2];
sx q[2];
rz(-2.1336485) q[2];
sx q[2];
rz(-0.70729296) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0551974) q[1];
sx q[1];
rz(-0.1493624) q[1];
sx q[1];
rz(-2.1131074) q[1];
x q[2];
rz(1.2779034) q[3];
sx q[3];
rz(-0.97145069) q[3];
sx q[3];
rz(2.1191611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1842492) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(0.59172612) q[2];
rz(2.7033778) q[3];
sx q[3];
rz(-0.44132909) q[3];
sx q[3];
rz(-2.2653968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.1205207) q[0];
sx q[0];
rz(-0.35728917) q[0];
sx q[0];
rz(2.0508118) q[0];
rz(2.4202994) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(-2.8968107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0185623) q[0];
sx q[0];
rz(-2.5726312) q[0];
sx q[0];
rz(0.91482343) q[0];
rz(1.1058979) q[2];
sx q[2];
rz(-0.80543488) q[2];
sx q[2];
rz(-1.9400846) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9823581) q[1];
sx q[1];
rz(-2.1042542) q[1];
sx q[1];
rz(0.70833556) q[1];
rz(-pi) q[2];
rz(1.4238722) q[3];
sx q[3];
rz(-1.8645446) q[3];
sx q[3];
rz(0.6361286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(2.2954693) q[2];
rz(2.7766679) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(-2.164771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1035128) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(0.14347759) q[0];
rz(-1.4682651) q[1];
sx q[1];
rz(-1.9915308) q[1];
sx q[1];
rz(3.0751244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9212355) q[0];
sx q[0];
rz(-2.4548303) q[0];
sx q[0];
rz(-3.0007017) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0280175) q[2];
sx q[2];
rz(-1.3508056) q[2];
sx q[2];
rz(0.56072158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0520774) q[1];
sx q[1];
rz(-1.4735392) q[1];
sx q[1];
rz(-0.11701028) q[1];
rz(-1.4044365) q[3];
sx q[3];
rz(-2.3041953) q[3];
sx q[3];
rz(2.211253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0299783) q[2];
sx q[2];
rz(-2.6252803) q[2];
sx q[2];
rz(-2.8288793) q[2];
rz(-0.2615658) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(2.3436782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0161491) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(3.0545767) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5785297) q[1];
sx q[1];
rz(-3.0932025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91735578) q[0];
sx q[0];
rz(-0.78460562) q[0];
sx q[0];
rz(2.9600701) q[0];
rz(-2.9913285) q[2];
sx q[2];
rz(-1.74904) q[2];
sx q[2];
rz(-1.5706289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.4573867) q[1];
sx q[1];
rz(-2.4438639) q[1];
sx q[1];
rz(-0.084650234) q[1];
x q[2];
rz(2.5278306) q[3];
sx q[3];
rz(-2.1832972) q[3];
sx q[3];
rz(-2.4341754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71063572) q[2];
sx q[2];
rz(-0.29971665) q[2];
sx q[2];
rz(2.7002913) q[2];
rz(-2.273061) q[3];
sx q[3];
rz(-1.7732311) q[3];
sx q[3];
rz(1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271673) q[0];
sx q[0];
rz(-1.4366356) q[0];
sx q[0];
rz(-2.4145678) q[0];
rz(-1.6983039) q[1];
sx q[1];
rz(-1.8836421) q[1];
sx q[1];
rz(1.4220994) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0682536) q[0];
sx q[0];
rz(-1.4311106) q[0];
sx q[0];
rz(0.54648593) q[0];
x q[1];
rz(2.9013394) q[2];
sx q[2];
rz(-1.7723871) q[2];
sx q[2];
rz(-0.83490419) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77249563) q[1];
sx q[1];
rz(-2.9926692) q[1];
sx q[1];
rz(2.5293674) q[1];
rz(-pi) q[2];
rz(-2.2660013) q[3];
sx q[3];
rz(-0.79447047) q[3];
sx q[3];
rz(-0.21847413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21165851) q[2];
sx q[2];
rz(-0.59017605) q[2];
sx q[2];
rz(-1.5809853) q[2];
rz(-1.8033146) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3703281) q[0];
sx q[0];
rz(-1.2900347) q[0];
sx q[0];
rz(-1.5465897) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9754378) q[2];
sx q[2];
rz(-2.2863467) q[2];
sx q[2];
rz(0.095794769) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64061058) q[1];
sx q[1];
rz(-2.3965008) q[1];
sx q[1];
rz(0.19265811) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1518549) q[3];
sx q[3];
rz(-0.92760689) q[3];
sx q[3];
rz(2.1758428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4973732) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(-2.3419044) q[2];
rz(-2.6175446) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(0.016949765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2972357) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(0.921184) q[0];
rz(-1.4211897) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(-0.99501077) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5667082) q[0];
sx q[0];
rz(-2.6972983) q[0];
sx q[0];
rz(1.0309451) q[0];
x q[1];
rz(-1.424355) q[2];
sx q[2];
rz(-0.42420039) q[2];
sx q[2];
rz(-2.7643124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2240963) q[1];
sx q[1];
rz(-0.69643785) q[1];
sx q[1];
rz(-2.4516979) q[1];
rz(-pi) q[2];
rz(-0.41982502) q[3];
sx q[3];
rz(-1.9732765) q[3];
sx q[3];
rz(2.4703006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93817389) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(2.7040238) q[2];
rz(0.82018745) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.1831128) q[0];
sx q[0];
rz(-0.11976972) q[0];
sx q[0];
rz(2.130765) q[0];
rz(-3.0567567) q[1];
sx q[1];
rz(-1.1654221) q[1];
sx q[1];
rz(2.5929677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6609782) q[0];
sx q[0];
rz(-1.3950431) q[0];
sx q[0];
rz(0.032175933) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43253501) q[2];
sx q[2];
rz(-1.8709595) q[2];
sx q[2];
rz(-1.6520239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1774166) q[1];
sx q[1];
rz(-2.2488222) q[1];
sx q[1];
rz(-1.3864338) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0579257) q[3];
sx q[3];
rz(-1.4691969) q[3];
sx q[3];
rz(-1.6179832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7944472) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(0.53317201) q[2];
rz(2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(-2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4248764) q[0];
sx q[0];
rz(-0.098669395) q[0];
sx q[0];
rz(2.3697321) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18780577) q[2];
sx q[2];
rz(-0.90131288) q[2];
sx q[2];
rz(-1.737843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6739686) q[1];
sx q[1];
rz(-1.6804763) q[1];
sx q[1];
rz(0.53916559) q[1];
x q[2];
rz(1.8217279) q[3];
sx q[3];
rz(-1.5925515) q[3];
sx q[3];
rz(-1.3298498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.7314343) q[2];
rz(-0.26257026) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(3.0098651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866289) q[0];
sx q[0];
rz(-0.2121191) q[0];
sx q[0];
rz(-2.9578399) q[0];
rz(2.9495268) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(-2.5180838) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.469891) q[0];
sx q[0];
rz(-0.56607038) q[0];
sx q[0];
rz(0.59622391) q[0];
rz(-pi) q[1];
rz(-0.76304014) q[2];
sx q[2];
rz(-2.3614778) q[2];
sx q[2];
rz(0.060465079) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3721613) q[1];
sx q[1];
rz(-1.4648533) q[1];
sx q[1];
rz(-1.547936) q[1];
x q[2];
rz(-3.1106408) q[3];
sx q[3];
rz(-1.9452458) q[3];
sx q[3];
rz(0.46377814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(-0.037671063) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(0.17879626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2035718) q[0];
sx q[0];
rz(-2.2311214) q[0];
sx q[0];
rz(2.9537383) q[0];
rz(2.6899295) q[1];
sx q[1];
rz(-1.3126806) q[1];
sx q[1];
rz(-1.5246593) q[1];
rz(2.6653566) q[2];
sx q[2];
rz(-0.58544896) q[2];
sx q[2];
rz(0.61665012) q[2];
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

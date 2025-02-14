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
rz(-2.63201) q[1];
sx q[1];
rz(-0.086960763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8243884) q[0];
sx q[0];
rz(-1.9476711) q[0];
sx q[0];
rz(-2.5254842) q[0];
rz(-pi) q[1];
rz(-1.7638788) q[2];
sx q[2];
rz(-2.1336485) q[2];
sx q[2];
rz(2.4342997) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46087825) q[1];
sx q[1];
rz(-1.4429922) q[1];
sx q[1];
rz(-3.0640814) q[1];
rz(-pi) q[2];
rz(2.7417408) q[3];
sx q[3];
rz(-0.65910554) q[3];
sx q[3];
rz(-2.6100998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.95734346) q[2];
sx q[2];
rz(-0.49746305) q[2];
sx q[2];
rz(2.5498665) q[2];
rz(-2.7033778) q[3];
sx q[3];
rz(-2.7002636) q[3];
sx q[3];
rz(0.87619585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1205207) q[0];
sx q[0];
rz(-2.7843035) q[0];
sx q[0];
rz(1.0907809) q[0];
rz(2.4202994) q[1];
sx q[1];
rz(-1.6585766) q[1];
sx q[1];
rz(-0.24478197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0140959) q[0];
sx q[0];
rz(-1.9056221) q[0];
sx q[0];
rz(-2.0398519) q[0];
x q[1];
rz(-2.0356947) q[2];
sx q[2];
rz(-0.80543488) q[2];
sx q[2];
rz(-1.9400846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2657425) q[1];
sx q[1];
rz(-0.85803723) q[1];
sx q[1];
rz(2.4044988) q[1];
x q[2];
rz(0.45072933) q[3];
sx q[3];
rz(-0.32748947) q[3];
sx q[3];
rz(2.9779676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.327534) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(-0.84612334) q[2];
rz(0.36492473) q[3];
sx q[3];
rz(-2.7155184) q[3];
sx q[3];
rz(-2.164771) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038079809) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(2.9981151) q[0];
rz(-1.6733276) q[1];
sx q[1];
rz(-1.9915308) q[1];
sx q[1];
rz(-3.0751244) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398427) q[0];
sx q[0];
rz(-0.89213138) q[0];
sx q[0];
rz(1.6854273) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25543737) q[2];
sx q[2];
rz(-1.0424926) q[2];
sx q[2];
rz(0.87919368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.089515226) q[1];
sx q[1];
rz(-1.4735392) q[1];
sx q[1];
rz(-0.11701028) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7371561) q[3];
sx q[3];
rz(-2.3041953) q[3];
sx q[3];
rz(0.93033965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0299783) q[2];
sx q[2];
rz(-0.51631236) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0161491) q[0];
sx q[0];
rz(-2.8479072) q[0];
sx q[0];
rz(-3.0545767) q[0];
rz(-1.3548939) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(-0.048390128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91735578) q[0];
sx q[0];
rz(-2.356987) q[0];
sx q[0];
rz(-2.9600701) q[0];
rz(-pi) q[1];
rz(-1.751028) q[2];
sx q[2];
rz(-1.7186621) q[2];
sx q[2];
rz(0.027005349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1783414) q[1];
sx q[1];
rz(-1.6251441) q[1];
sx q[1];
rz(-2.4456294) q[1];
x q[2];
rz(2.5278306) q[3];
sx q[3];
rz(-2.1832972) q[3];
sx q[3];
rz(-2.4341754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71063572) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(0.44130138) q[2];
rz(0.86853164) q[3];
sx q[3];
rz(-1.7732311) q[3];
sx q[3];
rz(1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.5271673) q[0];
sx q[0];
rz(-1.704957) q[0];
sx q[0];
rz(-0.72702485) q[0];
rz(1.4432888) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(1.7194933) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.559645) q[0];
sx q[0];
rz(-2.1113681) q[0];
sx q[0];
rz(1.7339043) q[0];
rz(1.3634113) q[2];
sx q[2];
rz(-1.3355012) q[2];
sx q[2];
rz(0.6868785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.369097) q[1];
sx q[1];
rz(-0.14892347) q[1];
sx q[1];
rz(-2.5293674) q[1];
rz(0.90713769) q[3];
sx q[3];
rz(-2.0454413) q[3];
sx q[3];
rz(1.2603708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9299341) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(-1.5809853) q[2];
rz(1.8033146) q[3];
sx q[3];
rz(-0.18733297) q[3];
sx q[3];
rz(-0.54792255) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
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
rz(-1.980282) q[1];
sx q[1];
rz(-0.24615157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7712646) q[0];
sx q[0];
rz(-1.2900347) q[0];
sx q[0];
rz(1.5465897) q[0];
rz(-1.1661548) q[2];
sx q[2];
rz(-2.2863467) q[2];
sx q[2];
rz(-3.0457979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64061058) q[1];
sx q[1];
rz(-0.74509186) q[1];
sx q[1];
rz(-0.19265811) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68709685) q[3];
sx q[3];
rz(-1.9023484) q[3];
sx q[3];
rz(-0.86602634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4973732) q[2];
sx q[2];
rz(-2.5770598) q[2];
sx q[2];
rz(0.79968828) q[2];
rz(2.6175446) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(3.1246429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84435695) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(0.921184) q[0];
rz(-1.720403) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(-0.99501077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5667082) q[0];
sx q[0];
rz(-0.44429438) q[0];
sx q[0];
rz(-1.0309451) q[0];
rz(-pi) q[1];
rz(-1.7172377) q[2];
sx q[2];
rz(-2.7173923) q[2];
sx q[2];
rz(0.37728024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2240963) q[1];
sx q[1];
rz(-0.69643785) q[1];
sx q[1];
rz(2.4516979) q[1];
rz(-2.7217676) q[3];
sx q[3];
rz(-1.1683162) q[3];
sx q[3];
rz(-0.67129204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93817389) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.1831128) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(1.0108277) q[0];
rz(-0.084835947) q[1];
sx q[1];
rz(-1.1654221) q[1];
sx q[1];
rz(-2.5929677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47893229) q[0];
sx q[0];
rz(-0.17864431) q[0];
sx q[0];
rz(-1.7500072) q[0];
x q[1];
rz(-1.2422416) q[2];
sx q[2];
rz(-1.1588084) q[2];
sx q[2];
rz(-3.0871473) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8884436) q[1];
sx q[1];
rz(-0.69880077) q[1];
sx q[1];
rz(-0.22380016) q[1];
x q[2];
rz(-1.4688427) q[3];
sx q[3];
rz(-1.6540308) q[3];
sx q[3];
rz(0.038681313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34714547) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(-0.53317201) q[2];
rz(-2.063607) q[3];
sx q[3];
rz(-2.2205133) q[3];
sx q[3];
rz(-2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086640373) q[0];
sx q[0];
rz(-2.7821879) q[0];
sx q[0];
rz(-2.3593498) q[0];
rz(3.0714463) q[1];
sx q[1];
rz(-2.6633496) q[1];
sx q[1];
rz(-0.22629647) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0840112) q[0];
sx q[0];
rz(-1.6414483) q[0];
sx q[0];
rz(1.6397301) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89266291) q[2];
sx q[2];
rz(-1.7177267) q[2];
sx q[2];
rz(0.28444296) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8574684) q[1];
sx q[1];
rz(-0.54912607) q[1];
sx q[1];
rz(-0.21128564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3198648) q[3];
sx q[3];
rz(-1.5490412) q[3];
sx q[3];
rz(1.8117429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.4101583) q[2];
rz(-0.26257026) q[3];
sx q[3];
rz(-1.5841443) q[3];
sx q[3];
rz(0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054963741) q[0];
sx q[0];
rz(-0.2121191) q[0];
sx q[0];
rz(-0.18375272) q[0];
rz(-2.9495268) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(0.62350887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62105084) q[0];
sx q[0];
rz(-1.8767002) q[0];
sx q[0];
rz(0.48407475) q[0];
x q[1];
rz(-0.76304014) q[2];
sx q[2];
rz(-0.78011489) q[2];
sx q[2];
rz(-0.060465079) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5564881) q[1];
sx q[1];
rz(-3.0332203) q[1];
sx q[1];
rz(-2.9298615) q[1];
rz(-1.6493813) q[3];
sx q[3];
rz(-2.7659263) q[3];
sx q[3];
rz(-2.5933655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3925675) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(-0.037671063) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-0.28067121) q[3];
sx q[3];
rz(-0.17879626) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9380209) q[0];
sx q[0];
rz(-2.2311214) q[0];
sx q[0];
rz(2.9537383) q[0];
rz(-2.6899295) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(-1.2757318) q[2];
sx q[2];
rz(-1.0574592) q[2];
sx q[2];
rz(-1.9707373) q[2];
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

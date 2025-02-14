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
rz(-1.085936) q[1];
sx q[1];
rz(-0.50958264) q[1];
sx q[1];
rz(0.086960763) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6330948) q[0];
sx q[0];
rz(-2.1380391) q[0];
sx q[0];
rz(2.022341) q[0];
x q[1];
rz(-0.29524191) q[2];
sx q[2];
rz(-2.549941) q[2];
sx q[2];
rz(1.0585143) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1198169) q[1];
sx q[1];
rz(-1.6476742) q[1];
sx q[1];
rz(-1.6989811) q[1];
x q[2];
rz(-2.5218202) q[3];
sx q[3];
rz(-1.3300782) q[3];
sx q[3];
rz(0.71686577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1842492) q[2];
sx q[2];
rz(-0.49746305) q[2];
sx q[2];
rz(0.59172612) q[2];
rz(-0.43821487) q[3];
sx q[3];
rz(-0.44132909) q[3];
sx q[3];
rz(-2.2653968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.483016) q[1];
sx q[1];
rz(-2.8968107) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2782841) q[0];
sx q[0];
rz(-1.1296891) q[0];
sx q[0];
rz(0.37190227) q[0];
x q[1];
rz(0.82142395) q[2];
sx q[2];
rz(-1.2415747) q[2];
sx q[2];
rz(3.1067348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87585014) q[1];
sx q[1];
rz(-2.2835554) q[1];
sx q[1];
rz(-0.73709388) q[1];
rz(-1.4238722) q[3];
sx q[3];
rz(-1.8645446) q[3];
sx q[3];
rz(-0.6361286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.327534) q[2];
sx q[2];
rz(-0.56698292) q[2];
sx q[2];
rz(-0.84612334) q[2];
rz(-2.7766679) q[3];
sx q[3];
rz(-2.7155184) q[3];
sx q[3];
rz(0.97682166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1035128) q[0];
sx q[0];
rz(-2.3338023) q[0];
sx q[0];
rz(2.9981151) q[0];
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
rz(-0.22035711) q[0];
sx q[0];
rz(-2.4548303) q[0];
sx q[0];
rz(3.0007017) q[0];
x q[1];
rz(-1.9793545) q[2];
sx q[2];
rz(-0.58149946) q[2];
sx q[2];
rz(1.7844019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6488978) q[1];
sx q[1];
rz(-1.6872511) q[1];
sx q[1];
rz(1.4728738) q[1];
x q[2];
rz(2.9598438) q[3];
sx q[3];
rz(-0.74859339) q[3];
sx q[3];
rz(-2.4570217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1116144) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(2.8288793) q[2];
rz(0.2615658) q[3];
sx q[3];
rz(-1.692619) q[3];
sx q[3];
rz(-2.3436782) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0161491) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(3.0545767) q[0];
rz(1.3548939) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(0.048390128) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91735578) q[0];
sx q[0];
rz(-2.356987) q[0];
sx q[0];
rz(-2.9600701) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87746967) q[2];
sx q[2];
rz(-2.9089768) q[2];
sx q[2];
rz(2.2777429) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.794487) q[1];
sx q[1];
rz(-0.87606591) q[1];
sx q[1];
rz(-1.5000276) q[1];
rz(-pi) q[2];
rz(0.88416962) q[3];
sx q[3];
rz(-2.3034952) q[3];
sx q[3];
rz(0.17894408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71063572) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(0.44130138) q[2];
rz(2.273061) q[3];
sx q[3];
rz(-1.3683616) q[3];
sx q[3];
rz(1.3335479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
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
rz(0.72702485) q[0];
rz(-1.4432888) q[1];
sx q[1];
rz(-1.8836421) q[1];
sx q[1];
rz(-1.4220994) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8691532) q[0];
sx q[0];
rz(-2.5792988) q[0];
sx q[0];
rz(0.2642239) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70962064) q[2];
sx q[2];
rz(-0.31236744) q[2];
sx q[2];
rz(-1.4210977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.369097) q[1];
sx q[1];
rz(-2.9926692) q[1];
sx q[1];
rz(2.5293674) q[1];
rz(-0.90713769) q[3];
sx q[3];
rz(-1.0961514) q[3];
sx q[3];
rz(-1.8812219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9299341) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(1.5809853) q[2];
rz(1.8033146) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(-2.5936701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018589858) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(2.4216968) q[0];
rz(0.24049354) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(-2.8954411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45748487) q[0];
sx q[0];
rz(-0.2817758) q[0];
sx q[0];
rz(-3.0578567) q[0];
rz(-0.42527898) q[2];
sx q[2];
rz(-0.80406791) q[2];
sx q[2];
rz(-0.67415392) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5009821) q[1];
sx q[1];
rz(-0.74509186) q[1];
sx q[1];
rz(-0.19265811) q[1];
x q[2];
rz(1.1518549) q[3];
sx q[3];
rz(-2.2139858) q[3];
sx q[3];
rz(0.96574984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4973732) q[2];
sx q[2];
rz(-0.56453288) q[2];
sx q[2];
rz(-2.3419044) q[2];
rz(-0.52404809) q[3];
sx q[3];
rz(-2.7553813) q[3];
sx q[3];
rz(0.016949765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2972357) q[0];
sx q[0];
rz(-0.083426282) q[0];
sx q[0];
rz(0.921184) q[0];
rz(-1.4211897) q[1];
sx q[1];
rz(-0.68266308) q[1];
sx q[1];
rz(-2.1465819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667082) q[0];
sx q[0];
rz(-2.6972983) q[0];
sx q[0];
rz(-1.0309451) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.065804577) q[2];
sx q[2];
rz(-1.9901681) q[2];
sx q[2];
rz(-0.21682993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2176358) q[1];
sx q[1];
rz(-1.1502277) q[1];
sx q[1];
rz(0.57284184) q[1];
rz(-pi) q[2];
rz(-2.0070442) q[3];
sx q[3];
rz(-1.1863669) q[3];
sx q[3];
rz(-2.4151797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2034188) q[2];
sx q[2];
rz(-2.9439681) q[2];
sx q[2];
rz(2.7040238) q[2];
rz(0.82018745) q[3];
sx q[3];
rz(-1.5929675) q[3];
sx q[3];
rz(-3.0062413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.54862499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4806145) q[0];
sx q[0];
rz(-1.7465495) q[0];
sx q[0];
rz(0.032175933) q[0];
rz(1.2422416) q[2];
sx q[2];
rz(-1.9827843) q[2];
sx q[2];
rz(-3.0871473) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6514142) q[1];
sx q[1];
rz(-1.7140577) q[1];
sx q[1];
rz(0.68639042) q[1];
rz(1.6727499) q[3];
sx q[3];
rz(-1.4875618) q[3];
sx q[3];
rz(-0.038681313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34714547) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(-0.53317201) q[2];
rz(2.063607) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(-2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.6633496) q[1];
sx q[1];
rz(-2.9152962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057581456) q[0];
sx q[0];
rz(-1.6414483) q[0];
sx q[0];
rz(1.5018626) q[0];
rz(-pi) q[1];
rz(-2.2489297) q[2];
sx q[2];
rz(-1.423866) q[2];
sx q[2];
rz(-0.28444296) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28412425) q[1];
sx q[1];
rz(-0.54912607) q[1];
sx q[1];
rz(-2.930307) q[1];
rz(-1.4833916) q[3];
sx q[3];
rz(-0.25185302) q[3];
sx q[3];
rz(2.9853068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(-1.4101583) q[2];
rz(-2.8790224) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054963741) q[0];
sx q[0];
rz(-0.2121191) q[0];
sx q[0];
rz(-0.18375272) q[0];
rz(-0.19206583) q[1];
sx q[1];
rz(-1.686325) q[1];
sx q[1];
rz(2.5180838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5205418) q[0];
sx q[0];
rz(-1.8767002) q[0];
sx q[0];
rz(2.6575179) q[0];
rz(-pi) q[1];
rz(-0.76304014) q[2];
sx q[2];
rz(-2.3614778) q[2];
sx q[2];
rz(0.060465079) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5564881) q[1];
sx q[1];
rz(-3.0332203) q[1];
sx q[1];
rz(2.9298615) q[1];
rz(-pi) q[2];
rz(1.9454089) q[3];
sx q[3];
rz(-1.5419898) q[3];
sx q[3];
rz(-1.0956941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7490251) q[2];
sx q[2];
rz(-2.5291269) q[2];
sx q[2];
rz(-3.1039216) q[2];
rz(-2.7231349) q[3];
sx q[3];
rz(-2.8609214) q[3];
sx q[3];
rz(-2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9380209) q[0];
sx q[0];
rz(-0.9104712) q[0];
sx q[0];
rz(-0.18785432) q[0];
rz(0.45166311) q[1];
sx q[1];
rz(-1.8289121) q[1];
sx q[1];
rz(1.6169333) q[1];
rz(-1.2757318) q[2];
sx q[2];
rz(-1.0574592) q[2];
sx q[2];
rz(-1.9707373) q[2];
rz(2.5463811) q[3];
sx q[3];
rz(-2.8122936) q[3];
sx q[3];
rz(-0.82174792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

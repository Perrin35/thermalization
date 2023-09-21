OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(2.8770652) q[0];
sx q[0];
rz(9.8192083) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.278468) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(1.6186884) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94304396) q[2];
sx q[2];
rz(-2.606751) q[2];
sx q[2];
rz(-0.56378555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4760121) q[1];
sx q[1];
rz(-1.1886547) q[1];
sx q[1];
rz(-1.0784472) q[1];
x q[2];
rz(0.37813152) q[3];
sx q[3];
rz(-1.6117125) q[3];
sx q[3];
rz(0.10842987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(2.8519894) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-2.7976024) q[0];
rz(-3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.3551691) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8665168) q[0];
sx q[0];
rz(-1.6356902) q[0];
sx q[0];
rz(2.1823723) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6927035) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(0.09300692) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.720571) q[1];
sx q[1];
rz(-0.45383006) q[1];
sx q[1];
rz(1.717091) q[1];
x q[2];
rz(1.387272) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(0.21218382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.4780592) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0707866) q[0];
sx q[0];
rz(-0.35937989) q[0];
sx q[0];
rz(1.552836) q[0];
x q[1];
rz(2.5550585) q[2];
sx q[2];
rz(-0.9622935) q[2];
sx q[2];
rz(1.3003132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0703735) q[1];
sx q[1];
rz(-0.95446903) q[1];
sx q[1];
rz(3.0973869) q[1];
x q[2];
rz(2.8912656) q[3];
sx q[3];
rz(-0.86391376) q[3];
sx q[3];
rz(1.3644497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-0.24969077) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-0.011118523) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3664368) q[0];
sx q[0];
rz(-1.9754793) q[0];
sx q[0];
rz(2.0648271) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3035994) q[2];
sx q[2];
rz(-1.0655155) q[2];
sx q[2];
rz(-0.2620286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8324082) q[1];
sx q[1];
rz(-2.5384181) q[1];
sx q[1];
rz(0.96997728) q[1];
x q[2];
rz(-2.6075881) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(-0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(-0.49452531) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.3269075) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715695) q[0];
sx q[0];
rz(-1.4172535) q[0];
sx q[0];
rz(0.75385401) q[0];
rz(-pi) q[1];
rz(0.7977428) q[2];
sx q[2];
rz(-1.1290871) q[2];
sx q[2];
rz(1.9517348) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2111645) q[1];
sx q[1];
rz(-1.2439562) q[1];
sx q[1];
rz(1.8153166) q[1];
rz(0.29019659) q[3];
sx q[3];
rz(-1.1786596) q[3];
sx q[3];
rz(-0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(1.2456606) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(2.025827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67416699) q[0];
sx q[0];
rz(-0.90303991) q[0];
sx q[0];
rz(0.52533124) q[0];
x q[1];
rz(0.43004604) q[2];
sx q[2];
rz(-0.94656813) q[2];
sx q[2];
rz(-0.66832322) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9876154) q[1];
sx q[1];
rz(-0.9824285) q[1];
sx q[1];
rz(1.5920559) q[1];
rz(0.66949087) q[3];
sx q[3];
rz(-1.7117501) q[3];
sx q[3];
rz(1.6657176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(1.1122423) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(2.5792714) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.651197) q[0];
sx q[0];
rz(-2.2389452) q[0];
sx q[0];
rz(-3.1227123) q[0];
x q[1];
rz(-1.6717031) q[2];
sx q[2];
rz(-1.2404053) q[2];
sx q[2];
rz(2.908387) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6866236) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(1.7186233) q[1];
x q[2];
rz(-0.087248487) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(2.0283386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(0.39644077) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.6202392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0884468) q[0];
sx q[0];
rz(-1.7674812) q[0];
sx q[0];
rz(1.6403857) q[0];
x q[1];
rz(1.0993768) q[2];
sx q[2];
rz(-1.2087012) q[2];
sx q[2];
rz(-0.21381703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.13008598) q[1];
sx q[1];
rz(-0.79400051) q[1];
sx q[1];
rz(-1.6542875) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11407125) q[3];
sx q[3];
rz(-0.89469203) q[3];
sx q[3];
rz(-1.8166208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-0.29433027) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(2.8709581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198869) q[0];
sx q[0];
rz(-1.6019551) q[0];
sx q[0];
rz(3.1057538) q[0];
rz(-pi) q[1];
rz(-1.0893906) q[2];
sx q[2];
rz(-2.319616) q[2];
sx q[2];
rz(-0.35441986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86079019) q[1];
sx q[1];
rz(-1.6165015) q[1];
sx q[1];
rz(1.6153107) q[1];
rz(2.4828033) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(2.7487019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(0.45544004) q[2];
rz(-2.3296302) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(-2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-2.4023138) q[0];
rz(2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(0.49490067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56787581) q[0];
sx q[0];
rz(-1.084869) q[0];
sx q[0];
rz(3.0987415) q[0];
rz(-pi) q[1];
rz(1.5636744) q[2];
sx q[2];
rz(-1.8986423) q[2];
sx q[2];
rz(0.42051007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6634076) q[1];
sx q[1];
rz(-2.2067398) q[1];
sx q[1];
rz(0.138476) q[1];
x q[2];
rz(-1.4570191) q[3];
sx q[3];
rz(-0.69996951) q[3];
sx q[3];
rz(-1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0080002) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(0.83256759) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-0.20559786) q[2];
sx q[2];
rz(-2.1046706) q[2];
sx q[2];
rz(-0.54866366) q[2];
rz(1.8264063) q[3];
sx q[3];
rz(-1.3127463) q[3];
sx q[3];
rz(-1.128872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

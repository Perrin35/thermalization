OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(4.5259024) q[0];
sx q[0];
rz(10.685267) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0215065) q[0];
sx q[0];
rz(-1.2354295) q[0];
sx q[0];
rz(1.0943227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14416868) q[2];
sx q[2];
rz(-1.2906133) q[2];
sx q[2];
rz(0.35137128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0337861) q[1];
sx q[1];
rz(-2.4596655) q[1];
sx q[1];
rz(-2.4297907) q[1];
x q[2];
rz(-1.190891) q[3];
sx q[3];
rz(-1.1067179) q[3];
sx q[3];
rz(-0.78117785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(-2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.686036) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3128132) q[0];
sx q[0];
rz(-1.4038741) q[0];
sx q[0];
rz(0.98915999) q[0];
x q[1];
rz(-0.37462072) q[2];
sx q[2];
rz(-1.6472367) q[2];
sx q[2];
rz(-0.74707109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62555186) q[1];
sx q[1];
rz(-1.5572773) q[1];
sx q[1];
rz(2.0936172) q[1];
rz(-pi) q[2];
rz(-3.0807642) q[3];
sx q[3];
rz(-2.7292477) q[3];
sx q[3];
rz(2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7130647) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-0.3119719) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.172539) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55484178) q[0];
sx q[0];
rz(-1.2625853) q[0];
sx q[0];
rz(-0.59535938) q[0];
rz(-1.8981947) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(-2.3936405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0769656) q[1];
sx q[1];
rz(-0.58276999) q[1];
sx q[1];
rz(1.1175734) q[1];
rz(-pi) q[2];
rz(-1.3200687) q[3];
sx q[3];
rz(-1.513895) q[3];
sx q[3];
rz(1.0292605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0671493) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(-0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610385) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(-3.0134841) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(2.6180843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2142221) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(-0.58332304) q[0];
rz(-0.24387118) q[2];
sx q[2];
rz(-0.3393617) q[2];
sx q[2];
rz(-1.3432168) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47141155) q[1];
sx q[1];
rz(-1.8243316) q[1];
sx q[1];
rz(-2.2815435) q[1];
rz(-pi) q[2];
rz(3.0047699) q[3];
sx q[3];
rz(-0.98383437) q[3];
sx q[3];
rz(-2.6329991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(2.5775487) q[2];
rz(0.28856746) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6600835) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-2.1496444) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49682) q[0];
sx q[0];
rz(-0.20475514) q[0];
sx q[0];
rz(0.55069189) q[0];
rz(1.8736585) q[2];
sx q[2];
rz(-1.5793243) q[2];
sx q[2];
rz(2.0451343) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79210287) q[1];
sx q[1];
rz(-1.2621242) q[1];
sx q[1];
rz(-1.5285138) q[1];
x q[2];
rz(1.4335853) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(1.2587794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(-0.22578421) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.7165002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544918) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(-1.8255193) q[0];
rz(3.0503057) q[2];
sx q[2];
rz(-1.316615) q[2];
sx q[2];
rz(-2.9448178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7879494) q[1];
sx q[1];
rz(-1.1218346) q[1];
sx q[1];
rz(3.06762) q[1];
x q[2];
rz(0.96958843) q[3];
sx q[3];
rz(-2.0242656) q[3];
sx q[3];
rz(-0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0075334) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(-2.5777204) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-0.82180506) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99188995) q[0];
sx q[0];
rz(-2.2737962) q[0];
sx q[0];
rz(2.4292612) q[0];
rz(2.2981811) q[2];
sx q[2];
rz(-1.9551829) q[2];
sx q[2];
rz(-0.20993983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51470876) q[1];
sx q[1];
rz(-0.85299546) q[1];
sx q[1];
rz(-1.2674598) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9415226) q[3];
sx q[3];
rz(-2.126174) q[3];
sx q[3];
rz(-2.6851418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3372779) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(2.896893) q[2];
rz(-0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.3051422) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96852036) q[0];
sx q[0];
rz(-0.99181306) q[0];
sx q[0];
rz(1.1770244) q[0];
rz(-pi) q[1];
rz(-1.201198) q[2];
sx q[2];
rz(-2.2410789) q[2];
sx q[2];
rz(-3.0073462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0933989) q[1];
sx q[1];
rz(-1.1953925) q[1];
sx q[1];
rz(-0.6742538) q[1];
x q[2];
rz(-0.023530258) q[3];
sx q[3];
rz(-1.2315893) q[3];
sx q[3];
rz(2.6587405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(1.0772609) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7548783) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(2.4386491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020333175) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(2.3938177) q[0];
rz(-0.039673294) q[2];
sx q[2];
rz(-2.0797605) q[2];
sx q[2];
rz(0.85658011) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1785537) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(-1.2600684) q[1];
rz(-pi) q[2];
rz(-1.1144936) q[3];
sx q[3];
rz(-0.66702402) q[3];
sx q[3];
rz(-1.7175355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(0.46863619) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43384957) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(1.3462523) q[0];
x q[1];
rz(0.87601985) q[2];
sx q[2];
rz(-0.48601905) q[2];
sx q[2];
rz(-1.8234058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84506449) q[1];
sx q[1];
rz(-1.5025286) q[1];
sx q[1];
rz(1.8370085) q[1];
rz(-2.7077984) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(-1.0292366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.998385) q[2];
rz(0.11463595) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(-1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951915) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(0.8016349) q[2];
sx q[2];
rz(-0.87052204) q[2];
sx q[2];
rz(-1.082765) q[2];
rz(2.2223496) q[3];
sx q[3];
rz(-1.5083434) q[3];
sx q[3];
rz(-1.2283243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

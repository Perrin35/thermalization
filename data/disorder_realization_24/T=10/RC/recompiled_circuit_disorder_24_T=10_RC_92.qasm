OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5819117) q[0];
sx q[0];
rz(-2.0547325) q[0];
sx q[0];
rz(1.342919) q[0];
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(2.4989541) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896048) q[0];
sx q[0];
rz(-1.5986686) q[0];
sx q[0];
rz(0.53919381) q[0];
rz(2.8060993) q[2];
sx q[2];
rz(-0.53288424) q[2];
sx q[2];
rz(-2.3032308) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5583615) q[1];
sx q[1];
rz(-1.4500631) q[1];
sx q[1];
rz(1.6912778) q[1];
rz(2.9075165) q[3];
sx q[3];
rz(-1.5353068) q[3];
sx q[3];
rz(-1.2897829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(-2.4228418) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(-0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46368018) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(-2.1610778) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(-0.78871361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0448488) q[0];
sx q[0];
rz(-1.6941438) q[0];
sx q[0];
rz(-0.11476536) q[0];
rz(-pi) q[1];
rz(2.8854495) q[2];
sx q[2];
rz(-1.5400585) q[2];
sx q[2];
rz(-2.5275633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9849557) q[1];
sx q[1];
rz(-0.64132323) q[1];
sx q[1];
rz(-0.5248458) q[1];
x q[2];
rz(2.2867145) q[3];
sx q[3];
rz(-0.93920556) q[3];
sx q[3];
rz(0.43254334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.366189) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(2.9555087) q[2];
rz(2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(-3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9300951) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(0.63013664) q[0];
rz(3.0139626) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(-0.72174597) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33357027) q[0];
sx q[0];
rz(-2.434242) q[0];
sx q[0];
rz(-2.7471514) q[0];
x q[1];
rz(2.0314991) q[2];
sx q[2];
rz(-0.93020541) q[2];
sx q[2];
rz(3.1415423) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0112146) q[1];
sx q[1];
rz(-1.5323258) q[1];
sx q[1];
rz(-0.4442261) q[1];
rz(-pi) q[2];
rz(-0.33629041) q[3];
sx q[3];
rz(-2.170917) q[3];
sx q[3];
rz(-2.1093413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7599941) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(2.2325113) q[2];
rz(2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(-2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5575314) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(0.042536143) q[0];
rz(2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(0.76400486) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1799058) q[0];
sx q[0];
rz(-0.67877239) q[0];
sx q[0];
rz(-1.5508482) q[0];
rz(-pi) q[1];
rz(0.64455428) q[2];
sx q[2];
rz(-1.4696215) q[2];
sx q[2];
rz(0.25601706) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2149787) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(0.39919969) q[1];
rz(-2.571991) q[3];
sx q[3];
rz(-0.74865018) q[3];
sx q[3];
rz(-2.1514055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2970695) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(0.26238966) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(-0.46245241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2247291) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(-1.2438783) q[0];
rz(-0.22661701) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(0.40333834) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4606844) q[0];
sx q[0];
rz(-1.7366689) q[0];
sx q[0];
rz(-0.18780639) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4558805) q[2];
sx q[2];
rz(-1.8572241) q[2];
sx q[2];
rz(0.098766947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4407318) q[1];
sx q[1];
rz(-1.7248389) q[1];
sx q[1];
rz(-0.53149077) q[1];
x q[2];
rz(2.1020528) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(1.2951375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8205745) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(2.3590951) q[2];
rz(-2.0292422) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(-1.0085683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.430442) q[0];
sx q[0];
rz(-2.7395881) q[0];
sx q[0];
rz(2.6859786) q[0];
rz(0.09952155) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(-3.1351556) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.54654) q[0];
sx q[0];
rz(-1.3188625) q[0];
sx q[0];
rz(-2.3810054) q[0];
x q[1];
rz(1.1941031) q[2];
sx q[2];
rz(-1.4720535) q[2];
sx q[2];
rz(1.8284947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7240119) q[1];
sx q[1];
rz(-1.15861) q[1];
sx q[1];
rz(-0.25556232) q[1];
x q[2];
rz(-0.75565831) q[3];
sx q[3];
rz(-2.14058) q[3];
sx q[3];
rz(-2.8017686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.5051179) q[2];
sx q[2];
rz(-2.2951365) q[2];
rz(-0.99772292) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098175123) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(3.085882) q[0];
rz(-0.78272351) q[1];
sx q[1];
rz(-2.0109773) q[1];
sx q[1];
rz(-1.9810716) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8502064) q[0];
sx q[0];
rz(-1.123748) q[0];
sx q[0];
rz(-0.38711754) q[0];
rz(-pi) q[1];
rz(-2.6703228) q[2];
sx q[2];
rz(-2.0435213) q[2];
sx q[2];
rz(-1.8007235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22735587) q[1];
sx q[1];
rz(-2.1319234) q[1];
sx q[1];
rz(2.6941872) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37961827) q[3];
sx q[3];
rz(-1.2109204) q[3];
sx q[3];
rz(0.66756638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.065585) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-0.78545061) q[2];
rz(-2.3857332) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(-2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128368) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(1.998741) q[0];
rz(-1.3061334) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(0.41608861) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04836719) q[0];
sx q[0];
rz(-1.0972293) q[0];
sx q[0];
rz(0.22305365) q[0];
rz(-0.14256723) q[2];
sx q[2];
rz(-2.0349742) q[2];
sx q[2];
rz(0.67205059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21129747) q[1];
sx q[1];
rz(-0.43356178) q[1];
sx q[1];
rz(-2.7525206) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3564524) q[3];
sx q[3];
rz(-2.6594901) q[3];
sx q[3];
rz(1.4516423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(0.32315928) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(-1.9966104) q[0];
rz(-1.1960944) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(-0.51913613) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0275092) q[0];
sx q[0];
rz(-2.7402096) q[0];
sx q[0];
rz(-1.9163577) q[0];
rz(-1.3349418) q[2];
sx q[2];
rz(-1.354885) q[2];
sx q[2];
rz(2.4049135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12431006) q[1];
sx q[1];
rz(-0.78622422) q[1];
sx q[1];
rz(3.0148274) q[1];
rz(1.395123) q[3];
sx q[3];
rz(-1.3364949) q[3];
sx q[3];
rz(2.7845886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(0.040977565) q[2];
rz(0.86769062) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(-1.7136259) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(-1.6428927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1536627) q[0];
sx q[0];
rz(-1.4953574) q[0];
sx q[0];
rz(3.1213785) q[0];
rz(-pi) q[1];
rz(0.047949009) q[2];
sx q[2];
rz(-1.7257294) q[2];
sx q[2];
rz(-0.63873728) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5012706) q[1];
sx q[1];
rz(-0.81201279) q[1];
sx q[1];
rz(-2.0444319) q[1];
rz(1.6937709) q[3];
sx q[3];
rz(-1.3696559) q[3];
sx q[3];
rz(0.5826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75858086) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(1.2659484) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(0.11309359) q[2];
sx q[2];
rz(-1.332915) q[2];
sx q[2];
rz(-1.1168196) q[2];
rz(-2.416715) q[3];
sx q[3];
rz(-2.1028044) q[3];
sx q[3];
rz(0.067306991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
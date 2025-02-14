OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71830463) q[0];
sx q[0];
rz(-0.77646065) q[0];
sx q[0];
rz(2.7756696) q[0];
rz(-2.794682) q[1];
sx q[1];
rz(3.4263098) q[1];
sx q[1];
rz(10.813536) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58908868) q[0];
sx q[0];
rz(-2.9285746) q[0];
sx q[0];
rz(2.8211843) q[0];
rz(1.8687145) q[2];
sx q[2];
rz(-0.90014825) q[2];
sx q[2];
rz(2.9811695) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60290135) q[1];
sx q[1];
rz(-2.5316892) q[1];
sx q[1];
rz(-2.2944488) q[1];
rz(-pi) q[2];
rz(1.0424869) q[3];
sx q[3];
rz(-1.0431223) q[3];
sx q[3];
rz(1.2880402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95691386) q[2];
sx q[2];
rz(-2.2214486) q[2];
sx q[2];
rz(0.93519768) q[2];
rz(2.8397371) q[3];
sx q[3];
rz(-1.5422834) q[3];
sx q[3];
rz(2.7479318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1644156) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(-2.6432977) q[0];
rz(1.4277108) q[1];
sx q[1];
rz(-1.9352813) q[1];
sx q[1];
rz(-2.0548342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55288494) q[0];
sx q[0];
rz(-2.4037139) q[0];
sx q[0];
rz(-0.35564519) q[0];
x q[1];
rz(-2.1939414) q[2];
sx q[2];
rz(-0.32440475) q[2];
sx q[2];
rz(1.0783006) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1732498) q[1];
sx q[1];
rz(-1.1361215) q[1];
sx q[1];
rz(-1.7477075) q[1];
x q[2];
rz(-0.03647904) q[3];
sx q[3];
rz(-0.31841296) q[3];
sx q[3];
rz(-2.1977001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9429417) q[2];
sx q[2];
rz(-0.31060878) q[2];
sx q[2];
rz(-1.3322213) q[2];
rz(-1.1474991) q[3];
sx q[3];
rz(-1.5787326) q[3];
sx q[3];
rz(-1.8089627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20951095) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(-0.15723666) q[0];
rz(1.2278185) q[1];
sx q[1];
rz(-2.9324052) q[1];
sx q[1];
rz(-0.42207178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7868766) q[0];
sx q[0];
rz(-1.7484957) q[0];
sx q[0];
rz(0.091376809) q[0];
rz(0.27760013) q[2];
sx q[2];
rz(-1.4034082) q[2];
sx q[2];
rz(0.041942747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2873309) q[1];
sx q[1];
rz(-2.7497083) q[1];
sx q[1];
rz(-1.3381097) q[1];
x q[2];
rz(1.8881599) q[3];
sx q[3];
rz(-2.3471213) q[3];
sx q[3];
rz(2.1411856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32855365) q[2];
sx q[2];
rz(-0.96031323) q[2];
sx q[2];
rz(2.4075244) q[2];
rz(-1.0775393) q[3];
sx q[3];
rz(-2.4534093) q[3];
sx q[3];
rz(-3.0663826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828736) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(0.016481312) q[0];
rz(-3.0670498) q[1];
sx q[1];
rz(-2.6838979) q[1];
sx q[1];
rz(-0.25513908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1265565) q[0];
sx q[0];
rz(-1.714596) q[0];
sx q[0];
rz(0.48499996) q[0];
rz(-pi) q[1];
rz(-1.0048149) q[2];
sx q[2];
rz(-1.2850956) q[2];
sx q[2];
rz(-0.20488747) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29683477) q[1];
sx q[1];
rz(-0.95176178) q[1];
sx q[1];
rz(2.1059787) q[1];
rz(0.11394088) q[3];
sx q[3];
rz(-1.9181532) q[3];
sx q[3];
rz(1.243227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5269346) q[2];
sx q[2];
rz(-2.4444828) q[2];
sx q[2];
rz(-0.69980168) q[2];
rz(0.091863306) q[3];
sx q[3];
rz(-1.2173165) q[3];
sx q[3];
rz(0.45472586) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62234539) q[0];
sx q[0];
rz(-0.82038251) q[0];
sx q[0];
rz(-1.1118332) q[0];
rz(2.9513997) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(-1.1598738) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96944189) q[0];
sx q[0];
rz(-1.5363201) q[0];
sx q[0];
rz(1.9146754) q[0];
x q[1];
rz(2.2091523) q[2];
sx q[2];
rz(-1.4529422) q[2];
sx q[2];
rz(1.7913851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76471699) q[1];
sx q[1];
rz(-2.5333166) q[1];
sx q[1];
rz(2.6981955) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0926863) q[3];
sx q[3];
rz(-1.5318222) q[3];
sx q[3];
rz(-3.1119973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5707034) q[2];
sx q[2];
rz(-2.3184226) q[2];
sx q[2];
rz(-2.9111351) q[2];
rz(-2.0795836) q[3];
sx q[3];
rz(-1.2758723) q[3];
sx q[3];
rz(-1.6331858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45288169) q[0];
sx q[0];
rz(-2.0827561) q[0];
sx q[0];
rz(1.3558615) q[0];
rz(-0.52976766) q[1];
sx q[1];
rz(-2.3593088) q[1];
sx q[1];
rz(-1.1518325) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760107) q[0];
sx q[0];
rz(-1.9261203) q[0];
sx q[0];
rz(1.4655345) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1124338) q[2];
sx q[2];
rz(-1.6525606) q[2];
sx q[2];
rz(2.0643016) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85688389) q[1];
sx q[1];
rz(-1.207106) q[1];
sx q[1];
rz(-1.3631352) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21422503) q[3];
sx q[3];
rz(-1.6856632) q[3];
sx q[3];
rz(1.8095995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4749703) q[2];
sx q[2];
rz(-0.64133659) q[2];
sx q[2];
rz(0.46710157) q[2];
rz(-2.1304255) q[3];
sx q[3];
rz(-0.34365383) q[3];
sx q[3];
rz(1.3694192) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8491299) q[0];
sx q[0];
rz(-2.9865773) q[0];
sx q[0];
rz(-2.9169061) q[0];
rz(2.977773) q[1];
sx q[1];
rz(-0.33912173) q[1];
sx q[1];
rz(-2.4952369) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689679) q[0];
sx q[0];
rz(-0.59212055) q[0];
sx q[0];
rz(2.0298241) q[0];
rz(-pi) q[1];
rz(0.21804131) q[2];
sx q[2];
rz(-1.2550809) q[2];
sx q[2];
rz(1.9670847) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9440981) q[1];
sx q[1];
rz(-1.0284327) q[1];
sx q[1];
rz(-1.4828862) q[1];
x q[2];
rz(1.7067753) q[3];
sx q[3];
rz(-1.1866015) q[3];
sx q[3];
rz(-2.3649491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15567638) q[2];
sx q[2];
rz(-1.6935657) q[2];
sx q[2];
rz(-2.7316366) q[2];
rz(-2.3112467) q[3];
sx q[3];
rz(-0.94873077) q[3];
sx q[3];
rz(1.693694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.603867) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(2.895288) q[0];
rz(2.314997) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(0.26396096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1136341) q[0];
sx q[0];
rz(-1.6662681) q[0];
sx q[0];
rz(-0.0072633538) q[0];
rz(2.0461778) q[2];
sx q[2];
rz(-0.84566294) q[2];
sx q[2];
rz(2.1755168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5689478) q[1];
sx q[1];
rz(-1.2897451) q[1];
sx q[1];
rz(0.87139178) q[1];
rz(-pi) q[2];
rz(2.554515) q[3];
sx q[3];
rz(-1.2012606) q[3];
sx q[3];
rz(1.0472421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.10494122) q[2];
sx q[2];
rz(-1.2678601) q[2];
sx q[2];
rz(2.4264753) q[2];
rz(-3.0702843) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(-2.3947072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0095373) q[0];
sx q[0];
rz(-1.6522464) q[0];
sx q[0];
rz(0.43553964) q[0];
rz(-2.9938193) q[1];
sx q[1];
rz(-1.7099893) q[1];
sx q[1];
rz(-1.6730283) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3776079) q[0];
sx q[0];
rz(-1.1398291) q[0];
sx q[0];
rz(1.3727643) q[0];
x q[1];
rz(-2.4442102) q[2];
sx q[2];
rz(-1.5101523) q[2];
sx q[2];
rz(1.9549119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.33959162) q[1];
sx q[1];
rz(-2.7164384) q[1];
sx q[1];
rz(-1.9849066) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7870258) q[3];
sx q[3];
rz(-1.9659967) q[3];
sx q[3];
rz(0.54602269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74464166) q[2];
sx q[2];
rz(-1.7640742) q[2];
sx q[2];
rz(-0.89293876) q[2];
rz(1.8630155) q[3];
sx q[3];
rz(-1.1568926) q[3];
sx q[3];
rz(-1.6781835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8466012) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(-2.5355329) q[0];
rz(3.1312969) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(0.079924718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91658501) q[0];
sx q[0];
rz(-1.3189303) q[0];
sx q[0];
rz(-2.1504137) q[0];
x q[1];
rz(-1.8907919) q[2];
sx q[2];
rz(-2.3872132) q[2];
sx q[2];
rz(-0.88448521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8680758) q[1];
sx q[1];
rz(-2.7339327) q[1];
sx q[1];
rz(0.056053921) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8013838) q[3];
sx q[3];
rz(-1.7713431) q[3];
sx q[3];
rz(1.5179364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0004878) q[2];
sx q[2];
rz(-2.2147369) q[2];
sx q[2];
rz(2.1263988) q[2];
rz(-0.60123932) q[3];
sx q[3];
rz(-2.4398118) q[3];
sx q[3];
rz(-2.0950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64423185) q[0];
sx q[0];
rz(-1.5174706) q[0];
sx q[0];
rz(0.68751412) q[0];
rz(-2.5524706) q[1];
sx q[1];
rz(-0.67710572) q[1];
sx q[1];
rz(2.6893375) q[1];
rz(-2.6964006) q[2];
sx q[2];
rz(-1.4839076) q[2];
sx q[2];
rz(-2.1640726) q[2];
rz(1.4429055) q[3];
sx q[3];
rz(-0.75027942) q[3];
sx q[3];
rz(1.1269509) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

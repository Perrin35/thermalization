OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(4.2545118) q[1];
sx q[1];
rz(1.7634044) q[1];
sx q[1];
rz(7.4982285) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62984798) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(1.3826136) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8589852) q[2];
sx q[2];
rz(-2.211314) q[2];
sx q[2];
rz(-0.033601947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79599586) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(0.3791581) q[1];
x q[2];
rz(-1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9248283) q[0];
sx q[0];
rz(-0.98470682) q[0];
sx q[0];
rz(3.0490962) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9039731) q[2];
sx q[2];
rz(-0.75287205) q[2];
sx q[2];
rz(1.876229) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7669249) q[1];
sx q[1];
rz(-0.51480773) q[1];
sx q[1];
rz(1.2109846) q[1];
rz(-0.48861309) q[3];
sx q[3];
rz(-2.6693137) q[3];
sx q[3];
rz(-0.74913914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-3.113134) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(-1.1921047) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(0.55535299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12268513) q[0];
sx q[0];
rz(-0.92120954) q[0];
sx q[0];
rz(-2.5193549) q[0];
x q[1];
rz(-2.136134) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(-0.29758673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66545031) q[1];
sx q[1];
rz(-0.94201554) q[1];
sx q[1];
rz(-0.20035845) q[1];
rz(-pi) q[2];
rz(-2.9402296) q[3];
sx q[3];
rz(-1.0346197) q[3];
sx q[3];
rz(0.45504967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(-2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92514738) q[0];
sx q[0];
rz(-1.2108004) q[0];
sx q[0];
rz(1.2544592) q[0];
x q[1];
rz(-0.63919477) q[2];
sx q[2];
rz(-1.8393469) q[2];
sx q[2];
rz(0.07721363) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7123588) q[1];
sx q[1];
rz(-2.7731967) q[1];
sx q[1];
rz(2.1894987) q[1];
x q[2];
rz(-1.6595483) q[3];
sx q[3];
rz(-0.52270652) q[3];
sx q[3];
rz(-0.50597092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(0.056093562) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984021) q[0];
sx q[0];
rz(-1.0023596) q[0];
sx q[0];
rz(2.0060904) q[0];
x q[1];
rz(2.202583) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(2.3914571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6009112) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(-1.5555698) q[1];
rz(1.7124743) q[3];
sx q[3];
rz(-2.5767527) q[3];
sx q[3];
rz(-0.39839881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1405979) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(-0.8941935) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(1.4661219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(2.1870959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4650824) q[0];
sx q[0];
rz(-2.0970779) q[0];
sx q[0];
rz(0.17980534) q[0];
rz(3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(0.45331732) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33659014) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(-2.5617983) q[1];
x q[2];
rz(-1.0682085) q[3];
sx q[3];
rz(-2.8121901) q[3];
sx q[3];
rz(-2.9941032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59297562) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(-2.3983811) q[0];
rz(1.6339533) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(-0.61002237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726495) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(-0.48604301) q[0];
rz(-pi) q[1];
rz(2.2248613) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(-2.9031861) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2513189) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(3.1097079) q[1];
x q[2];
rz(2.9355572) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(0.73658529) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026222762) q[0];
sx q[0];
rz(-1.6008458) q[0];
sx q[0];
rz(0.010859246) q[0];
rz(-pi) q[1];
rz(1.2543711) q[2];
sx q[2];
rz(-2.8036615) q[2];
sx q[2];
rz(-2.8527609) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8229586) q[1];
sx q[1];
rz(-1.4689323) q[1];
sx q[1];
rz(-2.676079) q[1];
rz(-pi) q[2];
rz(-0.96111091) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(1.5101658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(-2.4576808) q[2];
rz(1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(0.7448147) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72070044) q[0];
sx q[0];
rz(-1.2487131) q[0];
sx q[0];
rz(2.3548404) q[0];
x q[1];
rz(0.3955598) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(-1.2566483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71036584) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(2.9037895) q[1];
rz(-1.9635779) q[3];
sx q[3];
rz(-0.82735705) q[3];
sx q[3];
rz(2.3208997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334082) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(1.1709447) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084328018) q[0];
sx q[0];
rz(-1.143976) q[0];
sx q[0];
rz(1.5945934) q[0];
rz(-pi) q[1];
rz(-0.77565907) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(1.0722216) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78632894) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(-1.5934056) q[1];
rz(0.42697866) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(-0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(-0.56636089) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-0.70384937) q[2];
sx q[2];
rz(-2.2675632) q[2];
sx q[2];
rz(2.8013196) q[2];
rz(1.2373274) q[3];
sx q[3];
rz(-1.5860535) q[3];
sx q[3];
rz(0.64762583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

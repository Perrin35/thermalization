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
rz(2.365132) q[0];
sx q[0];
rz(9.790701) q[0];
rz(0.34691063) q[1];
sx q[1];
rz(-0.28471714) q[1];
sx q[1];
rz(-1.3887583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4735218) q[0];
sx q[0];
rz(-1.5041623) q[0];
sx q[0];
rz(0.20247831) q[0];
rz(-pi) q[1];
rz(-1.2728782) q[2];
sx q[2];
rz(-0.90014825) q[2];
sx q[2];
rz(2.9811695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5947111) q[1];
sx q[1];
rz(-1.9597907) q[1];
sx q[1];
rz(2.0531897) q[1];
rz(-pi) q[2];
rz(-2.5479814) q[3];
sx q[3];
rz(-1.1201914) q[3];
sx q[3];
rz(0.0030779989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95691386) q[2];
sx q[2];
rz(-0.92014402) q[2];
sx q[2];
rz(-0.93519768) q[2];
rz(-0.30185559) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(0.39366084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.97717706) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(2.6432977) q[0];
rz(-1.7138819) q[1];
sx q[1];
rz(-1.2063113) q[1];
sx q[1];
rz(2.0548342) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3918864) q[0];
sx q[0];
rz(-1.3343617) q[0];
sx q[0];
rz(0.70583418) q[0];
rz(-1.837376) q[2];
sx q[2];
rz(-1.757903) q[2];
sx q[2];
rz(-0.10548909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1732498) q[1];
sx q[1];
rz(-1.1361215) q[1];
sx q[1];
rz(-1.3938851) q[1];
rz(-pi) q[2];
x q[2];
rz(0.03647904) q[3];
sx q[3];
rz(-0.31841296) q[3];
sx q[3];
rz(-0.94389254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19865092) q[2];
sx q[2];
rz(-0.31060878) q[2];
sx q[2];
rz(-1.3322213) q[2];
rz(-1.9940935) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(1.3326299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20951095) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(0.15723666) q[0];
rz(-1.9137742) q[1];
sx q[1];
rz(-2.9324052) q[1];
sx q[1];
rz(2.7195209) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3738689) q[0];
sx q[0];
rz(-1.6607303) q[0];
sx q[0];
rz(-1.3923682) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8639925) q[2];
sx q[2];
rz(-1.4034082) q[2];
sx q[2];
rz(-3.0996499) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8542618) q[1];
sx q[1];
rz(-0.39188436) q[1];
sx q[1];
rz(1.803483) q[1];
rz(0.8019358) q[3];
sx q[3];
rz(-1.7953331) q[3];
sx q[3];
rz(2.3449986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32855365) q[2];
sx q[2];
rz(-0.96031323) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(1.0775393) q[3];
sx q[3];
rz(-2.4534093) q[3];
sx q[3];
rz(-0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587191) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(0.016481312) q[0];
rz(-3.0670498) q[1];
sx q[1];
rz(-2.6838979) q[1];
sx q[1];
rz(2.8864536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0150361) q[0];
sx q[0];
rz(-1.4269967) q[0];
sx q[0];
rz(2.6565927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33489494) q[2];
sx q[2];
rz(-1.0303632) q[2];
sx q[2];
rz(-1.9528567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8447579) q[1];
sx q[1];
rz(-2.1898309) q[1];
sx q[1];
rz(-1.0356139) q[1];
rz(-pi) q[2];
rz(-1.2665073) q[3];
sx q[3];
rz(-2.7767468) q[3];
sx q[3];
rz(-1.5675275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61465803) q[2];
sx q[2];
rz(-2.4444828) q[2];
sx q[2];
rz(0.69980168) q[2];
rz(3.0497293) q[3];
sx q[3];
rz(-1.2173165) q[3];
sx q[3];
rz(2.6868668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62234539) q[0];
sx q[0];
rz(-2.3212101) q[0];
sx q[0];
rz(2.0297594) q[0];
rz(-2.9513997) q[1];
sx q[1];
rz(-2.2564502) q[1];
sx q[1];
rz(-1.1598738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69731312) q[0];
sx q[0];
rz(-0.34553465) q[0];
sx q[0];
rz(1.6727425) q[0];
rz(0.93244035) q[2];
sx q[2];
rz(-1.6886504) q[2];
sx q[2];
rz(1.7913851) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3768757) q[1];
sx q[1];
rz(-0.60827601) q[1];
sx q[1];
rz(2.6981955) q[1];
rz(-1.0489063) q[3];
sx q[3];
rz(-1.6097704) q[3];
sx q[3];
rz(3.1119973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5707034) q[2];
sx q[2];
rz(-0.8231701) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(-1.0620091) q[3];
sx q[3];
rz(-1.8657203) q[3];
sx q[3];
rz(1.5084069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45288169) q[0];
sx q[0];
rz(-2.0827561) q[0];
sx q[0];
rz(1.7857312) q[0];
rz(-0.52976766) q[1];
sx q[1];
rz(-0.7822839) q[1];
sx q[1];
rz(1.1518325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6731176) q[0];
sx q[0];
rz(-1.4721319) q[0];
sx q[0];
rz(2.7844546) q[0];
rz(-1.0291589) q[2];
sx q[2];
rz(-1.6525606) q[2];
sx q[2];
rz(-2.0643016) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85688389) q[1];
sx q[1];
rz(-1.207106) q[1];
sx q[1];
rz(1.3631352) q[1];
rz(1.6883259) q[3];
sx q[3];
rz(-1.3580048) q[3];
sx q[3];
rz(0.26373395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6666224) q[2];
sx q[2];
rz(-0.64133659) q[2];
sx q[2];
rz(-2.6744911) q[2];
rz(1.0111672) q[3];
sx q[3];
rz(-0.34365383) q[3];
sx q[3];
rz(-1.7721734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2924627) q[0];
sx q[0];
rz(-0.15501538) q[0];
sx q[0];
rz(-2.9169061) q[0];
rz(-0.16381964) q[1];
sx q[1];
rz(-2.8024709) q[1];
sx q[1];
rz(-0.64635578) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3098329) q[0];
sx q[0];
rz(-1.0467967) q[0];
sx q[0];
rz(-2.8519467) q[0];
rz(-pi) q[1];
rz(1.8936526) q[2];
sx q[2];
rz(-1.7778998) q[2];
sx q[2];
rz(-2.8139909) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19749459) q[1];
sx q[1];
rz(-1.0284327) q[1];
sx q[1];
rz(1.4828862) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8180653) q[3];
sx q[3];
rz(-0.40641847) q[3];
sx q[3];
rz(0.42662963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15567638) q[2];
sx q[2];
rz(-1.4480269) q[2];
sx q[2];
rz(-0.4099561) q[2];
rz(-2.3112467) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.603867) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(-0.24630462) q[0];
rz(0.82659563) q[1];
sx q[1];
rz(-1.3984503) q[1];
sx q[1];
rz(0.26396096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1896828) q[0];
sx q[0];
rz(-3.0458458) q[0];
sx q[0];
rz(1.646498) q[0];
x q[1];
rz(0.7837495) q[2];
sx q[2];
rz(-1.2211868) q[2];
sx q[2];
rz(2.8658681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57264489) q[1];
sx q[1];
rz(-1.8518475) q[1];
sx q[1];
rz(2.2702009) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5313782) q[3];
sx q[3];
rz(-2.4596678) q[3];
sx q[3];
rz(-2.1206926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10494122) q[2];
sx q[2];
rz(-1.8737326) q[2];
sx q[2];
rz(2.4264753) q[2];
rz(3.0702843) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(2.3947072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(2.9938193) q[1];
sx q[1];
rz(-1.4316033) q[1];
sx q[1];
rz(1.4685644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298115) q[0];
sx q[0];
rz(-2.6699097) q[0];
sx q[0];
rz(-0.40431542) q[0];
rz(-pi) q[1];
rz(-0.69738241) q[2];
sx q[2];
rz(-1.6314403) q[2];
sx q[2];
rz(-1.1866807) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.802001) q[1];
sx q[1];
rz(-0.42515426) q[1];
sx q[1];
rz(-1.9849066) q[1];
rz(0.40364175) q[3];
sx q[3];
rz(-1.7701245) q[3];
sx q[3];
rz(0.9404054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74464166) q[2];
sx q[2];
rz(-1.7640742) q[2];
sx q[2];
rz(2.2486539) q[2];
rz(1.2785771) q[3];
sx q[3];
rz(-1.1568926) q[3];
sx q[3];
rz(1.6781835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29499149) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(2.5355329) q[0];
rz(-3.1312969) q[1];
sx q[1];
rz(-2.1538487) q[1];
sx q[1];
rz(0.079924718) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2250076) q[0];
sx q[0];
rz(-1.8226624) q[0];
sx q[0];
rz(0.99117898) q[0];
rz(-pi) q[1];
rz(-0.28744185) q[2];
sx q[2];
rz(-0.86311695) q[2];
sx q[2];
rz(-1.3112932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8680758) q[1];
sx q[1];
rz(-2.7339327) q[1];
sx q[1];
rz(-0.056053921) q[1];
rz(1.3402088) q[3];
sx q[3];
rz(-1.3702495) q[3];
sx q[3];
rz(-1.6236562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0004878) q[2];
sx q[2];
rz(-0.92685574) q[2];
sx q[2];
rz(-1.0151939) q[2];
rz(2.5403533) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(-1.0465485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4973608) q[0];
sx q[0];
rz(-1.5174706) q[0];
sx q[0];
rz(0.68751412) q[0];
rz(-0.58912206) q[1];
sx q[1];
rz(-2.4644869) q[1];
sx q[1];
rz(-0.45225515) q[1];
rz(-2.6964006) q[2];
sx q[2];
rz(-1.4839076) q[2];
sx q[2];
rz(-2.1640726) q[2];
rz(-0.82460578) q[3];
sx q[3];
rz(-1.6578703) q[3];
sx q[3];
rz(-0.35005611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

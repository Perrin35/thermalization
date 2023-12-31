OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(-1.0153271) q[0];
sx q[0];
rz(0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(-0.68038124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7079948) q[0];
sx q[0];
rz(-0.85277075) q[0];
sx q[0];
rz(-2.9684767) q[0];
rz(-1.8664076) q[2];
sx q[2];
rz(-1.3829074) q[2];
sx q[2];
rz(1.0653898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80402741) q[1];
sx q[1];
rz(-1.0684039) q[1];
sx q[1];
rz(-0.2663836) q[1];
x q[2];
rz(-2.8075571) q[3];
sx q[3];
rz(-1.4035657) q[3];
sx q[3];
rz(2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(-0.087466784) q[2];
rz(-2.4123689) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(-0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(2.9691866) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(-0.52406812) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755718) q[0];
sx q[0];
rz(-2.4040939) q[0];
sx q[0];
rz(-2.0185508) q[0];
rz(-pi) q[1];
rz(1.0753724) q[2];
sx q[2];
rz(-0.76250695) q[2];
sx q[2];
rz(-0.60053315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1986875) q[1];
sx q[1];
rz(-0.44829475) q[1];
sx q[1];
rz(0.56221902) q[1];
rz(1.4161795) q[3];
sx q[3];
rz(-2.1025476) q[3];
sx q[3];
rz(2.1173409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9841763) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(-1.8015507) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-3.1047344) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(2.5168193) q[0];
rz(-2.1773188) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-2.952081) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0441168) q[0];
sx q[0];
rz(-0.86320816) q[0];
sx q[0];
rz(1.9301027) q[0];
x q[1];
rz(1.6763716) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(2.2764652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.0057919766) q[1];
sx q[1];
rz(-1.6264919) q[1];
sx q[1];
rz(-2.9018351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1720445) q[3];
sx q[3];
rz(-1.4253972) q[3];
sx q[3];
rz(1.6325412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0345962) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(1.3872046) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4819734) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(1.5959928) q[0];
rz(-1.1384456) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(2.0514354) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7074993) q[0];
sx q[0];
rz(-0.543569) q[0];
sx q[0];
rz(-0.71732934) q[0];
x q[1];
rz(2.3218669) q[2];
sx q[2];
rz(-1.1851386) q[2];
sx q[2];
rz(-1.0939329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.412147) q[1];
sx q[1];
rz(-2.4033961) q[1];
sx q[1];
rz(-2.5022238) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7100536) q[3];
sx q[3];
rz(-2.6772237) q[3];
sx q[3];
rz(0.73016703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(1.1936197) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(-0.91317552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(2.9918616) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.9715462) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2802551) q[0];
sx q[0];
rz(-2.4147408) q[0];
sx q[0];
rz(1.0759541) q[0];
x q[1];
rz(-1.9636376) q[2];
sx q[2];
rz(-0.71937865) q[2];
sx q[2];
rz(-0.17578416) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3586732) q[1];
sx q[1];
rz(-1.9649319) q[1];
sx q[1];
rz(2.6637117) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7647469) q[3];
sx q[3];
rz(-1.0411106) q[3];
sx q[3];
rz(-2.291631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(0.13723792) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(-0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(0.67289105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56129365) q[0];
sx q[0];
rz(-1.467388) q[0];
sx q[0];
rz(-0.99594492) q[0];
x q[1];
rz(-0.024725155) q[2];
sx q[2];
rz(-2.284986) q[2];
sx q[2];
rz(1.6753472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91025483) q[1];
sx q[1];
rz(-0.64097039) q[1];
sx q[1];
rz(0.35229589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9665885) q[3];
sx q[3];
rz(-2.724218) q[3];
sx q[3];
rz(2.008703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(-0.43506452) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.172794) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(-2.0299714) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(1.7395082) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3081449) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(-1.863198) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76485302) q[2];
sx q[2];
rz(-1.3991038) q[2];
sx q[2];
rz(0.21966759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84343108) q[1];
sx q[1];
rz(-1.2097675) q[1];
sx q[1];
rz(-1.4069188) q[1];
rz(-pi) q[2];
rz(3.0038463) q[3];
sx q[3];
rz(-2.1546954) q[3];
sx q[3];
rz(-0.13970845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(2.7133572) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(-1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.5493786) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(-2.6002398) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4562948) q[0];
sx q[0];
rz(-1.0066427) q[0];
sx q[0];
rz(-1.5840522) q[0];
x q[1];
rz(-1.9114242) q[2];
sx q[2];
rz(-1.1414141) q[2];
sx q[2];
rz(0.18061772) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15335042) q[1];
sx q[1];
rz(-1.1594979) q[1];
sx q[1];
rz(2.3565355) q[1];
x q[2];
rz(1.5726611) q[3];
sx q[3];
rz(-1.9511173) q[3];
sx q[3];
rz(-2.5170381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(2.3642335) q[2];
rz(-0.85123953) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(-3.1316277) q[0];
rz(-2.126157) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.4452971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8953067) q[0];
sx q[0];
rz(-2.200897) q[0];
sx q[0];
rz(0.90755983) q[0];
rz(2.5601013) q[2];
sx q[2];
rz(-2.1610689) q[2];
sx q[2];
rz(-1.7828538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5053133) q[1];
sx q[1];
rz(-1.5863824) q[1];
sx q[1];
rz(2.511335) q[1];
rz(-pi) q[2];
rz(-0.61049283) q[3];
sx q[3];
rz(-0.99184147) q[3];
sx q[3];
rz(-2.0905153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69892591) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(0.6138531) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-0.39224958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948254) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(2.4368068) q[0];
x q[1];
rz(2.5160772) q[2];
sx q[2];
rz(-2.054347) q[2];
sx q[2];
rz(3.0160883) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85234648) q[1];
sx q[1];
rz(-2.1834071) q[1];
sx q[1];
rz(-0.5457408) q[1];
rz(-pi) q[2];
rz(-0.20499968) q[3];
sx q[3];
rz(-1.5991755) q[3];
sx q[3];
rz(1.320822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0739416) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7077211) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(2.1451163) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(-2.375013) q[2];
sx q[2];
rz(-1.3091514) q[2];
sx q[2];
rz(-1.8717629) q[2];
rz(2.2157833) q[3];
sx q[3];
rz(-2.3498597) q[3];
sx q[3];
rz(1.3510977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

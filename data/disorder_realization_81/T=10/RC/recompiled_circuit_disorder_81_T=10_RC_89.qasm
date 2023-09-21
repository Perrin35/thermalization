OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(4.9217304) q[0];
sx q[0];
rz(11.187727) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(4.6255914) q[1];
sx q[1];
rz(8.9738823) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8060018) q[0];
sx q[0];
rz(-1.5561034) q[0];
sx q[0];
rz(3.0528085) q[0];
rz(-pi) q[1];
rz(1.0916753) q[2];
sx q[2];
rz(-1.6699104) q[2];
sx q[2];
rz(-1.5168158) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2921819) q[1];
sx q[1];
rz(-1.0725478) q[1];
sx q[1];
rz(-1.8241747) q[1];
rz(0.095590683) q[3];
sx q[3];
rz(-0.88926892) q[3];
sx q[3];
rz(-2.3694627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(2.5355693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(2.8785008) q[0];
rz(0.94353765) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.9553604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3620136) q[0];
sx q[0];
rz(-1.5148666) q[0];
sx q[0];
rz(-1.5520142) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19940168) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(2.8859438) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4071719) q[1];
sx q[1];
rz(-2.2312575) q[1];
sx q[1];
rz(-2.9829626) q[1];
rz(1.7632742) q[3];
sx q[3];
rz(-0.88497439) q[3];
sx q[3];
rz(-1.6845077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(-0.31093591) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7611258) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(2.9280248) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504282) q[0];
sx q[0];
rz(-2.4463852) q[0];
sx q[0];
rz(-1.7466963) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0772466) q[2];
sx q[2];
rz(-1.764467) q[2];
sx q[2];
rz(2.1950302) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0606544) q[1];
sx q[1];
rz(-1.9229691) q[1];
sx q[1];
rz(-2.0209795) q[1];
rz(-pi) q[2];
rz(-1.3385653) q[3];
sx q[3];
rz(-2.665822) q[3];
sx q[3];
rz(2.8867718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(-0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61313066) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(-2.2303175) q[0];
rz(0.43831929) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(1.8211676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1962122) q[0];
sx q[0];
rz(-1.5604661) q[0];
sx q[0];
rz(1.9804079) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7337012) q[2];
sx q[2];
rz(-0.48569277) q[2];
sx q[2];
rz(2.8913468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1530694) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(0.70446976) q[1];
x q[2];
rz(-0.60453316) q[3];
sx q[3];
rz(-2.2570838) q[3];
sx q[3];
rz(-1.7741007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.8300366) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(-2.2763021) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.8409761) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0375992) q[0];
sx q[0];
rz(-2.7866057) q[0];
sx q[0];
rz(-0.07261891) q[0];
x q[1];
rz(-1.6426716) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(-0.14771151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2758267) q[1];
sx q[1];
rz(-1.2424801) q[1];
sx q[1];
rz(-2.5715716) q[1];
rz(-pi) q[2];
rz(1.424765) q[3];
sx q[3];
rz(-1.8304123) q[3];
sx q[3];
rz(-0.29823333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(-0.47719964) q[2];
rz(2.9495083) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(-0.55039644) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(-1.5884429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71205157) q[0];
sx q[0];
rz(-2.6459604) q[0];
sx q[0];
rz(-0.80722157) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1197238) q[2];
sx q[2];
rz(-1.936603) q[2];
sx q[2];
rz(-2.0621698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99332422) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(3.0118914) q[1];
rz(3.0670777) q[3];
sx q[3];
rz(-2.2722368) q[3];
sx q[3];
rz(0.45267347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(1.8590415) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5605374) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(0.67725956) q[0];
rz(2.989785) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(-0.97704926) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86977406) q[0];
sx q[0];
rz(-2.5442113) q[0];
sx q[0];
rz(0.13418829) q[0];
rz(-pi) q[1];
rz(-1.5708959) q[2];
sx q[2];
rz(-1.7028371) q[2];
sx q[2];
rz(3.0299203) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3358826) q[1];
sx q[1];
rz(-1.4937595) q[1];
sx q[1];
rz(1.9105934) q[1];
rz(-pi) q[2];
rz(-0.75120039) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(0.62869149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-2.4429328) q[0];
rz(-2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.2493856) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65803618) q[0];
sx q[0];
rz(-1.617096) q[0];
sx q[0];
rz(2.9306843) q[0];
rz(-pi) q[1];
rz(-2.2970389) q[2];
sx q[2];
rz(-0.65659467) q[2];
sx q[2];
rz(-0.086364634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.318995) q[1];
sx q[1];
rz(-1.5929475) q[1];
sx q[1];
rz(-1.5652565) q[1];
rz(-pi) q[2];
rz(0.93805712) q[3];
sx q[3];
rz(-1.3360099) q[3];
sx q[3];
rz(-0.42890047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(-1.9630986) q[2];
rz(-1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(-0.59757772) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1406527) q[0];
sx q[0];
rz(-1.6093328) q[0];
sx q[0];
rz(-1.4356104) q[0];
x q[1];
rz(-0.75002807) q[2];
sx q[2];
rz(-0.70493297) q[2];
sx q[2];
rz(-2.7761369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9012007) q[1];
sx q[1];
rz(-1.2482093) q[1];
sx q[1];
rz(-2.8278973) q[1];
x q[2];
rz(-1.647801) q[3];
sx q[3];
rz(-0.98620755) q[3];
sx q[3];
rz(-1.0494174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.2333599) q[2];
rz(-0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(-1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(0.67217174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32165747) q[0];
sx q[0];
rz(-1.5760839) q[0];
sx q[0];
rz(-1.5810285) q[0];
x q[1];
rz(0.01654458) q[2];
sx q[2];
rz(-0.51157727) q[2];
sx q[2];
rz(1.7553154) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9384267) q[1];
sx q[1];
rz(-0.90183479) q[1];
sx q[1];
rz(1.1346362) q[1];
rz(-2.8176114) q[3];
sx q[3];
rz(-0.45352648) q[3];
sx q[3];
rz(1.4676263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(-2.771634) q[2];
rz(-1.6379179) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(-1.2009719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(2.4178986) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(1.1719218) q[2];
sx q[2];
rz(-1.4234067) q[2];
sx q[2];
rz(1.1248551) q[2];
rz(0.29851144) q[3];
sx q[3];
rz(-1.1805503) q[3];
sx q[3];
rz(-1.3211484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

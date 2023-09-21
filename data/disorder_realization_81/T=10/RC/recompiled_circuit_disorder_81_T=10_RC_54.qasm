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
rz(-1.3614549) q[0];
sx q[0];
rz(1.7629495) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8060018) q[0];
sx q[0];
rz(-1.5854892) q[0];
sx q[0];
rz(-3.0528085) q[0];
x q[1];
rz(1.7832463) q[2];
sx q[2];
rz(-0.48848402) q[2];
sx q[2];
rz(-2.8993895) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8494107) q[1];
sx q[1];
rz(-2.0690448) q[1];
sx q[1];
rz(1.317418) q[1];
rz(-pi) q[2];
rz(2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(-2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1575872) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(-2.700581) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-2.5355693) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(-2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.9553604) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3620136) q[0];
sx q[0];
rz(-1.626726) q[0];
sx q[0];
rz(1.5895784) q[0];
x q[1];
rz(0.19940168) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(2.8859438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7344208) q[1];
sx q[1];
rz(-2.2312575) q[1];
sx q[1];
rz(0.15863005) q[1];
rz(1.3783185) q[3];
sx q[3];
rz(-2.2566183) q[3];
sx q[3];
rz(-1.6845077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(1.1594695) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(-0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7611258) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(0.80672112) q[0];
rz(-0.21356788) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(-0.82021964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0637174) q[0];
sx q[0];
rz(-0.88839196) q[0];
sx q[0];
rz(2.9966485) q[0];
rz(-1.9633031) q[2];
sx q[2];
rz(-0.5272534) q[2];
sx q[2];
rz(-0.96780992) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1097475) q[1];
sx q[1];
rz(-2.577563) q[1];
sx q[1];
rz(-0.86947039) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0355755) q[3];
sx q[3];
rz(-1.4651863) q[3];
sx q[3];
rz(-1.5231903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.31072581) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(-2.2303175) q[0];
rz(-0.43831929) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(1.320425) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62990084) q[0];
sx q[0];
rz(-1.9803847) q[0];
sx q[0];
rz(3.1303309) q[0];
rz(-pi) q[1];
rz(-1.3643866) q[2];
sx q[2];
rz(-2.0136535) q[2];
sx q[2];
rz(2.4368311) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.28543132) q[1];
sx q[1];
rz(-2.397575) q[1];
sx q[1];
rz(2.7475949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60453316) q[3];
sx q[3];
rz(-0.88450888) q[3];
sx q[3];
rz(-1.7741007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3115561) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(-2.2763021) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.3006166) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1039935) q[0];
sx q[0];
rz(-2.7866057) q[0];
sx q[0];
rz(0.07261891) q[0];
x q[1];
rz(-1.498921) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(-2.9938811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3805441) q[1];
sx q[1];
rz(-0.64861464) q[1];
sx q[1];
rz(0.56306871) q[1];
x q[2];
rz(-1.424765) q[3];
sx q[3];
rz(-1.8304123) q[3];
sx q[3];
rz(0.29823333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(-0.47719964) q[2];
rz(-2.9495083) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(-0.011750301) q[0];
rz(-0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5884429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4295411) q[0];
sx q[0];
rz(-2.6459604) q[0];
sx q[0];
rz(-0.80722157) q[0];
rz(-pi) q[1];
rz(-3.1197238) q[2];
sx q[2];
rz(-1.936603) q[2];
sx q[2];
rz(2.0621698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8223871) q[1];
sx q[1];
rz(-2.2777595) q[1];
sx q[1];
rz(-1.4591401) q[1];
rz(-pi) q[2];
rz(-0.074514975) q[3];
sx q[3];
rz(-2.2722368) q[3];
sx q[3];
rz(0.45267347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(1.8590415) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(-2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58105528) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(-2.1645434) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718186) q[0];
sx q[0];
rz(-2.5442113) q[0];
sx q[0];
rz(-0.13418829) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1408429) q[2];
sx q[2];
rz(-0.13204083) q[2];
sx q[2];
rz(-3.0291639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.44927412) q[1];
sx q[1];
rz(-0.34808967) q[1];
sx q[1];
rz(1.3432137) q[1];
rz(1.8924176) q[3];
sx q[3];
rz(-1.8971895) q[3];
sx q[3];
rz(-1.4332989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7523505) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(0.69865984) q[0];
rz(2.0195122) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4417393) q[0];
sx q[0];
rz(-0.2158567) q[0];
sx q[0];
rz(-2.9237843) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0935358) q[2];
sx q[2];
rz(-1.153423) q[2];
sx q[2];
rz(2.0975031) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.57751209) q[1];
sx q[1];
rz(-0.022833303) q[1];
sx q[1];
rz(-0.24502416) q[1];
x q[2];
rz(-0.28835339) q[3];
sx q[3];
rz(-0.95803146) q[3];
sx q[3];
rz(-1.310865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(-1.1784941) q[2];
rz(-1.684749) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(1.2930124) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(-2.5440149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1406527) q[0];
sx q[0];
rz(-1.6093328) q[0];
sx q[0];
rz(-1.4356104) q[0];
rz(-pi) q[1];
rz(1.045268) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(-1.8906821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24039195) q[1];
sx q[1];
rz(-1.2482093) q[1];
sx q[1];
rz(0.31369536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5556373) q[3];
sx q[3];
rz(-1.6349941) q[3];
sx q[3];
rz(2.5776598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22275816) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.2333599) q[2];
rz(0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(0.023660252) q[0];
rz(2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-2.4694209) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4155054) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(-2.0477717) q[0];
rz(1.5615084) q[2];
sx q[2];
rz(-2.0822968) q[2];
sx q[2];
rz(1.4052504) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0553531) q[1];
sx q[1];
rz(-1.9085911) q[1];
sx q[1];
rz(-0.71725459) q[1];
rz(-pi) q[2];
rz(-0.43283312) q[3];
sx q[3];
rz(-1.4308617) q[3];
sx q[3];
rz(0.39633745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(1.6379179) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621915) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(2.4178986) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(-1.9696708) q[2];
sx q[2];
rz(-1.4234067) q[2];
sx q[2];
rz(1.1248551) q[2];
rz(1.9772114) q[3];
sx q[3];
rz(-1.8462528) q[3];
sx q[3];
rz(0.13312199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
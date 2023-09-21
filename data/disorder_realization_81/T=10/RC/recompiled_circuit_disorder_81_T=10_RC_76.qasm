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
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(0.4508957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2338976) q[0];
sx q[0];
rz(-1.4820218) q[0];
sx q[0];
rz(-1.5855473) q[0];
x q[1];
rz(-2.0499174) q[2];
sx q[2];
rz(-1.6699104) q[2];
sx q[2];
rz(-1.5168158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7887468) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(0.43177859) q[1];
rz(-1.6879184) q[3];
sx q[3];
rz(-0.68713596) q[3];
sx q[3];
rz(-2.5205034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(2.297304) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-2.5355693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59250295) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(0.26309183) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(-1.1862322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037558) q[0];
sx q[0];
rz(-3.0825966) q[0];
sx q[0];
rz(-0.32365139) q[0];
rz(-pi) q[1];
rz(2.942191) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(0.25564889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1520878) q[1];
sx q[1];
rz(-0.67645914) q[1];
sx q[1];
rz(-1.3701887) q[1];
rz(-2.4466189) q[3];
sx q[3];
rz(-1.7193828) q[3];
sx q[3];
rz(2.9050764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(-0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.7611258) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-2.3348715) q[0];
rz(2.9280248) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(0.82021964) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7263111) q[0];
sx q[0];
rz(-1.4584686) q[0];
sx q[0];
rz(-2.2583654) q[0];
rz(-0.21913146) q[2];
sx q[2];
rz(-2.0543155) q[2];
sx q[2];
rz(-0.52106524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1097475) q[1];
sx q[1];
rz(-2.577563) q[1];
sx q[1];
rz(-2.2721223) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0235602) q[3];
sx q[3];
rz(-2.0327838) q[3];
sx q[3];
rz(3.1363917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(0.93079981) q[2];
rz(0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(-0.91127515) q[0];
rz(-0.43831929) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(-1.320425) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60162773) q[0];
sx q[0];
rz(-0.40973445) q[0];
sx q[0];
rz(-1.5967303) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45122066) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(2.1860683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.28543132) q[1];
sx q[1];
rz(-0.74401765) q[1];
sx q[1];
rz(2.7475949) q[1];
rz(2.5370595) q[3];
sx q[3];
rz(-0.88450888) q[3];
sx q[3];
rz(1.7741007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300366) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(-2.2763021) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.3006166) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39869719) q[0];
sx q[0];
rz(-1.5455751) q[0];
sx q[0];
rz(2.7874649) q[0];
rz(-1.498921) q[2];
sx q[2];
rz(-1.2051029) q[2];
sx q[2];
rz(2.9938811) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2758267) q[1];
sx q[1];
rz(-1.2424801) q[1];
sx q[1];
rz(-2.5715716) q[1];
x q[2];
rz(1.424765) q[3];
sx q[3];
rz(-1.3111804) q[3];
sx q[3];
rz(-2.8433593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(0.47719964) q[2];
rz(-2.9495083) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3451097) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(-0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5531497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9827305) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(1.1984675) q[0];
x q[1];
rz(-3.1197238) q[2];
sx q[2];
rz(-1.2049897) q[2];
sx q[2];
rz(-2.0621698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67571083) q[1];
sx q[1];
rz(-1.4859745) q[1];
sx q[1];
rz(0.71004962) q[1];
rz(0.074514975) q[3];
sx q[3];
rz(-2.2722368) q[3];
sx q[3];
rz(-0.45267347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(-1.2825512) q[2];
rz(1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.5605374) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(-0.67725956) q[0];
rz(2.989785) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517294) q[0];
sx q[0];
rz(-1.4954733) q[0];
sx q[0];
rz(0.59318869) q[0];
rz(-pi) q[1];
rz(0.13204079) q[2];
sx q[2];
rz(-1.5706976) q[2];
sx q[2];
rz(-1.4591109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3358826) q[1];
sx q[1];
rz(-1.4937595) q[1];
sx q[1];
rz(-1.9105934) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3903923) q[3];
sx q[3];
rz(-2.6874472) q[3];
sx q[3];
rz(-2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(2.1288669) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(-0.48721203) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-0.69865984) q[0];
rz(1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.8922071) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4835565) q[0];
sx q[0];
rz(-1.617096) q[0];
sx q[0];
rz(-2.9306843) q[0];
rz(-0.47301936) q[2];
sx q[2];
rz(-2.0447391) q[2];
sx q[2];
rz(0.75616403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74832143) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(0.022151532) q[1];
rz(-0.93805712) q[3];
sx q[3];
rz(-1.8055827) q[3];
sx q[3];
rz(-0.42890047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(-1.1784941) q[2];
rz(1.684749) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(1.8485803) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(0.59757772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4353838) q[0];
sx q[0];
rz(-1.4357114) q[0];
sx q[0];
rz(3.1027017) q[0];
rz(1.045268) q[2];
sx q[2];
rz(-1.0768441) q[2];
sx q[2];
rz(1.8906821) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4328879) q[1];
sx q[1];
rz(-1.8678027) q[1];
sx q[1];
rz(-1.9087285) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4937917) q[3];
sx q[3];
rz(-2.1553851) q[3];
sx q[3];
rz(-1.0494174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22275816) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.2333599) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(-1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
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
rz(-1.3095983) q[1];
sx q[1];
rz(2.4694209) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7260872) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(-1.0938209) q[0];
x q[1];
rz(0.01654458) q[2];
sx q[2];
rz(-2.6300154) q[2];
sx q[2];
rz(1.3862773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.086239554) q[1];
sx q[1];
rz(-1.2330016) q[1];
sx q[1];
rz(0.71725459) q[1];
x q[2];
rz(0.32398128) q[3];
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
rz(-1.2269292) q[2];
sx q[2];
rz(2.771634) q[2];
rz(-1.6379179) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.5794012) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(2.9818515) q[2];
sx q[2];
rz(-1.9651056) q[2];
sx q[2];
rz(-0.38412487) q[2];
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

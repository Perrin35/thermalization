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
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.742813) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(-0.16422693) q[0];
rz(1.3583463) q[2];
sx q[2];
rz(-0.48848402) q[2];
sx q[2];
rz(2.8993895) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3528459) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(-2.7098141) q[1];
rz(-pi) q[2];
rz(-2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1575872) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(2.8785008) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(-1.1862322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3514254) q[0];
sx q[0];
rz(-1.5520436) q[0];
sx q[0];
rz(3.0856531) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8430014) q[2];
sx q[2];
rz(-0.208374) q[2];
sx q[2];
rz(1.5339472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9895049) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(1.771404) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69497377) q[3];
sx q[3];
rz(-1.4222099) q[3];
sx q[3];
rz(-2.9050764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0120323) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(1.9821232) q[2];
rz(-0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7611258) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(-0.82021964) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0778753) q[0];
sx q[0];
rz(-2.2532007) q[0];
sx q[0];
rz(2.9966485) q[0];
x q[1];
rz(-2.9224612) q[2];
sx q[2];
rz(-2.0543155) q[2];
sx q[2];
rz(0.52106524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.324675) q[1];
sx q[1];
rz(-1.1500689) q[1];
sx q[1];
rz(0.38751985) q[1];
x q[2];
rz(1.8030274) q[3];
sx q[3];
rz(-2.665822) q[3];
sx q[3];
rz(-0.25482086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(2.8500407) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61313066) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(1.8211676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62990084) q[0];
sx q[0];
rz(-1.1612079) q[0];
sx q[0];
rz(0.011261777) q[0];
rz(2.7337012) q[2];
sx q[2];
rz(-0.48569277) q[2];
sx q[2];
rz(-2.8913468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98852324) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(-2.4371229) q[1];
rz(-pi) q[2];
rz(2.3539691) q[3];
sx q[3];
rz(-1.1155323) q[3];
sx q[3];
rz(0.61592197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1057672) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300366) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(1.8409761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1039935) q[0];
sx q[0];
rz(-0.35498699) q[0];
sx q[0];
rz(-0.07261891) q[0];
rz(-pi) q[1];
rz(0.36655764) q[2];
sx q[2];
rz(-1.5036811) q[2];
sx q[2];
rz(-1.3973438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.091150065) q[1];
sx q[1];
rz(-1.0346518) q[1];
sx q[1];
rz(1.9552783) q[1];
x q[2];
rz(2.8793094) q[3];
sx q[3];
rz(-1.7119006) q[3];
sx q[3];
rz(1.9067681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(-0.47719964) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(1.5884429) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71205157) q[0];
sx q[0];
rz(-0.49563227) q[0];
sx q[0];
rz(2.3343711) q[0];
x q[1];
rz(0.021868869) q[2];
sx q[2];
rz(-1.2049897) q[2];
sx q[2];
rz(1.0794229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3192056) q[1];
sx q[1];
rz(-2.2777595) q[1];
sx q[1];
rz(1.6824526) q[1];
x q[2];
rz(1.6586967) q[3];
sx q[3];
rz(-2.4368736) q[3];
sx q[3];
rz(0.3375012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(1.2825512) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(-2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5605374) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86977406) q[0];
sx q[0];
rz(-0.5973814) q[0];
sx q[0];
rz(3.0074044) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5706967) q[2];
sx q[2];
rz(-1.7028371) q[2];
sx q[2];
rz(-3.0299203) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2077142) q[1];
sx q[1];
rz(-1.9095451) q[1];
sx q[1];
rz(-3.0599041) q[1];
rz(-pi) q[2];
rz(-1.8924176) q[3];
sx q[3];
rz(-1.8971895) q[3];
sx q[3];
rz(-1.7082937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(0.48721203) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174719) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(-2.4429328) q[0];
rz(1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92266868) q[0];
sx q[0];
rz(-1.3601174) q[0];
sx q[0];
rz(-1.6181437) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6685733) q[2];
sx q[2];
rz(-1.0968536) q[2];
sx q[2];
rz(2.3854286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3932712) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(-3.1194411) q[1];
rz(2.8532393) q[3];
sx q[3];
rz(-2.1835612) q[3];
sx q[3];
rz(1.310865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8119048) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(-1.1784941) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(-1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(0.59757772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9874728) q[0];
sx q[0];
rz(-3.0010536) q[0];
sx q[0];
rz(-1.8494291) q[0];
x q[1];
rz(2.5848128) q[2];
sx q[2];
rz(-2.0282929) q[2];
sx q[2];
rz(-2.5533887) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.55652009) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(0.82533605) q[1];
x q[2];
rz(1.647801) q[3];
sx q[3];
rz(-2.1553851) q[3];
sx q[3];
rz(-1.0494174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22275816) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(-1.2333599) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(-1.4982769) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-0.023660252) q[0];
rz(-2.1854782) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(0.67217174) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8199352) q[0];
sx q[0];
rz(-1.5655087) q[0];
sx q[0];
rz(-1.5810285) q[0];
rz(3.1250481) q[2];
sx q[2];
rz(-0.51157727) q[2];
sx q[2];
rz(-1.7553154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2939261) q[1];
sx q[1];
rz(-2.3617509) q[1];
sx q[1];
rz(2.6508209) q[1];
x q[2];
rz(2.8176114) q[3];
sx q[3];
rz(-2.6880662) q[3];
sx q[3];
rz(1.4676263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6293634) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(2.771634) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.9406208) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621915) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(1.205668) q[2];
sx q[2];
rz(-0.42386133) q[2];
sx q[2];
rz(-0.78122666) q[2];
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
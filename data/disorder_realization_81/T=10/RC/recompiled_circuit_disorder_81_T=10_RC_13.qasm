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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.742813) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(2.9773657) q[0];
x q[1];
rz(0.11159201) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(-0.0026207844) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7887468) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(0.43177859) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4536742) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(-2.5205034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(2.297304) q[2];
rz(2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(0.94353765) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(-1.9553604) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79016722) q[0];
sx q[0];
rz(-1.5895491) q[0];
sx q[0];
rz(-3.0856531) q[0];
rz(-pi) q[1];
rz(-0.19940168) q[2];
sx q[2];
rz(-1.6316895) q[2];
sx q[2];
rz(2.8859438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9895049) q[1];
sx q[1];
rz(-0.67645914) q[1];
sx q[1];
rz(1.771404) q[1];
rz(1.7632742) q[3];
sx q[3];
rz(-2.2566183) q[3];
sx q[3];
rz(1.6845077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(-0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.3804669) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(-0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41528156) q[0];
sx q[0];
rz(-1.4584686) q[0];
sx q[0];
rz(-0.88322722) q[0];
rz(-pi) q[1];
rz(-2.9224612) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-0.52106524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.324675) q[1];
sx q[1];
rz(-1.9915238) q[1];
sx q[1];
rz(-2.7540728) q[1];
x q[2];
rz(-1.8030274) q[3];
sx q[3];
rz(-0.47577061) q[3];
sx q[3];
rz(2.8867718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(2.2107928) q[2];
rz(-0.15549774) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.3221909) q[1];
sx q[1];
rz(1.8211676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5399649) q[0];
sx q[0];
rz(-0.40973445) q[0];
sx q[0];
rz(1.5967303) q[0];
x q[1];
rz(-0.45122066) q[2];
sx q[2];
rz(-1.3845452) q[2];
sx q[2];
rz(2.1860683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79975407) q[1];
sx q[1];
rz(-2.2463887) q[1];
sx q[1];
rz(-1.9104596) q[1];
rz(-pi) q[2];
rz(-0.96418013) q[3];
sx q[3];
rz(-0.88084953) q[3];
sx q[3];
rz(2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1057672) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(0.17677447) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.86529055) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.8409761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39869719) q[0];
sx q[0];
rz(-1.5455751) q[0];
sx q[0];
rz(0.35412776) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.498921) q[2];
sx q[2];
rz(-1.2051029) q[2];
sx q[2];
rz(-0.14771151) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.865766) q[1];
sx q[1];
rz(-1.8991125) q[1];
sx q[1];
rz(-0.57002108) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50118581) q[3];
sx q[3];
rz(-2.8445344) q[3];
sx q[3];
rz(-2.3230769) q[3];
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
rz(0.19208433) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(-2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79648298) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(-1.5531497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1588622) q[0];
sx q[0];
rz(-1.2356865) q[0];
sx q[0];
rz(-1.1984675) q[0];
rz(1.2049098) q[2];
sx q[2];
rz(-1.591218) q[2];
sx q[2];
rz(0.4991971) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1482684) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(0.12970129) q[1];
x q[2];
rz(-1.6586967) q[3];
sx q[3];
rz(-0.70471901) q[3];
sx q[3];
rz(0.3375012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.50756303) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(1.7717308) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(-2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58986321) q[0];
sx q[0];
rz(-1.6461194) q[0];
sx q[0];
rz(-2.548404) q[0];
rz(-pi) q[1];
rz(-3.0095519) q[2];
sx q[2];
rz(-1.570895) q[2];
sx q[2];
rz(-1.6824818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9338785) q[1];
sx q[1];
rz(-1.9095451) q[1];
sx q[1];
rz(-3.0599041) q[1];
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
x q[1];
rz(-1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(1.1879454) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(0.69865984) q[0];
rz(1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.8922071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4835565) q[0];
sx q[0];
rz(-1.5244966) q[0];
sx q[0];
rz(-2.9306843) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84455372) q[2];
sx q[2];
rz(-2.484998) q[2];
sx q[2];
rz(-3.055228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.318995) q[1];
sx q[1];
rz(-1.5486451) q[1];
sx q[1];
rz(-1.5763361) q[1];
x q[2];
rz(-1.1864248) q[3];
sx q[3];
rz(-0.66925183) q[3];
sx q[3];
rz(-0.8347019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8119048) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(-1.684749) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(-0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6417398) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(-1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(0.59757772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15411988) q[0];
sx q[0];
rz(-3.0010536) q[0];
sx q[0];
rz(-1.8494291) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0963247) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(-1.2509105) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4328879) q[1];
sx q[1];
rz(-1.8678027) q[1];
sx q[1];
rz(1.9087285) q[1];
rz(-pi) q[2];
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
sx q[1];
rz(-pi/2) q[1];
rz(-2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.4982769) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047886588) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-0.023660252) q[0];
rz(2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-2.4694209) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2490847) q[0];
sx q[0];
rz(-1.5605643) q[0];
sx q[0];
rz(0.0052878629) q[0];
rz(-2.6300738) q[2];
sx q[2];
rz(-1.5626972) q[2];
sx q[2];
rz(0.17009232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.086239554) q[1];
sx q[1];
rz(-1.9085911) q[1];
sx q[1];
rz(-2.4243381) q[1];
x q[2];
rz(-1.4168596) q[3];
sx q[3];
rz(-1.9991176) q[3];
sx q[3];
rz(1.1101013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(2.4178986) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(1.9359246) q[2];
sx q[2];
rz(-2.7177313) q[2];
sx q[2];
rz(2.360366) q[2];
rz(-0.29851144) q[3];
sx q[3];
rz(-1.9610423) q[3];
sx q[3];
rz(1.8204443) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4616315) q[0];
sx q[0];
rz(-2.2357219) q[0];
sx q[0];
rz(1.0933956) q[0];
rz(-1.6361341) q[1];
sx q[1];
rz(-1.1727762) q[1];
sx q[1];
rz(2.7170031) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6494006) q[0];
sx q[0];
rz(-1.369242) q[0];
sx q[0];
rz(-0.038945065) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4167673) q[2];
sx q[2];
rz(-1.6209148) q[2];
sx q[2];
rz(2.199204) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72480768) q[1];
sx q[1];
rz(-2.8832286) q[1];
sx q[1];
rz(-1.2577673) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2810117) q[3];
sx q[3];
rz(-1.672555) q[3];
sx q[3];
rz(-1.4836131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7666185) q[2];
sx q[2];
rz(-1.582229) q[2];
sx q[2];
rz(-1.9826822) q[2];
rz(0.22289395) q[3];
sx q[3];
rz(-1.4054207) q[3];
sx q[3];
rz(-2.8515653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463773) q[0];
sx q[0];
rz(-1.615849) q[0];
sx q[0];
rz(2.0624397) q[0];
rz(1.4214628) q[1];
sx q[1];
rz(-2.454897) q[1];
sx q[1];
rz(-3.0770643) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9345029) q[0];
sx q[0];
rz(-0.10182589) q[0];
sx q[0];
rz(2.1125211) q[0];
rz(-pi) q[1];
rz(1.5428434) q[2];
sx q[2];
rz(-2.0468901) q[2];
sx q[2];
rz(-0.7715438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5777305) q[1];
sx q[1];
rz(-2.3268493) q[1];
sx q[1];
rz(-1.2352863) q[1];
rz(-1.7713624) q[3];
sx q[3];
rz(-2.3820357) q[3];
sx q[3];
rz(2.650039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9885538) q[2];
sx q[2];
rz(-1.3554074) q[2];
sx q[2];
rz(-1.1174196) q[2];
rz(0.85876632) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(-2.7045238) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3808463) q[0];
sx q[0];
rz(-2.8570211) q[0];
sx q[0];
rz(0.17307702) q[0];
rz(-2.1848047) q[1];
sx q[1];
rz(-1.0280949) q[1];
sx q[1];
rz(1.0622567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5231402) q[0];
sx q[0];
rz(-1.5644531) q[0];
sx q[0];
rz(2.2940318) q[0];
rz(-pi) q[1];
rz(-0.64086242) q[2];
sx q[2];
rz(-2.0483565) q[2];
sx q[2];
rz(-1.8030082) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77327496) q[1];
sx q[1];
rz(-2.2466806) q[1];
sx q[1];
rz(-2.2818517) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88835399) q[3];
sx q[3];
rz(-2.0815947) q[3];
sx q[3];
rz(-0.79456431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.7303077) q[2];
sx q[2];
rz(-1.9399119) q[2];
sx q[2];
rz(0.74472204) q[2];
rz(0.64106411) q[3];
sx q[3];
rz(-1.791626) q[3];
sx q[3];
rz(-2.6509269) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915801) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(2.2963754) q[0];
rz(-0.40329626) q[1];
sx q[1];
rz(-1.6019628) q[1];
sx q[1];
rz(-0.32803112) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.302264) q[0];
sx q[0];
rz(-1.3899046) q[0];
sx q[0];
rz(2.9995287) q[0];
rz(-pi) q[1];
rz(2.8745266) q[2];
sx q[2];
rz(-0.24644463) q[2];
sx q[2];
rz(1.3738969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5535132) q[1];
sx q[1];
rz(-2.5633286) q[1];
sx q[1];
rz(0.2167313) q[1];
rz(-0.27356903) q[3];
sx q[3];
rz(-1.0427339) q[3];
sx q[3];
rz(-2.7298965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8959117) q[2];
sx q[2];
rz(-2.3978105) q[2];
sx q[2];
rz(-0.15920676) q[2];
rz(-3.1253452) q[3];
sx q[3];
rz(-2.1135606) q[3];
sx q[3];
rz(-0.73133674) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74723393) q[0];
sx q[0];
rz(-1.1927698) q[0];
sx q[0];
rz(2.0514945) q[0];
rz(-3.0116426) q[1];
sx q[1];
rz(-2.3268685) q[1];
sx q[1];
rz(1.6158339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65924683) q[0];
sx q[0];
rz(-2.1123943) q[0];
sx q[0];
rz(-2.0603176) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1431057) q[2];
sx q[2];
rz(-0.18015716) q[2];
sx q[2];
rz(-0.13979543) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33541778) q[1];
sx q[1];
rz(-1.5433784) q[1];
sx q[1];
rz(-2.5078678) q[1];
x q[2];
rz(0.035338621) q[3];
sx q[3];
rz(-1.3900847) q[3];
sx q[3];
rz(-0.23739761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7502363) q[2];
sx q[2];
rz(-0.82166925) q[2];
sx q[2];
rz(2.6643122) q[2];
rz(2.9376049) q[3];
sx q[3];
rz(-0.19652772) q[3];
sx q[3];
rz(2.5175214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9452962) q[0];
sx q[0];
rz(-2.3260703) q[0];
sx q[0];
rz(-0.77769172) q[0];
rz(-0.6174736) q[1];
sx q[1];
rz(-1.4910881) q[1];
sx q[1];
rz(-0.9309353) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8645267) q[0];
sx q[0];
rz(-1.1204989) q[0];
sx q[0];
rz(-0.45876512) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7226376) q[2];
sx q[2];
rz(-1.1280931) q[2];
sx q[2];
rz(2.1328164) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1738064) q[1];
sx q[1];
rz(-1.5266982) q[1];
sx q[1];
rz(1.2323061) q[1];
rz(-pi) q[2];
x q[2];
rz(2.191675) q[3];
sx q[3];
rz(-1.2395879) q[3];
sx q[3];
rz(3.1049867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4882539) q[2];
sx q[2];
rz(-1.3335756) q[2];
sx q[2];
rz(-0.2674357) q[2];
rz(2.2087162) q[3];
sx q[3];
rz(-1.2837912) q[3];
sx q[3];
rz(0.4944087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6701508) q[0];
sx q[0];
rz(-1.9981367) q[0];
sx q[0];
rz(0.66584051) q[0];
rz(0.2039856) q[1];
sx q[1];
rz(-0.98122707) q[1];
sx q[1];
rz(-1.9542255) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4543633) q[0];
sx q[0];
rz(-2.6626514) q[0];
sx q[0];
rz(2.9402551) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6186278) q[2];
sx q[2];
rz(-2.209216) q[2];
sx q[2];
rz(2.0948727) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.07908917) q[1];
sx q[1];
rz(-2.3983208) q[1];
sx q[1];
rz(2.7447002) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3755619) q[3];
sx q[3];
rz(-2.8831006) q[3];
sx q[3];
rz(3.0743161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5953956) q[2];
sx q[2];
rz(-0.51837817) q[2];
sx q[2];
rz(-0.67016822) q[2];
rz(-2.8403122) q[3];
sx q[3];
rz(-1.8471085) q[3];
sx q[3];
rz(2.7231351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30348521) q[0];
sx q[0];
rz(-1.0024339) q[0];
sx q[0];
rz(2.0759034) q[0];
rz(2.4844555) q[1];
sx q[1];
rz(-2.4333351) q[1];
sx q[1];
rz(1.4297952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0428061) q[0];
sx q[0];
rz(-1.0142066) q[0];
sx q[0];
rz(2.7419529) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3892608) q[2];
sx q[2];
rz(-1.8971271) q[2];
sx q[2];
rz(-1.8472613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.011466) q[1];
sx q[1];
rz(-1.3523755) q[1];
sx q[1];
rz(2.7232724) q[1];
rz(-pi) q[2];
rz(1.4265238) q[3];
sx q[3];
rz(-1.7585837) q[3];
sx q[3];
rz(-0.71545631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8810205) q[2];
sx q[2];
rz(-1.2994956) q[2];
sx q[2];
rz(-1.0183498) q[2];
rz(-1.5492505) q[3];
sx q[3];
rz(-1.6377623) q[3];
sx q[3];
rz(2.3840267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26579478) q[0];
sx q[0];
rz(-2.6261411) q[0];
sx q[0];
rz(0.96762586) q[0];
rz(2.8014917) q[1];
sx q[1];
rz(-0.83931559) q[1];
sx q[1];
rz(-2.1077154) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.114703) q[0];
sx q[0];
rz(-1.4324766) q[0];
sx q[0];
rz(-1.0113495) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62693779) q[2];
sx q[2];
rz(-1.9276918) q[2];
sx q[2];
rz(-0.42742768) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.25701538) q[1];
sx q[1];
rz(-2.3126763) q[1];
sx q[1];
rz(-2.1121426) q[1];
x q[2];
rz(1.4236064) q[3];
sx q[3];
rz(-1.7924252) q[3];
sx q[3];
rz(-2.2972884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10244441) q[2];
sx q[2];
rz(-1.4463964) q[2];
sx q[2];
rz(-0.038912494) q[2];
rz(0.14032042) q[3];
sx q[3];
rz(-0.084429927) q[3];
sx q[3];
rz(-1.3154359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6572396) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(-1.2580309) q[0];
rz(2.925442) q[1];
sx q[1];
rz(-0.70544568) q[1];
sx q[1];
rz(0.36453882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98436873) q[0];
sx q[0];
rz(-1.8175392) q[0];
sx q[0];
rz(1.1086199) q[0];
x q[1];
rz(2.0272354) q[2];
sx q[2];
rz(-0.82122148) q[2];
sx q[2];
rz(2.8783893) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3033402) q[1];
sx q[1];
rz(-1.5152182) q[1];
sx q[1];
rz(1.5838632) q[1];
rz(-pi) q[2];
x q[2];
rz(2.632439) q[3];
sx q[3];
rz(-3.0753071) q[3];
sx q[3];
rz(-1.8788415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2898499) q[2];
sx q[2];
rz(-1.936329) q[2];
sx q[2];
rz(0.64129788) q[2];
rz(1.4801721) q[3];
sx q[3];
rz(-1.8931754) q[3];
sx q[3];
rz(-2.4680468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4079473) q[0];
sx q[0];
rz(-2.6483364) q[0];
sx q[0];
rz(-0.37089621) q[0];
rz(0.98450487) q[1];
sx q[1];
rz(-1.6592818) q[1];
sx q[1];
rz(2.0936113) q[1];
rz(3.0251236) q[2];
sx q[2];
rz(-1.3745728) q[2];
sx q[2];
rz(-0.41043888) q[2];
rz(1.0897286) q[3];
sx q[3];
rz(-1.0952598) q[3];
sx q[3];
rz(2.918603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

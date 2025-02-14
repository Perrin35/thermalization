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
rz(-0.42458951) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49219201) q[0];
sx q[0];
rz(-1.7723506) q[0];
sx q[0];
rz(0.038945065) q[0];
x q[1];
rz(1.7248254) q[2];
sx q[2];
rz(-1.5206778) q[2];
sx q[2];
rz(-0.94238867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0478617) q[1];
sx q[1];
rz(-1.3252531) q[1];
sx q[1];
rz(0.081201038) q[1];
rz(-1.4154451) q[3];
sx q[3];
rz(-2.4253824) q[3];
sx q[3];
rz(2.9368212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37497416) q[2];
sx q[2];
rz(-1.582229) q[2];
sx q[2];
rz(1.1589104) q[2];
rz(2.9186987) q[3];
sx q[3];
rz(-1.736172) q[3];
sx q[3];
rz(-2.8515653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7952154) q[0];
sx q[0];
rz(-1.615849) q[0];
sx q[0];
rz(1.079153) q[0];
rz(1.4214628) q[1];
sx q[1];
rz(-0.68669569) q[1];
sx q[1];
rz(3.0770643) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8046605) q[0];
sx q[0];
rz(-1.6580026) q[0];
sx q[0];
rz(0.052636458) q[0];
rz(-pi) q[1];
rz(3.0874443) q[2];
sx q[2];
rz(-2.6647419) q[2];
sx q[2];
rz(-2.3091174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5777305) q[1];
sx q[1];
rz(-2.3268493) q[1];
sx q[1];
rz(-1.2352863) q[1];
x q[2];
rz(2.9546176) q[3];
sx q[3];
rz(-0.83010736) q[3];
sx q[3];
rz(2.376698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1530389) q[2];
sx q[2];
rz(-1.3554074) q[2];
sx q[2];
rz(2.0241731) q[2];
rz(0.85876632) q[3];
sx q[3];
rz(-0.78341165) q[3];
sx q[3];
rz(2.7045238) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7607464) q[0];
sx q[0];
rz(-2.8570211) q[0];
sx q[0];
rz(-2.9685156) q[0];
rz(2.1848047) q[1];
sx q[1];
rz(-2.1134977) q[1];
sx q[1];
rz(-2.079336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040471023) q[0];
sx q[0];
rz(-0.72325828) q[0];
sx q[0];
rz(-1.5612119) q[0];
rz(-pi) q[1];
rz(-0.99750869) q[2];
sx q[2];
rz(-2.1305704) q[2];
sx q[2];
rz(2.5791175) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7159285) q[1];
sx q[1];
rz(-2.2032208) q[1];
sx q[1];
rz(0.683149) q[1];
x q[2];
rz(-2.5161602) q[3];
sx q[3];
rz(-2.1533826) q[3];
sx q[3];
rz(-2.7435477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(0.49066576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500126) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(-2.2963754) q[0];
rz(-2.7382964) q[1];
sx q[1];
rz(-1.5396298) q[1];
sx q[1];
rz(-0.32803112) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.302264) q[0];
sx q[0];
rz(-1.3899046) q[0];
sx q[0];
rz(0.14206391) q[0];
x q[1];
rz(0.26706605) q[2];
sx q[2];
rz(-2.895148) q[2];
sx q[2];
rz(-1.7676958) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58807948) q[1];
sx q[1];
rz(-0.57826406) q[1];
sx q[1];
rz(0.2167313) q[1];
x q[2];
rz(2.115504) q[3];
sx q[3];
rz(-1.8063365) q[3];
sx q[3];
rz(-2.1229471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8959117) q[2];
sx q[2];
rz(-0.74378219) q[2];
sx q[2];
rz(0.15920676) q[2];
rz(-0.016247449) q[3];
sx q[3];
rz(-2.1135606) q[3];
sx q[3];
rz(-2.4102559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74723393) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(-2.0514945) q[0];
rz(-3.0116426) q[1];
sx q[1];
rz(-2.3268685) q[1];
sx q[1];
rz(-1.5257588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14262202) q[0];
sx q[0];
rz(-2.4283613) q[0];
sx q[0];
rz(-0.66340982) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0661929) q[2];
sx q[2];
rz(-1.4070208) q[2];
sx q[2];
rz(-2.5679156) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.255521) q[1];
sx q[1];
rz(-0.93734765) q[1];
sx q[1];
rz(-1.5367763) q[1];
rz(-pi) q[2];
rz(0.035338621) q[3];
sx q[3];
rz(-1.7515079) q[3];
sx q[3];
rz(0.23739761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7502363) q[2];
sx q[2];
rz(-0.82166925) q[2];
sx q[2];
rz(0.4772805) q[2];
rz(-2.9376049) q[3];
sx q[3];
rz(-0.19652772) q[3];
sx q[3];
rz(0.6240713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1962965) q[0];
sx q[0];
rz(-0.81552234) q[0];
sx q[0];
rz(0.77769172) q[0];
rz(-2.5241191) q[1];
sx q[1];
rz(-1.4910881) q[1];
sx q[1];
rz(-2.2106574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57140152) q[0];
sx q[0];
rz(-0.63136111) q[0];
sx q[0];
rz(-0.82918138) q[0];
x q[1];
rz(1.0920877) q[2];
sx q[2];
rz(-1.9472113) q[2];
sx q[2];
rz(-0.37351721) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.41250788) q[1];
sx q[1];
rz(-1.9089444) q[1];
sx q[1];
rz(-0.04674712) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.10465) q[3];
sx q[3];
rz(-2.4483238) q[3];
sx q[3];
rz(1.9609089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4882539) q[2];
sx q[2];
rz(-1.3335756) q[2];
sx q[2];
rz(0.2674357) q[2];
rz(-0.93287647) q[3];
sx q[3];
rz(-1.2837912) q[3];
sx q[3];
rz(-2.647184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47144181) q[0];
sx q[0];
rz(-1.9981367) q[0];
sx q[0];
rz(0.66584051) q[0];
rz(-2.937607) q[1];
sx q[1];
rz(-0.98122707) q[1];
sx q[1];
rz(1.1873672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(0.86249528) q[2];
sx q[2];
rz(-1.1581026) q[2];
sx q[2];
rz(0.19315456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59661181) q[1];
sx q[1];
rz(-2.2447733) q[1];
sx q[1];
rz(1.9121714) q[1];
x q[2];
rz(1.7660308) q[3];
sx q[3];
rz(-2.8831006) q[3];
sx q[3];
rz(0.067276567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54619706) q[2];
sx q[2];
rz(-2.6232145) q[2];
sx q[2];
rz(2.4714244) q[2];
rz(-0.3012805) q[3];
sx q[3];
rz(-1.2944841) q[3];
sx q[3];
rz(-0.41845751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30348521) q[0];
sx q[0];
rz(-2.1391588) q[0];
sx q[0];
rz(-1.0656892) q[0];
rz(0.65713716) q[1];
sx q[1];
rz(-2.4333351) q[1];
sx q[1];
rz(1.7117975) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0428061) q[0];
sx q[0];
rz(-2.1273861) q[0];
sx q[0];
rz(2.7419529) q[0];
rz(-pi) q[1];
rz(-2.6515342) q[2];
sx q[2];
rz(-0.37184162) q[2];
sx q[2];
rz(-1.3272367) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.34463528) q[1];
sx q[1];
rz(-1.1630139) q[1];
sx q[1];
rz(1.8090882) q[1];
rz(-pi) q[2];
rz(2.9518806) q[3];
sx q[3];
rz(-1.7125152) q[3];
sx q[3];
rz(-0.82822463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26057217) q[2];
sx q[2];
rz(-1.2994956) q[2];
sx q[2];
rz(-2.1232429) q[2];
rz(1.5492505) q[3];
sx q[3];
rz(-1.5038303) q[3];
sx q[3];
rz(2.3840267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757979) q[0];
sx q[0];
rz(-0.51545155) q[0];
sx q[0];
rz(-2.1739668) q[0];
rz(2.8014917) q[1];
sx q[1];
rz(-2.3022771) q[1];
sx q[1];
rz(-1.0338773) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5115625) q[0];
sx q[0];
rz(-2.1242737) q[0];
sx q[0];
rz(-0.16279499) q[0];
x q[1];
rz(0.62693779) q[2];
sx q[2];
rz(-1.2139008) q[2];
sx q[2];
rz(0.42742768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.441923) q[1];
sx q[1];
rz(-1.1811387) q[1];
sx q[1];
rz(-0.81894919) q[1];
x q[2];
rz(-1.4236064) q[3];
sx q[3];
rz(-1.3491674) q[3];
sx q[3];
rz(0.8443043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10244441) q[2];
sx q[2];
rz(-1.6951963) q[2];
sx q[2];
rz(-0.038912494) q[2];
rz(0.14032042) q[3];
sx q[3];
rz(-0.084429927) q[3];
sx q[3];
rz(1.8261568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572396) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(-1.8835618) q[0];
rz(2.925442) q[1];
sx q[1];
rz(-2.436147) q[1];
sx q[1];
rz(2.7770538) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4340857) q[0];
sx q[0];
rz(-1.123652) q[0];
sx q[0];
rz(-0.27430496) q[0];
rz(-2.3380525) q[2];
sx q[2];
rz(-1.8992956) q[2];
sx q[2];
rz(0.9847509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8748624) q[1];
sx q[1];
rz(-1.583843) q[1];
sx q[1];
rz(-3.0860098) q[1];
rz(2.632439) q[3];
sx q[3];
rz(-0.066285523) q[3];
sx q[3];
rz(-1.2627511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2898499) q[2];
sx q[2];
rz(-1.2052636) q[2];
sx q[2];
rz(2.5002948) q[2];
rz(1.6614206) q[3];
sx q[3];
rz(-1.8931754) q[3];
sx q[3];
rz(-0.67354584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73364532) q[0];
sx q[0];
rz(-2.6483364) q[0];
sx q[0];
rz(-0.37089621) q[0];
rz(-2.1570878) q[1];
sx q[1];
rz(-1.6592818) q[1];
sx q[1];
rz(2.0936113) q[1];
rz(-2.099809) q[2];
sx q[2];
rz(-2.9137901) q[2];
sx q[2];
rz(0.13005039) q[2];
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

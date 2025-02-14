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
rz(-0.32831353) q[0];
sx q[0];
rz(5.1483122) q[0];
sx q[0];
rz(6.8490646) q[0];
rz(0.24428754) q[1];
sx q[1];
rz(-1.3574418) q[1];
sx q[1];
rz(0.53425962) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0097144011) q[0];
sx q[0];
rz(-1.8410972) q[0];
sx q[0];
rz(0.86608315) q[0];
rz(-pi) q[1];
rz(-0.59492971) q[2];
sx q[2];
rz(-1.8244566) q[2];
sx q[2];
rz(0.67753917) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45785832) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(0.079816743) q[1];
rz(-pi) q[2];
rz(0.45716469) q[3];
sx q[3];
rz(-1.7948397) q[3];
sx q[3];
rz(0.30686298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.013022097) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(2.8924083) q[2];
rz(-1.7976044) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(-3.0715166) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4942112) q[0];
sx q[0];
rz(-1.2263612) q[0];
sx q[0];
rz(-0.98980728) q[0];
rz(2.3509332) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(0.31165037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383331) q[0];
sx q[0];
rz(-2.0796392) q[0];
sx q[0];
rz(-2.3340387) q[0];
rz(-pi) q[1];
rz(1.7575532) q[2];
sx q[2];
rz(-1.4644564) q[2];
sx q[2];
rz(-0.76923907) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6271882) q[1];
sx q[1];
rz(-1.7429105) q[1];
sx q[1];
rz(1.8132416) q[1];
rz(-pi) q[2];
rz(2.1856863) q[3];
sx q[3];
rz(-1.2887738) q[3];
sx q[3];
rz(2.7789555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1404169) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(0.95138335) q[2];
rz(-0.22826711) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15762873) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(0.11446318) q[0];
rz(2.139367) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(0.1618298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9660258) q[0];
sx q[0];
rz(-1.1822261) q[0];
sx q[0];
rz(-2.2184664) q[0];
x q[1];
rz(1.6587018) q[2];
sx q[2];
rz(-1.958333) q[2];
sx q[2];
rz(-1.0532925) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.075167478) q[1];
sx q[1];
rz(-2.4801765) q[1];
sx q[1];
rz(0.39638806) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80064818) q[3];
sx q[3];
rz(-1.8654239) q[3];
sx q[3];
rz(-1.9484474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3942922) q[2];
sx q[2];
rz(-2.6960399) q[2];
sx q[2];
rz(-0.14075819) q[2];
rz(2.5898139) q[3];
sx q[3];
rz(-1.2740302) q[3];
sx q[3];
rz(2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.20525876) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(-2.4667013) q[0];
rz(2.9871125) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-0.45796576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7688528) q[0];
sx q[0];
rz(-0.64592664) q[0];
sx q[0];
rz(-0.67966184) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0003597) q[2];
sx q[2];
rz(-1.8515966) q[2];
sx q[2];
rz(-2.8605638) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6667216) q[1];
sx q[1];
rz(-2.7242766) q[1];
sx q[1];
rz(0.46837193) q[1];
rz(-pi) q[2];
rz(1.4639888) q[3];
sx q[3];
rz(-2.0764917) q[3];
sx q[3];
rz(3.0042574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13350479) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(2.8221455) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(-2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9820246) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(-1.2215479) q[0];
rz(2.9087861) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(-2.1585042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1637835) q[0];
sx q[0];
rz(-2.0931232) q[0];
sx q[0];
rz(0.99951007) q[0];
rz(2.1914761) q[2];
sx q[2];
rz(-1.4779556) q[2];
sx q[2];
rz(0.29618057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5447588) q[1];
sx q[1];
rz(-1.7167517) q[1];
sx q[1];
rz(-0.18827855) q[1];
rz(1.2125254) q[3];
sx q[3];
rz(-2.6942109) q[3];
sx q[3];
rz(3.0606396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7588707) q[2];
sx q[2];
rz(-1.657234) q[2];
sx q[2];
rz(2.7663084) q[2];
rz(-2.0659857) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4009092) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(1.4373454) q[0];
rz(-3.0628693) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9288612) q[0];
sx q[0];
rz(-0.91758801) q[0];
sx q[0];
rz(-1.3034588) q[0];
x q[1];
rz(0.55193837) q[2];
sx q[2];
rz(-1.9935529) q[2];
sx q[2];
rz(1.1814337) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.50997616) q[1];
sx q[1];
rz(-2.2215034) q[1];
sx q[1];
rz(-1.0751616) q[1];
x q[2];
rz(0.6491361) q[3];
sx q[3];
rz(-0.86926354) q[3];
sx q[3];
rz(2.8809406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6846201) q[2];
sx q[2];
rz(-2.520884) q[2];
sx q[2];
rz(0.6753298) q[2];
rz(1.5572549) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(2.5244782) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5439344) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-2.6455998) q[0];
rz(1.687382) q[1];
sx q[1];
rz(-1.0554689) q[1];
sx q[1];
rz(-1.1167663) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1710715) q[0];
sx q[0];
rz(-1.9200033) q[0];
sx q[0];
rz(1.1961351) q[0];
rz(-pi) q[1];
rz(2.1082868) q[2];
sx q[2];
rz(-1.8040207) q[2];
sx q[2];
rz(1.6526412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45124829) q[1];
sx q[1];
rz(-2.6670229) q[1];
sx q[1];
rz(0.75466538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12567606) q[3];
sx q[3];
rz(-0.45914859) q[3];
sx q[3];
rz(2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3108959) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(-2.7867219) q[2];
rz(-2.3762459) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(3.0520458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68882051) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(-0.76559693) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(-1.3227468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046503456) q[0];
sx q[0];
rz(-2.8757767) q[0];
sx q[0];
rz(2.9429716) q[0];
rz(-pi) q[1];
rz(2.2501588) q[2];
sx q[2];
rz(-1.7560894) q[2];
sx q[2];
rz(-2.7012555) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8812528) q[1];
sx q[1];
rz(-1.5872721) q[1];
sx q[1];
rz(-0.50838281) q[1];
rz(-pi) q[2];
rz(-2.3115308) q[3];
sx q[3];
rz(-0.30127159) q[3];
sx q[3];
rz(0.64727441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65779296) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(2.9198809) q[2];
rz(2.0486369) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(-2.4634585) q[3];
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
rz(-pi/2) q[0];
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
rz(-3.0677277) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(2.8644323) q[0];
rz(0.74835888) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(1.998273) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7171913) q[0];
sx q[0];
rz(-0.36866185) q[0];
sx q[0];
rz(-0.98510806) q[0];
x q[1];
rz(1.3242993) q[2];
sx q[2];
rz(-0.72028226) q[2];
sx q[2];
rz(1.5697073) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4797693) q[1];
sx q[1];
rz(-2.7854438) q[1];
sx q[1];
rz(-1.1060064) q[1];
rz(-pi) q[2];
rz(-2.7128814) q[3];
sx q[3];
rz(-2.747251) q[3];
sx q[3];
rz(-2.0286146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9144168) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(0.048967036) q[2];
rz(1.0660727) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(-2.9858885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0186036) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(-3.1220806) q[0];
rz(-0.77876577) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(1.5787517) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9492968) q[0];
sx q[0];
rz(-1.4817294) q[0];
sx q[0];
rz(1.7013447) q[0];
x q[1];
rz(0.18264211) q[2];
sx q[2];
rz(-2.172875) q[2];
sx q[2];
rz(-1.6537567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.409118) q[1];
sx q[1];
rz(-0.98608398) q[1];
sx q[1];
rz(0.52701108) q[1];
x q[2];
rz(0.14871493) q[3];
sx q[3];
rz(-0.84245719) q[3];
sx q[3];
rz(-3.0535798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.19266263) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(0.16555244) q[2];
rz(-2.5340269) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7623357) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(1.634585) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(1.4277369) q[2];
sx q[2];
rz(-2.3410627) q[2];
sx q[2];
rz(-1.9973199) q[2];
rz(2.8589917) q[3];
sx q[3];
rz(-0.62487124) q[3];
sx q[3];
rz(1.665653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

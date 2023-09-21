OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4484654) q[0];
sx q[0];
rz(-2.6187596) q[0];
sx q[0];
rz(-2.5180106) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(2.6410988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60808027) q[0];
sx q[0];
rz(-2.170616) q[0];
sx q[0];
rz(1.5754726) q[0];
x q[1];
rz(3.0481553) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(1.1364394) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1653633) q[1];
sx q[1];
rz(-2.2803218) q[1];
sx q[1];
rz(-0.43433365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1327098) q[3];
sx q[3];
rz(-1.4287018) q[3];
sx q[3];
rz(-2.4349468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3502675) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(1.2228489) q[2];
rz(-1.6932999) q[3];
sx q[3];
rz(-0.99213123) q[3];
sx q[3];
rz(2.1712415) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7704849) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(2.0626542) q[0];
rz(1.3868015) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-2.4761377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5011713) q[0];
sx q[0];
rz(-1.8542395) q[0];
sx q[0];
rz(-0.47407504) q[0];
rz(0.47768728) q[2];
sx q[2];
rz(-0.3590695) q[2];
sx q[2];
rz(-1.8600841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14109719) q[1];
sx q[1];
rz(-2.5971203) q[1];
sx q[1];
rz(-0.46846868) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1075675) q[3];
sx q[3];
rz(-0.861654) q[3];
sx q[3];
rz(2.2505086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60454303) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(-3.0664505) q[2];
rz(-1.6710619) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(-2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(1.8485908) q[1];
sx q[1];
rz(-1.8483775) q[1];
sx q[1];
rz(-1.4000777) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132719) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(-0.4011641) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4351575) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(-0.26091012) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2565508) q[1];
sx q[1];
rz(-2.1532144) q[1];
sx q[1];
rz(-2.7529653) q[1];
rz(-2.9754144) q[3];
sx q[3];
rz(-1.2697392) q[3];
sx q[3];
rz(-2.500246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.489958) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(1.970132) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(1.4452274) q[0];
rz(-1.6943278) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(-2.7935374) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8349583) q[0];
sx q[0];
rz(-1.8638896) q[0];
sx q[0];
rz(-2.4819863) q[0];
x q[1];
rz(1.939417) q[2];
sx q[2];
rz(-2.7432132) q[2];
sx q[2];
rz(-0.57952651) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0192249) q[1];
sx q[1];
rz(-0.61673635) q[1];
sx q[1];
rz(1.4160181) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9799558) q[3];
sx q[3];
rz(-1.4378387) q[3];
sx q[3];
rz(1.1383575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7248914) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(0.42281881) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(0.53043956) q[0];
rz(0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.8431429) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2753678) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(2.2989681) q[0];
x q[1];
rz(-0.41579397) q[2];
sx q[2];
rz(-2.4928164) q[2];
sx q[2];
rz(-0.4817889) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7309155) q[1];
sx q[1];
rz(-1.2123322) q[1];
sx q[1];
rz(1.4645542) q[1];
rz(2.2523746) q[3];
sx q[3];
rz(-0.71435706) q[3];
sx q[3];
rz(2.0444972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(-0.47362622) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(2.1437058) q[0];
rz(2.2672794) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-0.46674892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7346863) q[0];
sx q[0];
rz(-2.2285301) q[0];
sx q[0];
rz(-1.2975733) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7919962) q[2];
sx q[2];
rz(-1.9890519) q[2];
sx q[2];
rz(-2.7115371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14086831) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(2.9904757) q[1];
rz(-pi) q[2];
rz(-1.0420226) q[3];
sx q[3];
rz(-0.18897945) q[3];
sx q[3];
rz(0.46940645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85990396) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(-1.5647282) q[2];
rz(0.62670296) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(-0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3307813) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(-0.41982857) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5931169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67328582) q[0];
sx q[0];
rz(-1.562403) q[0];
sx q[0];
rz(-0.087455672) q[0];
rz(-pi) q[1];
rz(1.8553472) q[2];
sx q[2];
rz(-1.3853405) q[2];
sx q[2];
rz(-2.4539349) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0737887) q[1];
sx q[1];
rz(-1.5556591) q[1];
sx q[1];
rz(-0.10555663) q[1];
rz(-pi) q[2];
rz(-0.94957955) q[3];
sx q[3];
rz(-0.55286828) q[3];
sx q[3];
rz(2.051193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5856813) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(0.79706556) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(-1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109011) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(0.28924334) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(1.8274868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.700456) q[0];
sx q[0];
rz(-1.4608129) q[0];
sx q[0];
rz(-2.3814047) q[0];
rz(-pi) q[1];
rz(2.3938789) q[2];
sx q[2];
rz(-1.8166102) q[2];
sx q[2];
rz(2.4112548) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5319547) q[1];
sx q[1];
rz(-1.4502118) q[1];
sx q[1];
rz(3.0871255) q[1];
x q[2];
rz(-2.833509) q[3];
sx q[3];
rz(-0.83762533) q[3];
sx q[3];
rz(1.777491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(-0.96380487) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(2.5949809) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8209936) q[0];
sx q[0];
rz(-0.59251596) q[0];
sx q[0];
rz(2.2792363) q[0];
rz(-pi) q[1];
rz(-2.390929) q[2];
sx q[2];
rz(-0.51807907) q[2];
sx q[2];
rz(2.2028365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4709028) q[1];
sx q[1];
rz(-2.7109475) q[1];
sx q[1];
rz(1.1501269) q[1];
rz(0.48874493) q[3];
sx q[3];
rz(-0.97526032) q[3];
sx q[3];
rz(2.3084156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(-1.6938422) q[2];
rz(2.0041806) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(2.494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(0.71722537) q[0];
rz(-1.2099129) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(-0.70770121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9428064) q[0];
sx q[0];
rz(-1.2278623) q[0];
sx q[0];
rz(0.85533157) q[0];
x q[1];
rz(-2.131358) q[2];
sx q[2];
rz(-2.8786504) q[2];
sx q[2];
rz(1.5514785) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51703875) q[1];
sx q[1];
rz(-2.3561764) q[1];
sx q[1];
rz(-2.024827) q[1];
rz(-1.0912861) q[3];
sx q[3];
rz(-0.67688739) q[3];
sx q[3];
rz(-0.53938473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5132961) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(-0.90325242) q[2];
rz(1.5385657) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(-0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(2.3256336) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(1.1076526) q[2];
sx q[2];
rz(-1.9911498) q[2];
sx q[2];
rz(-0.1291612) q[2];
rz(-2.49414) q[3];
sx q[3];
rz(-3.0734607) q[3];
sx q[3];
rz(-1.4510696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
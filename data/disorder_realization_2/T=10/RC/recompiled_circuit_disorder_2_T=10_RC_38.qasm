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
rz(-0.50049385) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252285) q[0];
sx q[0];
rz(-2.541757) q[0];
sx q[0];
rz(3.1347549) q[0];
rz(-pi) q[1];
rz(-1.4209461) q[2];
sx q[2];
rz(-1.4784001) q[2];
sx q[2];
rz(-0.42042755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1653633) q[1];
sx q[1];
rz(-0.86127087) q[1];
sx q[1];
rz(-0.43433365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1327098) q[3];
sx q[3];
rz(-1.7128908) q[3];
sx q[3];
rz(-0.70664584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3502675) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(1.9187437) q[2];
rz(1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-0.97035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37110776) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(1.0789385) q[0];
rz(1.7547912) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-0.66545495) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5011713) q[0];
sx q[0];
rz(-1.8542395) q[0];
sx q[0];
rz(2.6675176) q[0];
x q[1];
rz(-2.8198492) q[2];
sx q[2];
rz(-1.7330568) q[2];
sx q[2];
rz(-2.4010047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7484819) q[1];
sx q[1];
rz(-1.0903653) q[1];
sx q[1];
rz(-1.3039116) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5311702) q[3];
sx q[3];
rz(-0.70981662) q[3];
sx q[3];
rz(-0.94330793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60454303) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(0.075142168) q[2];
rz(-1.6710619) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(-2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(1.2930019) q[1];
sx q[1];
rz(-1.8483775) q[1];
sx q[1];
rz(-1.741515) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62832075) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(-0.4011641) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18685762) q[2];
sx q[2];
rz(-1.7041022) q[2];
sx q[2];
rz(1.8568298) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2565508) q[1];
sx q[1];
rz(-0.98837822) q[1];
sx q[1];
rz(-0.38862733) q[1];
x q[2];
rz(1.0812976) q[3];
sx q[3];
rz(-0.34265095) q[3];
sx q[3];
rz(0.12658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65163461) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(-1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.3477429) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(1.6963652) q[0];
rz(1.6943278) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(2.7935374) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2340654) q[0];
sx q[0];
rz(-2.4287927) q[0];
sx q[0];
rz(0.45760052) q[0];
rz(-pi) q[1];
rz(1.939417) q[2];
sx q[2];
rz(-0.39837948) q[2];
sx q[2];
rz(0.57952651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0749803) q[1];
sx q[1];
rz(-2.1790824) q[1];
sx q[1];
rz(0.10886701) q[1];
x q[2];
rz(-0.16163687) q[3];
sx q[3];
rz(-1.703754) q[3];
sx q[3];
rz(-2.0032351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41670123) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-0.42281881) q[2];
rz(-2.4041798) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(-0.92492217) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(-1.2984498) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2753678) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(0.84262459) q[0];
x q[1];
rz(-2.5351296) q[2];
sx q[2];
rz(-1.3242553) q[2];
sx q[2];
rz(-2.390887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4358208) q[1];
sx q[1];
rz(-0.37322361) q[1];
sx q[1];
rz(0.2758287) q[1];
x q[2];
rz(2.1634444) q[3];
sx q[3];
rz(-1.1453298) q[3];
sx q[3];
rz(-0.076171906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(2.6679664) q[2];
rz(0.099362699) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(-2.1437058) q[0];
rz(-2.2672794) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(-0.46674892) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3335553) q[0];
sx q[0];
rz(-1.3555962) q[0];
sx q[0];
rz(-0.67610418) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7919962) q[2];
sx q[2];
rz(-1.9890519) q[2];
sx q[2];
rz(-2.7115371) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14086831) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(-0.15111698) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7344597) q[3];
sx q[3];
rz(-1.4758849) q[3];
sx q[3];
rz(-1.6223736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2816887) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(-1.5768645) q[2];
rz(0.62670296) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108114) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(0.41982857) q[0];
rz(2.9201674) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(1.5931169) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80207434) q[0];
sx q[0];
rz(-3.0537362) q[0];
sx q[0];
rz(0.095803424) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9485547) q[2];
sx q[2];
rz(-1.8503354) q[2];
sx q[2];
rz(-2.2045731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.067804) q[1];
sx q[1];
rz(-1.5859335) q[1];
sx q[1];
rz(3.036036) q[1];
rz(-1.1057304) q[3];
sx q[3];
rz(-1.8814058) q[3];
sx q[3];
rz(0.066699337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(-1.7377724) q[2];
rz(-0.79706556) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.9246624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109011) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(2.8523493) q[0];
rz(0.62492433) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(1.8274868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24459141) q[0];
sx q[0];
rz(-2.3750711) q[0];
sx q[0];
rz(-2.9826829) q[0];
rz(-0.74771379) q[2];
sx q[2];
rz(-1.3249825) q[2];
sx q[2];
rz(-2.4112548) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.96771679) q[1];
sx q[1];
rz(-1.6248676) q[1];
sx q[1];
rz(1.6915583) q[1];
rz(-pi) q[2];
rz(-1.2460327) q[3];
sx q[3];
rz(-2.357558) q[3];
sx q[3];
rz(-1.8079545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(-1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(-1.8956986) q[0];
rz(-0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(-2.5949809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7738757) q[0];
sx q[0];
rz(-1.9426632) q[0];
sx q[0];
rz(1.0982151) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.941628) q[2];
sx q[2];
rz(-1.9413345) q[2];
sx q[2];
rz(1.759699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48601549) q[1];
sx q[1];
rz(-1.7421107) q[1];
sx q[1];
rz(-1.1737215) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96484465) q[3];
sx q[3];
rz(-0.75111872) q[3];
sx q[3];
rz(1.5918819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0299915) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(-1.4477504) q[2];
rz(-2.0041806) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(-2.494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(0.71722537) q[0];
rz(-1.9316797) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(2.4338914) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087698547) q[0];
sx q[0];
rz(-0.90488926) q[0];
sx q[0];
rz(-2.6997487) q[0];
rz(-pi) q[1];
rz(2.131358) q[2];
sx q[2];
rz(-0.26294225) q[2];
sx q[2];
rz(-1.5901142) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6245539) q[1];
sx q[1];
rz(-0.78541628) q[1];
sx q[1];
rz(-2.024827) q[1];
x q[2];
rz(2.7865949) q[3];
sx q[3];
rz(-0.98155752) q[3];
sx q[3];
rz(-1.1276576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(-2.2383402) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(-2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.306504) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(-2.6782398) q[2];
sx q[2];
rz(-1.9909161) q[2];
sx q[2];
rz(-1.9009895) q[2];
rz(0.054374183) q[3];
sx q[3];
rz(-1.6118703) q[3];
sx q[3];
rz(0.76606228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

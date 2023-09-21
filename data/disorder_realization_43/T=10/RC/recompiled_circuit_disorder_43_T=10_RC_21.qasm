OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(1.8068846) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061493) q[0];
sx q[0];
rz(-0.9978928) q[0];
sx q[0];
rz(1.2826305) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.072233) q[2];
sx q[2];
rz(-2.517759) q[2];
sx q[2];
rz(-1.2944348) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18260278) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3372041) q[3];
sx q[3];
rz(-1.2803004) q[3];
sx q[3];
rz(0.94798541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41539899) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(-2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(0.17856199) q[0];
rz(-1.3372955) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(-3.1352502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7550678) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(0.13271876) q[0];
rz(-pi) q[1];
rz(-0.96887529) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(2.7768163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9176863) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(-0.19885893) q[1];
rz(-pi) q[2];
rz(2.7553495) q[3];
sx q[3];
rz(-2.5442903) q[3];
sx q[3];
rz(-2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-0.87810278) q[2];
rz(-2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5851615) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(0.095741622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75671065) q[0];
sx q[0];
rz(-1.4178935) q[0];
sx q[0];
rz(-1.3923313) q[0];
rz(2.6355866) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(-2.4207052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84903753) q[1];
sx q[1];
rz(-2.0325066) q[1];
sx q[1];
rz(2.1828116) q[1];
rz(-pi) q[2];
rz(-2.5903969) q[3];
sx q[3];
rz(-1.6008198) q[3];
sx q[3];
rz(-0.72202819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0720955) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(-0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9469706) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(1.0338763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-2.3405502) q[0];
sx q[0];
rz(-0.69181504) q[0];
x q[1];
rz(-2.1448574) q[2];
sx q[2];
rz(-1.1166995) q[2];
sx q[2];
rz(0.10276375) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32767195) q[1];
sx q[1];
rz(-1.3797626) q[1];
sx q[1];
rz(-0.066285985) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28169607) q[3];
sx q[3];
rz(-1.8857737) q[3];
sx q[3];
rz(1.4178993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088257) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(2.72686) q[0];
rz(-1.746009) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(0.57410747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2643471) q[0];
sx q[0];
rz(-2.1640722) q[0];
sx q[0];
rz(0.1546774) q[0];
x q[1];
rz(0.026002361) q[2];
sx q[2];
rz(-2.4261195) q[2];
sx q[2];
rz(-2.0040087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2327323) q[1];
sx q[1];
rz(-2.6870011) q[1];
sx q[1];
rz(2.4652387) q[1];
rz(-0.90548924) q[3];
sx q[3];
rz(-1.4169766) q[3];
sx q[3];
rz(-2.820435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.3621618) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(-2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(2.4093157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1351654) q[0];
sx q[0];
rz(-0.26627243) q[0];
sx q[0];
rz(2.4977495) q[0];
x q[1];
rz(-0.90768355) q[2];
sx q[2];
rz(-2.216202) q[2];
sx q[2];
rz(-1.4484608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9975035) q[1];
sx q[1];
rz(-2.2170076) q[1];
sx q[1];
rz(1.8540107) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7888277) q[3];
sx q[3];
rz(-2.1302967) q[3];
sx q[3];
rz(0.41527173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7531062) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.1614655) q[2];
rz(2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-0.77004534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2609445) q[0];
sx q[0];
rz(-2.5136607) q[0];
sx q[0];
rz(2.8448366) q[0];
x q[1];
rz(0.40600834) q[2];
sx q[2];
rz(-2.1210665) q[2];
sx q[2];
rz(-1.5232435) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60939497) q[1];
sx q[1];
rz(-2.1708793) q[1];
sx q[1];
rz(0.9432015) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50290147) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(-2.6754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(-2.1264123) q[2];
rz(1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5148233) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-0.36639211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.217427) q[0];
sx q[0];
rz(-2.1806742) q[0];
sx q[0];
rz(-2.2141371) q[0];
rz(-pi) q[1];
rz(0.56717746) q[2];
sx q[2];
rz(-0.71225538) q[2];
sx q[2];
rz(-2.8380307) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.061325039) q[1];
sx q[1];
rz(-2.419201) q[1];
sx q[1];
rz(-0.25218833) q[1];
x q[2];
rz(1.4922769) q[3];
sx q[3];
rz(-2.588387) q[3];
sx q[3];
rz(1.6012524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1239132) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(-2.4813095) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(0.79137897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521097) q[0];
sx q[0];
rz(-1.2533422) q[0];
sx q[0];
rz(2.2788458) q[0];
x q[1];
rz(-2.0558526) q[2];
sx q[2];
rz(-1.1527449) q[2];
sx q[2];
rz(0.56318356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3190805) q[1];
sx q[1];
rz(-1.1573536) q[1];
sx q[1];
rz(-0.49088571) q[1];
x q[2];
rz(1.295624) q[3];
sx q[3];
rz(-2.5857946) q[3];
sx q[3];
rz(2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(0.3785454) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1191694) q[0];
sx q[0];
rz(-1.4403617) q[0];
sx q[0];
rz(2.2020257) q[0];
rz(-pi) q[1];
rz(-1.1353178) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.6053111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.32523649) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(-0.34646323) q[1];
x q[2];
rz(2.4139666) q[3];
sx q[3];
rz(-2.8907223) q[3];
sx q[3];
rz(0.48654702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-0.33622462) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(-2.5972988) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(0.33967321) q[2];
sx q[2];
rz(-0.85396955) q[2];
sx q[2];
rz(0.11243482) q[2];
rz(-2.1401134) q[3];
sx q[3];
rz(-1.3712728) q[3];
sx q[3];
rz(-1.2029592) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

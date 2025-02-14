OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71830463) q[0];
sx q[0];
rz(-0.77646065) q[0];
sx q[0];
rz(2.7756696) q[0];
rz(0.34691063) q[1];
sx q[1];
rz(-0.28471714) q[1];
sx q[1];
rz(-1.3887583) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66807085) q[0];
sx q[0];
rz(-1.6374303) q[0];
sx q[0];
rz(-2.9391143) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69268815) q[2];
sx q[2];
rz(-1.8028304) q[2];
sx q[2];
rz(-1.9197861) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9216278) q[1];
sx q[1];
rz(-2.014451) q[1];
sx q[1];
rz(0.43334623) q[1];
rz(-pi) q[2];
rz(-2.4284989) q[3];
sx q[3];
rz(-2.4131107) q[3];
sx q[3];
rz(-2.147004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1846788) q[2];
sx q[2];
rz(-2.2214486) q[2];
sx q[2];
rz(-0.93519768) q[2];
rz(0.30185559) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(-0.39366084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1644156) q[0];
sx q[0];
rz(-1.875832) q[0];
sx q[0];
rz(2.6432977) q[0];
rz(1.4277108) q[1];
sx q[1];
rz(-1.9352813) q[1];
sx q[1];
rz(1.0867585) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0181429) q[0];
sx q[0];
rz(-0.88844127) q[0];
sx q[0];
rz(-1.2641973) q[0];
rz(2.9478023) q[2];
sx q[2];
rz(-1.3089798) q[2];
sx q[2];
rz(1.4145535) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9683428) q[1];
sx q[1];
rz(-2.0054711) q[1];
sx q[1];
rz(1.7477075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.558775) q[3];
sx q[3];
rz(-1.88899) q[3];
sx q[3];
rz(2.1592922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9429417) q[2];
sx q[2];
rz(-0.31060878) q[2];
sx q[2];
rz(-1.3322213) q[2];
rz(-1.9940935) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(1.3326299) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20951095) q[0];
sx q[0];
rz(-1.7387583) q[0];
sx q[0];
rz(-2.984356) q[0];
rz(-1.2278185) q[1];
sx q[1];
rz(-0.20918748) q[1];
sx q[1];
rz(2.7195209) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76772375) q[0];
sx q[0];
rz(-1.4808624) q[0];
sx q[0];
rz(1.3923682) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8639925) q[2];
sx q[2];
rz(-1.7381845) q[2];
sx q[2];
rz(3.0996499) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0737368) q[1];
sx q[1];
rz(-1.658981) q[1];
sx q[1];
rz(-1.1884618) q[1];
rz(1.2534327) q[3];
sx q[3];
rz(-2.3471213) q[3];
sx q[3];
rz(1.0004071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.813039) q[2];
sx q[2];
rz(-2.1812794) q[2];
sx q[2];
rz(2.4075244) q[2];
rz(-2.0640533) q[3];
sx q[3];
rz(-0.68818337) q[3];
sx q[3];
rz(0.075210007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.5828736) q[0];
sx q[0];
rz(-1.0105157) q[0];
sx q[0];
rz(-3.1251113) q[0];
rz(3.0670498) q[1];
sx q[1];
rz(-0.45769474) q[1];
sx q[1];
rz(-0.25513908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82127386) q[0];
sx q[0];
rz(-2.6373568) q[0];
sx q[0];
rz(-2.8404499) q[0];
x q[1];
rz(-2.071923) q[2];
sx q[2];
rz(-0.6269031) q[2];
sx q[2];
rz(1.3581585) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.093016) q[1];
sx q[1];
rz(-0.79461351) q[1];
sx q[1];
rz(-0.62127415) q[1];
x q[2];
rz(-1.2665073) q[3];
sx q[3];
rz(-2.7767468) q[3];
sx q[3];
rz(-1.5675275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5269346) q[2];
sx q[2];
rz(-2.4444828) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(-0.091863306) q[3];
sx q[3];
rz(-1.2173165) q[3];
sx q[3];
rz(-0.45472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62234539) q[0];
sx q[0];
rz(-2.3212101) q[0];
sx q[0];
rz(-1.1118332) q[0];
rz(-0.19019292) q[1];
sx q[1];
rz(-0.88514248) q[1];
sx q[1];
rz(1.9817188) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4442795) q[0];
sx q[0];
rz(-0.34553465) q[0];
sx q[0];
rz(-1.4688501) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2091523) q[2];
sx q[2];
rz(-1.4529422) q[2];
sx q[2];
rz(1.3502075) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3768757) q[1];
sx q[1];
rz(-0.60827601) q[1];
sx q[1];
rz(2.6981955) q[1];
rz(-1.0489063) q[3];
sx q[3];
rz(-1.5318222) q[3];
sx q[3];
rz(-3.1119973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5708892) q[2];
sx q[2];
rz(-0.8231701) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(-2.0795836) q[3];
sx q[3];
rz(-1.2758723) q[3];
sx q[3];
rz(1.5084069) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.688711) q[0];
sx q[0];
rz(-1.0588366) q[0];
sx q[0];
rz(-1.3558615) q[0];
rz(0.52976766) q[1];
sx q[1];
rz(-2.3593088) q[1];
sx q[1];
rz(1.1518325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065581948) q[0];
sx q[0];
rz(-1.2154723) q[0];
sx q[0];
rz(1.6760582) q[0];
x q[1];
rz(-3.0462469) q[2];
sx q[2];
rz(-1.0311677) q[2];
sx q[2];
rz(-0.44440545) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7499654) q[1];
sx q[1];
rz(-0.4164975) q[1];
sx q[1];
rz(-0.49642621) q[1];
rz(-pi) q[2];
rz(-1.4532667) q[3];
sx q[3];
rz(-1.7835878) q[3];
sx q[3];
rz(-0.26373395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4749703) q[2];
sx q[2];
rz(-2.5002561) q[2];
sx q[2];
rz(0.46710157) q[2];
rz(2.1304255) q[3];
sx q[3];
rz(-2.7979388) q[3];
sx q[3];
rz(1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2924627) q[0];
sx q[0];
rz(-0.15501538) q[0];
sx q[0];
rz(2.9169061) q[0];
rz(2.977773) q[1];
sx q[1];
rz(-2.8024709) q[1];
sx q[1];
rz(-0.64635578) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3098329) q[0];
sx q[0];
rz(-1.0467967) q[0];
sx q[0];
rz(2.8519467) q[0];
rz(-pi) q[1];
rz(-1.24794) q[2];
sx q[2];
rz(-1.3636929) q[2];
sx q[2];
rz(2.8139909) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.028370628) q[1];
sx q[1];
rz(-2.592855) q[1];
sx q[1];
rz(0.1446677) q[1];
x q[2];
rz(0.32352738) q[3];
sx q[3];
rz(-0.40641847) q[3];
sx q[3];
rz(2.714963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15567638) q[2];
sx q[2];
rz(-1.4480269) q[2];
sx q[2];
rz(-0.4099561) q[2];
rz(-0.83034596) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(-1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.603867) q[0];
sx q[0];
rz(-0.69320885) q[0];
sx q[0];
rz(-0.24630462) q[0];
rz(-2.314997) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(-0.26396096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1136341) q[0];
sx q[0];
rz(-1.6662681) q[0];
sx q[0];
rz(-3.1343293) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6648681) q[2];
sx q[2];
rz(-0.84273224) q[2];
sx q[2];
rz(1.6260894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3726793) q[1];
sx q[1];
rz(-2.2376334) q[1];
sx q[1];
rz(0.36075488) q[1];
rz(-pi) q[2];
x q[2];
rz(2.554515) q[3];
sx q[3];
rz(-1.9403321) q[3];
sx q[3];
rz(2.0943506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10494122) q[2];
sx q[2];
rz(-1.8737326) q[2];
sx q[2];
rz(2.4264753) q[2];
rz(-0.07130833) q[3];
sx q[3];
rz(-1.5846059) q[3];
sx q[3];
rz(0.74688545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.0095373) q[0];
sx q[0];
rz(-1.4893463) q[0];
sx q[0];
rz(-0.43553964) q[0];
rz(-2.9938193) q[1];
sx q[1];
rz(-1.7099893) q[1];
sx q[1];
rz(1.4685644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72318059) q[0];
sx q[0];
rz(-1.7505129) q[0];
sx q[0];
rz(0.43850684) q[0];
rz(-pi) q[1];
rz(0.69738241) q[2];
sx q[2];
rz(-1.5101523) q[2];
sx q[2];
rz(-1.1866807) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10994153) q[1];
sx q[1];
rz(-1.9579971) q[1];
sx q[1];
rz(-0.18020682) q[1];
x q[2];
rz(-0.40364175) q[3];
sx q[3];
rz(-1.3714681) q[3];
sx q[3];
rz(-2.2011873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74464166) q[2];
sx q[2];
rz(-1.3775185) q[2];
sx q[2];
rz(2.2486539) q[2];
rz(-1.8630155) q[3];
sx q[3];
rz(-1.9847001) q[3];
sx q[3];
rz(1.4634092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.29499149) q[0];
sx q[0];
rz(-0.33116594) q[0];
sx q[0];
rz(-2.5355329) q[0];
rz(0.010295708) q[1];
sx q[1];
rz(-2.1538487) q[1];
sx q[1];
rz(-3.0616679) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.851136) q[0];
sx q[0];
rz(-2.5154167) q[0];
sx q[0];
rz(-1.1315704) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84236161) q[2];
sx q[2];
rz(-1.7879221) q[2];
sx q[2];
rz(2.6922243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8680758) q[1];
sx q[1];
rz(-0.40765992) q[1];
sx q[1];
rz(-3.0855387) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3402088) q[3];
sx q[3];
rz(-1.7713431) q[3];
sx q[3];
rz(1.6236562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0004878) q[2];
sx q[2];
rz(-0.92685574) q[2];
sx q[2];
rz(1.0151939) q[2];
rz(-2.5403533) q[3];
sx q[3];
rz(-0.70178086) q[3];
sx q[3];
rz(-2.0950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4973608) q[0];
sx q[0];
rz(-1.624122) q[0];
sx q[0];
rz(-2.4540785) q[0];
rz(-0.58912206) q[1];
sx q[1];
rz(-2.4644869) q[1];
sx q[1];
rz(-0.45225515) q[1];
rz(0.44519201) q[2];
sx q[2];
rz(-1.4839076) q[2];
sx q[2];
rz(-2.1640726) q[2];
rz(-1.6986871) q[3];
sx q[3];
rz(-0.75027942) q[3];
sx q[3];
rz(1.1269509) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

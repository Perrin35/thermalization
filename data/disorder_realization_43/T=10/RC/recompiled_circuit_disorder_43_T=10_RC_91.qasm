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
rz(-1.3347081) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061493) q[0];
sx q[0];
rz(-2.1436999) q[0];
sx q[0];
rz(1.8589622) q[0];
rz(2.072233) q[2];
sx q[2];
rz(-2.517759) q[2];
sx q[2];
rz(1.2944348) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5648956) q[1];
sx q[1];
rz(-1.4370059) q[1];
sx q[1];
rz(2.1810075) q[1];
rz(-0.8043886) q[3];
sx q[3];
rz(-1.2803004) q[3];
sx q[3];
rz(0.94798541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41539899) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(-0.78380084) q[2];
rz(0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(1.3403085) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587104) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(-3.1352502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3621688) q[0];
sx q[0];
rz(-0.2340901) q[0];
sx q[0];
rz(-0.97658821) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2750521) q[2];
sx q[2];
rz(-2.4079977) q[2];
sx q[2];
rz(-1.3726335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9176863) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(-0.19885893) q[1];
rz(-0.56224058) q[3];
sx q[3];
rz(-1.3573109) q[3];
sx q[3];
rz(-1.436304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(-3.1306144) q[0];
rz(-2.7745461) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(-3.045851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.028713) q[0];
sx q[0];
rz(-2.9071147) q[0];
sx q[0];
rz(2.2857091) q[0];
rz(-0.7272561) q[2];
sx q[2];
rz(-2.5033592) q[2];
sx q[2];
rz(-0.22932316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0248191) q[1];
sx q[1];
rz(-1.0305335) q[1];
sx q[1];
rz(-2.5953672) q[1];
x q[2];
rz(0.057283244) q[3];
sx q[3];
rz(-2.5896642) q[3];
sx q[3];
rz(0.79997593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(2.3068008) q[2];
rz(0.21162027) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(1.0338763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945275) q[0];
sx q[0];
rz(-0.80104242) q[0];
sx q[0];
rz(0.69181504) q[0];
rz(-0.83861645) q[2];
sx q[2];
rz(-0.71574434) q[2];
sx q[2];
rz(2.0640304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8139207) q[1];
sx q[1];
rz(-1.3797626) q[1];
sx q[1];
rz(-0.066285985) q[1];
rz(-2.8598966) q[3];
sx q[3];
rz(-1.8857737) q[3];
sx q[3];
rz(-1.7236934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(1.3467849) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(-0.41473266) q[0];
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
rz(-0.87724553) q[0];
sx q[0];
rz(-2.1640722) q[0];
sx q[0];
rz(0.1546774) q[0];
rz(3.1155903) q[2];
sx q[2];
rz(-0.71547316) q[2];
sx q[2];
rz(-2.0040087) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9088604) q[1];
sx q[1];
rz(-0.4545916) q[1];
sx q[1];
rz(2.4652387) q[1];
x q[2];
rz(0.19458171) q[3];
sx q[3];
rz(-2.2268725) q[3];
sx q[3];
rz(-1.3692828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(0.73227698) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4740144) q[0];
sx q[0];
rz(-1.7828373) q[0];
sx q[0];
rz(1.408512) q[0];
rz(-pi) q[1];
rz(-2.3789669) q[2];
sx q[2];
rz(-1.0566933) q[2];
sx q[2];
rz(-2.5800173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9975035) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(-1.8540107) q[1];
rz(-pi) q[2];
rz(-2.5712588) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(-2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-0.86507559) q[1];
sx q[1];
rz(-0.77004534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89966398) q[0];
sx q[0];
rz(-0.97424346) q[0];
sx q[0];
rz(-1.7799737) q[0];
rz(-pi) q[1];
rz(2.1595575) q[2];
sx q[2];
rz(-1.9141478) q[2];
sx q[2];
rz(2.8729168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60939497) q[1];
sx q[1];
rz(-2.1708793) q[1];
sx q[1];
rz(-2.1983912) q[1];
rz(-pi) q[2];
rz(-1.2392063) q[3];
sx q[3];
rz(-2.050403) q[3];
sx q[3];
rz(0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4975171) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62676936) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(0.24169895) q[0];
rz(-0.73879755) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(-0.36639211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05224932) q[0];
sx q[0];
rz(-1.0567259) q[0];
sx q[0];
rz(-0.71787562) q[0];
rz(1.1364469) q[2];
sx q[2];
rz(-0.98698101) q[2];
sx q[2];
rz(2.7455612) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7494292) q[1];
sx q[1];
rz(-2.2656419) q[1];
sx q[1];
rz(1.7872582) q[1];
x q[2];
rz(3.0931926) q[3];
sx q[3];
rz(-2.1221005) q[3];
sx q[3];
rz(-1.5090404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1239132) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(0.29385847) q[2];
rz(-3.1270694) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(-2.4709539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83207399) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-0.88395399) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(2.3502137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732849) q[0];
sx q[0];
rz(-2.377016) q[0];
sx q[0];
rz(-2.0386049) q[0];
x q[1];
rz(-1.0857401) q[2];
sx q[2];
rz(-1.1527449) q[2];
sx q[2];
rz(-0.56318356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8225122) q[1];
sx q[1];
rz(-1.984239) q[1];
sx q[1];
rz(-2.6507069) q[1];
rz(-pi) q[2];
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
rz(-3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(1.0212612) q[2];
rz(-2.7630473) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-2.5436201) q[3];
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
rz(-3.1148949) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(0.60780203) q[0];
rz(2.9027477) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35683435) q[0];
sx q[0];
rz(-0.94576242) q[0];
sx q[0];
rz(-2.9805095) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1353178) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.5362816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8163562) q[1];
sx q[1];
rz(-1.163985) q[1];
sx q[1];
rz(2.7951294) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4139666) q[3];
sx q[3];
rz(-0.25087038) q[3];
sx q[3];
rz(2.6550456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-0.33622462) q[2];
rz(-2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(-2.5972988) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(-0.33967321) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(1.2119157) q[3];
sx q[3];
rz(-0.59960312) q[3];
sx q[3];
rz(-3.0740769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
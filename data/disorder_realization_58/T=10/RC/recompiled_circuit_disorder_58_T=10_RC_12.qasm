OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(2.7352754) q[0];
sx q[0];
rz(8.6046594) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93847371) q[0];
sx q[0];
rz(-1.9460558) q[0];
sx q[0];
rz(1.9012326) q[0];
rz(-pi) q[1];
rz(2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(0.29104656) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1931397) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(-0.83955168) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9740231) q[3];
sx q[3];
rz(-1.1551757) q[3];
sx q[3];
rz(2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(-3.0214156) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08081089) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90601901) q[0];
sx q[0];
rz(-0.86750194) q[0];
sx q[0];
rz(-2.2554382) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15906449) q[2];
sx q[2];
rz(-1.1317562) q[2];
sx q[2];
rz(-0.95869267) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8596526) q[1];
sx q[1];
rz(-1.167206) q[1];
sx q[1];
rz(-0.90359009) q[1];
rz(2.2585906) q[3];
sx q[3];
rz(-1.6371173) q[3];
sx q[3];
rz(-1.4955213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-0.32989311) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(2.0842016) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43305106) q[0];
sx q[0];
rz(-1.9110702) q[0];
sx q[0];
rz(0.021854594) q[0];
x q[1];
rz(2.6234987) q[2];
sx q[2];
rz(-2.7535541) q[2];
sx q[2];
rz(2.5990017) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96506572) q[1];
sx q[1];
rz(-0.9170734) q[1];
sx q[1];
rz(-0.42018004) q[1];
rz(-0.40773817) q[3];
sx q[3];
rz(-0.45555112) q[3];
sx q[3];
rz(1.7689442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.7827591) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(-0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7968984) q[0];
sx q[0];
rz(-1.7555153) q[0];
sx q[0];
rz(3.0152507) q[0];
rz(-pi) q[1];
rz(-0.16343127) q[2];
sx q[2];
rz(-1.391841) q[2];
sx q[2];
rz(0.93726678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.0488102) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(0.55744967) q[1];
x q[2];
rz(-1.0159675) q[3];
sx q[3];
rz(-1.8053683) q[3];
sx q[3];
rz(-2.6382584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-0.20544927) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.8849461) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(-0.35476312) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1629099) q[0];
sx q[0];
rz(-1.9918348) q[0];
sx q[0];
rz(2.1091503) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9734128) q[2];
sx q[2];
rz(-2.4390704) q[2];
sx q[2];
rz(1.0934193) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9310589) q[1];
sx q[1];
rz(-1.7909044) q[1];
sx q[1];
rz(1.348043) q[1];
rz(-pi) q[2];
rz(-0.050458126) q[3];
sx q[3];
rz(-2.7726463) q[3];
sx q[3];
rz(-1.7794533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(2.1441377) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31842445) q[0];
sx q[0];
rz(-2.423375) q[0];
sx q[0];
rz(1.3658001) q[0];
rz(0.41001292) q[2];
sx q[2];
rz(-0.7747246) q[2];
sx q[2];
rz(1.1020401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4221229) q[1];
sx q[1];
rz(-2.649253) q[1];
sx q[1];
rz(-1.0455529) q[1];
x q[2];
rz(-0.84647471) q[3];
sx q[3];
rz(-1.9123565) q[3];
sx q[3];
rz(-2.4404756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(2.690199) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.9394978) q[3];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-0.62414449) q[0];
rz(1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543024) q[0];
sx q[0];
rz(-0.76061941) q[0];
sx q[0];
rz(1.8421696) q[0];
x q[1];
rz(2.9250547) q[2];
sx q[2];
rz(-0.29007402) q[2];
sx q[2];
rz(1.1494344) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2144199) q[1];
sx q[1];
rz(-1.629429) q[1];
sx q[1];
rz(0.10417948) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11717637) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43341407) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(0.60837778) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(2.1571295) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0122089) q[0];
sx q[0];
rz(-2.1166271) q[0];
sx q[0];
rz(1.6145541) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7266399) q[2];
sx q[2];
rz(-2.3256362) q[2];
sx q[2];
rz(-2.6331537) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22673785) q[1];
sx q[1];
rz(-1.2054772) q[1];
sx q[1];
rz(1.5772596) q[1];
rz(-pi) q[2];
rz(1.2582614) q[3];
sx q[3];
rz(-1.1141889) q[3];
sx q[3];
rz(-1.4708335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-0.53517503) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64131367) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(2.2156782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029322421) q[0];
sx q[0];
rz(-1.7824714) q[0];
sx q[0];
rz(-2.9979343) q[0];
x q[1];
rz(-0.7943031) q[2];
sx q[2];
rz(-2.1199193) q[2];
sx q[2];
rz(-1.6910451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.029286413) q[1];
sx q[1];
rz(-2.6767113) q[1];
sx q[1];
rz(-1.2390562) q[1];
x q[2];
rz(1.4521452) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(-0.53393902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(-1.013247) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17908827) q[0];
sx q[0];
rz(-2.8694186) q[0];
sx q[0];
rz(1.0286691) q[0];
rz(0.11864885) q[2];
sx q[2];
rz(-0.35602202) q[2];
sx q[2];
rz(-2.7729386) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80876795) q[1];
sx q[1];
rz(-1.9908394) q[1];
sx q[1];
rz(-1.6767477) q[1];
rz(-2.3771044) q[3];
sx q[3];
rz(-2.0937243) q[3];
sx q[3];
rz(0.58512277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(-1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(0.57327523) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(-1.7309932) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
rz(-3.0782386) q[3];
sx q[3];
rz(-0.8739211) q[3];
sx q[3];
rz(-2.7779761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
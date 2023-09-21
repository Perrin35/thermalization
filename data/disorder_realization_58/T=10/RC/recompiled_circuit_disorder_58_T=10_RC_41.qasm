OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18618628) q[0];
sx q[0];
rz(-0.49476981) q[0];
sx q[0];
rz(2.4525989) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4029998) q[2];
sx q[2];
rz(-1.9361708) q[2];
sx q[2];
rz(1.5942758) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0349866) q[1];
sx q[1];
rz(-2.258746) q[1];
sx q[1];
rz(0.41253849) q[1];
x q[2];
rz(1.1675695) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.43710199) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-2.1531838) q[1];
sx q[1];
rz(2.8149014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17900285) q[0];
sx q[0];
rz(-1.0674745) q[0];
sx q[0];
rz(2.3110564) q[0];
x q[1];
rz(1.8961043) q[2];
sx q[2];
rz(-2.6763958) q[2];
sx q[2];
rz(-2.543769) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5888728) q[1];
sx q[1];
rz(-0.96541222) q[1];
sx q[1];
rz(-0.49792109) q[1];
x q[2];
rz(3.0558415) q[3];
sx q[3];
rz(-2.2567856) q[3];
sx q[3];
rz(-3.0119197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5102753) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(2.8116995) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(-2.0842016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1304504) q[0];
sx q[0];
rz(-1.550195) q[0];
sx q[0];
rz(-1.9111454) q[0];
rz(-1.3710652) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(-2.0470326) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2705546) q[1];
sx q[1];
rz(-1.241031) q[1];
sx q[1];
rz(-2.2689181) q[1];
rz(-pi) q[2];
rz(2.7189414) q[3];
sx q[3];
rz(-1.395441) q[3];
sx q[3];
rz(0.56817504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3588336) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(-1.3624297) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(2.9761956) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1918614) q[0];
sx q[0];
rz(-0.22338578) q[0];
sx q[0];
rz(0.97747691) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(2.2043259) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71388984) q[1];
sx q[1];
rz(-2.1117003) q[1];
sx q[1];
rz(1.1869903) q[1];
x q[2];
rz(2.1256251) q[3];
sx q[3];
rz(-1.8053683) q[3];
sx q[3];
rz(0.50333422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(-1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(0.23194557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35271586) q[0];
sx q[0];
rz(-2.0577601) q[0];
sx q[0];
rz(-2.6608174) q[0];
rz(-pi) q[1];
rz(-2.4460692) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(-0.34851375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9310589) q[1];
sx q[1];
rz(-1.3506883) q[1];
sx q[1];
rz(1.7935497) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5902953) q[3];
sx q[3];
rz(-1.9392506) q[3];
sx q[3];
rz(-1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(-1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0971138) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(2.1441377) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0925804) q[0];
sx q[0];
rz(-0.870734) q[0];
sx q[0];
rz(2.9655365) q[0];
rz(-pi) q[1];
rz(1.1987655) q[2];
sx q[2];
rz(-2.2673006) q[2];
sx q[2];
rz(-0.55559413) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84038489) q[1];
sx q[1];
rz(-1.1494698) q[1];
sx q[1];
rz(-0.2627443) q[1];
rz(-pi) q[2];
rz(0.44316767) q[3];
sx q[3];
rz(-2.245095) q[3];
sx q[3];
rz(-2.5600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.30615) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(0.61378941) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543024) q[0];
sx q[0];
rz(-2.3809732) q[0];
sx q[0];
rz(-1.2994231) q[0];
x q[1];
rz(2.8579312) q[2];
sx q[2];
rz(-1.6322871) q[2];
sx q[2];
rz(0.21360699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1544513) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(-2.6277072) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0244163) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(1.5069435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(-1.6237367) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(-2.1571295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9280633) q[0];
sx q[0];
rz(-0.54740471) q[0];
sx q[0];
rz(3.0696966) q[0];
x q[1];
rz(0.41495277) q[2];
sx q[2];
rz(-2.3256362) q[2];
sx q[2];
rz(0.50843898) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8967646) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(-3.1246964) q[1];
rz(-pi) q[2];
rz(-1.2582614) q[3];
sx q[3];
rz(-2.0274037) q[3];
sx q[3];
rz(1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(-0.92591441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5092963) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(-2.1584828) q[0];
rz(-2.2885867) q[2];
sx q[2];
rz(-2.2249613) q[2];
sx q[2];
rz(-2.5329563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8402108) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(-1.128003) q[1];
x q[2];
rz(2.9374398) q[3];
sx q[3];
rz(-1.4545868) q[3];
sx q[3];
rz(-2.1287363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(0.53608981) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(-2.1283456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17908827) q[0];
sx q[0];
rz(-2.8694186) q[0];
sx q[0];
rz(-1.0286691) q[0];
x q[1];
rz(3.0229438) q[2];
sx q[2];
rz(-2.7855706) q[2];
sx q[2];
rz(0.3686541) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0638748) q[1];
sx q[1];
rz(-2.7091654) q[1];
sx q[1];
rz(-2.9090911) q[1];
x q[2];
rz(-0.76448828) q[3];
sx q[3];
rz(-1.0478684) q[3];
sx q[3];
rz(0.58512277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41291819) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(-1.7309932) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
rz(0.87293252) q[3];
sx q[3];
rz(-1.5222266) q[3];
sx q[3];
rz(1.8937187) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

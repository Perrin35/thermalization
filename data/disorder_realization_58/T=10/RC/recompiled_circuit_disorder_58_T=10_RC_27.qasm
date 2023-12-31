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
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75738534) q[0];
sx q[0];
rz(-1.2641347) q[0];
sx q[0];
rz(2.7469809) q[0];
rz(0.51672207) q[2];
sx q[2];
rz(-0.80846723) q[2];
sx q[2];
rz(0.39743039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1931397) q[1];
sx q[1];
rz(-1.2558736) q[1];
sx q[1];
rz(2.302041) q[1];
rz(-0.72676267) q[3];
sx q[3];
rz(-2.5708377) q[3];
sx q[3];
rz(1.6970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(-3.0214156) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(-0.91896287) q[3];
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
rz(0.08081089) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(0.3266913) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17900285) q[0];
sx q[0];
rz(-1.0674745) q[0];
sx q[0];
rz(0.83053629) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9825282) q[2];
sx q[2];
rz(-2.0098364) q[2];
sx q[2];
rz(-2.1829) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82636278) q[1];
sx q[1];
rz(-2.3781207) q[1];
sx q[1];
rz(2.1748494) q[1];
rz(-pi) q[2];
rz(1.6750402) q[3];
sx q[3];
rz(-2.4511271) q[3];
sx q[3];
rz(-0.0052099293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-2.0842016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43305106) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(-0.021854594) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6234987) q[2];
sx q[2];
rz(-2.7535541) q[2];
sx q[2];
rz(-2.5990017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8101465) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(2.0600832) q[1];
x q[2];
rz(-0.40773817) q[3];
sx q[3];
rz(-0.45555112) q[3];
sx q[3];
rz(-1.3726485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3588336) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(2.6065361) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(2.9761956) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3446942) q[0];
sx q[0];
rz(-1.3860774) q[0];
sx q[0];
rz(-0.12634191) q[0];
rz(1.7521162) q[2];
sx q[2];
rz(-1.7315947) q[2];
sx q[2];
rz(0.60418512) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0796697) q[1];
sx q[1];
rz(-1.2440145) q[1];
sx q[1];
rz(0.57481874) q[1];
x q[2];
rz(1.1449279) q[3];
sx q[3];
rz(-2.5440359) q[3];
sx q[3];
rz(1.4262517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(-2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797392) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-2.9096471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922798) q[0];
sx q[0];
rz(-2.4711907) q[0];
sx q[0];
rz(0.85286661) q[0];
rz(-0.69552341) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(-2.7930789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0143118) q[1];
sx q[1];
rz(-2.8297272) q[1];
sx q[1];
rz(2.36253) q[1];
x q[2];
rz(-1.5902953) q[3];
sx q[3];
rz(-1.202342) q[3];
sx q[3];
rz(-1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.3034472) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(2.1525106) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049012262) q[0];
sx q[0];
rz(-0.870734) q[0];
sx q[0];
rz(0.17605619) q[0];
rz(-0.73156725) q[2];
sx q[2];
rz(-1.288207) q[2];
sx q[2];
rz(0.76991316) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62085405) q[1];
sx q[1];
rz(-1.3315017) q[1];
sx q[1];
rz(1.1362856) q[1];
rz(-pi) q[2];
rz(-2.0632083) q[3];
sx q[3];
rz(-2.3541854) q[3];
sx q[3];
rz(-1.9099964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(-0.4513936) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.9394978) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928771) q[0];
sx q[0];
rz(-2.2971417) q[0];
sx q[0];
rz(2.8918299) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5067528) q[2];
sx q[2];
rz(-1.2876858) q[2];
sx q[2];
rz(-1.7664906) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92717273) q[1];
sx q[1];
rz(-1.629429) q[1];
sx q[1];
rz(0.10417948) q[1];
rz(-pi) q[2];
rz(1.8011439) q[3];
sx q[3];
rz(-0.47603546) q[3];
sx q[3];
rz(1.7649094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(1.9667352) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(1.5178559) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(-0.98446313) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280633) q[0];
sx q[0];
rz(-2.5941879) q[0];
sx q[0];
rz(3.0696966) q[0];
rz(-pi) q[1];
rz(0.41495277) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(2.6331537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8967646) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(-pi) q[2];
rz(-0.55926178) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(0.53517503) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(-2.2156782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5092963) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(2.1584828) q[0];
x q[1];
rz(-2.4326153) q[2];
sx q[2];
rz(-2.211494) q[2];
sx q[2];
rz(-0.35352732) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8402108) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(-1.128003) q[1];
x q[2];
rz(2.6191606) q[3];
sx q[3];
rz(-2.9070832) q[3];
sx q[3];
rz(1.0684551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(1.4887811) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(-1.6292054) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(2.1283456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(-0.14302111) q[0];
x q[1];
rz(-0.35372325) q[2];
sx q[2];
rz(-1.612066) q[2];
sx q[2];
rz(1.8281787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0638748) q[1];
sx q[1];
rz(-0.43242726) q[1];
sx q[1];
rz(-2.9090911) q[1];
x q[2];
rz(0.69443955) q[3];
sx q[3];
rz(-2.246292) q[3];
sx q[3];
rz(-1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-2.3821793) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(-0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(2.5683174) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(1.4105994) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
rz(2.2686601) q[3];
sx q[3];
rz(-1.6193661) q[3];
sx q[3];
rz(-1.247874) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

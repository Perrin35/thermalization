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
rz(-2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18618628) q[0];
sx q[0];
rz(-2.6468228) q[0];
sx q[0];
rz(-0.68899378) q[0];
x q[1];
rz(-1.0933502) q[2];
sx q[2];
rz(-2.2508143) q[2];
sx q[2];
rz(-0.29104656) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9484529) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(0.83955168) q[1];
rz(-pi) q[2];
rz(-0.72676267) q[3];
sx q[3];
rz(-0.57075497) q[3];
sx q[3];
rz(-1.6970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.08081089) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-2.2297915) q[0];
rz(-0.78951019) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3347496) q[0];
sx q[0];
rz(-0.93870367) q[0];
sx q[0];
rz(-0.64081162) q[0];
x q[1];
rz(-1.2454883) q[2];
sx q[2];
rz(-2.6763958) q[2];
sx q[2];
rz(0.59782366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8596526) q[1];
sx q[1];
rz(-1.9743866) q[1];
sx q[1];
rz(2.2380026) q[1];
x q[2];
rz(-1.6750402) q[3];
sx q[3];
rz(-2.4511271) q[3];
sx q[3];
rz(-3.1363827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085416) q[0];
sx q[0];
rz(-1.9110702) q[0];
sx q[0];
rz(3.1197381) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6234987) q[2];
sx q[2];
rz(-0.38803852) q[2];
sx q[2];
rz(0.54259091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33144618) q[1];
sx q[1];
rz(-2.3815037) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
rz(1.7626761) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9388158) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(1.3846272) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7521162) q[2];
sx q[2];
rz(-1.7315947) q[2];
sx q[2];
rz(2.5374075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0927825) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(-2.584143) q[1];
rz(2.8675219) q[3];
sx q[3];
rz(-1.032853) q[3];
sx q[3];
rz(-2.2171973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(0.20544927) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(-1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(0.35476312) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(0.23194557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1922798) q[0];
sx q[0];
rz(-2.4711907) q[0];
sx q[0];
rz(-0.85286661) q[0];
x q[1];
rz(1.4300214) q[2];
sx q[2];
rz(-0.88015926) q[2];
sx q[2];
rz(1.8292793) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8307454) q[1];
sx q[1];
rz(-1.3535045) q[1];
sx q[1];
rz(2.9160935) q[1];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(-2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(0.99745497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049012262) q[0];
sx q[0];
rz(-2.2708587) q[0];
sx q[0];
rz(0.17605619) q[0];
x q[1];
rz(-0.41001292) q[2];
sx q[2];
rz(-0.7747246) q[2];
sx q[2];
rz(-1.1020401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5207386) q[1];
sx q[1];
rz(-1.3315017) q[1];
sx q[1];
rz(-1.1362856) q[1];
x q[2];
rz(1.0783844) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(1.9099964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(-0.4513936) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-2.5278032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1873338) q[0];
sx q[0];
rz(-1.7566534) q[0];
sx q[0];
rz(2.3128187) q[0];
rz(-pi) q[1];
rz(-2.8579312) q[2];
sx q[2];
rz(-1.5093056) q[2];
sx q[2];
rz(0.21360699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.504096) q[1];
sx q[1];
rz(-1.6747961) q[1];
sx q[1];
rz(-1.5118447) q[1];
rz(1.1055787) q[3];
sx q[3];
rz(-1.4659766) q[3];
sx q[3];
rz(-3.1302111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(-2.5332149) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-2.6532145) q[0];
rz(1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0122089) q[0];
sx q[0];
rz(-1.0249656) q[0];
sx q[0];
rz(1.5270385) q[0];
rz(2.7266399) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(-2.6331537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3417495) q[1];
sx q[1];
rz(-1.5768331) q[1];
sx q[1];
rz(0.36532613) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47655388) q[3];
sx q[3];
rz(-1.8503975) q[3];
sx q[3];
rz(-0.24147803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-0.53517503) q[2];
rz(-1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(2.1323269) q[3];
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
rz(-0.64131367) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(-2.7507239) q[1];
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
rz(0.029322421) q[0];
sx q[0];
rz(-1.3591213) q[0];
sx q[0];
rz(-0.14365833) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3472896) q[2];
sx q[2];
rz(-2.1199193) q[2];
sx q[2];
rz(1.4505475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.029286413) q[1];
sx q[1];
rz(-2.6767113) q[1];
sx q[1];
rz(1.9025365) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52243201) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(-1.0684551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9459076) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(2.1283456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.917406) q[0];
sx q[0];
rz(-1.4316443) q[0];
sx q[0];
rz(1.3361206) q[0];
rz(0.11864885) q[2];
sx q[2];
rz(-2.7855706) q[2];
sx q[2];
rz(2.7729386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0777178) q[1];
sx q[1];
rz(-2.7091654) q[1];
sx q[1];
rz(-0.23250154) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89684422) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(2.602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.2686601) q[3];
sx q[3];
rz(-1.5222266) q[3];
sx q[3];
rz(1.8937187) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

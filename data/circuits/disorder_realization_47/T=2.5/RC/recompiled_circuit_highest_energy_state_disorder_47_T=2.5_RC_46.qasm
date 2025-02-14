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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(-0.43098488) q[0];
rz(0.35297901) q[1];
sx q[1];
rz(-1.9184435) q[1];
sx q[1];
rz(-2.7117742) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5777464) q[0];
sx q[0];
rz(-1.0809076) q[0];
sx q[0];
rz(-1.6728982) q[0];
x q[1];
rz(1.5360898) q[2];
sx q[2];
rz(-0.43259753) q[2];
sx q[2];
rz(-2.2207867) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4861075) q[1];
sx q[1];
rz(-1.7799417) q[1];
sx q[1];
rz(0.36693348) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6453703) q[3];
sx q[3];
rz(-1.9093976) q[3];
sx q[3];
rz(-0.34211516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0777883) q[2];
sx q[2];
rz(-0.52738515) q[2];
sx q[2];
rz(-2.1303614) q[2];
rz(-2.1003335) q[3];
sx q[3];
rz(-2.6285089) q[3];
sx q[3];
rz(0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1083199) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(2.394115) q[0];
rz(-0.29391089) q[1];
sx q[1];
rz(-2.0103318) q[1];
sx q[1];
rz(-0.25516137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9440112) q[0];
sx q[0];
rz(-1.7367878) q[0];
sx q[0];
rz(-0.98910467) q[0];
rz(0.5396217) q[2];
sx q[2];
rz(-2.789342) q[2];
sx q[2];
rz(-2.7045111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16431674) q[1];
sx q[1];
rz(-0.9292047) q[1];
sx q[1];
rz(-0.56820984) q[1];
x q[2];
rz(0.94631715) q[3];
sx q[3];
rz(-1.7229652) q[3];
sx q[3];
rz(0.63380235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.030563844) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(2.2500706) q[2];
rz(2.8920178) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(2.7454624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7340649) q[0];
sx q[0];
rz(-1.1363109) q[0];
sx q[0];
rz(-0.044128142) q[0];
rz(1.2491501) q[1];
sx q[1];
rz(-2.7154778) q[1];
sx q[1];
rz(2.6557907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51903782) q[0];
sx q[0];
rz(-1.7960894) q[0];
sx q[0];
rz(-1.5908498) q[0];
rz(-pi) q[1];
rz(-3.0032512) q[2];
sx q[2];
rz(-1.7560343) q[2];
sx q[2];
rz(-2.0092162) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0016077) q[1];
sx q[1];
rz(-2.1473653) q[1];
sx q[1];
rz(-0.92553161) q[1];
rz(-pi) q[2];
rz(-2.7445469) q[3];
sx q[3];
rz(-1.9730803) q[3];
sx q[3];
rz(-2.8232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32758632) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(0.69474727) q[2];
rz(2.5658549) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(-1.4076788) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95076743) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(2.8488391) q[0];
rz(-2.5726908) q[1];
sx q[1];
rz(-2.3141373) q[1];
sx q[1];
rz(-2.2946766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0305188) q[0];
sx q[0];
rz(-1.1630217) q[0];
sx q[0];
rz(-3.0570583) q[0];
x q[1];
rz(-0.90130536) q[2];
sx q[2];
rz(-1.5141103) q[2];
sx q[2];
rz(2.0079119) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6509622) q[1];
sx q[1];
rz(-1.3329026) q[1];
sx q[1];
rz(-1.3722653) q[1];
x q[2];
rz(-1.2388703) q[3];
sx q[3];
rz(-0.20461018) q[3];
sx q[3];
rz(1.0910891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7920821) q[2];
sx q[2];
rz(-1.258472) q[2];
sx q[2];
rz(-0.2087896) q[2];
rz(-2.4236692) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(-2.2010402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0402886) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(-2.8029602) q[0];
rz(-0.70308095) q[1];
sx q[1];
rz(-0.73892006) q[1];
sx q[1];
rz(-2.3032761) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95264065) q[0];
sx q[0];
rz(-1.8243655) q[0];
sx q[0];
rz(-2.6978561) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82413582) q[2];
sx q[2];
rz(-1.7172684) q[2];
sx q[2];
rz(0.15268385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5456713) q[1];
sx q[1];
rz(-1.537408) q[1];
sx q[1];
rz(-0.8707134) q[1];
x q[2];
rz(-1.1456212) q[3];
sx q[3];
rz(-2.1347649) q[3];
sx q[3];
rz(-2.8759046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4579953) q[2];
sx q[2];
rz(-2.0532875) q[2];
sx q[2];
rz(2.1387157) q[2];
rz(2.9966127) q[3];
sx q[3];
rz(-0.093791157) q[3];
sx q[3];
rz(3.0785479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1231287) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(-0.96328324) q[0];
rz(-2.9604498) q[1];
sx q[1];
rz(-2.2434442) q[1];
sx q[1];
rz(-1.9021665) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9299663) q[0];
sx q[0];
rz(-0.74290448) q[0];
sx q[0];
rz(-2.5725098) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.507429) q[2];
sx q[2];
rz(-2.6118738) q[2];
sx q[2];
rz(1.5629753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9037348) q[1];
sx q[1];
rz(-1.9457726) q[1];
sx q[1];
rz(-1.9720248) q[1];
rz(-pi) q[2];
rz(0.12902) q[3];
sx q[3];
rz(-1.4203398) q[3];
sx q[3];
rz(2.4935088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88508254) q[2];
sx q[2];
rz(-0.91826597) q[2];
sx q[2];
rz(2.2086842) q[2];
rz(-1.4126623) q[3];
sx q[3];
rz(-1.4224854) q[3];
sx q[3];
rz(0.2093813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4466208) q[0];
sx q[0];
rz(-2.7734723) q[0];
sx q[0];
rz(2.3387261) q[0];
rz(2.39957) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-0.41025695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0502346) q[0];
sx q[0];
rz(-1.0438393) q[0];
sx q[0];
rz(1.9684674) q[0];
rz(-pi) q[1];
rz(0.76426498) q[2];
sx q[2];
rz(-1.4857123) q[2];
sx q[2];
rz(-2.6972023) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8150543) q[1];
sx q[1];
rz(-1.4966623) q[1];
sx q[1];
rz(1.733041) q[1];
rz(-pi) q[2];
rz(-1.0954082) q[3];
sx q[3];
rz(-2.1719526) q[3];
sx q[3];
rz(2.9132995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2731169) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(-2.9203019) q[2];
rz(-2.7224329) q[3];
sx q[3];
rz(-2.6594888) q[3];
sx q[3];
rz(-2.6656849) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23685037) q[0];
sx q[0];
rz(-0.99822799) q[0];
sx q[0];
rz(-1.5297484) q[0];
rz(1.8086241) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(1.6882247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8771031) q[0];
sx q[0];
rz(-0.96142381) q[0];
sx q[0];
rz(0.86848702) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6906934) q[2];
sx q[2];
rz(-1.7603619) q[2];
sx q[2];
rz(-1.8360863) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0537117) q[1];
sx q[1];
rz(-2.2715927) q[1];
sx q[1];
rz(1.3739963) q[1];
rz(1.0681719) q[3];
sx q[3];
rz(-0.31538559) q[3];
sx q[3];
rz(-0.50704765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2754485) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(3.1308543) q[2];
rz(0.3206611) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(2.5519154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084932) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(0.47472111) q[0];
rz(2.8482598) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(-1.6620103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5005762) q[0];
sx q[0];
rz(-1.5307284) q[0];
sx q[0];
rz(-3.1047719) q[0];
x q[1];
rz(2.4655645) q[2];
sx q[2];
rz(-1.3814298) q[2];
sx q[2];
rz(2.2631136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.868728) q[1];
sx q[1];
rz(-1.1050779) q[1];
sx q[1];
rz(1.8293423) q[1];
rz(2.9607356) q[3];
sx q[3];
rz(-1.3279727) q[3];
sx q[3];
rz(2.9202785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.044363) q[2];
sx q[2];
rz(-1.2791415) q[2];
sx q[2];
rz(0.43107671) q[2];
rz(0.19032446) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(-0.88461191) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1650319) q[0];
sx q[0];
rz(-0.29104069) q[0];
sx q[0];
rz(-0.2188368) q[0];
rz(-0.95895514) q[1];
sx q[1];
rz(-0.7364277) q[1];
sx q[1];
rz(-1.3921907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0420211) q[0];
sx q[0];
rz(-2.9416658) q[0];
sx q[0];
rz(2.8801092) q[0];
rz(-0.31453374) q[2];
sx q[2];
rz(-0.43930231) q[2];
sx q[2];
rz(0.9010992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2710378) q[1];
sx q[1];
rz(-2.1907594) q[1];
sx q[1];
rz(1.2024173) q[1];
rz(-2.1851319) q[3];
sx q[3];
rz(-1.2117361) q[3];
sx q[3];
rz(2.7484855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86768156) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(1.8217746) q[2];
rz(-2.6093318) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(-1.5925315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.559448) q[0];
sx q[0];
rz(-1.5503217) q[0];
sx q[0];
rz(0.4009608) q[0];
rz(-0.22790146) q[1];
sx q[1];
rz(-1.0453929) q[1];
sx q[1];
rz(1.4452404) q[1];
rz(0.54351844) q[2];
sx q[2];
rz(-0.47558161) q[2];
sx q[2];
rz(-2.4408792) q[2];
rz(-2.1501272) q[3];
sx q[3];
rz(-1.4351063) q[3];
sx q[3];
rz(-1.2096005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

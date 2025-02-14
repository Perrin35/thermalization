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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0864705) q[0];
sx q[0];
rz(-1.4807379) q[0];
sx q[0];
rz(0.49205972) q[0];
rz(-pi) q[1];
rz(1.6055029) q[2];
sx q[2];
rz(-0.43259753) q[2];
sx q[2];
rz(2.2207867) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.9949405) q[1];
sx q[1];
rz(-1.9293679) q[1];
sx q[1];
rz(1.3472122) q[1];
rz(-pi) q[2];
rz(-0.49622239) q[3];
sx q[3];
rz(-1.9093976) q[3];
sx q[3];
rz(0.34211516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.063804403) q[2];
sx q[2];
rz(-0.52738515) q[2];
sx q[2];
rz(-1.0112313) q[2];
rz(-1.0412591) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0332727) q[0];
sx q[0];
rz(-2.1450873) q[0];
sx q[0];
rz(-2.394115) q[0];
rz(2.8476818) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(0.25516137) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6601488) q[0];
sx q[0];
rz(-2.143476) q[0];
sx q[0];
rz(-0.19788589) q[0];
x q[1];
rz(2.601971) q[2];
sx q[2];
rz(-2.789342) q[2];
sx q[2];
rz(-0.43708153) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9772759) q[1];
sx q[1];
rz(-0.9292047) q[1];
sx q[1];
rz(0.56820984) q[1];
rz(-0.94631715) q[3];
sx q[3];
rz(-1.7229652) q[3];
sx q[3];
rz(2.5077903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1110288) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(0.89152208) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(-2.7454624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4075277) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(0.044128142) q[0];
rz(1.2491501) q[1];
sx q[1];
rz(-2.7154778) q[1];
sx q[1];
rz(2.6557907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7120948) q[0];
sx q[0];
rz(-0.22616874) q[0];
sx q[0];
rz(0.087271376) q[0];
x q[1];
rz(0.13834149) q[2];
sx q[2];
rz(-1.3855583) q[2];
sx q[2];
rz(2.0092162) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9456182) q[1];
sx q[1];
rz(-2.3047119) q[1];
sx q[1];
rz(-0.7463781) q[1];
rz(-pi) q[2];
rz(-2.3084675) q[3];
sx q[3];
rz(-0.55748788) q[3];
sx q[3];
rz(-1.1380029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8140063) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(-2.4468454) q[2];
rz(-2.5658549) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(1.4076788) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1908252) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(-2.8488391) q[0];
rz(-2.5726908) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(-0.84691602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6349061) q[0];
sx q[0];
rz(-1.6483848) q[0];
sx q[0];
rz(-1.1617178) q[0];
rz(2.2402873) q[2];
sx q[2];
rz(-1.6274823) q[2];
sx q[2];
rz(1.1336807) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6509622) q[1];
sx q[1];
rz(-1.8086901) q[1];
sx q[1];
rz(1.7693273) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7645232) q[3];
sx q[3];
rz(-1.5045369) q[3];
sx q[3];
rz(2.3364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7920821) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(-0.2087896) q[2];
rz(-2.4236692) q[3];
sx q[3];
rz(-2.8585298) q[3];
sx q[3];
rz(2.2010402) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1013041) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(2.8029602) q[0];
rz(-2.4385117) q[1];
sx q[1];
rz(-0.73892006) q[1];
sx q[1];
rz(-0.83831659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6421239) q[0];
sx q[0];
rz(-1.1422061) q[0];
sx q[0];
rz(-1.2913676) q[0];
x q[1];
rz(1.3569068) q[2];
sx q[2];
rz(-0.75816064) q[2];
sx q[2];
rz(1.2616829) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.053239659) q[1];
sx q[1];
rz(-0.87118282) q[1];
sx q[1];
rz(3.0979473) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1456212) q[3];
sx q[3];
rz(-1.0068277) q[3];
sx q[3];
rz(-2.8759046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68359739) q[2];
sx q[2];
rz(-2.0532875) q[2];
sx q[2];
rz(-2.1387157) q[2];
rz(-0.14497997) q[3];
sx q[3];
rz(-0.093791157) q[3];
sx q[3];
rz(3.0785479) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.018464) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(-2.1783094) q[0];
rz(2.9604498) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(-1.9021665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9299663) q[0];
sx q[0];
rz(-2.3986882) q[0];
sx q[0];
rz(0.56908281) q[0];
rz(2.6685817) q[2];
sx q[2];
rz(-1.8188698) q[2];
sx q[2];
rz(-0.45516994) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1787923) q[1];
sx q[1];
rz(-1.1988678) q[1];
sx q[1];
rz(0.40403251) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86706676) q[3];
sx q[3];
rz(-0.19788225) q[3];
sx q[3];
rz(1.3615511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2565101) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(0.93290848) q[2];
rz(-1.7289303) q[3];
sx q[3];
rz(-1.7191073) q[3];
sx q[3];
rz(-2.9322114) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6949718) q[0];
sx q[0];
rz(-2.7734723) q[0];
sx q[0];
rz(-0.80286655) q[0];
rz(-0.74202263) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-0.41025695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78722505) q[0];
sx q[0];
rz(-2.4929308) q[0];
sx q[0];
rz(0.58726585) q[0];
rz(-1.4531936) q[2];
sx q[2];
rz(-0.80999331) q[2];
sx q[2];
rz(1.2076898) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32653836) q[1];
sx q[1];
rz(-1.4966623) q[1];
sx q[1];
rz(-1.733041) q[1];
x q[2];
rz(2.5531261) q[3];
sx q[3];
rz(-0.7477254) q[3];
sx q[3];
rz(-0.51008734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2731169) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(-0.22129076) q[2];
rz(2.7224329) q[3];
sx q[3];
rz(-2.6594888) q[3];
sx q[3];
rz(-0.47590772) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047423) q[0];
sx q[0];
rz(-0.99822799) q[0];
sx q[0];
rz(-1.5297484) q[0];
rz(1.8086241) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(1.6882247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3842363) q[0];
sx q[0];
rz(-1.0125375) q[0];
sx q[0];
rz(0.7406969) q[0];
rz(-pi) q[1];
rz(-1.4508993) q[2];
sx q[2];
rz(-1.7603619) q[2];
sx q[2];
rz(1.8360863) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.087880922) q[1];
sx q[1];
rz(-0.86999991) q[1];
sx q[1];
rz(1.7675964) q[1];
rz(-pi) q[2];
rz(2.9856921) q[3];
sx q[3];
rz(-1.8460801) q[3];
sx q[3];
rz(0.017214765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86614418) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(0.010738372) q[2];
rz(0.3206611) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(2.5519154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084932) q[0];
sx q[0];
rz(-2.1864317) q[0];
sx q[0];
rz(2.6668715) q[0];
rz(-2.8482598) q[1];
sx q[1];
rz(-2.0359998) q[1];
sx q[1];
rz(-1.6620103) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5005762) q[0];
sx q[0];
rz(-1.5307284) q[0];
sx q[0];
rz(3.1047719) q[0];
rz(-pi) q[1];
rz(-2.4655645) q[2];
sx q[2];
rz(-1.7601628) q[2];
sx q[2];
rz(2.2631136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.868728) q[1];
sx q[1];
rz(-1.1050779) q[1];
sx q[1];
rz(-1.8293423) q[1];
x q[2];
rz(0.18085705) q[3];
sx q[3];
rz(-1.3279727) q[3];
sx q[3];
rz(-2.9202785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.097229615) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(-0.43107671) q[2];
rz(-0.19032446) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(-2.2569807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97656074) q[0];
sx q[0];
rz(-0.29104069) q[0];
sx q[0];
rz(0.2188368) q[0];
rz(-2.1826375) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(-1.3921907) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16698027) q[0];
sx q[0];
rz(-1.7638399) q[0];
sx q[0];
rz(1.5184605) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8270589) q[2];
sx q[2];
rz(-2.7022903) q[2];
sx q[2];
rz(0.9010992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6619685) q[1];
sx q[1];
rz(-1.2733409) q[1];
sx q[1];
rz(0.65315078) q[1];
rz(-2.1851319) q[3];
sx q[3];
rz(-1.9298565) q[3];
sx q[3];
rz(0.39310716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2739111) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(-1.8217746) q[2];
rz(-0.53226081) q[3];
sx q[3];
rz(-1.3479193) q[3];
sx q[3];
rz(-1.5925315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.559448) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(-2.9136912) q[1];
sx q[1];
rz(-2.0961998) q[1];
sx q[1];
rz(-1.6963522) q[1];
rz(-2.7264222) q[2];
sx q[2];
rz(-1.3317458) q[2];
sx q[2];
rz(1.7785704) q[2];
rz(-1.3263973) q[3];
sx q[3];
rz(-2.5483589) q[3];
sx q[3];
rz(0.15729558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

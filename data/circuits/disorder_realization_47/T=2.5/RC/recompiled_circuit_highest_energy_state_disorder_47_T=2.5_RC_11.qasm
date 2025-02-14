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
rz(-2.6055251) q[0];
sx q[0];
rz(-1.1626838) q[0];
sx q[0];
rz(-2.7106078) q[0];
rz(0.35297901) q[1];
sx q[1];
rz(-1.9184435) q[1];
sx q[1];
rz(0.42981848) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777464) q[0];
sx q[0];
rz(-1.0809076) q[0];
sx q[0];
rz(1.6728982) q[0];
x q[1];
rz(-1.138428) q[2];
sx q[2];
rz(-1.5562487) q[2];
sx q[2];
rz(0.61847875) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42001549) q[1];
sx q[1];
rz(-0.41999451) q[1];
sx q[1];
rz(0.53424044) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6453703) q[3];
sx q[3];
rz(-1.2321951) q[3];
sx q[3];
rz(2.7994775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.063804403) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(-2.1303614) q[2];
rz(1.0412591) q[3];
sx q[3];
rz(-2.6285089) q[3];
sx q[3];
rz(0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-1.1312609) q[1];
sx q[1];
rz(0.25516137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19758148) q[0];
sx q[0];
rz(-1.7367878) q[0];
sx q[0];
rz(0.98910467) q[0];
rz(-pi) q[1];
rz(2.8361144) q[2];
sx q[2];
rz(-1.7490088) q[2];
sx q[2];
rz(0.62159789) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1000864) q[1];
sx q[1];
rz(-2.0165245) q[1];
sx q[1];
rz(-2.2959832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8273081) q[3];
sx q[3];
rz(-0.64033901) q[3];
sx q[3];
rz(2.4119055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.030563844) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(0.89152208) q[2];
rz(2.8920178) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(-0.3961302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4075277) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(0.044128142) q[0];
rz(-1.8924425) q[1];
sx q[1];
rz(-2.7154778) q[1];
sx q[1];
rz(2.6557907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51903782) q[0];
sx q[0];
rz(-1.3455032) q[0];
sx q[0];
rz(1.5507429) q[0];
rz(-pi) q[1];
rz(3.0032512) q[2];
sx q[2];
rz(-1.3855583) q[2];
sx q[2];
rz(1.1323765) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9456182) q[1];
sx q[1];
rz(-0.83688078) q[1];
sx q[1];
rz(0.7463781) q[1];
rz(-pi) q[2];
rz(1.1385192) q[3];
sx q[3];
rz(-1.2069993) q[3];
sx q[3];
rz(1.0897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32758632) q[2];
sx q[2];
rz(-2.0531211) q[2];
sx q[2];
rz(0.69474727) q[2];
rz(0.57573777) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(-1.7339138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1908252) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(-2.8488391) q[0];
rz(0.5689019) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(-0.84691602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0305188) q[0];
sx q[0];
rz(-1.9785709) q[0];
sx q[0];
rz(3.0570583) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90130536) q[2];
sx q[2];
rz(-1.6274823) q[2];
sx q[2];
rz(-1.1336807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9255815) q[1];
sx q[1];
rz(-2.8329509) q[1];
sx q[1];
rz(2.4587548) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0740739) q[3];
sx q[3];
rz(-1.3775) q[3];
sx q[3];
rz(-0.75261469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34951052) q[2];
sx q[2];
rz(-1.258472) q[2];
sx q[2];
rz(-2.9328031) q[2];
rz(2.4236692) q[3];
sx q[3];
rz(-2.8585298) q[3];
sx q[3];
rz(0.94055241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(2.4385117) q[1];
sx q[1];
rz(-0.73892006) q[1];
sx q[1];
rz(0.83831659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.188952) q[0];
sx q[0];
rz(-1.3172272) q[0];
sx q[0];
rz(2.6978561) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82413582) q[2];
sx q[2];
rz(-1.4243243) q[2];
sx q[2];
rz(2.9889088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.014480249) q[1];
sx q[1];
rz(-2.4408484) q[1];
sx q[1];
rz(1.6225918) q[1];
rz(-2.5637066) q[3];
sx q[3];
rz(-0.69212038) q[3];
sx q[3];
rz(0.96847615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68359739) q[2];
sx q[2];
rz(-1.0883051) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.018464) q[0];
sx q[0];
rz(-2.5557684) q[0];
sx q[0];
rz(-2.1783094) q[0];
rz(0.18114289) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(-1.2394261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2148149) q[0];
sx q[0];
rz(-0.96450761) q[0];
sx q[0];
rz(-2.0303594) q[0];
rz(-2.6341637) q[2];
sx q[2];
rz(-2.6118738) q[2];
sx q[2];
rz(-1.5629753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1787923) q[1];
sx q[1];
rz(-1.9427248) q[1];
sx q[1];
rz(-0.40403251) q[1];
x q[2];
rz(1.7224947) q[3];
sx q[3];
rz(-1.4432419) q[3];
sx q[3];
rz(2.2383245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2565101) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(-0.93290848) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949718) q[0];
sx q[0];
rz(-2.7734723) q[0];
sx q[0];
rz(-2.3387261) q[0];
rz(-0.74202263) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-0.41025695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091358034) q[0];
sx q[0];
rz(-2.0977533) q[0];
sx q[0];
rz(1.9684674) q[0];
x q[1];
rz(-0.76426498) q[2];
sx q[2];
rz(-1.4857123) q[2];
sx q[2];
rz(2.6972023) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8852119) q[1];
sx q[1];
rz(-1.4090013) q[1];
sx q[1];
rz(3.0664758) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4845552) q[3];
sx q[3];
rz(-1.1837623) q[3];
sx q[3];
rz(-1.6258193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2731169) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(-0.22129076) q[2];
rz(2.7224329) q[3];
sx q[3];
rz(-0.48210382) q[3];
sx q[3];
rz(0.47590772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047423) q[0];
sx q[0];
rz(-0.99822799) q[0];
sx q[0];
rz(-1.6118443) q[0];
rz(-1.3329685) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(-1.453368) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4299592) q[0];
sx q[0];
rz(-0.89444133) q[0];
sx q[0];
rz(2.3948689) q[0];
rz(-1.4508993) q[2];
sx q[2];
rz(-1.7603619) q[2];
sx q[2];
rz(-1.3055064) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9295974) q[1];
sx q[1];
rz(-0.72337615) q[1];
sx q[1];
rz(2.9138447) q[1];
rz(-pi) q[2];
rz(-2.0734208) q[3];
sx q[3];
rz(-0.31538559) q[3];
sx q[3];
rz(-0.50704765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2754485) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(-3.1308543) q[2];
rz(0.3206611) q[3];
sx q[3];
rz(-1.843822) q[3];
sx q[3];
rz(0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084932) q[0];
sx q[0];
rz(-2.1864317) q[0];
sx q[0];
rz(-0.47472111) q[0];
rz(0.2933329) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(-1.4795823) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06874456) q[0];
sx q[0];
rz(-1.6075875) q[0];
sx q[0];
rz(1.6108914) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4655645) q[2];
sx q[2];
rz(-1.7601628) q[2];
sx q[2];
rz(-2.2631136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.868728) q[1];
sx q[1];
rz(-2.0365148) q[1];
sx q[1];
rz(-1.8293423) q[1];
x q[2];
rz(0.18085705) q[3];
sx q[3];
rz(-1.81362) q[3];
sx q[3];
rz(-0.22131418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.097229615) q[2];
sx q[2];
rz(-1.2791415) q[2];
sx q[2];
rz(-2.7105159) q[2];
rz(-2.9512682) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(2.2569807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1650319) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(2.9227559) q[0];
rz(-0.95895514) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(1.3921907) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099571596) q[0];
sx q[0];
rz(-2.9416658) q[0];
sx q[0];
rz(-2.8801092) q[0];
x q[1];
rz(0.31453374) q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2710378) q[1];
sx q[1];
rz(-2.1907594) q[1];
sx q[1];
rz(-1.9391754) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1851319) q[3];
sx q[3];
rz(-1.2117361) q[3];
sx q[3];
rz(0.39310716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.86768156) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(1.8217746) q[2];
rz(2.6093318) q[3];
sx q[3];
rz(-1.3479193) q[3];
sx q[3];
rz(-1.5925315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821447) q[0];
sx q[0];
rz(-1.5503217) q[0];
sx q[0];
rz(0.4009608) q[0];
rz(2.9136912) q[1];
sx q[1];
rz(-1.0453929) q[1];
sx q[1];
rz(1.4452404) q[1];
rz(1.3105021) q[2];
sx q[2];
rz(-1.1681265) q[2];
sx q[2];
rz(0.10377965) q[2];
rz(-0.16172541) q[3];
sx q[3];
rz(-2.1441256) q[3];
sx q[3];
rz(-2.6921289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

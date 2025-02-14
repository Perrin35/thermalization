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
rz(7.4458691) q[0];
sx q[0];
rz(9.8557628) q[0];
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
rz(2.0864705) q[0];
sx q[0];
rz(-1.6608547) q[0];
sx q[0];
rz(-2.6495329) q[0];
rz(-pi) q[1];
rz(-1.5360898) q[2];
sx q[2];
rz(-2.7089951) q[2];
sx q[2];
rz(0.92080599) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1466521) q[1];
sx q[1];
rz(-1.9293679) q[1];
sx q[1];
rz(1.7943804) q[1];
x q[2];
rz(2.6453703) q[3];
sx q[3];
rz(-1.2321951) q[3];
sx q[3];
rz(2.7994775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.063804403) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(-1.0112313) q[2];
rz(1.0412591) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(2.7019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1083199) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(-2.394115) q[0];
rz(-2.8476818) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(-0.25516137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19758148) q[0];
sx q[0];
rz(-1.4048049) q[0];
sx q[0];
rz(2.152488) q[0];
rz(-1.7574649) q[2];
sx q[2];
rz(-1.2703086) q[2];
sx q[2];
rz(2.1365502) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9772759) q[1];
sx q[1];
rz(-0.9292047) q[1];
sx q[1];
rz(-2.5733828) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1952755) q[3];
sx q[3];
rz(-1.7229652) q[3];
sx q[3];
rz(-2.5077903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.030563844) q[2];
sx q[2];
rz(-0.054628987) q[2];
sx q[2];
rz(-0.89152208) q[2];
rz(-0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(-0.3961302) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7340649) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(0.044128142) q[0];
rz(1.8924425) q[1];
sx q[1];
rz(-2.7154778) q[1];
sx q[1];
rz(-2.6557907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7120948) q[0];
sx q[0];
rz(-2.9154239) q[0];
sx q[0];
rz(-3.0543213) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.75778) q[2];
sx q[2];
rz(-1.4348363) q[2];
sx q[2];
rz(-2.7288108) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9618157) q[1];
sx q[1];
rz(-2.0992341) q[1];
sx q[1];
rz(-2.4584199) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3084675) q[3];
sx q[3];
rz(-2.5841048) q[3];
sx q[3];
rz(-2.0035898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8140063) q[2];
sx q[2];
rz(-2.0531211) q[2];
sx q[2];
rz(-0.69474727) q[2];
rz(-0.57573777) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(-1.4076788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.1908252) q[0];
sx q[0];
rz(-2.8468843) q[0];
sx q[0];
rz(-2.8488391) q[0];
rz(0.5689019) q[1];
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
rz(-1.0305188) q[0];
sx q[0];
rz(-1.1630217) q[0];
sx q[0];
rz(-3.0570583) q[0];
x q[1];
rz(-2.2402873) q[2];
sx q[2];
rz(-1.5141103) q[2];
sx q[2];
rz(-2.0079119) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0140527) q[1];
sx q[1];
rz(-1.7636646) q[1];
sx q[1];
rz(-2.8991153) q[1];
rz(1.9027223) q[3];
sx q[3];
rz(-2.9369825) q[3];
sx q[3];
rz(-1.0910891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7920821) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(-0.2087896) q[2];
rz(2.4236692) q[3];
sx q[3];
rz(-2.8585298) q[3];
sx q[3];
rz(0.94055241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402886) q[0];
sx q[0];
rz(-2.6191819) q[0];
sx q[0];
rz(-2.8029602) q[0];
rz(2.4385117) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(-0.83831659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037443) q[0];
sx q[0];
rz(-2.6347343) q[0];
sx q[0];
rz(2.5985107) q[0];
x q[1];
rz(1.3569068) q[2];
sx q[2];
rz(-2.383432) q[2];
sx q[2];
rz(1.8799097) q[2];
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
rz(-pi) q[2];
rz(0.57788604) q[3];
sx q[3];
rz(-2.4494723) q[3];
sx q[3];
rz(-0.96847615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4579953) q[2];
sx q[2];
rz(-2.0532875) q[2];
sx q[2];
rz(-1.0028769) q[2];
rz(0.14497997) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(-0.063044757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.018464) q[0];
sx q[0];
rz(-2.5557684) q[0];
sx q[0];
rz(-0.96328324) q[0];
rz(-0.18114289) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(1.2394261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604804) q[0];
sx q[0];
rz(-1.1976997) q[0];
sx q[0];
rz(-2.4831071) q[0];
rz(1.2935898) q[2];
sx q[2];
rz(-1.1133901) q[2];
sx q[2];
rz(0.99062571) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0964862) q[1];
sx q[1];
rz(-0.54212057) q[1];
sx q[1];
rz(0.78150918) q[1];
rz(-pi) q[2];
rz(-0.86706676) q[3];
sx q[3];
rz(-0.19788225) q[3];
sx q[3];
rz(-1.7800415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88508254) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(2.2086842) q[2];
rz(-1.7289303) q[3];
sx q[3];
rz(-1.4224854) q[3];
sx q[3];
rz(-0.2093813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4466208) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(0.80286655) q[0];
rz(-0.74202263) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-0.41025695) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8703394) q[0];
sx q[0];
rz(-1.9121207) q[0];
sx q[0];
rz(-0.56296157) q[0];
x q[1];
rz(-2.3773277) q[2];
sx q[2];
rz(-1.6558803) q[2];
sx q[2];
rz(-0.44439038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2563808) q[1];
sx q[1];
rz(-1.7325914) q[1];
sx q[1];
rz(0.075116861) q[1];
rz(-pi) q[2];
rz(2.5531261) q[3];
sx q[3];
rz(-0.7477254) q[3];
sx q[3];
rz(-0.51008734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86847574) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(0.22129076) q[2];
rz(-0.41915974) q[3];
sx q[3];
rz(-0.48210382) q[3];
sx q[3];
rz(-2.6656849) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047423) q[0];
sx q[0];
rz(-0.99822799) q[0];
sx q[0];
rz(1.5297484) q[0];
rz(1.8086241) q[1];
sx q[1];
rz(-2.1905441) q[1];
sx q[1];
rz(1.453368) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4299592) q[0];
sx q[0];
rz(-0.89444133) q[0];
sx q[0];
rz(-2.3948689) q[0];
rz(1.6906934) q[2];
sx q[2];
rz(-1.3812307) q[2];
sx q[2];
rz(1.3055064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5308127) q[1];
sx q[1];
rz(-1.4207834) q[1];
sx q[1];
rz(2.4311745) q[1];
rz(-pi) q[2];
rz(2.0734208) q[3];
sx q[3];
rz(-2.8262071) q[3];
sx q[3];
rz(2.634545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2754485) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(3.1308543) q[2];
rz(-2.8209316) q[3];
sx q[3];
rz(-1.843822) q[3];
sx q[3];
rz(-2.5519154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084932) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(-2.6668715) q[0];
rz(2.8482598) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(-1.6620103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89746633) q[0];
sx q[0];
rz(-3.0871824) q[0];
sx q[0];
rz(-0.82798402) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67602812) q[2];
sx q[2];
rz(-1.7601628) q[2];
sx q[2];
rz(0.87847906) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.41614) q[1];
sx q[1];
rz(-1.3403157) q[1];
sx q[1];
rz(-2.6621755) q[1];
x q[2];
rz(-0.94274272) q[3];
sx q[3];
rz(-0.30170479) q[3];
sx q[3];
rz(-2.7127271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.044363) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(2.7105159) q[2];
rz(-2.9512682) q[3];
sx q[3];
rz(-0.35076916) q[3];
sx q[3];
rz(-2.2569807) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1650319) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(2.9227559) q[0];
rz(0.95895514) q[1];
sx q[1];
rz(-0.7364277) q[1];
sx q[1];
rz(1.3921907) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0420211) q[0];
sx q[0];
rz(-0.19992689) q[0];
sx q[0];
rz(-0.26148343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31453374) q[2];
sx q[2];
rz(-2.7022903) q[2];
sx q[2];
rz(-2.2404935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.47962418) q[1];
sx q[1];
rz(-1.8682518) q[1];
sx q[1];
rz(0.65315078) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(2.2739111) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(1.8217746) q[2];
rz(2.6093318) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(-1.5490612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.22790146) q[1];
sx q[1];
rz(-2.0961998) q[1];
sx q[1];
rz(-1.6963522) q[1];
rz(1.3105021) q[2];
sx q[2];
rz(-1.1681265) q[2];
sx q[2];
rz(0.10377965) q[2];
rz(0.16172541) q[3];
sx q[3];
rz(-0.99746708) q[3];
sx q[3];
rz(0.4494638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

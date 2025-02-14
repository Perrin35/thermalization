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
rz(-2.7886136) q[1];
sx q[1];
rz(-1.2231491) q[1];
sx q[1];
rz(-0.42981848) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0864705) q[0];
sx q[0];
rz(-1.6608547) q[0];
sx q[0];
rz(-0.49205972) q[0];
rz(-pi) q[1];
rz(-1.6055029) q[2];
sx q[2];
rz(-0.43259753) q[2];
sx q[2];
rz(0.92080599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65548518) q[1];
sx q[1];
rz(-1.7799417) q[1];
sx q[1];
rz(-2.7746592) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63685959) q[3];
sx q[3];
rz(-0.59266337) q[3];
sx q[3];
rz(-0.67837472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0777883) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(-2.1303614) q[2];
rz(-1.0412591) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(-2.7019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.1083199) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(2.394115) q[0];
rz(-0.29391089) q[1];
sx q[1];
rz(-2.0103318) q[1];
sx q[1];
rz(2.8864313) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19758148) q[0];
sx q[0];
rz(-1.4048049) q[0];
sx q[0];
rz(-2.152488) q[0];
rz(-pi) q[1];
rz(0.30547826) q[2];
sx q[2];
rz(-1.3925838) q[2];
sx q[2];
rz(-2.5199948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0415062) q[1];
sx q[1];
rz(-1.1250682) q[1];
sx q[1];
rz(0.84560945) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18682602) q[3];
sx q[3];
rz(-0.95462026) q[3];
sx q[3];
rz(1.0458067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.030563844) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(-0.89152208) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(0.3961302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075277) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(-0.044128142) q[0];
rz(-1.2491501) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(2.6557907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0853538) q[0];
sx q[0];
rz(-1.590343) q[0];
sx q[0];
rz(-2.9162558) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2052231) q[2];
sx q[2];
rz(-0.23072019) q[2];
sx q[2];
rz(-1.7796734) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9618157) q[1];
sx q[1];
rz(-2.0992341) q[1];
sx q[1];
rz(-0.68317271) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83312513) q[3];
sx q[3];
rz(-0.55748788) q[3];
sx q[3];
rz(-1.1380029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8140063) q[2];
sx q[2];
rz(-1.0884716) q[2];
sx q[2];
rz(-0.69474727) q[2];
rz(0.57573777) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(1.4076788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95076743) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(0.29275352) q[0];
rz(-2.5726908) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(2.2946766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0305188) q[0];
sx q[0];
rz(-1.1630217) q[0];
sx q[0];
rz(0.084534377) q[0];
x q[1];
rz(-1.4796094) q[2];
sx q[2];
rz(-2.4700748) q[2];
sx q[2];
rz(2.6330122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.12754) q[1];
sx q[1];
rz(-1.7636646) q[1];
sx q[1];
rz(-2.8991153) q[1];
rz(1.2388703) q[3];
sx q[3];
rz(-2.9369825) q[3];
sx q[3];
rz(-2.0505035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34951052) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(-2.9328031) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402886) q[0];
sx q[0];
rz(-2.6191819) q[0];
sx q[0];
rz(-2.8029602) q[0];
rz(-2.4385117) q[1];
sx q[1];
rz(-0.73892006) q[1];
sx q[1];
rz(2.3032761) q[1];
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
rz(-0.75816064) q[2];
sx q[2];
rz(1.2616829) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.014480249) q[1];
sx q[1];
rz(-0.70074425) q[1];
sx q[1];
rz(1.6225918) q[1];
x q[2];
rz(1.1456212) q[3];
sx q[3];
rz(-1.0068277) q[3];
sx q[3];
rz(-2.8759046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4579953) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(2.1387157) q[2];
rz(-0.14497997) q[3];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1231287) q[0];
sx q[0];
rz(-2.5557684) q[0];
sx q[0];
rz(0.96328324) q[0];
rz(2.9604498) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(1.2394261) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9267777) q[0];
sx q[0];
rz(-2.177085) q[0];
sx q[0];
rz(-2.0303594) q[0];
x q[1];
rz(-0.507429) q[2];
sx q[2];
rz(-0.52971887) q[2];
sx q[2];
rz(-1.5629753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2378578) q[1];
sx q[1];
rz(-1.1958201) q[1];
sx q[1];
rz(1.1695678) q[1];
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
rz(-0.88508254) q[2];
sx q[2];
rz(-0.91826597) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6949718) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(2.3387261) q[0];
rz(-0.74202263) q[1];
sx q[1];
rz(-1.6525533) q[1];
sx q[1];
rz(-2.7313357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8703394) q[0];
sx q[0];
rz(-1.9121207) q[0];
sx q[0];
rz(-0.56296157) q[0];
x q[1];
rz(-0.76426498) q[2];
sx q[2];
rz(-1.6558803) q[2];
sx q[2];
rz(0.44439038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32653836) q[1];
sx q[1];
rz(-1.4966623) q[1];
sx q[1];
rz(-1.4085517) q[1];
rz(1.0954082) q[3];
sx q[3];
rz(-0.96964004) q[3];
sx q[3];
rz(2.9132995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2731169) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(0.22129076) q[2];
rz(0.41915974) q[3];
sx q[3];
rz(-2.6594888) q[3];
sx q[3];
rz(-2.6656849) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047423) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(1.5297484) q[0];
rz(-1.3329685) q[1];
sx q[1];
rz(-2.1905441) q[1];
sx q[1];
rz(-1.6882247) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3842363) q[0];
sx q[0];
rz(-1.0125375) q[0];
sx q[0];
rz(-2.4008958) q[0];
rz(-pi) q[1];
rz(-2.58415) q[2];
sx q[2];
rz(-2.9176783) q[2];
sx q[2];
rz(-2.4049408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9295974) q[1];
sx q[1];
rz(-0.72337615) q[1];
sx q[1];
rz(-0.22774793) q[1];
x q[2];
rz(-0.15590053) q[3];
sx q[3];
rz(-1.8460801) q[3];
sx q[3];
rz(-3.1243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2754485) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(0.010738372) q[2];
rz(0.3206611) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(-0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43309942) q[0];
sx q[0];
rz(-2.1864317) q[0];
sx q[0];
rz(-2.6668715) q[0];
rz(-2.8482598) q[1];
sx q[1];
rz(-2.0359998) q[1];
sx q[1];
rz(-1.6620103) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89746633) q[0];
sx q[0];
rz(-0.054410283) q[0];
sx q[0];
rz(-2.3136086) q[0];
x q[1];
rz(-1.8117254) q[2];
sx q[2];
rz(-0.90902599) q[2];
sx q[2];
rz(-0.84217254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8822215) q[1];
sx q[1];
rz(-2.6135635) q[1];
sx q[1];
rz(2.6709983) q[1];
rz(-pi) q[2];
rz(0.94274272) q[3];
sx q[3];
rz(-2.8398879) q[3];
sx q[3];
rz(-2.7127271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.044363) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(-2.7105159) q[2];
rz(0.19032446) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(-0.88461191) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97656074) q[0];
sx q[0];
rz(-0.29104069) q[0];
sx q[0];
rz(0.2188368) q[0];
rz(-0.95895514) q[1];
sx q[1];
rz(-0.7364277) q[1];
sx q[1];
rz(-1.3921907) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7277273) q[0];
sx q[0];
rz(-1.5194335) q[0];
sx q[0];
rz(2.9482909) q[0];
rz(-2.8270589) q[2];
sx q[2];
rz(-2.7022903) q[2];
sx q[2];
rz(0.9010992) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6619685) q[1];
sx q[1];
rz(-1.8682518) q[1];
sx q[1];
rz(0.65315078) q[1];
rz(-pi) q[2];
rz(-2.1851319) q[3];
sx q[3];
rz(-1.9298565) q[3];
sx q[3];
rz(0.39310716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2739111) q[2];
sx q[2];
rz(-0.76649222) q[2];
sx q[2];
rz(-1.319818) q[2];
rz(2.6093318) q[3];
sx q[3];
rz(-1.3479193) q[3];
sx q[3];
rz(1.5490612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5821447) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(2.9136912) q[1];
sx q[1];
rz(-1.0453929) q[1];
sx q[1];
rz(1.4452404) q[1];
rz(2.7264222) q[2];
sx q[2];
rz(-1.8098469) q[2];
sx q[2];
rz(-1.3630223) q[2];
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

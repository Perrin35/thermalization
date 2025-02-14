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
rz(-2.7886136) q[1];
sx q[1];
rz(-1.2231491) q[1];
sx q[1];
rz(-0.42981848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777464) q[0];
sx q[0];
rz(-2.060685) q[0];
sx q[0];
rz(-1.6728982) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6055029) q[2];
sx q[2];
rz(-0.43259753) q[2];
sx q[2];
rz(2.2207867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65548518) q[1];
sx q[1];
rz(-1.361651) q[1];
sx q[1];
rz(-0.36693348) q[1];
rz(-pi) q[2];
rz(0.63685959) q[3];
sx q[3];
rz(-2.5489293) q[3];
sx q[3];
rz(0.67837472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0777883) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(2.1303614) q[2];
rz(1.0412591) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(-0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0332727) q[0];
sx q[0];
rz(-2.1450873) q[0];
sx q[0];
rz(2.394115) q[0];
rz(2.8476818) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(-2.8864313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19758148) q[0];
sx q[0];
rz(-1.4048049) q[0];
sx q[0];
rz(2.152488) q[0];
x q[1];
rz(0.5396217) q[2];
sx q[2];
rz(-2.789342) q[2];
sx q[2];
rz(0.43708153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1595711) q[1];
sx q[1];
rz(-0.82948331) q[1];
sx q[1];
rz(2.1950566) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9547666) q[3];
sx q[3];
rz(-2.1869724) q[3];
sx q[3];
rz(-2.0957859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.030563844) q[2];
sx q[2];
rz(-0.054628987) q[2];
sx q[2];
rz(2.2500706) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(0.3961302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075277) q[0];
sx q[0];
rz(-1.1363109) q[0];
sx q[0];
rz(-3.0974645) q[0];
rz(-1.8924425) q[1];
sx q[1];
rz(-2.7154778) q[1];
sx q[1];
rz(-0.485802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42949781) q[0];
sx q[0];
rz(-2.9154239) q[0];
sx q[0];
rz(3.0543213) q[0];
rz(-pi) q[1];
rz(1.75778) q[2];
sx q[2];
rz(-1.4348363) q[2];
sx q[2];
rz(-0.41278186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9618157) q[1];
sx q[1];
rz(-1.0423585) q[1];
sx q[1];
rz(0.68317271) q[1];
x q[2];
rz(-2.0030735) q[3];
sx q[3];
rz(-1.2069993) q[3];
sx q[3];
rz(1.0897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32758632) q[2];
sx q[2];
rz(-2.0531211) q[2];
sx q[2];
rz(2.4468454) q[2];
rz(0.57573777) q[3];
sx q[3];
rz(-1.271194) q[3];
sx q[3];
rz(-1.4076788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1908252) q[0];
sx q[0];
rz(-2.8468843) q[0];
sx q[0];
rz(-0.29275352) q[0];
rz(0.5689019) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(2.2946766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0305188) q[0];
sx q[0];
rz(-1.9785709) q[0];
sx q[0];
rz(0.084534377) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.072242468) q[2];
sx q[2];
rz(-2.2390167) q[2];
sx q[2];
rz(2.7492858) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49063045) q[1];
sx q[1];
rz(-1.3329026) q[1];
sx q[1];
rz(-1.3722653) q[1];
x q[2];
rz(-1.3770695) q[3];
sx q[3];
rz(-1.5045369) q[3];
sx q[3];
rz(2.3364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34951052) q[2];
sx q[2];
rz(-1.258472) q[2];
sx q[2];
rz(0.2087896) q[2];
rz(0.71792349) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(0.94055241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402886) q[0];
sx q[0];
rz(-2.6191819) q[0];
sx q[0];
rz(0.33863246) q[0];
rz(2.4385117) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(-0.83831659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95264065) q[0];
sx q[0];
rz(-1.3172272) q[0];
sx q[0];
rz(0.44373657) q[0];
rz(-pi) q[1];
rz(-1.3569068) q[2];
sx q[2];
rz(-2.383432) q[2];
sx q[2];
rz(1.2616829) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5456713) q[1];
sx q[1];
rz(-1.537408) q[1];
sx q[1];
rz(0.8707134) q[1];
rz(-1.9959715) q[3];
sx q[3];
rz(-2.1347649) q[3];
sx q[3];
rz(2.8759046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4579953) q[2];
sx q[2];
rz(-2.0532875) q[2];
sx q[2];
rz(1.0028769) q[2];
rz(-2.9966127) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(3.0785479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.018464) q[0];
sx q[0];
rz(-2.5557684) q[0];
sx q[0];
rz(0.96328324) q[0];
rz(0.18114289) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(1.9021665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081112208) q[0];
sx q[0];
rz(-1.1976997) q[0];
sx q[0];
rz(-2.4831071) q[0];
x q[1];
rz(2.6685817) q[2];
sx q[2];
rz(-1.8188698) q[2];
sx q[2];
rz(-0.45516994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9037348) q[1];
sx q[1];
rz(-1.9457726) q[1];
sx q[1];
rz(-1.1695678) q[1];
x q[2];
rz(0.86706676) q[3];
sx q[3];
rz(-0.19788225) q[3];
sx q[3];
rz(1.7800415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88508254) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(2.2086842) q[2];
rz(1.7289303) q[3];
sx q[3];
rz(-1.4224854) q[3];
sx q[3];
rz(-2.9322114) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4466208) q[0];
sx q[0];
rz(-2.7734723) q[0];
sx q[0];
rz(2.3387261) q[0];
rz(-2.39957) q[1];
sx q[1];
rz(-1.6525533) q[1];
sx q[1];
rz(2.7313357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3543676) q[0];
sx q[0];
rz(-2.4929308) q[0];
sx q[0];
rz(0.58726585) q[0];
rz(-pi) q[1];
rz(-1.4531936) q[2];
sx q[2];
rz(-0.80999331) q[2];
sx q[2];
rz(1.2076898) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8852119) q[1];
sx q[1];
rz(-1.7325914) q[1];
sx q[1];
rz(3.0664758) q[1];
x q[2];
rz(-1.0954082) q[3];
sx q[3];
rz(-0.96964004) q[3];
sx q[3];
rz(-2.9132995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86847574) q[2];
sx q[2];
rz(-1.5831524) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23685037) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(-1.6118443) q[0];
rz(-1.3329685) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(1.6882247) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8771031) q[0];
sx q[0];
rz(-0.96142381) q[0];
sx q[0];
rz(2.2731056) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55744268) q[2];
sx q[2];
rz(-2.9176783) q[2];
sx q[2];
rz(-2.4049408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.087880922) q[1];
sx q[1];
rz(-2.2715927) q[1];
sx q[1];
rz(1.3739963) q[1];
rz(-pi) q[2];
rz(-1.8492886) q[3];
sx q[3];
rz(-1.4208111) q[3];
sx q[3];
rz(-1.5962792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86614418) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7084932) q[0];
sx q[0];
rz(-2.1864317) q[0];
sx q[0];
rz(-2.6668715) q[0];
rz(-2.8482598) q[1];
sx q[1];
rz(-2.0359998) q[1];
sx q[1];
rz(1.4795823) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89746633) q[0];
sx q[0];
rz(-0.054410283) q[0];
sx q[0];
rz(0.82798402) q[0];
rz(2.4655645) q[2];
sx q[2];
rz(-1.7601628) q[2];
sx q[2];
rz(-2.2631136) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8822215) q[1];
sx q[1];
rz(-2.6135635) q[1];
sx q[1];
rz(-2.6709983) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1988499) q[3];
sx q[3];
rz(-0.30170479) q[3];
sx q[3];
rz(0.42886558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.044363) q[2];
sx q[2];
rz(-1.8624511) q[2];
sx q[2];
rz(2.7105159) q[2];
rz(-2.9512682) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(-0.88461191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.7494019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0420211) q[0];
sx q[0];
rz(-2.9416658) q[0];
sx q[0];
rz(-2.8801092) q[0];
rz(0.31453374) q[2];
sx q[2];
rz(-2.7022903) q[2];
sx q[2];
rz(-2.2404935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6619685) q[1];
sx q[1];
rz(-1.2733409) q[1];
sx q[1];
rz(-0.65315078) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99361657) q[3];
sx q[3];
rz(-0.69968596) q[3];
sx q[3];
rz(-2.4260739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86768156) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(1.319818) q[2];
rz(0.53226081) q[3];
sx q[3];
rz(-1.3479193) q[3];
sx q[3];
rz(-1.5490612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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

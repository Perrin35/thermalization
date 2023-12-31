OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(2.8136301) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56337315) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(-2.5114775) q[0];
rz(-pi) q[1];
x q[1];
rz(2.882471) q[2];
sx q[2];
rz(-1.7012351) q[2];
sx q[2];
rz(-2.5082617) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5603197) q[1];
sx q[1];
rz(-2.8898507) q[1];
sx q[1];
rz(-1.9074349) q[1];
rz(-pi) q[2];
rz(-2.1655733) q[3];
sx q[3];
rz(-1.752749) q[3];
sx q[3];
rz(2.6973157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(-2.8664355) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317516) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(-0.69357187) q[0];
rz(2.0454848) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(0.19031659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41102558) q[0];
sx q[0];
rz(-1.9733631) q[0];
sx q[0];
rz(2.9257724) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1778189) q[2];
sx q[2];
rz(-1.0699546) q[2];
sx q[2];
rz(2.9289392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7029876) q[1];
sx q[1];
rz(-2.4191796) q[1];
sx q[1];
rz(-1.5370675) q[1];
x q[2];
rz(0.72782794) q[3];
sx q[3];
rz(-1.2065294) q[3];
sx q[3];
rz(-1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.50743121) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(2.9988585) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(-0.8202585) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.8935727) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67278636) q[0];
sx q[0];
rz(-0.60657036) q[0];
sx q[0];
rz(0.16288217) q[0];
rz(-pi) q[1];
rz(0.2072316) q[2];
sx q[2];
rz(-1.6336294) q[2];
sx q[2];
rz(0.40916967) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.9995211) q[1];
sx q[1];
rz(-1.6377249) q[1];
sx q[1];
rz(0.66758572) q[1];
rz(2.2218024) q[3];
sx q[3];
rz(-1.3295104) q[3];
sx q[3];
rz(1.1834708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.9281663) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(-0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(1.8444555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38940934) q[0];
sx q[0];
rz(-1.7011189) q[0];
sx q[0];
rz(-0.043686314) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0879891) q[2];
sx q[2];
rz(-0.32978168) q[2];
sx q[2];
rz(1.5151378) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12249494) q[1];
sx q[1];
rz(-2.2593453) q[1];
sx q[1];
rz(2.3282939) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6468871) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(-1.0877346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-0.91919351) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(-1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(-0.87819535) q[0];
rz(-1.325266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.426288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.250524) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(3.0380681) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23767383) q[2];
sx q[2];
rz(-1.2665247) q[2];
sx q[2];
rz(-1.301847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1713472) q[1];
sx q[1];
rz(-2.4023224) q[1];
sx q[1];
rz(1.7319748) q[1];
rz(-pi) q[2];
rz(-1.8910847) q[3];
sx q[3];
rz(-0.8257782) q[3];
sx q[3];
rz(-0.68009963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(-0.041794725) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(-1.0304931) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-0.57156634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540928) q[0];
sx q[0];
rz(-2.0797074) q[0];
sx q[0];
rz(-0.66977588) q[0];
rz(0.12708725) q[2];
sx q[2];
rz(-2.2327144) q[2];
sx q[2];
rz(1.454929) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7476269) q[1];
sx q[1];
rz(-2.6193301) q[1];
sx q[1];
rz(1.0737435) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84047517) q[3];
sx q[3];
rz(-2.2395036) q[3];
sx q[3];
rz(0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(-3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(0.71227658) q[0];
rz(-2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-0.66551048) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3411361) q[0];
sx q[0];
rz(-1.8130214) q[0];
sx q[0];
rz(-0.0010629396) q[0];
rz(2.268157) q[2];
sx q[2];
rz(-1.6849815) q[2];
sx q[2];
rz(2.9261677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9741386) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(-2.8970701) q[1];
rz(-pi) q[2];
rz(-2.9969278) q[3];
sx q[3];
rz(-1.8105227) q[3];
sx q[3];
rz(-0.02558115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15726382) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(1.4228014) q[2];
rz(-3.026399) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049906235) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-0.33871067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539464) q[0];
sx q[0];
rz(-0.37746143) q[0];
sx q[0];
rz(-1.1520555) q[0];
rz(-pi) q[1];
rz(-2.6892002) q[2];
sx q[2];
rz(-2.6466742) q[2];
sx q[2];
rz(0.53168833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45704493) q[1];
sx q[1];
rz(-1.0984048) q[1];
sx q[1];
rz(-0.12462516) q[1];
rz(-pi) q[2];
rz(1.8975774) q[3];
sx q[3];
rz(-1.3076926) q[3];
sx q[3];
rz(-0.14740482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64615858) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.609628) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(-0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(3.0019965) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1698648) q[0];
sx q[0];
rz(-1.5507878) q[0];
sx q[0];
rz(-1.6318897) q[0];
x q[1];
rz(-0.058768674) q[2];
sx q[2];
rz(-0.81548703) q[2];
sx q[2];
rz(1.2440484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86233625) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(-0.52211296) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.052501909) q[3];
sx q[3];
rz(-0.92696654) q[3];
sx q[3];
rz(-2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6167986) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(-0.33111462) q[2];
rz(0.75774276) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(3.0537135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452633) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(-0.18558003) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508704) q[0];
sx q[0];
rz(-1.8403247) q[0];
sx q[0];
rz(2.0275293) q[0];
x q[1];
rz(0.22557232) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(-3.0849506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(1.9468716) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2757169) q[3];
sx q[3];
rz(-1.5745592) q[3];
sx q[3];
rz(2.3072484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(0.55220848) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-2.6681343) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(-2.3579303) q[3];
sx q[3];
rz(-2.2829934) q[3];
sx q[3];
rz(-0.34096277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

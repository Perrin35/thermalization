OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(-0.8015269) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(-0.094223082) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4163001) q[0];
sx q[0];
rz(-0.8612809) q[0];
sx q[0];
rz(-0.50523357) q[0];
x q[1];
rz(-1.6802189) q[2];
sx q[2];
rz(-0.27421194) q[2];
sx q[2];
rz(-2.4502129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6871145) q[1];
sx q[1];
rz(-1.2839437) q[1];
sx q[1];
rz(2.0583365) q[1];
rz(-2.9889936) q[3];
sx q[3];
rz(-1.5729701) q[3];
sx q[3];
rz(1.847037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5376544) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(1.630265) q[2];
rz(0.18621914) q[3];
sx q[3];
rz(-2.7112466) q[3];
sx q[3];
rz(0.65096861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179825) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(-3.004916) q[0];
rz(-0.094820529) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(-0.76900855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.626132) q[0];
sx q[0];
rz(-2.7410586) q[0];
sx q[0];
rz(-0.35450165) q[0];
rz(-pi) q[1];
rz(-0.42119512) q[2];
sx q[2];
rz(-0.73891089) q[2];
sx q[2];
rz(0.74904672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0817854) q[1];
sx q[1];
rz(-1.4084653) q[1];
sx q[1];
rz(-2.3947706) q[1];
rz(0.034763408) q[3];
sx q[3];
rz(-1.4357908) q[3];
sx q[3];
rz(-0.95187675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7424221) q[2];
sx q[2];
rz(-2.6593282) q[2];
sx q[2];
rz(1.4313618) q[2];
rz(0.4906022) q[3];
sx q[3];
rz(-2.6503745) q[3];
sx q[3];
rz(-2.365999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(1.8785777) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(2.0468792) q[0];
rz(-1.8393983) q[1];
sx q[1];
rz(-0.98919386) q[1];
sx q[1];
rz(-3.1375695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4152609) q[0];
sx q[0];
rz(-1.7021966) q[0];
sx q[0];
rz(2.2982928) q[0];
rz(-0.74538576) q[2];
sx q[2];
rz(-2.2507387) q[2];
sx q[2];
rz(2.5542575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5044671) q[1];
sx q[1];
rz(-1.5404623) q[1];
sx q[1];
rz(0.25719579) q[1];
x q[2];
rz(1.324292) q[3];
sx q[3];
rz(-2.26943) q[3];
sx q[3];
rz(-0.43432954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5502988) q[2];
sx q[2];
rz(-2.55547) q[2];
sx q[2];
rz(-2.6934521) q[2];
rz(-0.48629931) q[3];
sx q[3];
rz(-0.86217642) q[3];
sx q[3];
rz(1.1264616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124741) q[0];
sx q[0];
rz(-1.7121226) q[0];
sx q[0];
rz(-2.4718156) q[0];
rz(-1.9122596) q[1];
sx q[1];
rz(-2.6826617) q[1];
sx q[1];
rz(0.545048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9733423) q[0];
sx q[0];
rz(-2.0297883) q[0];
sx q[0];
rz(2.5516872) q[0];
x q[1];
rz(2.3178905) q[2];
sx q[2];
rz(-0.81497008) q[2];
sx q[2];
rz(-2.4832249) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5601756) q[1];
sx q[1];
rz(-0.22318527) q[1];
sx q[1];
rz(1.2944004) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3432076) q[3];
sx q[3];
rz(-1.9163864) q[3];
sx q[3];
rz(-1.2452728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5547319) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(-0.16793212) q[2];
rz(-1.933291) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(-1.4471853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55030441) q[0];
sx q[0];
rz(-1.2606324) q[0];
sx q[0];
rz(-3.1377129) q[0];
rz(-0.30715352) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(-0.53877962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8545495) q[0];
sx q[0];
rz(-1.3417305) q[0];
sx q[0];
rz(1.4255217) q[0];
x q[1];
rz(0.42947763) q[2];
sx q[2];
rz(-1.4444077) q[2];
sx q[2];
rz(-0.9867368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6314124) q[1];
sx q[1];
rz(-2.4459527) q[1];
sx q[1];
rz(2.5744777) q[1];
rz(-pi) q[2];
rz(2.0224376) q[3];
sx q[3];
rz(-2.2502101) q[3];
sx q[3];
rz(-1.7070626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6946081) q[2];
sx q[2];
rz(-1.0993404) q[2];
sx q[2];
rz(1.2896607) q[2];
rz(-1.6192216) q[3];
sx q[3];
rz(-2.3963942) q[3];
sx q[3];
rz(1.4815909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30037844) q[0];
sx q[0];
rz(-0.9391681) q[0];
sx q[0];
rz(3.0389431) q[0];
rz(-0.73219055) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(0.57714677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10023663) q[0];
sx q[0];
rz(-1.5994342) q[0];
sx q[0];
rz(3.1379382) q[0];
rz(1.4056751) q[2];
sx q[2];
rz(-0.99011865) q[2];
sx q[2];
rz(2.6613147) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53458339) q[1];
sx q[1];
rz(-1.680169) q[1];
sx q[1];
rz(1.6499728) q[1];
x q[2];
rz(1.243337) q[3];
sx q[3];
rz(-0.84210448) q[3];
sx q[3];
rz(-0.00086833894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.010043667) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(1.5232167) q[2];
rz(0.067954436) q[3];
sx q[3];
rz(-0.32502919) q[3];
sx q[3];
rz(2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286569) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(-0.77142429) q[0];
rz(-0.5386638) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(-2.0753986) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5911884) q[0];
sx q[0];
rz(-0.26445358) q[0];
sx q[0];
rz(-1.2729086) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9087139) q[2];
sx q[2];
rz(-1.0693197) q[2];
sx q[2];
rz(-1.5053144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2354338) q[1];
sx q[1];
rz(-1.3536436) q[1];
sx q[1];
rz(1.557334) q[1];
rz(-pi) q[2];
rz(-0.33447075) q[3];
sx q[3];
rz(-1.3190509) q[3];
sx q[3];
rz(1.2125963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6923339) q[2];
sx q[2];
rz(-0.80654311) q[2];
sx q[2];
rz(-2.6991357) q[2];
rz(2.9509406) q[3];
sx q[3];
rz(-0.72398829) q[3];
sx q[3];
rz(0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2451179) q[0];
sx q[0];
rz(-0.18558311) q[0];
sx q[0];
rz(1.8316487) q[0];
rz(2.7438121) q[1];
sx q[1];
rz(-0.82292992) q[1];
sx q[1];
rz(-1.3936874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0017477) q[0];
sx q[0];
rz(-2.4187517) q[0];
sx q[0];
rz(-0.87840898) q[0];
rz(-0.23598052) q[2];
sx q[2];
rz(-1.8644635) q[2];
sx q[2];
rz(0.27494173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2893148) q[1];
sx q[1];
rz(-0.56486579) q[1];
sx q[1];
rz(2.5603065) q[1];
x q[2];
rz(1.3721136) q[3];
sx q[3];
rz(-1.4666838) q[3];
sx q[3];
rz(0.35752359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1132249) q[2];
sx q[2];
rz(-2.4312879) q[2];
sx q[2];
rz(0.073471546) q[2];
rz(-2.4469589) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(-1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6756814) q[0];
sx q[0];
rz(-0.2008734) q[0];
sx q[0];
rz(2.1821816) q[0];
rz(-2.361946) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(-2.8693105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1355541) q[0];
sx q[0];
rz(-1.557543) q[0];
sx q[0];
rz(-0.0051965836) q[0];
x q[1];
rz(2.831976) q[2];
sx q[2];
rz(-0.88954347) q[2];
sx q[2];
rz(-0.75971425) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12676621) q[1];
sx q[1];
rz(-2.0527168) q[1];
sx q[1];
rz(-2.6792106) q[1];
rz(-2.5078154) q[3];
sx q[3];
rz(-0.85481516) q[3];
sx q[3];
rz(-0.26308003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4043364) q[2];
sx q[2];
rz(-1.6360838) q[2];
sx q[2];
rz(-0.20757248) q[2];
rz(0.10457822) q[3];
sx q[3];
rz(-2.6360377) q[3];
sx q[3];
rz(-1.526621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.680147) q[0];
sx q[0];
rz(-1.2274281) q[0];
sx q[0];
rz(1.0567868) q[0];
rz(2.5959004) q[1];
sx q[1];
rz(-2.1634407) q[1];
sx q[1];
rz(-2.812885) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0820887) q[0];
sx q[0];
rz(-0.57841136) q[0];
sx q[0];
rz(2.2057981) q[0];
rz(-0.07241197) q[2];
sx q[2];
rz(-1.6254404) q[2];
sx q[2];
rz(1.8397728) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67519611) q[1];
sx q[1];
rz(-1.5541142) q[1];
sx q[1];
rz(-1.9128602) q[1];
rz(-pi) q[2];
rz(1.2966818) q[3];
sx q[3];
rz(-2.3233827) q[3];
sx q[3];
rz(-2.3880098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1591961) q[2];
sx q[2];
rz(-3.0089162) q[2];
sx q[2];
rz(1.1245493) q[2];
rz(0.076796181) q[3];
sx q[3];
rz(-2.7087637) q[3];
sx q[3];
rz(-0.92397773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9818253) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(-2.3253597) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(-0.53173595) q[2];
sx q[2];
rz(-0.81960631) q[2];
sx q[2];
rz(-1.3255957) q[2];
rz(-2.7362551) q[3];
sx q[3];
rz(-1.5673076) q[3];
sx q[3];
rz(0.81893541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

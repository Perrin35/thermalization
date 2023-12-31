OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(3.8656524) q[0];
sx q[0];
rz(11.081628) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(-2.2533921) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0892544) q[0];
sx q[0];
rz(-1.5262145) q[0];
sx q[0];
rz(2.9705439) q[0];
x q[1];
rz(-1.2970096) q[2];
sx q[2];
rz(-2.5878776) q[2];
sx q[2];
rz(-2.6562429) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8623212) q[1];
sx q[1];
rz(-2.9746378) q[1];
sx q[1];
rz(-0.0032940666) q[1];
x q[2];
rz(2.4429951) q[3];
sx q[3];
rz(-1.5454096) q[3];
sx q[3];
rz(-0.7198402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28329864) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-2.0236012) q[2];
rz(-2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(-3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730826) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(-2.4867687) q[0];
rz(1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(0.12589802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038371041) q[0];
sx q[0];
rz(-1.201259) q[0];
sx q[0];
rz(2.080337) q[0];
rz(-0.89181487) q[2];
sx q[2];
rz(-2.0795155) q[2];
sx q[2];
rz(1.5430792) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0071348) q[1];
sx q[1];
rz(-1.4440698) q[1];
sx q[1];
rz(3.111372) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1308261) q[3];
sx q[3];
rz(-2.2942703) q[3];
sx q[3];
rz(0.41124287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.048432365) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(-1.7175425) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(0.079332381) q[0];
rz(0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(-1.1598587) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9281611) q[0];
sx q[0];
rz(-1.0614938) q[0];
sx q[0];
rz(3.0470949) q[0];
rz(-pi) q[1];
rz(-2.2247505) q[2];
sx q[2];
rz(-2.4008304) q[2];
sx q[2];
rz(0.11292085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0283454) q[1];
sx q[1];
rz(-2.5302437) q[1];
sx q[1];
rz(0.69797413) q[1];
rz(-pi) q[2];
rz(1.6634116) q[3];
sx q[3];
rz(-2.5941656) q[3];
sx q[3];
rz(-0.80143354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(0.1082871) q[2];
rz(2.4978499) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.9765587) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(-1.4820341) q[0];
rz(0.7011134) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(-0.28465095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4581504) q[0];
sx q[0];
rz(-2.0476641) q[0];
sx q[0];
rz(0.083806888) q[0];
rz(-0.20547159) q[2];
sx q[2];
rz(-1.2887508) q[2];
sx q[2];
rz(-0.13946433) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9173911) q[1];
sx q[1];
rz(-2.9926139) q[1];
sx q[1];
rz(1.0148744) q[1];
x q[2];
rz(-0.48216362) q[3];
sx q[3];
rz(-2.1820222) q[3];
sx q[3];
rz(1.5507444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.231679) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(-2.5740734) q[2];
rz(-2.7275758) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48859566) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(-1.1524221) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(2.2036536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65813488) q[0];
sx q[0];
rz(-1.5926653) q[0];
sx q[0];
rz(1.4786426) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3100501) q[2];
sx q[2];
rz(-2.4387896) q[2];
sx q[2];
rz(1.0520832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2069619) q[1];
sx q[1];
rz(-3.0983371) q[1];
sx q[1];
rz(-2.2224777) q[1];
rz(0.35477562) q[3];
sx q[3];
rz(-1.6924904) q[3];
sx q[3];
rz(-2.0718758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-2.1002634) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(-2.5842216) q[0];
rz(-0.56466651) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-0.55647892) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59069955) q[0];
sx q[0];
rz(-1.7771582) q[0];
sx q[0];
rz(2.7569689) q[0];
rz(-pi) q[1];
rz(-0.85530497) q[2];
sx q[2];
rz(-1.7024634) q[2];
sx q[2];
rz(1.5065187) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52292697) q[1];
sx q[1];
rz(-1.7523125) q[1];
sx q[1];
rz(1.8261441) q[1];
x q[2];
rz(-0.052080215) q[3];
sx q[3];
rz(-2.3651603) q[3];
sx q[3];
rz(-1.4610425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.53211987) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.2247941) q[2];
rz(1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0790134) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(-3.0986837) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(2.6409805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6734877) q[0];
sx q[0];
rz(-1.3930495) q[0];
sx q[0];
rz(0.034981473) q[0];
rz(2.826564) q[2];
sx q[2];
rz(-2.0909967) q[2];
sx q[2];
rz(2.5556285) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40655207) q[1];
sx q[1];
rz(-2.7222689) q[1];
sx q[1];
rz(-0.29962824) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75617366) q[3];
sx q[3];
rz(-1.2740967) q[3];
sx q[3];
rz(2.365436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(0.67561692) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(2.6202257) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(-1.0674397) q[0];
rz(2.7087129) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(-1.4656461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059793652) q[0];
sx q[0];
rz(-2.0311211) q[0];
sx q[0];
rz(2.9914809) q[0];
x q[1];
rz(-0.18657121) q[2];
sx q[2];
rz(-1.6779043) q[2];
sx q[2];
rz(-2.4466116) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17897478) q[1];
sx q[1];
rz(-1.0777506) q[1];
sx q[1];
rz(-0.38260539) q[1];
rz(0.74387868) q[3];
sx q[3];
rz(-0.79259593) q[3];
sx q[3];
rz(2.511123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33891588) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(1.6652997) q[2];
rz(-0.87351292) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840238) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(-2.4043758) q[0];
rz(-3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-0.92528701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4243471) q[0];
sx q[0];
rz(-2.2204917) q[0];
sx q[0];
rz(-2.2109277) q[0];
rz(-pi) q[1];
rz(1.4383573) q[2];
sx q[2];
rz(-2.2297511) q[2];
sx q[2];
rz(0.027564136) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8268938) q[1];
sx q[1];
rz(-1.7539382) q[1];
sx q[1];
rz(-1.4383184) q[1];
x q[2];
rz(0.85261811) q[3];
sx q[3];
rz(-1.2468306) q[3];
sx q[3];
rz(-0.94911239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3725738) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(2.1378689) q[2];
rz(0.090099661) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039624778) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(-1.2596624) q[0];
rz(-2.8758077) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(-0.75751799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82523358) q[0];
sx q[0];
rz(-2.9037583) q[0];
sx q[0];
rz(2.1912078) q[0];
rz(1.3390113) q[2];
sx q[2];
rz(-0.60283649) q[2];
sx q[2];
rz(1.9581118) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.65749189) q[1];
sx q[1];
rz(-1.3490281) q[1];
sx q[1];
rz(1.1360672) q[1];
rz(-pi) q[2];
rz(-2.431589) q[3];
sx q[3];
rz(-1.5922976) q[3];
sx q[3];
rz(1.8332675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99047986) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-0.63888597) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4929852) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(0.63207788) q[2];
sx q[2];
rz(-0.52543228) q[2];
sx q[2];
rz(0.57731522) q[2];
rz(-2.0414447) q[3];
sx q[3];
rz(-0.63334076) q[3];
sx q[3];
rz(0.28950194) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

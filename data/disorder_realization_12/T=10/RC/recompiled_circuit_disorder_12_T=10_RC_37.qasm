OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.84500399) q[0];
sx q[0];
rz(-0.72405976) q[0];
sx q[0];
rz(1.4847423) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(-2.2533921) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22904299) q[0];
sx q[0];
rz(-0.17670822) q[0];
sx q[0];
rz(0.25632174) q[0];
rz(1.2970096) q[2];
sx q[2];
rz(-0.55371504) q[2];
sx q[2];
rz(0.48534976) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29477316) q[1];
sx q[1];
rz(-1.5713437) q[1];
sx q[1];
rz(-2.9746387) q[1];
x q[2];
rz(-0.69859759) q[3];
sx q[3];
rz(-1.5454096) q[3];
sx q[3];
rz(-0.7198402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.858294) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(2.0236012) q[2];
rz(-2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685101) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(-0.65482393) q[0];
rz(-1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-0.12589802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038371041) q[0];
sx q[0];
rz(-1.201259) q[0];
sx q[0];
rz(-1.0612556) q[0];
rz(-pi) q[1];
rz(-2.5198031) q[2];
sx q[2];
rz(-0.99025531) q[2];
sx q[2];
rz(0.34678005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44015917) q[1];
sx q[1];
rz(-1.5408181) q[1];
sx q[1];
rz(1.4440126) q[1];
x q[2];
rz(0.84729362) q[3];
sx q[3];
rz(-1.5788659) q[3];
sx q[3];
rz(-1.1524259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0931603) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(0.38802567) q[2];
rz(1.4240501) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5858784) q[0];
sx q[0];
rz(-0.782574) q[0];
sx q[0];
rz(0.079332381) q[0];
rz(3.0575867) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(1.9817339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0214329) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(-1.73818) q[0];
rz(-0.91684219) q[2];
sx q[2];
rz(-2.4008304) q[2];
sx q[2];
rz(3.0286718) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0821655) q[1];
sx q[1];
rz(-1.1929999) q[1];
sx q[1];
rz(0.49281812) q[1];
x q[2];
rz(0.056315259) q[3];
sx q[3];
rz(-1.0259797) q[3];
sx q[3];
rz(2.44851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(2.4978499) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.165034) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(1.6595586) q[0];
rz(0.7011134) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(0.28465095) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50243044) q[0];
sx q[0];
rz(-2.6579755) q[0];
sx q[0];
rz(1.731427) q[0];
rz(-pi) q[1];
rz(2.1842723) q[2];
sx q[2];
rz(-2.7942604) q[2];
sx q[2];
rz(2.6383102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33672562) q[1];
sx q[1];
rz(-1.6972099) q[1];
sx q[1];
rz(-3.0625507) q[1];
x q[2];
rz(-0.48216362) q[3];
sx q[3];
rz(-2.1820222) q[3];
sx q[3];
rz(1.5507444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9099137) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(2.5740734) q[2];
rz(0.41401687) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48859566) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.8150785) q[0];
rz(1.1524221) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(0.93793905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6803166) q[0];
sx q[0];
rz(-0.094705908) q[0];
sx q[0];
rz(-1.3374431) q[0];
rz(-pi) q[1];
rz(0.51860923) q[2];
sx q[2];
rz(-2.0687639) q[2];
sx q[2];
rz(1.2155611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2873958) q[1];
sx q[1];
rz(-1.5970267) q[1];
sx q[1];
rz(1.6051952) q[1];
rz(-pi) q[2];
rz(-1.7004847) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(2.6854533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1753297) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-1.0413292) q[2];
rz(-1.8390309) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(1.3180102) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(2.5842216) q[0];
rz(-2.5769261) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(0.55647892) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0787449) q[0];
sx q[0];
rz(-1.9468465) q[0];
sx q[0];
rz(1.7929121) q[0];
rz(-2.9678992) q[2];
sx q[2];
rz(-2.2787893) q[2];
sx q[2];
rz(3.0922572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4887052) q[1];
sx q[1];
rz(-2.8294551) q[1];
sx q[1];
rz(-0.94241347) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5197093) q[3];
sx q[3];
rz(-2.3458977) q[3];
sx q[3];
rz(1.6076128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53211987) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(-1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0625793) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(-0.042908948) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(2.6409805) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096503784) q[0];
sx q[0];
rz(-1.6052264) q[0];
sx q[0];
rz(1.3929429) q[0];
x q[1];
rz(2.826564) q[2];
sx q[2];
rz(-1.050596) q[2];
sx q[2];
rz(-2.5556285) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2523633) q[1];
sx q[1];
rz(-1.6912618) q[1];
sx q[1];
rz(2.7389588) q[1];
rz(1.9686437) q[3];
sx q[3];
rz(-2.2865191) q[3];
sx q[3];
rz(-1.0637103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(0.44089857) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(0.52136695) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0764517) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(1.0674397) q[0];
rz(0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(-1.4656461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059793652) q[0];
sx q[0];
rz(-2.0311211) q[0];
sx q[0];
rz(-2.9914809) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52532105) q[2];
sx q[2];
rz(-0.21481951) q[2];
sx q[2];
rz(0.36052442) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2035985) q[1];
sx q[1];
rz(-1.9059056) q[1];
sx q[1];
rz(-1.0458228) q[1];
x q[2];
rz(0.64126863) q[3];
sx q[3];
rz(-1.0675758) q[3];
sx q[3];
rz(0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(2.2680797) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(-2.6749558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2575689) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(-0.018741477) q[1];
sx q[1];
rz(-0.32967162) q[1];
sx q[1];
rz(-0.92528701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7172456) q[0];
sx q[0];
rz(-2.2204917) q[0];
sx q[0];
rz(2.2109277) q[0];
x q[1];
rz(-1.7032353) q[2];
sx q[2];
rz(-0.91184154) q[2];
sx q[2];
rz(3.1140285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9097594) q[1];
sx q[1];
rz(-1.4405466) q[1];
sx q[1];
rz(0.18472437) q[1];
rz(-pi) q[2];
rz(-0.41947375) q[3];
sx q[3];
rz(-2.2443716) q[3];
sx q[3];
rz(0.89299612) q[3];
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
rz(-3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(-1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(3.1019679) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(0.75751799) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82523358) q[0];
sx q[0];
rz(-2.9037583) q[0];
sx q[0];
rz(-0.95038484) q[0];
x q[1];
rz(-1.8025814) q[2];
sx q[2];
rz(-2.5387562) q[2];
sx q[2];
rz(1.1834809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4841008) q[1];
sx q[1];
rz(-1.7925646) q[1];
sx q[1];
rz(-2.0055254) q[1];
rz(-0.032978756) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1511128) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.5007301) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(-0.43754229) q[2];
sx q[2];
rz(-1.8716639) q[2];
sx q[2];
rz(1.5834091) q[2];
rz(2.0414447) q[3];
sx q[3];
rz(-2.5082519) q[3];
sx q[3];
rz(-2.8520907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
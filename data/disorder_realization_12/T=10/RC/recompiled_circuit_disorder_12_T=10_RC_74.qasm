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
rz(-1.6568503) q[0];
rz(1.2031263) q[1];
sx q[1];
rz(-0.523518) q[1];
sx q[1];
rz(2.2533921) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48924016) q[0];
sx q[0];
rz(-1.7416735) q[0];
sx q[0];
rz(-1.6160374) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1076803) q[2];
sx q[2];
rz(-1.7134588) q[2];
sx q[2];
rz(-2.2906274) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2759309) q[1];
sx q[1];
rz(-1.7377503) q[1];
sx q[1];
rz(-1.5702412) q[1];
x q[2];
rz(-0.69859759) q[3];
sx q[3];
rz(-1.5961831) q[3];
sx q[3];
rz(0.7198402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28329864) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(1.1179914) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685101) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(0.65482393) q[0];
rz(-1.2163935) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(3.0156946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9585796) q[0];
sx q[0];
rz(-2.5218681) q[0];
sx q[0];
rz(2.241914) q[0];
rz(-2.5198031) q[2];
sx q[2];
rz(-2.1513373) q[2];
sx q[2];
rz(-0.34678005) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8996846) q[1];
sx q[1];
rz(-3.0113314) q[1];
sx q[1];
rz(-1.8036519) q[1];
rz(0.010766518) q[3];
sx q[3];
rz(-0.84732238) q[3];
sx q[3];
rz(2.7303498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.048432365) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(-1.4240501) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(0.58732906) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5858784) q[0];
sx q[0];
rz(-0.782574) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(-0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(1.1598587) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0214329) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(1.73818) q[0];
x q[1];
rz(-0.91684219) q[2];
sx q[2];
rz(-2.4008304) q[2];
sx q[2];
rz(-0.11292085) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0283454) q[1];
sx q[1];
rz(-0.61134895) q[1];
sx q[1];
rz(-2.4436185) q[1];
rz(-pi) q[2];
rz(-1.6634116) q[3];
sx q[3];
rz(-2.5941656) q[3];
sx q[3];
rz(0.80143354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(-2.4978499) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.9765587) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.6595586) q[0];
rz(-0.7011134) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(0.28465095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4581504) q[0];
sx q[0];
rz(-1.0939286) q[0];
sx q[0];
rz(-0.083806888) q[0];
x q[1];
rz(0.95732032) q[2];
sx q[2];
rz(-0.34733221) q[2];
sx q[2];
rz(-0.50328244) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8975361) q[1];
sx q[1];
rz(-1.6492062) q[1];
sx q[1];
rz(-1.4439911) q[1];
rz(-pi) q[2];
rz(-2.2399726) q[3];
sx q[3];
rz(-1.1812783) q[3];
sx q[3];
rz(-0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9099137) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.652997) q[0];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91468231) q[0];
sx q[0];
rz(-1.662928) q[0];
sx q[0];
rz(0.021962086) q[0];
rz(-1.0114952) q[2];
sx q[2];
rz(-2.0213631) q[2];
sx q[2];
rz(-3.0524888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2873958) q[1];
sx q[1];
rz(-1.5445659) q[1];
sx q[1];
rz(1.5363974) q[1];
x q[2];
rz(-1.4411079) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(-2.6854533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1753297) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(-1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.702521) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(2.5842216) q[0];
rz(-0.56466651) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-0.55647892) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508931) q[0];
sx q[0];
rz(-1.7771582) q[0];
sx q[0];
rz(0.38462374) q[0];
x q[1];
rz(0.17369341) q[2];
sx q[2];
rz(-0.86280338) q[2];
sx q[2];
rz(0.049335418) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0007799) q[1];
sx q[1];
rz(-1.3197348) q[1];
sx q[1];
rz(2.954133) q[1];
x q[2];
rz(3.0895124) q[3];
sx q[3];
rz(-0.77643231) q[3];
sx q[3];
rz(1.4610425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(-1.9167985) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.640655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0625793) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(-0.042908948) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(-0.50061217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4780873) q[0];
sx q[0];
rz(-2.9604719) q[0];
sx q[0];
rz(1.3785133) q[0];
x q[1];
rz(-2.0666276) q[2];
sx q[2];
rz(-2.5410286) q[2];
sx q[2];
rz(1.1662837) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88922933) q[1];
sx q[1];
rz(-1.4503308) q[1];
sx q[1];
rz(-2.7389588) q[1];
rz(-pi) q[2];
rz(-2.385419) q[3];
sx q[3];
rz(-1.2740967) q[3];
sx q[3];
rz(-2.365436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5381955) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(0.67561692) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(2.7087129) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(-1.6759466) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.081799) q[0];
sx q[0];
rz(-2.0311211) q[0];
sx q[0];
rz(2.9914809) q[0];
rz(-pi) q[1];
rz(0.52532105) q[2];
sx q[2];
rz(-0.21481951) q[2];
sx q[2];
rz(2.7810682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2035985) q[1];
sx q[1];
rz(-1.2356871) q[1];
sx q[1];
rz(-1.0458228) q[1];
x q[2];
rz(-2.500324) q[3];
sx q[3];
rz(-2.0740168) q[3];
sx q[3];
rz(2.7748231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(-2.2680797) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(2.6749558) q[3];
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
rz(-1.8840238) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(0.73721686) q[0];
rz(0.018741477) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(2.2163056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46426526) q[0];
sx q[0];
rz(-0.87809169) q[0];
sx q[0];
rz(2.4753184) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7032353) q[2];
sx q[2];
rz(-0.91184154) q[2];
sx q[2];
rz(-3.1140285) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9097594) q[1];
sx q[1];
rz(-1.4405466) q[1];
sx q[1];
rz(-0.18472437) q[1];
rz(-pi) q[2];
rz(-2.7221189) q[3];
sx q[3];
rz(-0.89722108) q[3];
sx q[3];
rz(0.89299612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(1.0037237) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039624778) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(1.2596624) q[0];
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
rz(-1.3525317) q[0];
sx q[0];
rz(-1.4333945) q[0];
sx q[0];
rz(-1.3760516) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9847758) q[2];
sx q[2];
rz(-0.98625253) q[2];
sx q[2];
rz(-2.2371694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0150891) q[1];
sx q[1];
rz(-1.9941829) q[1];
sx q[1];
rz(2.8979315) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.032978756) q[3];
sx q[3];
rz(-0.71027256) q[3];
sx q[3];
rz(0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1511128) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(2.347836) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(2.5095148) q[2];
sx q[2];
rz(-2.6161604) q[2];
sx q[2];
rz(-2.5642774) q[2];
rz(-0.32140857) q[3];
sx q[3];
rz(-2.1264429) q[3];
sx q[3];
rz(0.85254729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(-0.40894142) q[0];
sx q[0];
rz(2.1483108) q[0];
rz(2.9945057) q[1];
sx q[1];
rz(-2.1451201) q[1];
sx q[1];
rz(1.7239404) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2679358) q[0];
sx q[0];
rz(-1.2136154) q[0];
sx q[0];
rz(2.3620679) q[0];
rz(-pi) q[1];
rz(1.0101914) q[2];
sx q[2];
rz(-0.94807887) q[2];
sx q[2];
rz(-1.6499008) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65517814) q[1];
sx q[1];
rz(-1.298212) q[1];
sx q[1];
rz(-2.4943939) q[1];
rz(-pi) q[2];
rz(-1.8090114) q[3];
sx q[3];
rz(-3.0034667) q[3];
sx q[3];
rz(-2.3254243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.950497) q[2];
sx q[2];
rz(-1.3206626) q[2];
sx q[2];
rz(2.0868059) q[2];
rz(-2.6705006) q[3];
sx q[3];
rz(-1.6867009) q[3];
sx q[3];
rz(-0.86827046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3926113) q[0];
sx q[0];
rz(-1.457021) q[0];
sx q[0];
rz(2.9066322) q[0];
rz(-1.2745534) q[1];
sx q[1];
rz(-1.3095368) q[1];
sx q[1];
rz(-2.4539006) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0318555) q[0];
sx q[0];
rz(-2.8068455) q[0];
sx q[0];
rz(2.9358923) q[0];
rz(-0.85854395) q[2];
sx q[2];
rz(-1.1919824) q[2];
sx q[2];
rz(1.0114947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5357757) q[1];
sx q[1];
rz(-1.5763842) q[1];
sx q[1];
rz(-1.5833012) q[1];
rz(3.0491203) q[3];
sx q[3];
rz(-0.95917976) q[3];
sx q[3];
rz(-0.05449748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4253652) q[2];
sx q[2];
rz(-1.3024412) q[2];
sx q[2];
rz(2.9322374) q[2];
rz(-0.9730722) q[3];
sx q[3];
rz(-0.89997411) q[3];
sx q[3];
rz(0.59276855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3000325) q[0];
sx q[0];
rz(-2.7356) q[0];
sx q[0];
rz(2.9969065) q[0];
rz(1.2380098) q[1];
sx q[1];
rz(-2.369945) q[1];
sx q[1];
rz(3.0498116) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84052575) q[0];
sx q[0];
rz(-0.65945461) q[0];
sx q[0];
rz(-0.98697885) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94987671) q[2];
sx q[2];
rz(-0.11813049) q[2];
sx q[2];
rz(1.0850414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8615177) q[1];
sx q[1];
rz(-2.8535378) q[1];
sx q[1];
rz(-2.3029598) q[1];
rz(1.1385659) q[3];
sx q[3];
rz(-1.5097268) q[3];
sx q[3];
rz(1.7342126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.38516513) q[2];
sx q[2];
rz(-1.5084718) q[2];
sx q[2];
rz(-0.37919322) q[2];
rz(2.3181629) q[3];
sx q[3];
rz(-0.40612602) q[3];
sx q[3];
rz(2.8520975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8150811) q[0];
sx q[0];
rz(-0.87711763) q[0];
sx q[0];
rz(2.1429578) q[0];
rz(-2.6594992) q[1];
sx q[1];
rz(-0.95078743) q[1];
sx q[1];
rz(-1.3179717) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40967746) q[0];
sx q[0];
rz(-1.9008667) q[0];
sx q[0];
rz(2.5797352) q[0];
x q[1];
rz(-2.5561294) q[2];
sx q[2];
rz(-1.778086) q[2];
sx q[2];
rz(-2.6475865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.98612477) q[1];
sx q[1];
rz(-1.7805702) q[1];
sx q[1];
rz(-2.5986157) q[1];
rz(-pi) q[2];
x q[2];
rz(2.34818) q[3];
sx q[3];
rz(-1.9559091) q[3];
sx q[3];
rz(-0.18369964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92129293) q[2];
sx q[2];
rz(-0.90196323) q[2];
sx q[2];
rz(-2.656142) q[2];
rz(-2.7576533) q[3];
sx q[3];
rz(-1.5860312) q[3];
sx q[3];
rz(0.78401047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601198) q[0];
sx q[0];
rz(-3.0373242) q[0];
sx q[0];
rz(0.01165788) q[0];
rz(-3.0043789) q[1];
sx q[1];
rz(-1.5882746) q[1];
sx q[1];
rz(-1.8395909) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6377651) q[0];
sx q[0];
rz(-1.3877467) q[0];
sx q[0];
rz(-1.4603015) q[0];
x q[1];
rz(-0.51568337) q[2];
sx q[2];
rz(-2.1123675) q[2];
sx q[2];
rz(2.6560419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4123789) q[1];
sx q[1];
rz(-1.679031) q[1];
sx q[1];
rz(1.8844876) q[1];
x q[2];
rz(2.7632159) q[3];
sx q[3];
rz(-1.9073745) q[3];
sx q[3];
rz(2.5314296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4181218) q[2];
sx q[2];
rz(-1.301441) q[2];
sx q[2];
rz(-0.29329014) q[2];
rz(0.063118525) q[3];
sx q[3];
rz(-1.2226353) q[3];
sx q[3];
rz(-2.2466834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1191331) q[0];
sx q[0];
rz(-0.28894153) q[0];
sx q[0];
rz(-0.038507842) q[0];
rz(2.0571845) q[1];
sx q[1];
rz(-0.42250982) q[1];
sx q[1];
rz(0.96673059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16822505) q[0];
sx q[0];
rz(-1.7666139) q[0];
sx q[0];
rz(1.7652579) q[0];
x q[1];
rz(-1.3036846) q[2];
sx q[2];
rz(-1.583364) q[2];
sx q[2];
rz(-1.723701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.048592) q[1];
sx q[1];
rz(-1.9680259) q[1];
sx q[1];
rz(-2.7723958) q[1];
rz(-pi) q[2];
rz(-0.24992515) q[3];
sx q[3];
rz(-0.93533134) q[3];
sx q[3];
rz(2.9912586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66095573) q[2];
sx q[2];
rz(-2.9065242) q[2];
sx q[2];
rz(0.98141518) q[2];
rz(1.4977945) q[3];
sx q[3];
rz(-1.8794329) q[3];
sx q[3];
rz(1.5952544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82866955) q[0];
sx q[0];
rz(-2.4607615) q[0];
sx q[0];
rz(2.8269826) q[0];
rz(-0.18868748) q[1];
sx q[1];
rz(-0.43470946) q[1];
sx q[1];
rz(-0.46868086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4026827) q[0];
sx q[0];
rz(-0.27640585) q[0];
sx q[0];
rz(-2.6457647) q[0];
rz(-pi) q[1];
rz(-1.5596703) q[2];
sx q[2];
rz(-0.16777786) q[2];
sx q[2];
rz(-2.0027225) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0464161) q[1];
sx q[1];
rz(-0.79907387) q[1];
sx q[1];
rz(-2.0002736) q[1];
rz(-2.5545337) q[3];
sx q[3];
rz(-0.54577845) q[3];
sx q[3];
rz(3.0445638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9357052) q[2];
sx q[2];
rz(-2.2998655) q[2];
sx q[2];
rz(-0.97071281) q[2];
rz(0.5361706) q[3];
sx q[3];
rz(-1.2090809) q[3];
sx q[3];
rz(2.4255588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4416874) q[0];
sx q[0];
rz(-0.55460414) q[0];
sx q[0];
rz(1.6957138) q[0];
rz(2.1031759) q[1];
sx q[1];
rz(-2.0380135) q[1];
sx q[1];
rz(-3.0605002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38302416) q[0];
sx q[0];
rz(-1.283535) q[0];
sx q[0];
rz(-0.69604915) q[0];
x q[1];
rz(2.7620188) q[2];
sx q[2];
rz(-0.69131572) q[2];
sx q[2];
rz(-2.8065681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1782185) q[1];
sx q[1];
rz(-1.4556754) q[1];
sx q[1];
rz(-1.2769775) q[1];
rz(-1.1015973) q[3];
sx q[3];
rz(-2.0939671) q[3];
sx q[3];
rz(-1.2468456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1991835) q[2];
sx q[2];
rz(-2.3641391) q[2];
sx q[2];
rz(2.8988163) q[2];
rz(-1.5593922) q[3];
sx q[3];
rz(-0.42583164) q[3];
sx q[3];
rz(1.2529713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2389857) q[0];
sx q[0];
rz(-2.1098397) q[0];
sx q[0];
rz(1.8865939) q[0];
rz(-1.3764508) q[1];
sx q[1];
rz(-1.0245198) q[1];
sx q[1];
rz(2.4741516) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8534626) q[0];
sx q[0];
rz(-1.9024182) q[0];
sx q[0];
rz(1.0378855) q[0];
rz(-1.6133283) q[2];
sx q[2];
rz(-1.0746403) q[2];
sx q[2];
rz(-0.68788487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0170198) q[1];
sx q[1];
rz(-1.8609398) q[1];
sx q[1];
rz(0.97495074) q[1];
rz(-2.7605704) q[3];
sx q[3];
rz(-1.2667315) q[3];
sx q[3];
rz(-1.9532494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5218375) q[2];
sx q[2];
rz(-2.4027368) q[2];
sx q[2];
rz(-0.45836207) q[2];
rz(-1.6058263) q[3];
sx q[3];
rz(-1.0777487) q[3];
sx q[3];
rz(-2.2569236) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89852029) q[0];
sx q[0];
rz(-2.5694818) q[0];
sx q[0];
rz(2.7761053) q[0];
rz(-2.1416523) q[1];
sx q[1];
rz(-1.4391856) q[1];
sx q[1];
rz(1.684729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3709049) q[0];
sx q[0];
rz(-1.7720776) q[0];
sx q[0];
rz(1.007249) q[0];
x q[1];
rz(1.1816977) q[2];
sx q[2];
rz(-0.28841296) q[2];
sx q[2];
rz(-0.696495) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36023203) q[1];
sx q[1];
rz(-1.7443027) q[1];
sx q[1];
rz(0.75299112) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54014628) q[3];
sx q[3];
rz(-0.77278256) q[3];
sx q[3];
rz(2.5908822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0008056) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(-1.0055536) q[2];
rz(-0.65070659) q[3];
sx q[3];
rz(-1.7020099) q[3];
sx q[3];
rz(-2.2910291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477796) q[0];
sx q[0];
rz(-0.78421264) q[0];
sx q[0];
rz(-1.0017851) q[0];
rz(1.7306937) q[1];
sx q[1];
rz(-1.7041364) q[1];
sx q[1];
rz(-2.0028353) q[1];
rz(-1.9024814) q[2];
sx q[2];
rz(-2.7948236) q[2];
sx q[2];
rz(0.75955392) q[2];
rz(-1.6639573) q[3];
sx q[3];
rz(-1.6747083) q[3];
sx q[3];
rz(-1.2887521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

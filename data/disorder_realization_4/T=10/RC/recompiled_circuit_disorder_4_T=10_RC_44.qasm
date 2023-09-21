OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(2.7058869) q[0];
sx q[0];
rz(11.640179) q[0];
rz(1.9594833) q[1];
sx q[1];
rz(-0.73298454) q[1];
sx q[1];
rz(-2.7690673) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6719138) q[0];
sx q[0];
rz(-0.9979453) q[0];
sx q[0];
rz(-1.3290622) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4650605) q[2];
sx q[2];
rz(-2.1991962) q[2];
sx q[2];
rz(0.65199967) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8718611) q[1];
sx q[1];
rz(-1.6850123) q[1];
sx q[1];
rz(1.7608587) q[1];
rz(1.094369) q[3];
sx q[3];
rz(-0.24260394) q[3];
sx q[3];
rz(2.0403595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42840502) q[2];
sx q[2];
rz(-1.5822072) q[2];
sx q[2];
rz(-0.92450809) q[2];
rz(-1.472578) q[3];
sx q[3];
rz(-1.2481097) q[3];
sx q[3];
rz(1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18838841) q[0];
sx q[0];
rz(-0.062347118) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(-2.9470782) q[1];
sx q[1];
rz(-1.3214) q[1];
sx q[1];
rz(0.054873437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48557278) q[0];
sx q[0];
rz(-1.4783854) q[0];
sx q[0];
rz(-1.4569605) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61972159) q[2];
sx q[2];
rz(-2.0267068) q[2];
sx q[2];
rz(-1.0085307) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3000852) q[1];
sx q[1];
rz(-1.9133139) q[1];
sx q[1];
rz(3.1373346) q[1];
rz(-0.26344928) q[3];
sx q[3];
rz(-2.4512495) q[3];
sx q[3];
rz(3.1171947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7028971) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(1.4820209) q[2];
rz(-2.1510018) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(-1.3668485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(2.7666132) q[0];
rz(0.24770501) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(1.8050271) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5216615) q[0];
sx q[0];
rz(-1.7921899) q[0];
sx q[0];
rz(-3.0822166) q[0];
rz(-pi) q[1];
rz(1.1281563) q[2];
sx q[2];
rz(-2.1146362) q[2];
sx q[2];
rz(-1.7511055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6429236) q[1];
sx q[1];
rz(-1.8734697) q[1];
sx q[1];
rz(-2.32294) q[1];
rz(-pi) q[2];
rz(-2.7258354) q[3];
sx q[3];
rz(-2.3299179) q[3];
sx q[3];
rz(-0.31879253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(-1.0245727) q[2];
rz(1.2767977) q[3];
sx q[3];
rz(-1.5494616) q[3];
sx q[3];
rz(-1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17722002) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(-1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(1.4039325) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6934024) q[0];
sx q[0];
rz(-1.9172137) q[0];
sx q[0];
rz(1.3643826) q[0];
rz(0.83909859) q[2];
sx q[2];
rz(-2.1144146) q[2];
sx q[2];
rz(-2.3630138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2851706) q[1];
sx q[1];
rz(-0.31177786) q[1];
sx q[1];
rz(-1.9429893) q[1];
x q[2];
rz(1.674026) q[3];
sx q[3];
rz(-2.4475615) q[3];
sx q[3];
rz(0.3895143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6491062) q[2];
sx q[2];
rz(-0.84473574) q[2];
sx q[2];
rz(-0.21326324) q[2];
rz(3.1212741) q[3];
sx q[3];
rz(-1.3213986) q[3];
sx q[3];
rz(-0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7810818) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(-1.0158585) q[0];
rz(1.5083183) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(0.034084592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8071472) q[0];
sx q[0];
rz(-2.7140518) q[0];
sx q[0];
rz(2.7464675) q[0];
rz(-pi) q[1];
rz(2.7427865) q[2];
sx q[2];
rz(-0.6119234) q[2];
sx q[2];
rz(1.2955701) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7923813) q[1];
sx q[1];
rz(-1.5776909) q[1];
sx q[1];
rz(1.2864119) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50877737) q[3];
sx q[3];
rz(-2.1366589) q[3];
sx q[3];
rz(0.98154991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4013227) q[2];
sx q[2];
rz(-2.5677742) q[2];
sx q[2];
rz(-1.9704698) q[2];
rz(1.8088388) q[3];
sx q[3];
rz(-1.4274024) q[3];
sx q[3];
rz(-3.0122421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68833441) q[0];
sx q[0];
rz(-1.9535221) q[0];
sx q[0];
rz(-2.7057498) q[0];
rz(2.0360937) q[1];
sx q[1];
rz(-1.2970129) q[1];
sx q[1];
rz(1.2058535) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379369) q[0];
sx q[0];
rz(-2.3559542) q[0];
sx q[0];
rz(1.2334137) q[0];
x q[1];
rz(-1.4044936) q[2];
sx q[2];
rz(-3.0450833) q[2];
sx q[2];
rz(-1.4788747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44991325) q[1];
sx q[1];
rz(-2.0813585) q[1];
sx q[1];
rz(-1.6010067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.033127012) q[3];
sx q[3];
rz(-1.9702385) q[3];
sx q[3];
rz(1.2101733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.897052) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(-3.0899866) q[2];
rz(-2.7339325) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(-2.6962962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62149858) q[0];
sx q[0];
rz(-0.61722732) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(2.3576221) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(-0.53692445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2818031) q[0];
sx q[0];
rz(-1.6543232) q[0];
sx q[0];
rz(-3.0755755) q[0];
x q[1];
rz(2.0836673) q[2];
sx q[2];
rz(-0.29209902) q[2];
sx q[2];
rz(0.7823173) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9432224) q[1];
sx q[1];
rz(-2.3343625) q[1];
sx q[1];
rz(3.0109349) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7241931) q[3];
sx q[3];
rz(-0.985257) q[3];
sx q[3];
rz(-0.034754001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9817104) q[2];
sx q[2];
rz(-1.2572224) q[2];
sx q[2];
rz(0.19101492) q[2];
rz(-2.8602709) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087588) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(-0.38761815) q[0];
rz(-0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96466366) q[0];
sx q[0];
rz(-1.7690072) q[0];
sx q[0];
rz(-2.6508509) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6667716) q[2];
sx q[2];
rz(-0.64138597) q[2];
sx q[2];
rz(-2.1542187) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2355289) q[1];
sx q[1];
rz(-1.7019094) q[1];
sx q[1];
rz(0.99793418) q[1];
x q[2];
rz(-2.3861305) q[3];
sx q[3];
rz(-1.5243693) q[3];
sx q[3];
rz(1.69343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9541624) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(1.8743275) q[2];
rz(-0.77504843) q[3];
sx q[3];
rz(-1.6036443) q[3];
sx q[3];
rz(2.3601941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.18701126) q[0];
sx q[0];
rz(-2.1033852) q[0];
sx q[0];
rz(-2.9206081) q[0];
rz(2.2194608) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(2.1386713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98175883) q[0];
sx q[0];
rz(-2.8528385) q[0];
sx q[0];
rz(-0.50555484) q[0];
x q[1];
rz(-2.2641364) q[2];
sx q[2];
rz(-0.71091953) q[2];
sx q[2];
rz(0.029821776) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50833118) q[1];
sx q[1];
rz(-2.0610626) q[1];
sx q[1];
rz(2.0983178) q[1];
x q[2];
rz(-2.8691611) q[3];
sx q[3];
rz(-0.90297943) q[3];
sx q[3];
rz(-1.9200793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(2.9830902) q[2];
rz(-0.45378271) q[3];
sx q[3];
rz(-0.78275371) q[3];
sx q[3];
rz(-2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52699387) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(2.9504543) q[0];
rz(-0.29516164) q[1];
sx q[1];
rz(-0.8894397) q[1];
sx q[1];
rz(2.2492762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.17642) q[0];
sx q[0];
rz(-1.3942379) q[0];
sx q[0];
rz(-0.82133349) q[0];
rz(0.92992444) q[2];
sx q[2];
rz(-1.6460437) q[2];
sx q[2];
rz(1.6844695) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82980624) q[1];
sx q[1];
rz(-1.3327507) q[1];
sx q[1];
rz(1.5029552) q[1];
rz(-2.7097706) q[3];
sx q[3];
rz(-1.8027657) q[3];
sx q[3];
rz(0.98917978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7717379) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(-0.85956335) q[2];
rz(-1.2236979) q[3];
sx q[3];
rz(-2.2151291) q[3];
sx q[3];
rz(-0.74444509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201465) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(-1.6620811) q[1];
sx q[1];
rz(-0.48304396) q[1];
sx q[1];
rz(1.2190291) q[1];
rz(-0.027364846) q[2];
sx q[2];
rz(-1.27925) q[2];
sx q[2];
rz(-1.7640511) q[2];
rz(-1.5020694) q[3];
sx q[3];
rz(-1.6106265) q[3];
sx q[3];
rz(2.0911218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

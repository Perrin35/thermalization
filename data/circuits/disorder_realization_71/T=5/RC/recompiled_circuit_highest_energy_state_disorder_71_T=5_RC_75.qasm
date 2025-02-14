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
rz(1.1310391) q[0];
sx q[0];
rz(-0.57589632) q[0];
sx q[0];
rz(1.1269215) q[0];
rz(-0.95070401) q[1];
sx q[1];
rz(-0.60125142) q[1];
sx q[1];
rz(-0.33369219) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93984825) q[0];
sx q[0];
rz(-0.44838312) q[0];
sx q[0];
rz(-1.9223619) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0804685) q[2];
sx q[2];
rz(-2.6303419) q[2];
sx q[2];
rz(1.0488269) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0036949) q[1];
sx q[1];
rz(-2.7337791) q[1];
sx q[1];
rz(-0.56038709) q[1];
rz(-pi) q[2];
rz(0.80218704) q[3];
sx q[3];
rz(-1.7628551) q[3];
sx q[3];
rz(-3.1046228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8814964) q[2];
sx q[2];
rz(-2.0659955) q[2];
sx q[2];
rz(2.0913701) q[2];
rz(0.29551926) q[3];
sx q[3];
rz(-0.76885709) q[3];
sx q[3];
rz(-1.3692921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722026) q[0];
sx q[0];
rz(-2.1222332) q[0];
sx q[0];
rz(2.4743359) q[0];
rz(-0.15070209) q[1];
sx q[1];
rz(-1.0548016) q[1];
sx q[1];
rz(-2.8299455) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0829179) q[0];
sx q[0];
rz(-2.6150515) q[0];
sx q[0];
rz(-1.7114) q[0];
rz(-pi) q[1];
rz(0.30782757) q[2];
sx q[2];
rz(-1.1500949) q[2];
sx q[2];
rz(-1.6604648) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2978486) q[1];
sx q[1];
rz(-2.9103855) q[1];
sx q[1];
rz(-0.65194269) q[1];
rz(-pi) q[2];
rz(0.39763173) q[3];
sx q[3];
rz(-1.5036229) q[3];
sx q[3];
rz(-2.6844418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.47505891) q[2];
sx q[2];
rz(-0.89343137) q[2];
sx q[2];
rz(0.7676355) q[2];
rz(2.660699) q[3];
sx q[3];
rz(-2.3700263) q[3];
sx q[3];
rz(-1.5546999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84592205) q[0];
sx q[0];
rz(-1.6039811) q[0];
sx q[0];
rz(2.3620102) q[0];
rz(-2.7698611) q[1];
sx q[1];
rz(-2.1340243) q[1];
sx q[1];
rz(-2.380611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26139382) q[0];
sx q[0];
rz(-0.43873271) q[0];
sx q[0];
rz(-2.3068271) q[0];
rz(-pi) q[1];
x q[1];
rz(2.089431) q[2];
sx q[2];
rz(-2.385879) q[2];
sx q[2];
rz(1.99988) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8345388) q[1];
sx q[1];
rz(-2.4014856) q[1];
sx q[1];
rz(-1.8174875) q[1];
x q[2];
rz(-0.85270057) q[3];
sx q[3];
rz(-1.5670683) q[3];
sx q[3];
rz(1.119233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3373105) q[2];
sx q[2];
rz(-2.2411942) q[2];
sx q[2];
rz(2.6864181) q[2];
rz(0.69194397) q[3];
sx q[3];
rz(-2.4271836) q[3];
sx q[3];
rz(0.80879319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96875018) q[0];
sx q[0];
rz(-1.7946365) q[0];
sx q[0];
rz(2.8853048) q[0];
rz(1.7068656) q[1];
sx q[1];
rz(-1.4204357) q[1];
sx q[1];
rz(-0.11875471) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0066322346) q[0];
sx q[0];
rz(-1.7259571) q[0];
sx q[0];
rz(-0.4520024) q[0];
x q[1];
rz(-1.8639455) q[2];
sx q[2];
rz(-1.4481232) q[2];
sx q[2];
rz(-0.49201595) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12417158) q[1];
sx q[1];
rz(-1.3391437) q[1];
sx q[1];
rz(0.79504658) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3975735) q[3];
sx q[3];
rz(-2.3459917) q[3];
sx q[3];
rz(0.52594409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8541096) q[2];
sx q[2];
rz(-0.27366769) q[2];
sx q[2];
rz(-1.0556833) q[2];
rz(-1.4500729) q[3];
sx q[3];
rz(-1.6943211) q[3];
sx q[3];
rz(0.84901866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5797822) q[0];
sx q[0];
rz(-0.17794839) q[0];
sx q[0];
rz(1.6135038) q[0];
rz(-1.1771419) q[1];
sx q[1];
rz(-1.2700932) q[1];
sx q[1];
rz(1.5783763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10229853) q[0];
sx q[0];
rz(-1.6400518) q[0];
sx q[0];
rz(-2.5152492) q[0];
x q[1];
rz(-1.2424395) q[2];
sx q[2];
rz(-1.9617726) q[2];
sx q[2];
rz(2.7308488) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0381752) q[1];
sx q[1];
rz(-0.6486054) q[1];
sx q[1];
rz(-0.23596779) q[1];
rz(1.568119) q[3];
sx q[3];
rz(-1.9802046) q[3];
sx q[3];
rz(-2.1684627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13581181) q[2];
sx q[2];
rz(-1.1567189) q[2];
sx q[2];
rz(-1.8394252) q[2];
rz(2.8016413) q[3];
sx q[3];
rz(-0.8067185) q[3];
sx q[3];
rz(2.4792041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.0078916773) q[0];
sx q[0];
rz(-1.5971203) q[0];
sx q[0];
rz(1.1997461) q[0];
rz(-1.3613191) q[1];
sx q[1];
rz(-0.97313762) q[1];
sx q[1];
rz(2.5637085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21450522) q[0];
sx q[0];
rz(-1.5772235) q[0];
sx q[0];
rz(-1.584573) q[0];
x q[1];
rz(1.4125713) q[2];
sx q[2];
rz(-1.998868) q[2];
sx q[2];
rz(-0.23382631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.84915811) q[1];
sx q[1];
rz(-2.6657988) q[1];
sx q[1];
rz(-2.8413296) q[1];
rz(-pi) q[2];
rz(-3.1158524) q[3];
sx q[3];
rz(-1.8696068) q[3];
sx q[3];
rz(2.4853277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5100688) q[2];
sx q[2];
rz(-0.57454595) q[2];
sx q[2];
rz(2.6681382) q[2];
rz(1.8691285) q[3];
sx q[3];
rz(-1.748964) q[3];
sx q[3];
rz(-0.22245358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94721395) q[0];
sx q[0];
rz(-2.4838303) q[0];
sx q[0];
rz(2.8084602) q[0];
rz(-0.22854742) q[1];
sx q[1];
rz(-1.3401778) q[1];
sx q[1];
rz(0.57565912) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7299907) q[0];
sx q[0];
rz(-1.4648998) q[0];
sx q[0];
rz(-0.81231711) q[0];
x q[1];
rz(0.47708738) q[2];
sx q[2];
rz(-1.2735575) q[2];
sx q[2];
rz(-1.5444259) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.043283431) q[1];
sx q[1];
rz(-0.38615882) q[1];
sx q[1];
rz(-1.015806) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6603465) q[3];
sx q[3];
rz(-1.3433045) q[3];
sx q[3];
rz(-0.15531048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8754742) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(1.4761338) q[2];
rz(-2.0406145) q[3];
sx q[3];
rz(-1.063238) q[3];
sx q[3];
rz(0.5180009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5975033) q[0];
sx q[0];
rz(-2.2552555) q[0];
sx q[0];
rz(0.30712095) q[0];
rz(-0.33044526) q[1];
sx q[1];
rz(-1.5978866) q[1];
sx q[1];
rz(1.365136) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86548818) q[0];
sx q[0];
rz(-1.5886602) q[0];
sx q[0];
rz(-0.081327036) q[0];
rz(-pi) q[1];
rz(-1.3826755) q[2];
sx q[2];
rz(-1.9131919) q[2];
sx q[2];
rz(-1.2255648) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4933984) q[1];
sx q[1];
rz(-1.1432363) q[1];
sx q[1];
rz(-2.8570064) q[1];
rz(-pi) q[2];
rz(-1.4151814) q[3];
sx q[3];
rz(-2.3978516) q[3];
sx q[3];
rz(-1.9187601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.27552989) q[2];
sx q[2];
rz(-2.5538462) q[2];
sx q[2];
rz(-0.24869871) q[2];
rz(1.9281049) q[3];
sx q[3];
rz(-1.4444193) q[3];
sx q[3];
rz(3.0799589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.98951775) q[0];
sx q[0];
rz(-2.2353421) q[0];
sx q[0];
rz(2.4549947) q[0];
rz(1.1129414) q[1];
sx q[1];
rz(-2.3390892) q[1];
sx q[1];
rz(-0.31563219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9292727) q[0];
sx q[0];
rz(-1.2480191) q[0];
sx q[0];
rz(-3.0930711) q[0];
x q[1];
rz(-1.8191843) q[2];
sx q[2];
rz(-0.82652107) q[2];
sx q[2];
rz(-2.7374817) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13309578) q[1];
sx q[1];
rz(-1.6649455) q[1];
sx q[1];
rz(-0.90742438) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6417362) q[3];
sx q[3];
rz(-0.8345064) q[3];
sx q[3];
rz(0.77674538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.687872) q[2];
sx q[2];
rz(-2.5747955) q[2];
sx q[2];
rz(-2.4570214) q[2];
rz(-1.941393) q[3];
sx q[3];
rz(-2.0122416) q[3];
sx q[3];
rz(2.9620192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0803364) q[0];
sx q[0];
rz(-0.48114023) q[0];
sx q[0];
rz(-3.1367593) q[0];
rz(-1.6239411) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(0.25371107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2462354) q[0];
sx q[0];
rz(-2.9140436) q[0];
sx q[0];
rz(-1.1832817) q[0];
rz(-2.0433304) q[2];
sx q[2];
rz(-0.8923549) q[2];
sx q[2];
rz(-1.1393762) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1349134) q[1];
sx q[1];
rz(-0.90613922) q[1];
sx q[1];
rz(2.8921739) q[1];
x q[2];
rz(1.7697479) q[3];
sx q[3];
rz(-1.0291489) q[3];
sx q[3];
rz(-1.0197848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97810513) q[2];
sx q[2];
rz(-1.2349962) q[2];
sx q[2];
rz(1.9592436) q[2];
rz(-0.67665082) q[3];
sx q[3];
rz(-0.38656056) q[3];
sx q[3];
rz(1.0138018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5273298) q[0];
sx q[0];
rz(-2.3557721) q[0];
sx q[0];
rz(-2.8489805) q[0];
rz(1.0489427) q[1];
sx q[1];
rz(-1.4194149) q[1];
sx q[1];
rz(0.70379757) q[1];
rz(0.47164698) q[2];
sx q[2];
rz(-1.2670965) q[2];
sx q[2];
rz(-0.22128174) q[2];
rz(1.9378035) q[3];
sx q[3];
rz(-1.6631775) q[3];
sx q[3];
rz(0.1803995) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

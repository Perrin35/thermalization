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
rz(-2.3075624) q[0];
sx q[0];
rz(-0.40606719) q[0];
sx q[0];
rz(-1.9598444) q[0];
rz(-1.4015247) q[1];
sx q[1];
rz(-2.6670246) q[1];
sx q[1];
rz(-0.13303703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0358462) q[0];
sx q[0];
rz(-1.8324569) q[0];
sx q[0];
rz(-0.91069551) q[0];
x q[1];
rz(0.70212097) q[2];
sx q[2];
rz(-2.5062541) q[2];
sx q[2];
rz(-1.8566657) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5488447) q[1];
sx q[1];
rz(-1.8858375) q[1];
sx q[1];
rz(2.9058019) q[1];
rz(-pi) q[2];
rz(-2.7772589) q[3];
sx q[3];
rz(-1.0291463) q[3];
sx q[3];
rz(2.2026103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0066321) q[2];
sx q[2];
rz(-1.0198318) q[2];
sx q[2];
rz(2.4572065) q[2];
rz(2.0002174) q[3];
sx q[3];
rz(-0.26641521) q[3];
sx q[3];
rz(-0.086708955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619277) q[0];
sx q[0];
rz(-2.4517038) q[0];
sx q[0];
rz(2.6806114) q[0];
rz(2.0224679) q[1];
sx q[1];
rz(-2.7179167) q[1];
sx q[1];
rz(0.31203312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49544612) q[0];
sx q[0];
rz(-1.5428041) q[0];
sx q[0];
rz(-1.5585446) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4909541) q[2];
sx q[2];
rz(-2.3791358) q[2];
sx q[2];
rz(1.2079639) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65822433) q[1];
sx q[1];
rz(-1.3951571) q[1];
sx q[1];
rz(3.1356922) q[1];
rz(2.2226238) q[3];
sx q[3];
rz(-1.7953835) q[3];
sx q[3];
rz(-2.74204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95388874) q[2];
sx q[2];
rz(-2.0362594) q[2];
sx q[2];
rz(0.40404955) q[2];
rz(-0.36014253) q[3];
sx q[3];
rz(-1.4934243) q[3];
sx q[3];
rz(2.4586316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78524154) q[0];
sx q[0];
rz(-2.1389565) q[0];
sx q[0];
rz(2.8175765) q[0];
rz(2.0815381) q[1];
sx q[1];
rz(-2.8476604) q[1];
sx q[1];
rz(2.4411328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4840317) q[0];
sx q[0];
rz(-0.23398384) q[0];
sx q[0];
rz(-2.3768539) q[0];
x q[1];
rz(2.926658) q[2];
sx q[2];
rz(-1.2329458) q[2];
sx q[2];
rz(-0.95906729) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30875242) q[1];
sx q[1];
rz(-2.3165417) q[1];
sx q[1];
rz(1.8925335) q[1];
rz(-pi) q[2];
rz(-0.19566124) q[3];
sx q[3];
rz(-1.1563468) q[3];
sx q[3];
rz(-2.5959542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5141653) q[2];
sx q[2];
rz(-0.24719396) q[2];
sx q[2];
rz(2.3449281) q[2];
rz(1.1967777) q[3];
sx q[3];
rz(-1.6007042) q[3];
sx q[3];
rz(-2.5483907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6881123) q[0];
sx q[0];
rz(-2.1402833) q[0];
sx q[0];
rz(1.8185115) q[0];
rz(-2.6755013) q[1];
sx q[1];
rz(-2.2345462) q[1];
sx q[1];
rz(1.1493433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93724429) q[0];
sx q[0];
rz(-2.4597628) q[0];
sx q[0];
rz(2.3005062) q[0];
x q[1];
rz(-0.39300275) q[2];
sx q[2];
rz(-0.98969029) q[2];
sx q[2];
rz(-2.4266134) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2298849) q[1];
sx q[1];
rz(-2.9344892) q[1];
sx q[1];
rz(2.9207499) q[1];
x q[2];
rz(1.2364616) q[3];
sx q[3];
rz(-1.4282246) q[3];
sx q[3];
rz(2.4041374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65679067) q[2];
sx q[2];
rz(-0.66978407) q[2];
sx q[2];
rz(-0.82526866) q[2];
rz(1.2800062) q[3];
sx q[3];
rz(-1.620159) q[3];
sx q[3];
rz(0.72117225) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.698302) q[0];
sx q[0];
rz(-1.8809141) q[0];
sx q[0];
rz(-1.3662421) q[0];
rz(2.0514936) q[1];
sx q[1];
rz(-1.6366199) q[1];
sx q[1];
rz(-1.3390138) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96025673) q[0];
sx q[0];
rz(-0.68013817) q[0];
sx q[0];
rz(0.28016092) q[0];
x q[1];
rz(2.2750762) q[2];
sx q[2];
rz(-1.0703868) q[2];
sx q[2];
rz(1.4943708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8112534) q[1];
sx q[1];
rz(-2.898138) q[1];
sx q[1];
rz(-2.6324243) q[1];
rz(-pi) q[2];
rz(2.3424266) q[3];
sx q[3];
rz(-0.46119565) q[3];
sx q[3];
rz(-0.14652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.94685405) q[2];
sx q[2];
rz(-1.7297435) q[2];
sx q[2];
rz(-2.8080688) q[2];
rz(-0.58878318) q[3];
sx q[3];
rz(-0.47313658) q[3];
sx q[3];
rz(-2.3401882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3280585) q[0];
sx q[0];
rz(-1.3575587) q[0];
sx q[0];
rz(1.300977) q[0];
rz(-2.1133555) q[1];
sx q[1];
rz(-2.6515549) q[1];
sx q[1];
rz(-2.6992544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.418405) q[0];
sx q[0];
rz(-0.53779477) q[0];
sx q[0];
rz(-2.6502953) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47719854) q[2];
sx q[2];
rz(-2.0522567) q[2];
sx q[2];
rz(-1.19966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79696733) q[1];
sx q[1];
rz(-1.7360002) q[1];
sx q[1];
rz(-0.43032531) q[1];
x q[2];
rz(1.6746927) q[3];
sx q[3];
rz(-1.8891617) q[3];
sx q[3];
rz(1.1863134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77171317) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(-0.4717353) q[2];
rz(-1.6996023) q[3];
sx q[3];
rz(-2.580018) q[3];
sx q[3];
rz(2.877318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0278552) q[0];
sx q[0];
rz(-1.6738482) q[0];
sx q[0];
rz(-0.16673985) q[0];
rz(0.79404229) q[1];
sx q[1];
rz(-1.9414976) q[1];
sx q[1];
rz(0.099743191) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82631153) q[0];
sx q[0];
rz(-2.4874685) q[0];
sx q[0];
rz(1.0423866) q[0];
rz(0.95364596) q[2];
sx q[2];
rz(-1.2600449) q[2];
sx q[2];
rz(-1.0342395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7078898) q[1];
sx q[1];
rz(-1.6968445) q[1];
sx q[1];
rz(0.57278244) q[1];
rz(-pi) q[2];
rz(-0.41856237) q[3];
sx q[3];
rz(-2.1893614) q[3];
sx q[3];
rz(-1.5833861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2793067) q[2];
sx q[2];
rz(-1.0655094) q[2];
sx q[2];
rz(0.12969895) q[2];
rz(1.8779523) q[3];
sx q[3];
rz(-1.940515) q[3];
sx q[3];
rz(2.9960846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.170914) q[0];
sx q[0];
rz(-0.92249528) q[0];
sx q[0];
rz(0.42651919) q[0];
rz(1.0921987) q[1];
sx q[1];
rz(-2.009232) q[1];
sx q[1];
rz(1.0027286) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.730091) q[0];
sx q[0];
rz(-2.1637332) q[0];
sx q[0];
rz(2.8286295) q[0];
rz(-pi) q[1];
rz(-2.4187047) q[2];
sx q[2];
rz(-1.5010415) q[2];
sx q[2];
rz(2.7107014) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39053171) q[1];
sx q[1];
rz(-0.56274429) q[1];
sx q[1];
rz(-0.90969245) q[1];
rz(-pi) q[2];
rz(-2.0277983) q[3];
sx q[3];
rz(-0.49906956) q[3];
sx q[3];
rz(2.6036865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2670474) q[2];
sx q[2];
rz(-0.55570498) q[2];
sx q[2];
rz(-3.0181273) q[2];
rz(-2.5115783) q[3];
sx q[3];
rz(-1.29653) q[3];
sx q[3];
rz(2.4251078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89710871) q[0];
sx q[0];
rz(-2.4429595) q[0];
sx q[0];
rz(-2.3022292) q[0];
rz(0.15946236) q[1];
sx q[1];
rz(-0.99226743) q[1];
sx q[1];
rz(2.2119567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5073701) q[0];
sx q[0];
rz(-0.78222021) q[0];
sx q[0];
rz(-0.36391516) q[0];
x q[1];
rz(2.5505374) q[2];
sx q[2];
rz(-1.672571) q[2];
sx q[2];
rz(2.7391694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5413721) q[1];
sx q[1];
rz(-1.445679) q[1];
sx q[1];
rz(-0.88612963) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22945424) q[3];
sx q[3];
rz(-0.94564518) q[3];
sx q[3];
rz(1.8895686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.9440426) q[2];
sx q[2];
rz(-1.9746747) q[2];
sx q[2];
rz(0.76773947) q[2];
rz(-2.7770216) q[3];
sx q[3];
rz(-1.7578099) q[3];
sx q[3];
rz(3.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1345054) q[0];
sx q[0];
rz(-2.5275079) q[0];
sx q[0];
rz(-2.2579204) q[0];
rz(2.9332352) q[1];
sx q[1];
rz(-0.50269428) q[1];
sx q[1];
rz(-1.3349894) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3985838) q[0];
sx q[0];
rz(-1.5129876) q[0];
sx q[0];
rz(-0.065280838) q[0];
rz(0.52749525) q[2];
sx q[2];
rz(-1.732382) q[2];
sx q[2];
rz(1.8277663) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.523004) q[1];
sx q[1];
rz(-0.72806137) q[1];
sx q[1];
rz(-2.2656206) q[1];
rz(2.317071) q[3];
sx q[3];
rz(-2.093165) q[3];
sx q[3];
rz(-0.29545142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16365446) q[2];
sx q[2];
rz(-1.5074707) q[2];
sx q[2];
rz(3.0955691) q[2];
rz(-0.16511354) q[3];
sx q[3];
rz(-0.29296994) q[3];
sx q[3];
rz(-1.9273812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5788427) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(-0.1651172) q[1];
sx q[1];
rz(-1.4585635) q[1];
sx q[1];
rz(-0.30452902) q[1];
rz(-2.5223059) q[2];
sx q[2];
rz(-0.52792567) q[2];
sx q[2];
rz(-0.96326258) q[2];
rz(-1.0551999) q[3];
sx q[3];
rz(-1.1451117) q[3];
sx q[3];
rz(-2.3449863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

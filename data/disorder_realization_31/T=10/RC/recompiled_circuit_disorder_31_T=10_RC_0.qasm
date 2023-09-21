OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(-2.2990062) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(1.458118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320839) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(-0.30755933) q[0];
x q[1];
rz(3.1137142) q[2];
sx q[2];
rz(-0.61402245) q[2];
sx q[2];
rz(2.8369454) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9021437) q[1];
sx q[1];
rz(-1.4414756) q[1];
sx q[1];
rz(-2.7620478) q[1];
x q[2];
rz(-0.59431608) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(3.1241824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(2.9620985) q[2];
rz(1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(-2.3195482) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-0.32546145) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852859) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(0.73861648) q[0];
x q[1];
rz(2.1968368) q[2];
sx q[2];
rz(-1.2465887) q[2];
sx q[2];
rz(0.39694436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0305811) q[1];
sx q[1];
rz(-1.4861373) q[1];
sx q[1];
rz(-0.76534033) q[1];
rz(-2.9566544) q[3];
sx q[3];
rz(-1.3122845) q[3];
sx q[3];
rz(-2.1188494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.5154243) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-2.0498958) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9515581) q[0];
sx q[0];
rz(-1.0150195) q[0];
sx q[0];
rz(-0.65727289) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0902363) q[2];
sx q[2];
rz(-2.4161077) q[2];
sx q[2];
rz(-1.5067593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4746689) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(2.073642) q[1];
rz(-1.5063498) q[3];
sx q[3];
rz(-1.6985053) q[3];
sx q[3];
rz(-2.8116022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-2.2606405) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(-0.4367035) q[0];
rz(-2.9084335) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46582857) q[0];
sx q[0];
rz(-2.7390263) q[0];
sx q[0];
rz(-2.7990544) q[0];
x q[1];
rz(-0.15375806) q[2];
sx q[2];
rz(-0.69090828) q[2];
sx q[2];
rz(-1.5916057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51018184) q[1];
sx q[1];
rz(-1.5893755) q[1];
sx q[1];
rz(1.9287964) q[1];
rz(-1.6944828) q[3];
sx q[3];
rz(-1.6061155) q[3];
sx q[3];
rz(0.24231054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0115396) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524277) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4207626) q[0];
sx q[0];
rz(-1.2569068) q[0];
sx q[0];
rz(1.785196) q[0];
x q[1];
rz(-2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(-2.2774334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.00759) q[1];
sx q[1];
rz(-0.78467272) q[1];
sx q[1];
rz(2.6919634) q[1];
rz(1.3051885) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34981397) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.314635) q[0];
sx q[0];
rz(-1.1420113) q[0];
sx q[0];
rz(1.6786806) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1509402) q[2];
sx q[2];
rz(-0.93517762) q[2];
sx q[2];
rz(0.90204001) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4149474) q[1];
sx q[1];
rz(-0.84268314) q[1];
sx q[1];
rz(-0.77264087) q[1];
x q[2];
rz(0.10156472) q[3];
sx q[3];
rz(-1.933681) q[3];
sx q[3];
rz(0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(-2.2475524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(2.7303625) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(-0.033989865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2798529) q[0];
sx q[0];
rz(-1.732677) q[0];
sx q[0];
rz(0.78761657) q[0];
x q[1];
rz(-2.2529644) q[2];
sx q[2];
rz(-0.76105984) q[2];
sx q[2];
rz(-1.467848) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36044869) q[1];
sx q[1];
rz(-0.89372674) q[1];
sx q[1];
rz(-0.064868049) q[1];
rz(0.61492413) q[3];
sx q[3];
rz(-1.691754) q[3];
sx q[3];
rz(-1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(2.7476655) q[0];
rz(-2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.6961018) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5515585) q[0];
sx q[0];
rz(-0.43338767) q[0];
sx q[0];
rz(-1.5831468) q[0];
rz(0.94930737) q[2];
sx q[2];
rz(-2.0645803) q[2];
sx q[2];
rz(0.91044237) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8135012) q[1];
sx q[1];
rz(-0.5760759) q[1];
sx q[1];
rz(0.2920132) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66283488) q[3];
sx q[3];
rz(-1.0532866) q[3];
sx q[3];
rz(-0.55596065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(1.5173222) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.2517713) q[0];
rz(2.4720526) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(-2.8318185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944086) q[0];
sx q[0];
rz(-1.697288) q[0];
sx q[0];
rz(0.5267609) q[0];
x q[1];
rz(3.1321399) q[2];
sx q[2];
rz(-1.3938892) q[2];
sx q[2];
rz(2.0131907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0032776912) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(0.059271952) q[1];
rz(-0.59659776) q[3];
sx q[3];
rz(-1.784424) q[3];
sx q[3];
rz(-1.9006157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70242515) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(-1.9343728) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(-2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.6419798) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4831055) q[0];
sx q[0];
rz(-0.30880901) q[0];
sx q[0];
rz(1.2312066) q[0];
rz(-pi) q[1];
rz(-1.502938) q[2];
sx q[2];
rz(-2.3336764) q[2];
sx q[2];
rz(0.18005904) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3824532) q[1];
sx q[1];
rz(-1.1083974) q[1];
sx q[1];
rz(2.9049302) q[1];
rz(-pi) q[2];
rz(-1.5564735) q[3];
sx q[3];
rz(-0.17281547) q[3];
sx q[3];
rz(0.9848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(-0.45599109) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-3.070667) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(2.3161841) q[2];
sx q[2];
rz(-2.5054629) q[2];
sx q[2];
rz(2.7873743) q[2];
rz(1.0740888) q[3];
sx q[3];
rz(-2.2278193) q[3];
sx q[3];
rz(-0.57886119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

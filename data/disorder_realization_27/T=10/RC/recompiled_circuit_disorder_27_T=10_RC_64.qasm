OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(-1.385153) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18583488) q[0];
sx q[0];
rz(-1.251936) q[0];
sx q[0];
rz(-2.3715109) q[0];
rz(-pi) q[1];
rz(1.295747) q[2];
sx q[2];
rz(-1.0359456) q[2];
sx q[2];
rz(1.5799074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4298809) q[1];
sx q[1];
rz(-1.3379339) q[1];
sx q[1];
rz(1.4038605) q[1];
rz(-pi) q[2];
rz(0.6119421) q[3];
sx q[3];
rz(-0.7512593) q[3];
sx q[3];
rz(-0.9179759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2549071) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(2.936426) q[2];
rz(-0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.40760621) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(-2.6876887) q[0];
rz(-1.0247963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.227238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95943806) q[0];
sx q[0];
rz(-1.6998788) q[0];
sx q[0];
rz(3.0653619) q[0];
rz(-0.71364673) q[2];
sx q[2];
rz(-0.64672856) q[2];
sx q[2];
rz(1.5681859) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8933567) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(-0.30171079) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9407528) q[3];
sx q[3];
rz(-1.4174995) q[3];
sx q[3];
rz(2.4917345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.041302117) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(0.96810961) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48297468) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(2.2429402) q[0];
rz(0.99575106) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(2.8083037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5211398) q[0];
sx q[0];
rz(-0.95690173) q[0];
sx q[0];
rz(0.99434538) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3938445) q[2];
sx q[2];
rz(-1.4743917) q[2];
sx q[2];
rz(2.5071438) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7283199) q[1];
sx q[1];
rz(-1.4188758) q[1];
sx q[1];
rz(-2.1876213) q[1];
rz(-pi) q[2];
x q[2];
rz(0.962978) q[3];
sx q[3];
rz(-0.64219785) q[3];
sx q[3];
rz(-1.1841139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68625346) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(2.3006556) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44160098) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(-0.68471318) q[0];
rz(-2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(1.205014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3048153) q[0];
sx q[0];
rz(-2.4495821) q[0];
sx q[0];
rz(1.5185604) q[0];
x q[1];
rz(-1.8393458) q[2];
sx q[2];
rz(-2.4797202) q[2];
sx q[2];
rz(0.72999398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1531096) q[1];
sx q[1];
rz(-1.298269) q[1];
sx q[1];
rz(2.2793798) q[1];
rz(-pi) q[2];
rz(-2.2756696) q[3];
sx q[3];
rz(-1.62543) q[3];
sx q[3];
rz(1.5254525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(0.37115804) q[2];
rz(1.4012198) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(1.1192809) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.086833) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(0.99331028) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(-0.55508074) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9585882) q[0];
sx q[0];
rz(-2.8259813) q[0];
sx q[0];
rz(1.5297699) q[0];
x q[1];
rz(3.073408) q[2];
sx q[2];
rz(-1.9848595) q[2];
sx q[2];
rz(-2.806864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86075393) q[1];
sx q[1];
rz(-2.6959531) q[1];
sx q[1];
rz(-1.6153567) q[1];
x q[2];
rz(2.5693232) q[3];
sx q[3];
rz(-0.096147691) q[3];
sx q[3];
rz(0.26086807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.30620265) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(-3.0026657) q[2];
rz(-0.94240087) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(0.55148235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979364) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(-2.561835) q[0];
rz(-3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.5396083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67028763) q[0];
sx q[0];
rz(-1.3677214) q[0];
sx q[0];
rz(-0.85104403) q[0];
rz(-3.0110353) q[2];
sx q[2];
rz(-1.5255543) q[2];
sx q[2];
rz(3.0322078) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3967267) q[1];
sx q[1];
rz(-0.52280451) q[1];
sx q[1];
rz(-2.1641939) q[1];
rz(1.9220011) q[3];
sx q[3];
rz(-2.425819) q[3];
sx q[3];
rz(0.90482611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5876028) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(0.26947752) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6948029) q[0];
sx q[0];
rz(-0.61820522) q[0];
sx q[0];
rz(-2.457298) q[0];
rz(-3.0220095) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(0.51876846) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8914688) q[0];
sx q[0];
rz(-2.2995673) q[0];
sx q[0];
rz(2.9617589) q[0];
rz(1.7905551) q[2];
sx q[2];
rz(-2.3292654) q[2];
sx q[2];
rz(2.1213558) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3541975) q[1];
sx q[1];
rz(-2.0320315) q[1];
sx q[1];
rz(2.4128777) q[1];
rz(-pi) q[2];
rz(1.1701059) q[3];
sx q[3];
rz(-1.5466994) q[3];
sx q[3];
rz(-1.3093684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8873022) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(0.56345144) q[2];
rz(0.051579483) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(-2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96034399) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(1.3990336) q[0];
rz(-2.3545806) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(-0.74434892) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8101013) q[0];
sx q[0];
rz(-1.5163251) q[0];
sx q[0];
rz(-2.9936552) q[0];
rz(-pi) q[1];
rz(2.0853945) q[2];
sx q[2];
rz(-1.4375763) q[2];
sx q[2];
rz(1.0366057) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0155322) q[1];
sx q[1];
rz(-1.0711526) q[1];
sx q[1];
rz(1.3152895) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.048150496) q[3];
sx q[3];
rz(-0.45504323) q[3];
sx q[3];
rz(-2.9330848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4259592) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(2.5320833) q[2];
rz(0.65731796) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(-1.2414744) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(-0.77493587) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0477714) q[0];
sx q[0];
rz(-1.9592966) q[0];
sx q[0];
rz(-2.2872778) q[0];
rz(2.0140892) q[2];
sx q[2];
rz(-1.7121676) q[2];
sx q[2];
rz(-3.1183426) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7271125) q[1];
sx q[1];
rz(-0.4948805) q[1];
sx q[1];
rz(0.40463573) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14258607) q[3];
sx q[3];
rz(-1.7877868) q[3];
sx q[3];
rz(-2.2246974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.187414) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(-0.15787086) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071844) q[0];
sx q[0];
rz(-0.97244111) q[0];
sx q[0];
rz(0.21433314) q[0];
rz(-2.4841323) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(1.0459895) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25046529) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(1.6842708) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5506053) q[2];
sx q[2];
rz(-1.5948442) q[2];
sx q[2];
rz(-1.9901333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18098772) q[1];
sx q[1];
rz(-2.1121896) q[1];
sx q[1];
rz(0.42591806) q[1];
rz(-pi) q[2];
rz(-2.7361715) q[3];
sx q[3];
rz(-1.5463721) q[3];
sx q[3];
rz(3.0693808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7252698) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(-1.0894758) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(-2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8626704) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(-1.6090341) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(-2.5074742) q[2];
sx q[2];
rz(-0.49370439) q[2];
sx q[2];
rz(1.6903071) q[2];
rz(2.5915626) q[3];
sx q[3];
rz(-1.6246272) q[3];
sx q[3];
rz(0.18679242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
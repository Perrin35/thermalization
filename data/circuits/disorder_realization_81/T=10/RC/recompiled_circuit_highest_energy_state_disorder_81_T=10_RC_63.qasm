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
rz(0.089224815) q[0];
sx q[0];
rz(2.3951946) q[0];
sx q[0];
rz(8.7328773) q[0];
rz(3.6045868) q[1];
sx q[1];
rz(4.7582518) q[1];
sx q[1];
rz(7.313348) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97447831) q[0];
sx q[0];
rz(-1.8016707) q[0];
sx q[0];
rz(2.5690325) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18333034) q[2];
sx q[2];
rz(-1.1784484) q[2];
sx q[2];
rz(0.13428706) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15265118) q[1];
sx q[1];
rz(-1.8957152) q[1];
sx q[1];
rz(-0.12293651) q[1];
x q[2];
rz(0.68891565) q[3];
sx q[3];
rz(-1.6647881) q[3];
sx q[3];
rz(1.9436364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2618711) q[2];
sx q[2];
rz(-0.77360952) q[2];
sx q[2];
rz(0.36641463) q[2];
rz(-1.0877747) q[3];
sx q[3];
rz(-2.352114) q[3];
sx q[3];
rz(2.4521949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3159897) q[0];
sx q[0];
rz(-3.0593384) q[0];
sx q[0];
rz(0.91973037) q[0];
rz(0.74291825) q[1];
sx q[1];
rz(-1.2580322) q[1];
sx q[1];
rz(1.119335) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9237083) q[0];
sx q[0];
rz(-2.9842434) q[0];
sx q[0];
rz(-1.1697392) q[0];
x q[1];
rz(0.7190312) q[2];
sx q[2];
rz(-1.708548) q[2];
sx q[2];
rz(0.097830398) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97174469) q[1];
sx q[1];
rz(-1.3339284) q[1];
sx q[1];
rz(-2.4011517) q[1];
rz(-2.1503937) q[3];
sx q[3];
rz(-2.0127735) q[3];
sx q[3];
rz(0.22721618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6433581) q[2];
sx q[2];
rz(-1.1946119) q[2];
sx q[2];
rz(-0.052113459) q[2];
rz(1.107996) q[3];
sx q[3];
rz(-0.85927695) q[3];
sx q[3];
rz(0.064854709) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4614918) q[0];
sx q[0];
rz(-1.9100186) q[0];
sx q[0];
rz(-1.8517866) q[0];
rz(-2.3884804) q[1];
sx q[1];
rz(-1.6878004) q[1];
sx q[1];
rz(2.6444816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4418934) q[0];
sx q[0];
rz(-1.9632247) q[0];
sx q[0];
rz(0.59458574) q[0];
rz(-0.75211033) q[2];
sx q[2];
rz(-0.83319887) q[2];
sx q[2];
rz(-1.9832265) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0785574) q[1];
sx q[1];
rz(-0.37606701) q[1];
sx q[1];
rz(-3.0381812) q[1];
rz(-pi) q[2];
rz(3.0158832) q[3];
sx q[3];
rz(-2.4851228) q[3];
sx q[3];
rz(-1.2510811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46811732) q[2];
sx q[2];
rz(-2.868729) q[2];
sx q[2];
rz(-2.4172778) q[2];
rz(-2.916548) q[3];
sx q[3];
rz(-1.8658172) q[3];
sx q[3];
rz(2.3465274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1555136) q[0];
sx q[0];
rz(-0.88548311) q[0];
sx q[0];
rz(-0.28050637) q[0];
rz(-1.6579423) q[1];
sx q[1];
rz(-1.9397441) q[1];
sx q[1];
rz(-2.6703506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4330632) q[0];
sx q[0];
rz(-1.9380599) q[0];
sx q[0];
rz(-3.1172453) q[0];
x q[1];
rz(2.0986841) q[2];
sx q[2];
rz(-3.0110741) q[2];
sx q[2];
rz(-0.077613398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4268643) q[1];
sx q[1];
rz(-2.8694723) q[1];
sx q[1];
rz(-1.4674027) q[1];
x q[2];
rz(-2.1625278) q[3];
sx q[3];
rz(-2.2869898) q[3];
sx q[3];
rz(-1.5997052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9435141) q[2];
sx q[2];
rz(-2.5835865) q[2];
sx q[2];
rz(1.0151218) q[2];
rz(-2.4070168) q[3];
sx q[3];
rz(-1.3646804) q[3];
sx q[3];
rz(-2.6034897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.9839142) q[0];
sx q[0];
rz(-1.9240802) q[0];
sx q[0];
rz(0.96018803) q[0];
rz(-3.0958815) q[1];
sx q[1];
rz(-2.263133) q[1];
sx q[1];
rz(-3.1083621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7901781) q[0];
sx q[0];
rz(-1.1006257) q[0];
sx q[0];
rz(2.1432502) q[0];
x q[1];
rz(1.8383547) q[2];
sx q[2];
rz(-1.55624) q[2];
sx q[2];
rz(0.47786682) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7864125) q[1];
sx q[1];
rz(-2.7828453) q[1];
sx q[1];
rz(-1.9500109) q[1];
x q[2];
rz(-2.4016693) q[3];
sx q[3];
rz(-2.5725274) q[3];
sx q[3];
rz(-0.25404938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9783322) q[2];
sx q[2];
rz(-2.1425118) q[2];
sx q[2];
rz(-2.8283289) q[2];
rz(2.7836109) q[3];
sx q[3];
rz(-1.0397747) q[3];
sx q[3];
rz(3.0430005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.1545496) q[0];
sx q[0];
rz(-1.6313666) q[0];
sx q[0];
rz(0.48598591) q[0];
rz(-1.9588574) q[1];
sx q[1];
rz(-2.4458838) q[1];
sx q[1];
rz(2.8114496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8517338) q[0];
sx q[0];
rz(-2.8427474) q[0];
sx q[0];
rz(1.4695704) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2171872) q[2];
sx q[2];
rz(-0.87714783) q[2];
sx q[2];
rz(0.93111699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.15974074) q[1];
sx q[1];
rz(-1.0969437) q[1];
sx q[1];
rz(0.57094806) q[1];
rz(-pi) q[2];
rz(1.347769) q[3];
sx q[3];
rz(-1.2320326) q[3];
sx q[3];
rz(1.2028002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5064454) q[2];
sx q[2];
rz(-0.61656419) q[2];
sx q[2];
rz(-0.35663024) q[2];
rz(2.5470274) q[3];
sx q[3];
rz(-1.4831355) q[3];
sx q[3];
rz(-0.73896343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67721382) q[0];
sx q[0];
rz(-0.76736275) q[0];
sx q[0];
rz(-2.5546524) q[0];
rz(2.7691973) q[1];
sx q[1];
rz(-1.3601235) q[1];
sx q[1];
rz(0.24758235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0816168) q[0];
sx q[0];
rz(-0.46591407) q[0];
sx q[0];
rz(1.640725) q[0];
x q[1];
rz(-1.3273986) q[2];
sx q[2];
rz(-2.1888424) q[2];
sx q[2];
rz(2.7400561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7786488) q[1];
sx q[1];
rz(-1.9041859) q[1];
sx q[1];
rz(-2.2466546) q[1];
x q[2];
rz(-2.7917737) q[3];
sx q[3];
rz(-1.3365776) q[3];
sx q[3];
rz(-1.2522335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8783012) q[2];
sx q[2];
rz(-1.9204488) q[2];
sx q[2];
rz(-0.2090052) q[2];
rz(2.5448997) q[3];
sx q[3];
rz(-2.7946819) q[3];
sx q[3];
rz(-1.5359115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7335032) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(-0.39304131) q[0];
rz(-2.4709002) q[1];
sx q[1];
rz(-1.1436983) q[1];
sx q[1];
rz(1.6938946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5590305) q[0];
sx q[0];
rz(-2.6295289) q[0];
sx q[0];
rz(1.3150231) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1850519) q[2];
sx q[2];
rz(-1.8280085) q[2];
sx q[2];
rz(-1.1957912) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.61673855) q[1];
sx q[1];
rz(-2.3890807) q[1];
sx q[1];
rz(0.77857154) q[1];
x q[2];
rz(2.3438934) q[3];
sx q[3];
rz(-1.2273066) q[3];
sx q[3];
rz(-2.2960312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4286917) q[2];
sx q[2];
rz(-2.28076) q[2];
sx q[2];
rz(0.11670308) q[2];
rz(1.0137001) q[3];
sx q[3];
rz(-1.2319535) q[3];
sx q[3];
rz(1.8359756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75334615) q[0];
sx q[0];
rz(-1.4660864) q[0];
sx q[0];
rz(1.7062794) q[0];
rz(-2.1460136) q[1];
sx q[1];
rz(-1.8228056) q[1];
sx q[1];
rz(0.54042712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30975809) q[0];
sx q[0];
rz(-2.0481526) q[0];
sx q[0];
rz(-1.8670797) q[0];
rz(2.5535271) q[2];
sx q[2];
rz(-2.1508475) q[2];
sx q[2];
rz(1.560245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1241972) q[1];
sx q[1];
rz(-0.81995839) q[1];
sx q[1];
rz(0.52738054) q[1];
rz(-0.53934259) q[3];
sx q[3];
rz(-2.4423222) q[3];
sx q[3];
rz(2.5756096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.871375) q[2];
sx q[2];
rz(-1.8367218) q[2];
sx q[2];
rz(0.77667856) q[2];
rz(-3.1089879) q[3];
sx q[3];
rz(-1.5332581) q[3];
sx q[3];
rz(-1.5743272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33083522) q[0];
sx q[0];
rz(-0.36892712) q[0];
sx q[0];
rz(2.708013) q[0];
rz(-2.4724204) q[1];
sx q[1];
rz(-1.9586261) q[1];
sx q[1];
rz(-0.42632595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2992914) q[0];
sx q[0];
rz(-0.29336818) q[0];
sx q[0];
rz(0.82404739) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74684116) q[2];
sx q[2];
rz(-1.3051118) q[2];
sx q[2];
rz(2.8248252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43526442) q[1];
sx q[1];
rz(-2.5245275) q[1];
sx q[1];
rz(0.63539482) q[1];
x q[2];
rz(2.4678585) q[3];
sx q[3];
rz(-2.4387283) q[3];
sx q[3];
rz(1.7029312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7508037) q[2];
sx q[2];
rz(-0.28826737) q[2];
sx q[2];
rz(-1.1061579) q[2];
rz(3.028051) q[3];
sx q[3];
rz(-2.2083486) q[3];
sx q[3];
rz(0.043702628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15585598) q[0];
sx q[0];
rz(-2.2423797) q[0];
sx q[0];
rz(2.5806497) q[0];
rz(-0.38324311) q[1];
sx q[1];
rz(-2.0189197) q[1];
sx q[1];
rz(-0.22519208) q[1];
rz(-1.2651934) q[2];
sx q[2];
rz(-2.1060747) q[2];
sx q[2];
rz(1.0155564) q[2];
rz(-1.0998955) q[3];
sx q[3];
rz(-2.3318137) q[3];
sx q[3];
rz(-1.6398738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

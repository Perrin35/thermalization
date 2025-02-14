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
rz(2.2284722) q[0];
sx q[0];
rz(-2.1051536) q[0];
sx q[0];
rz(1.5358465) q[0];
rz(3.8029382) q[1];
sx q[1];
rz(4.3284419) q[1];
sx q[1];
rz(8.9882589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3876037) q[0];
sx q[0];
rz(-1.4719098) q[0];
sx q[0];
rz(-1.8900699) q[0];
rz(2.5769325) q[2];
sx q[2];
rz(-1.1810978) q[2];
sx q[2];
rz(1.0124568) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7750262) q[1];
sx q[1];
rz(-2.4373551) q[1];
sx q[1];
rz(1.5077059) q[1];
rz(0.8869041) q[3];
sx q[3];
rz(-1.7085848) q[3];
sx q[3];
rz(-0.52007127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6785757) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(0.92301816) q[2];
rz(3.0758514) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(-1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(-0.03751066) q[0];
rz(0.58049479) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(2.4494749) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1624533) q[0];
sx q[0];
rz(-1.395749) q[0];
sx q[0];
rz(-1.2412423) q[0];
rz(0.28022556) q[2];
sx q[2];
rz(-1.0839274) q[2];
sx q[2];
rz(3.1161199) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2381535) q[1];
sx q[1];
rz(-0.80793706) q[1];
sx q[1];
rz(-2.417516) q[1];
rz(-pi) q[2];
rz(-1.7440657) q[3];
sx q[3];
rz(-0.82635802) q[3];
sx q[3];
rz(-0.92191523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47824255) q[2];
sx q[2];
rz(-1.1967836) q[2];
sx q[2];
rz(-1.7459858) q[2];
rz(0.1869525) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-0.54005867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.4726987) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(0.5589267) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5907954) q[1];
sx q[1];
rz(2.8899946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4959481) q[0];
sx q[0];
rz(-2.6584266) q[0];
sx q[0];
rz(-2.7208485) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0195658) q[2];
sx q[2];
rz(-1.2740268) q[2];
sx q[2];
rz(-0.084159764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14519599) q[1];
sx q[1];
rz(-0.72219488) q[1];
sx q[1];
rz(-1.9737555) q[1];
rz(-pi) q[2];
rz(0.91355998) q[3];
sx q[3];
rz(-1.0042448) q[3];
sx q[3];
rz(2.4853137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.805213) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(-1.7097998) q[2];
rz(2.924887) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8828204) q[0];
sx q[0];
rz(-0.50676218) q[0];
sx q[0];
rz(-1.850542) q[0];
rz(-0.82921118) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(-2.1270027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050796789) q[0];
sx q[0];
rz(-1.3279339) q[0];
sx q[0];
rz(-2.5953351) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5780771) q[2];
sx q[2];
rz(-1.2402099) q[2];
sx q[2];
rz(-0.78621582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9660849) q[1];
sx q[1];
rz(-2.1455742) q[1];
sx q[1];
rz(-0.0055343363) q[1];
x q[2];
rz(-1.4271554) q[3];
sx q[3];
rz(-0.82089409) q[3];
sx q[3];
rz(0.088975541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0097051) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(2.3809872) q[2];
rz(0.46918121) q[3];
sx q[3];
rz(-1.4696591) q[3];
sx q[3];
rz(-2.40707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483265) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(0.58116466) q[0];
rz(-1.5900853) q[1];
sx q[1];
rz(-2.4947512) q[1];
sx q[1];
rz(-0.69492984) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0804176) q[0];
sx q[0];
rz(-1.8616849) q[0];
sx q[0];
rz(-2.8982452) q[0];
x q[1];
rz(0.29232358) q[2];
sx q[2];
rz(-1.6609123) q[2];
sx q[2];
rz(-2.9499049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.018539704) q[1];
sx q[1];
rz(-2.596302) q[1];
sx q[1];
rz(-3.0569949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65443734) q[3];
sx q[3];
rz(-2.8469466) q[3];
sx q[3];
rz(-1.0279931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(-2.8387873) q[2];
rz(1.7152202) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(-2.5572131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(-0.67702174) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(-0.69497481) q[0];
rz(-2.8547844) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(1.8668176) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5329842) q[0];
sx q[0];
rz(-1.2600945) q[0];
sx q[0];
rz(2.0771189) q[0];
x q[1];
rz(0.54194258) q[2];
sx q[2];
rz(-1.9795609) q[2];
sx q[2];
rz(2.8084076) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66836897) q[1];
sx q[1];
rz(-0.76938564) q[1];
sx q[1];
rz(-2.0955057) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8096089) q[3];
sx q[3];
rz(-1.0019571) q[3];
sx q[3];
rz(-2.5425926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0810762) q[2];
sx q[2];
rz(-2.923107) q[2];
sx q[2];
rz(-2.1742353) q[2];
rz(-0.71990144) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(-1.8657743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8571781) q[0];
sx q[0];
rz(-0.030817742) q[0];
sx q[0];
rz(1.6946633) q[0];
rz(2.6743496) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(1.9460868) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1074813) q[0];
sx q[0];
rz(-1.1426289) q[0];
sx q[0];
rz(-1.0635682) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81584064) q[2];
sx q[2];
rz(-1.3357478) q[2];
sx q[2];
rz(-1.5790958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0250281) q[1];
sx q[1];
rz(-2.0691074) q[1];
sx q[1];
rz(-3.0710941) q[1];
x q[2];
rz(0.012717207) q[3];
sx q[3];
rz(-1.5348292) q[3];
sx q[3];
rz(-1.4096703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82470992) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(-0.52428025) q[2];
rz(2.8892062) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1687932) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(1.9148781) q[0];
rz(0.75617689) q[1];
sx q[1];
rz(-0.94804472) q[1];
sx q[1];
rz(-0.39047584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062138012) q[0];
sx q[0];
rz(-1.3083757) q[0];
sx q[0];
rz(-0.45824893) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6092714) q[2];
sx q[2];
rz(-2.0145825) q[2];
sx q[2];
rz(2.2930068) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2578963) q[1];
sx q[1];
rz(-1.9054277) q[1];
sx q[1];
rz(2.9154214) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49440439) q[3];
sx q[3];
rz(-2.1526511) q[3];
sx q[3];
rz(2.0928008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8172424) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(1.3520757) q[2];
rz(-2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(-1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7290203) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(-0.66226688) q[0];
rz(-0.63016713) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(2.9060649) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2199687) q[0];
sx q[0];
rz(-1.3597836) q[0];
sx q[0];
rz(-1.3930265) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56925242) q[2];
sx q[2];
rz(-1.4919089) q[2];
sx q[2];
rz(-1.1752664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.16727) q[1];
sx q[1];
rz(-0.065792699) q[1];
sx q[1];
rz(-1.2029103) q[1];
x q[2];
rz(-0.8742378) q[3];
sx q[3];
rz(-1.6739419) q[3];
sx q[3];
rz(2.646605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0887289) q[2];
sx q[2];
rz(-2.1910618) q[2];
sx q[2];
rz(-0.55727422) q[2];
rz(0.82734621) q[3];
sx q[3];
rz(-1.6297623) q[3];
sx q[3];
rz(1.6415589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92852229) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(0.56761566) q[0];
rz(2.7396743) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(-1.7793122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80776725) q[0];
sx q[0];
rz(-2.5596715) q[0];
sx q[0];
rz(-2.3533217) q[0];
rz(-pi) q[1];
rz(-0.95133199) q[2];
sx q[2];
rz(-1.3580322) q[2];
sx q[2];
rz(-1.484953) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6475923) q[1];
sx q[1];
rz(-2.5060182) q[1];
sx q[1];
rz(-1.4383491) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3506593) q[3];
sx q[3];
rz(-0.91703992) q[3];
sx q[3];
rz(-2.4733651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42084971) q[2];
sx q[2];
rz(-1.9385312) q[2];
sx q[2];
rz(2.8813349) q[2];
rz(-0.69089729) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33547587) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(-2.7083022) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(-0.92930195) q[2];
sx q[2];
rz(-0.96972307) q[2];
sx q[2];
rz(-2.0676421) q[2];
rz(1.7078441) q[3];
sx q[3];
rz(-2.3564586) q[3];
sx q[3];
rz(-2.3740507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

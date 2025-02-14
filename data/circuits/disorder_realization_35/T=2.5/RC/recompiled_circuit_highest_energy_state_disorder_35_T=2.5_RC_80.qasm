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
rz(-0.91312042) q[0];
sx q[0];
rz(-1.0364391) q[0];
sx q[0];
rz(1.6057462) q[0];
rz(0.66134557) q[1];
sx q[1];
rz(-1.1868492) q[1];
sx q[1];
rz(-2.7050736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8494328) q[0];
sx q[0];
rz(-1.2531373) q[0];
sx q[0];
rz(0.10411291) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0233193) q[2];
sx q[2];
rz(-1.0529537) q[2];
sx q[2];
rz(0.3221727) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4492717) q[1];
sx q[1];
rz(-0.8682478) q[1];
sx q[1];
rz(-3.0880804) q[1];
rz(-2.2546886) q[3];
sx q[3];
rz(-1.4330079) q[3];
sx q[3];
rz(0.52007127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46301699) q[2];
sx q[2];
rz(-0.43872321) q[2];
sx q[2];
rz(-0.92301816) q[2];
rz(-3.0758514) q[3];
sx q[3];
rz(-1.8140503) q[3];
sx q[3];
rz(-1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-1.0634402) q[0];
sx q[0];
rz(0.03751066) q[0];
rz(-2.5610979) q[1];
sx q[1];
rz(-2.4721804) q[1];
sx q[1];
rz(-2.4494749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1624533) q[0];
sx q[0];
rz(-1.395749) q[0];
sx q[0];
rz(-1.2412423) q[0];
rz(-2.0522444) q[2];
sx q[2];
rz(-0.55608597) q[2];
sx q[2];
rz(2.5646445) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.925732) q[1];
sx q[1];
rz(-2.070148) q[1];
sx q[1];
rz(-2.4768922) q[1];
rz(0.18499891) q[3];
sx q[3];
rz(-2.3810412) q[3];
sx q[3];
rz(0.669125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47824255) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.7459858) q[2];
rz(2.9546402) q[3];
sx q[3];
rz(-1.170265) q[3];
sx q[3];
rz(2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66889399) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(0.5589267) q[0];
rz(-0.82866296) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-0.25159803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54788816) q[0];
sx q[0];
rz(-1.761709) q[0];
sx q[0];
rz(0.44661354) q[0];
x q[1];
rz(-2.7969486) q[2];
sx q[2];
rz(-2.0953669) q[2];
sx q[2];
rz(1.3087496) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1160186) q[1];
sx q[1];
rz(-1.3085828) q[1];
sx q[1];
rz(0.88974726) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2280327) q[3];
sx q[3];
rz(-1.0042448) q[3];
sx q[3];
rz(-2.4853137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.805213) q[2];
sx q[2];
rz(-2.7153335) q[2];
sx q[2];
rz(-1.4317929) q[2];
rz(-2.924887) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(-1.7823035) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8828204) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(1.850542) q[0];
rz(0.82921118) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(-1.0145899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9982581) q[0];
sx q[0];
rz(-0.59276544) q[0];
sx q[0];
rz(0.44500764) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9564187) q[2];
sx q[2];
rz(-2.1004538) q[2];
sx q[2];
rz(-0.5822863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9559052) q[1];
sx q[1];
rz(-0.57480156) q[1];
sx q[1];
rz(-1.5793403) q[1];
rz(-pi) q[2];
rz(0.75506702) q[3];
sx q[3];
rz(-1.4658548) q[3];
sx q[3];
rz(1.5800832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0097051) q[2];
sx q[2];
rz(-1.2010801) q[2];
sx q[2];
rz(-0.76060549) q[2];
rz(2.6724114) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0932662) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(-0.58116466) q[0];
rz(1.5900853) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(-0.69492984) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.061175) q[0];
sx q[0];
rz(-1.8616849) q[0];
sx q[0];
rz(2.8982452) q[0];
rz(0.29232358) q[2];
sx q[2];
rz(-1.4806804) q[2];
sx q[2];
rz(2.9499049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1174005) q[1];
sx q[1];
rz(-2.1139189) q[1];
sx q[1];
rz(1.5195816) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3881237) q[3];
sx q[3];
rz(-1.338306) q[3];
sx q[3];
rz(0.35216613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27802262) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(0.30280534) q[2];
rz(1.4263724) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(2.5572131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.4645709) q[0];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5329842) q[0];
sx q[0];
rz(-1.2600945) q[0];
sx q[0];
rz(2.0771189) q[0];
rz(-pi) q[1];
rz(2.0389245) q[2];
sx q[2];
rz(-1.0777567) q[2];
sx q[2];
rz(-1.6690892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66836897) q[1];
sx q[1];
rz(-0.76938564) q[1];
sx q[1];
rz(-1.046087) q[1];
rz(2.5596095) q[3];
sx q[3];
rz(-1.3701539) q[3];
sx q[3];
rz(2.0394005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(-2.1742353) q[2];
rz(2.4216912) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(-1.8657743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8571781) q[0];
sx q[0];
rz(-0.030817742) q[0];
sx q[0];
rz(-1.6946633) q[0];
rz(-2.6743496) q[1];
sx q[1];
rz(-1.2930361) q[1];
sx q[1];
rz(-1.1955059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.104983) q[0];
sx q[0];
rz(-0.65147841) q[0];
sx q[0];
rz(-2.3250513) q[0];
rz(-pi) q[1];
rz(2.325752) q[2];
sx q[2];
rz(-1.3357478) q[2];
sx q[2];
rz(1.5624969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1165645) q[1];
sx q[1];
rz(-1.0724853) q[1];
sx q[1];
rz(-3.0710941) q[1];
rz(-3.1288754) q[3];
sx q[3];
rz(-1.5348292) q[3];
sx q[3];
rz(1.7319224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3168827) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(-2.6173124) q[2];
rz(2.8892062) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(-3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1687932) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(-1.2267145) q[0];
rz(-2.3854158) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(0.39047584) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1487685) q[0];
sx q[0];
rz(-2.618194) q[0];
sx q[0];
rz(-0.54570178) q[0];
rz(-pi) q[1];
rz(-0.4440733) q[2];
sx q[2];
rz(-1.5360498) q[2];
sx q[2];
rz(2.4028558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23747102) q[1];
sx q[1];
rz(-1.3573705) q[1];
sx q[1];
rz(1.9135114) q[1];
rz(2.195695) q[3];
sx q[3];
rz(-2.3970553) q[3];
sx q[3];
rz(-2.8685399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8172424) q[2];
sx q[2];
rz(-2.4411026) q[2];
sx q[2];
rz(-1.7895169) q[2];
rz(-0.17524854) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(-1.9372743) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4125724) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(-0.66226688) q[0];
rz(2.5114255) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(-0.23552775) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5109638) q[0];
sx q[0];
rz(-2.8665344) q[0];
sx q[0];
rz(-2.4514626) q[0];
rz(-pi) q[1];
rz(2.9959683) q[2];
sx q[2];
rz(-2.5674985) q[2];
sx q[2];
rz(0.51806322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79865658) q[1];
sx q[1];
rz(-1.6321811) q[1];
sx q[1];
rz(-0.023691698) q[1];
x q[2];
rz(-2.2673549) q[3];
sx q[3];
rz(-1.4676508) q[3];
sx q[3];
rz(-0.49498765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0887289) q[2];
sx q[2];
rz(-2.1910618) q[2];
sx q[2];
rz(-0.55727422) q[2];
rz(-2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(1.5000337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92852229) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(-0.56761566) q[0];
rz(0.40191832) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(1.7793122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0716996) q[0];
sx q[0];
rz(-1.1728468) q[0];
sx q[0];
rz(-1.1342628) q[0];
x q[1];
rz(-2.882233) q[2];
sx q[2];
rz(-2.1742714) q[2];
sx q[2];
rz(3.0779787) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1115845) q[1];
sx q[1];
rz(-1.6492731) q[1];
sx q[1];
rz(2.2021738) q[1];
rz(-pi) q[2];
rz(-0.79093334) q[3];
sx q[3];
rz(-2.2245527) q[3];
sx q[3];
rz(-0.66822754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42084971) q[2];
sx q[2];
rz(-1.9385312) q[2];
sx q[2];
rz(-0.26025772) q[2];
rz(-2.4506954) q[3];
sx q[3];
rz(-2.2862209) q[3];
sx q[3];
rz(1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33547587) q[0];
sx q[0];
rz(-1.5806883) q[0];
sx q[0];
rz(-1.6089532) q[0];
rz(-0.43329049) q[1];
sx q[1];
rz(-1.8186124) q[1];
sx q[1];
rz(-1.1846452) q[1];
rz(-2.4240906) q[2];
sx q[2];
rz(-2.2926471) q[2];
sx q[2];
rz(-2.9903464) q[2];
rz(-2.3512202) q[3];
sx q[3];
rz(-1.4740667) q[3];
sx q[3];
rz(2.4355751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

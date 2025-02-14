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
rz(4.1780317) q[0];
sx q[0];
rz(10.960624) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2921599) q[0];
sx q[0];
rz(-1.8884553) q[0];
sx q[0];
rz(3.0374797) q[0];
rz(2.0233193) q[2];
sx q[2];
rz(-1.0529537) q[2];
sx q[2];
rz(2.81942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4492717) q[1];
sx q[1];
rz(-2.2733449) q[1];
sx q[1];
rz(0.0535123) q[1];
rz(-pi) q[2];
rz(-2.9645676) q[3];
sx q[3];
rz(-0.89460556) q[3];
sx q[3];
rz(2.2023622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46301699) q[2];
sx q[2];
rz(-0.43872321) q[2];
sx q[2];
rz(0.92301816) q[2];
rz(-3.0758514) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(-1.340723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(3.104082) q[0];
rz(2.5610979) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(-2.4494749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0789675) q[0];
sx q[0];
rz(-0.37165549) q[0];
sx q[0];
rz(-2.0709447) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0672928) q[2];
sx q[2];
rz(-1.8177336) q[2];
sx q[2];
rz(1.7301138) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1375274) q[1];
sx q[1];
rz(-2.1430227) q[1];
sx q[1];
rz(0.96478421) q[1];
rz(-1.397527) q[3];
sx q[3];
rz(-0.82635802) q[3];
sx q[3];
rz(-2.2196774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.47824255) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.3956068) q[2];
rz(-0.1869525) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4726987) q[0];
sx q[0];
rz(-0.63303328) q[0];
sx q[0];
rz(-2.582666) q[0];
rz(-2.3129297) q[1];
sx q[1];
rz(-1.5907954) q[1];
sx q[1];
rz(-0.25159803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64564451) q[0];
sx q[0];
rz(-0.48316607) q[0];
sx q[0];
rz(2.7208485) q[0];
rz(-pi) q[1];
rz(2.0992804) q[2];
sx q[2];
rz(-2.5229075) q[2];
sx q[2];
rz(1.2109735) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9963967) q[1];
sx q[1];
rz(-0.72219488) q[1];
sx q[1];
rz(1.1678371) q[1];
x q[2];
rz(2.4647242) q[3];
sx q[3];
rz(-2.1123611) q[3];
sx q[3];
rz(1.3071909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.805213) q[2];
sx q[2];
rz(-2.7153335) q[2];
sx q[2];
rz(-1.4317929) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25877229) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(1.850542) q[0];
rz(2.3123815) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(-2.1270027) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0907959) q[0];
sx q[0];
rz(-1.3279339) q[0];
sx q[0];
rz(2.5953351) q[0];
rz(-1.185174) q[2];
sx q[2];
rz(-1.0411388) q[2];
sx q[2];
rz(-0.5822863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18568745) q[1];
sx q[1];
rz(-0.57480156) q[1];
sx q[1];
rz(1.5622524) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7144373) q[3];
sx q[3];
rz(-0.82089409) q[3];
sx q[3];
rz(0.088975541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0097051) q[2];
sx q[2];
rz(-1.2010801) q[2];
sx q[2];
rz(-0.76060549) q[2];
rz(-2.6724114) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(2.40707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0483265) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(-0.58116466) q[0];
rz(-1.5900853) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(-2.4466628) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7030554) q[0];
sx q[0];
rz(-1.8037272) q[0];
sx q[0];
rz(1.2715879) q[0];
rz(1.6648817) q[2];
sx q[2];
rz(-1.861899) q[2];
sx q[2];
rz(1.3520319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1174005) q[1];
sx q[1];
rz(-2.1139189) q[1];
sx q[1];
rz(-1.5195816) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7534689) q[3];
sx q[3];
rz(-1.8032866) q[3];
sx q[3];
rz(-2.7894265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.86357) q[2];
sx q[2];
rz(-1.6081955) q[2];
sx q[2];
rz(2.8387873) q[2];
rz(1.4263724) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(-0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702174) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(-2.4466178) q[0];
rz(-2.8547844) q[1];
sx q[1];
rz(-2.5359055) q[1];
sx q[1];
rz(1.274775) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5329842) q[0];
sx q[0];
rz(-1.2600945) q[0];
sx q[0];
rz(-1.0644738) q[0];
rz(1.1026682) q[2];
sx q[2];
rz(-1.0777567) q[2];
sx q[2];
rz(1.6690892) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4732237) q[1];
sx q[1];
rz(-2.372207) q[1];
sx q[1];
rz(1.046087) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7872179) q[3];
sx q[3];
rz(-2.5297926) q[3];
sx q[3];
rz(2.967088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0605165) q[2];
sx q[2];
rz(-2.923107) q[2];
sx q[2];
rz(-2.1742353) q[2];
rz(-2.4216912) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(-1.2758183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(1.6946633) q[0];
rz(-0.46724304) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(-1.1955059) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0366096) q[0];
sx q[0];
rz(-0.65147841) q[0];
sx q[0];
rz(2.3250513) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81584064) q[2];
sx q[2];
rz(-1.8058449) q[2];
sx q[2];
rz(1.5790958) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.171716) q[1];
sx q[1];
rz(-2.6387353) q[1];
sx q[1];
rz(-1.4420532) q[1];
rz(-1.9105114) q[3];
sx q[3];
rz(-0.038148316) q[3];
sx q[3];
rz(-1.7496141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3168827) q[2];
sx q[2];
rz(-0.83493817) q[2];
sx q[2];
rz(-2.6173124) q[2];
rz(-0.25238642) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(-0.061080385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97279945) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(-1.9148781) q[0];
rz(-2.3854158) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(0.39047584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3813864) q[0];
sx q[0];
rz(-1.1293656) q[0];
sx q[0];
rz(1.8618097) q[0];
x q[1];
rz(-1.5323213) q[2];
sx q[2];
rz(-1.1270102) q[2];
sx q[2];
rz(-2.2930068) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9041216) q[1];
sx q[1];
rz(-1.7842222) q[1];
sx q[1];
rz(-1.2280812) q[1];
rz(-pi) q[2];
rz(-2.195695) q[3];
sx q[3];
rz(-0.74453738) q[3];
sx q[3];
rz(0.27305279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8172424) q[2];
sx q[2];
rz(-2.4411026) q[2];
sx q[2];
rz(-1.7895169) q[2];
rz(-0.17524854) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.2043183) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7290203) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(2.4793258) q[0];
rz(0.63016713) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(-2.9060649) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92162392) q[0];
sx q[0];
rz(-1.7818091) q[0];
sx q[0];
rz(-1.3930265) q[0];
rz(-0.56925242) q[2];
sx q[2];
rz(-1.4919089) q[2];
sx q[2];
rz(-1.9663262) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9743226) q[1];
sx q[1];
rz(-3.0758) q[1];
sx q[1];
rz(-1.2029103) q[1];
rz(-pi) q[2];
rz(1.730758) q[3];
sx q[3];
rz(-2.4387038) q[3];
sx q[3];
rz(0.95332586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0528637) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(0.55727422) q[2];
rz(0.82734621) q[3];
sx q[3];
rz(-1.6297623) q[3];
sx q[3];
rz(1.6415589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2130704) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(-0.56761566) q[0];
rz(0.40191832) q[1];
sx q[1];
rz(-0.64359507) q[1];
sx q[1];
rz(-1.7793122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0716996) q[0];
sx q[0];
rz(-1.9687459) q[0];
sx q[0];
rz(2.0073298) q[0];
x q[1];
rz(0.2593597) q[2];
sx q[2];
rz(-0.96732124) q[2];
sx q[2];
rz(0.063613907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4835254) q[1];
sx q[1];
rz(-2.1999252) q[1];
sx q[1];
rz(3.0444798) q[1];
rz(0.79093334) q[3];
sx q[3];
rz(-0.91703992) q[3];
sx q[3];
rz(-0.66822754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7207429) q[2];
sx q[2];
rz(-1.9385312) q[2];
sx q[2];
rz(0.26025772) q[2];
rz(0.69089729) q[3];
sx q[3];
rz(-2.2862209) q[3];
sx q[3];
rz(1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33547587) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(0.43329049) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(-2.2122907) q[2];
sx q[2];
rz(-2.1718696) q[2];
sx q[2];
rz(1.0739506) q[2];
rz(-1.7078441) q[3];
sx q[3];
rz(-0.78513405) q[3];
sx q[3];
rz(0.76754192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

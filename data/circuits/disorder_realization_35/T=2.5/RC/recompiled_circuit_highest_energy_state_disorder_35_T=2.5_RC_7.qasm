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
rz(0.66134557) q[1];
sx q[1];
rz(-1.1868492) q[1];
sx q[1];
rz(-2.7050736) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8494328) q[0];
sx q[0];
rz(-1.8884553) q[0];
sx q[0];
rz(-3.0374797) q[0];
rz(2.4869957) q[2];
sx q[2];
rz(-0.67383781) q[2];
sx q[2];
rz(2.0430273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3665665) q[1];
sx q[1];
rz(-0.70423755) q[1];
sx q[1];
rz(-1.6338867) q[1];
rz(-pi) q[2];
rz(2.9645676) q[3];
sx q[3];
rz(-2.2469871) q[3];
sx q[3];
rz(-0.9392304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6785757) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(2.2185745) q[2];
rz(0.065741278) q[3];
sx q[3];
rz(-1.8140503) q[3];
sx q[3];
rz(1.340723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(3.104082) q[0];
rz(0.58049479) q[1];
sx q[1];
rz(-2.4721804) q[1];
sx q[1];
rz(-2.4494749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624533) q[0];
sx q[0];
rz(-1.395749) q[0];
sx q[0];
rz(1.2412423) q[0];
rz(-pi) q[1];
rz(1.0672928) q[2];
sx q[2];
rz(-1.8177336) q[2];
sx q[2];
rz(-1.4114789) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2158606) q[1];
sx q[1];
rz(-2.070148) q[1];
sx q[1];
rz(2.4768922) q[1];
rz(-pi) q[2];
rz(0.75196079) q[3];
sx q[3];
rz(-1.6979361) q[3];
sx q[3];
rz(-0.76691915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6633501) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.3956068) q[2];
rz(-0.1869525) q[3];
sx q[3];
rz(-1.170265) q[3];
sx q[3];
rz(2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.66889399) q[0];
sx q[0];
rz(-0.63303328) q[0];
sx q[0];
rz(-2.582666) q[0];
rz(-2.3129297) q[1];
sx q[1];
rz(-1.5907954) q[1];
sx q[1];
rz(2.8899946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54788816) q[0];
sx q[0];
rz(-1.761709) q[0];
sx q[0];
rz(-2.6949791) q[0];
x q[1];
rz(2.1220268) q[2];
sx q[2];
rz(-1.2740268) q[2];
sx q[2];
rz(3.0574329) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14519599) q[1];
sx q[1];
rz(-2.4193978) q[1];
sx q[1];
rz(-1.1678371) q[1];
rz(0.67686848) q[3];
sx q[3];
rz(-1.0292316) q[3];
sx q[3];
rz(1.3071909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33637968) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(1.4317929) q[2];
rz(2.924887) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(-1.7823035) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828204) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(1.2910507) q[0];
rz(0.82921118) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(2.1270027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1433346) q[0];
sx q[0];
rz(-0.59276544) q[0];
sx q[0];
rz(-0.44500764) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5710602) q[2];
sx q[2];
rz(-2.4974365) q[2];
sx q[2];
rz(-1.2591435) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9559052) q[1];
sx q[1];
rz(-2.5667911) q[1];
sx q[1];
rz(1.5793403) q[1];
rz(-pi) q[2];
rz(1.4271554) q[3];
sx q[3];
rz(-2.3206986) q[3];
sx q[3];
rz(-3.0526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0097051) q[2];
sx q[2];
rz(-1.9405126) q[2];
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
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483265) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(-0.58116466) q[0];
rz(1.5515074) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(0.69492984) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3669156) q[0];
sx q[0];
rz(-0.37702501) q[0];
sx q[0];
rz(2.24848) q[0];
x q[1];
rz(-2.8492691) q[2];
sx q[2];
rz(-1.6609123) q[2];
sx q[2];
rz(0.19168774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.018539704) q[1];
sx q[1];
rz(-0.54529069) q[1];
sx q[1];
rz(0.084597708) q[1];
x q[2];
rz(-2.9053123) q[3];
sx q[3];
rz(-1.3930916) q[3];
sx q[3];
rz(-1.1760933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27802262) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(2.8387873) q[2];
rz(-1.7152202) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(-2.5572131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.67702174) q[0];
sx q[0];
rz(-2.7320221) q[0];
sx q[0];
rz(0.69497481) q[0];
rz(0.28680828) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(1.8668176) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0114637) q[0];
sx q[0];
rz(-2.0507567) q[0];
sx q[0];
rz(-0.35188727) q[0];
rz(-pi) q[1];
rz(-2.5996501) q[2];
sx q[2];
rz(-1.1620318) q[2];
sx q[2];
rz(-2.8084076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7949603) q[1];
sx q[1];
rz(-2.2169276) q[1];
sx q[1];
rz(2.6898795) q[1];
x q[2];
rz(-1.8096089) q[3];
sx q[3];
rz(-1.0019571) q[3];
sx q[3];
rz(-0.59900008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(-0.9673574) q[2];
rz(-2.4216912) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(1.8657743) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(-1.6946633) q[0];
rz(0.46724304) q[1];
sx q[1];
rz(-1.2930361) q[1];
sx q[1];
rz(-1.1955059) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.104983) q[0];
sx q[0];
rz(-2.4901142) q[0];
sx q[0];
rz(2.3250513) q[0];
rz(-pi) q[1];
x q[1];
rz(2.823916) q[2];
sx q[2];
rz(-2.3001849) q[2];
sx q[2];
rz(0.22401545) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1165645) q[1];
sx q[1];
rz(-2.0691074) q[1];
sx q[1];
rz(-0.070498585) q[1];
rz(1.6067664) q[3];
sx q[3];
rz(-1.5580873) q[3];
sx q[3];
rz(2.980924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82470992) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(-2.6173124) q[2];
rz(-2.8892062) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97279945) q[0];
sx q[0];
rz(-1.1058829) q[0];
sx q[0];
rz(1.9148781) q[0];
rz(2.3854158) q[1];
sx q[1];
rz(-0.94804472) q[1];
sx q[1];
rz(0.39047584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9928241) q[0];
sx q[0];
rz(-0.52339866) q[0];
sx q[0];
rz(-2.5958909) q[0];
rz(-pi) q[1];
rz(-0.4440733) q[2];
sx q[2];
rz(-1.6055428) q[2];
sx q[2];
rz(-2.4028558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8690586) q[1];
sx q[1];
rz(-0.40149899) q[1];
sx q[1];
rz(-2.1436006) q[1];
rz(-pi) q[2];
x q[2];
rz(2.195695) q[3];
sx q[3];
rz(-0.74453738) q[3];
sx q[3];
rz(2.8685399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32435027) q[2];
sx q[2];
rz(-2.4411026) q[2];
sx q[2];
rz(1.7895169) q[2];
rz(0.17524854) q[3];
sx q[3];
rz(-1.3388355) q[3];
sx q[3];
rz(1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290203) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(0.66226688) q[0];
rz(-0.63016713) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(2.9060649) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5109638) q[0];
sx q[0];
rz(-2.8665344) q[0];
sx q[0];
rz(2.4514626) q[0];
x q[1];
rz(0.56925242) q[2];
sx q[2];
rz(-1.6496837) q[2];
sx q[2];
rz(-1.9663262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79865658) q[1];
sx q[1];
rz(-1.5094116) q[1];
sx q[1];
rz(-3.117901) q[1];
rz(-pi) q[2];
rz(3.0074545) q[3];
sx q[3];
rz(-0.87867498) q[3];
sx q[3];
rz(-1.1617171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0887289) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(2.5843184) q[2];
rz(2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(1.6415589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2130704) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(-2.573977) q[0];
rz(2.7396743) q[1];
sx q[1];
rz(-0.64359507) q[1];
sx q[1];
rz(1.7793122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3338254) q[0];
sx q[0];
rz(-0.58192116) q[0];
sx q[0];
rz(-2.3533217) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2145813) q[2];
sx q[2];
rz(-2.491174) q[2];
sx q[2];
rz(-2.7678571) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4835254) q[1];
sx q[1];
rz(-2.1999252) q[1];
sx q[1];
rz(3.0444798) q[1];
rz(-2.3506593) q[3];
sx q[3];
rz(-2.2245527) q[3];
sx q[3];
rz(0.66822754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7207429) q[2];
sx q[2];
rz(-1.9385312) q[2];
sx q[2];
rz(0.26025772) q[2];
rz(2.4506954) q[3];
sx q[3];
rz(-2.2862209) q[3];
sx q[3];
rz(-1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.70788371) q[2];
sx q[2];
rz(-2.08692) q[2];
sx q[2];
rz(-0.89649123) q[2];
rz(-3.0058849) q[3];
sx q[3];
rz(-0.79499034) q[3];
sx q[3];
rz(-2.181481) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(1.3833157) q[0];
sx q[0];
rz(-1.436469) q[0];
sx q[0];
rz(0.96487784) q[0];
rz(2.3896253) q[1];
sx q[1];
rz(-2.7120092) q[1];
sx q[1];
rz(1.0493976) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0239068) q[0];
sx q[0];
rz(-1.2313594) q[0];
sx q[0];
rz(2.1855559) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8951178) q[2];
sx q[2];
rz(-1.4359546) q[2];
sx q[2];
rz(-0.91712778) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.4360848) q[1];
sx q[1];
rz(-0.75172808) q[1];
sx q[1];
rz(-2.1421894) q[1];
x q[2];
rz(-1.6635632) q[3];
sx q[3];
rz(-2.1481107) q[3];
sx q[3];
rz(-2.4924459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6875978) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(0.018608658) q[2];
rz(-0.54443693) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(2.4936567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675156) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(-0.53502214) q[0];
rz(2.6016443) q[1];
sx q[1];
rz(-0.63553634) q[1];
sx q[1];
rz(1.3998869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2311418) q[0];
sx q[0];
rz(-0.60908356) q[0];
sx q[0];
rz(-0.85394359) q[0];
rz(-pi) q[1];
rz(2.3154357) q[2];
sx q[2];
rz(-2.3914876) q[2];
sx q[2];
rz(0.78924417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78838879) q[1];
sx q[1];
rz(-1.3223159) q[1];
sx q[1];
rz(-0.39239359) q[1];
rz(-pi) q[2];
rz(-2.1407645) q[3];
sx q[3];
rz(-2.1080351) q[3];
sx q[3];
rz(1.1972103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0743559) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(2.2155217) q[2];
rz(-0.98008424) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493018) q[0];
sx q[0];
rz(-1.8632357) q[0];
sx q[0];
rz(1.5128304) q[0];
rz(0.28383645) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(2.0294752) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6215206) q[0];
sx q[0];
rz(-0.010082897) q[0];
sx q[0];
rz(-2.3501189) q[0];
rz(2.1826964) q[2];
sx q[2];
rz(-3.009353) q[2];
sx q[2];
rz(-2.7041777) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92924547) q[1];
sx q[1];
rz(-2.622859) q[1];
sx q[1];
rz(-0.38350819) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9341957) q[3];
sx q[3];
rz(-1.3280686) q[3];
sx q[3];
rz(3.0128765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6185559) q[2];
sx q[2];
rz(-1.3448389) q[2];
sx q[2];
rz(2.0396566) q[2];
rz(-0.88895041) q[3];
sx q[3];
rz(-0.44822732) q[3];
sx q[3];
rz(-2.2811208) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76982826) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(0.68921047) q[0];
rz(-2.8269732) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(-1.7291732) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012216598) q[0];
sx q[0];
rz(-2.3123992) q[0];
sx q[0];
rz(1.1393121) q[0];
rz(-pi) q[1];
rz(-1.6666404) q[2];
sx q[2];
rz(-1.7828701) q[2];
sx q[2];
rz(2.1532358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3312348) q[1];
sx q[1];
rz(-0.38429934) q[1];
sx q[1];
rz(-2.9879301) q[1];
rz(-pi) q[2];
rz(-2.9474218) q[3];
sx q[3];
rz(-1.0946858) q[3];
sx q[3];
rz(0.74060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7249001) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(0.55595428) q[2];
rz(-1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(-2.8506193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.81412643) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(2.5471174) q[0];
rz(0.56198436) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(2.1655703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291479) q[0];
sx q[0];
rz(-1.9322104) q[0];
sx q[0];
rz(-2.0301719) q[0];
rz(1.541829) q[2];
sx q[2];
rz(-0.81695643) q[2];
sx q[2];
rz(-3.0532233) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0208066) q[1];
sx q[1];
rz(-2.401742) q[1];
sx q[1];
rz(2.6915361) q[1];
x q[2];
rz(-2.8401883) q[3];
sx q[3];
rz(-1.3632953) q[3];
sx q[3];
rz(2.1077771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6550265) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(0.61069926) q[2];
rz(-0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(-1.3293728) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139451) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(1.6754643) q[0];
rz(-2.3620391) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(-0.73878845) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2125428) q[0];
sx q[0];
rz(-1.7746141) q[0];
sx q[0];
rz(-1.2840367) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30679484) q[2];
sx q[2];
rz(-1.5606797) q[2];
sx q[2];
rz(-0.9900569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14911095) q[1];
sx q[1];
rz(-2.457239) q[1];
sx q[1];
rz(-2.1696287) q[1];
rz(-pi) q[2];
rz(0.8035369) q[3];
sx q[3];
rz(-0.49852405) q[3];
sx q[3];
rz(0.62832181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5256727) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(-0.8626779) q[2];
rz(0.45977965) q[3];
sx q[3];
rz(-1.6254057) q[3];
sx q[3];
rz(-1.1427243) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2343242) q[0];
sx q[0];
rz(-2.7643804) q[0];
sx q[0];
rz(2.1211076) q[0];
rz(-0.057295784) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(0.78757706) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60815281) q[0];
sx q[0];
rz(-0.69204563) q[0];
sx q[0];
rz(-0.37779053) q[0];
rz(-pi) q[1];
rz(-0.26489218) q[2];
sx q[2];
rz(-0.26703003) q[2];
sx q[2];
rz(1.1970113) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0021506) q[1];
sx q[1];
rz(-2.6137335) q[1];
sx q[1];
rz(-1.2107641) q[1];
rz(1.6224788) q[3];
sx q[3];
rz(-0.58850901) q[3];
sx q[3];
rz(1.9543598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3202177) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(-0.58464948) q[2];
rz(2.4892877) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(-1.8023796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.898191) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(0.32824326) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(-1.6627056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2770689) q[0];
sx q[0];
rz(-1.796699) q[0];
sx q[0];
rz(0.27033131) q[0];
rz(-pi) q[1];
rz(-0.72490643) q[2];
sx q[2];
rz(-2.4838243) q[2];
sx q[2];
rz(2.7810514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.11301576) q[1];
sx q[1];
rz(-1.8488819) q[1];
sx q[1];
rz(-0.098829513) q[1];
x q[2];
rz(-2.0746873) q[3];
sx q[3];
rz(-0.41620884) q[3];
sx q[3];
rz(1.4598626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0216003) q[2];
sx q[2];
rz(-1.768521) q[2];
sx q[2];
rz(-0.78869406) q[2];
rz(-2.9005519) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(-0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929844) q[0];
sx q[0];
rz(-2.6032175) q[0];
sx q[0];
rz(2.1797144) q[0];
rz(0.23421639) q[1];
sx q[1];
rz(-1.1385671) q[1];
sx q[1];
rz(-0.73807565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80958145) q[0];
sx q[0];
rz(-1.9652848) q[0];
sx q[0];
rz(-1.8277192) q[0];
rz(-0.57508399) q[2];
sx q[2];
rz(-2.421139) q[2];
sx q[2];
rz(-1.6937814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2710423) q[1];
sx q[1];
rz(-1.1525453) q[1];
sx q[1];
rz(-2.8224432) q[1];
x q[2];
rz(-1.8515737) q[3];
sx q[3];
rz(-2.2178136) q[3];
sx q[3];
rz(-0.10296497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48478475) q[2];
sx q[2];
rz(-2.1190376) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(1.3336746) q[3];
sx q[3];
rz(-2.1402054) q[3];
sx q[3];
rz(2.8792152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(1.242908) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(-2.2743478) q[0];
rz(2.1278837) q[1];
sx q[1];
rz(-1.3288682) q[1];
sx q[1];
rz(2.0713461) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27425048) q[0];
sx q[0];
rz(-0.19495067) q[0];
sx q[0];
rz(2.250953) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81999166) q[2];
sx q[2];
rz(-0.89897663) q[2];
sx q[2];
rz(-1.9209678) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6079191) q[1];
sx q[1];
rz(-1.9551139) q[1];
sx q[1];
rz(-2.584105) q[1];
x q[2];
rz(0.13954969) q[3];
sx q[3];
rz(-1.680239) q[3];
sx q[3];
rz(-2.8982996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9762207) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(1.2316068) q[2];
rz(1.6407137) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-2.4098082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3043542) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(2.7267743) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(1.5157221) q[2];
sx q[2];
rz(-2.016042) q[2];
sx q[2];
rz(2.945937) q[2];
rz(-2.5255193) q[3];
sx q[3];
rz(-0.85312927) q[3];
sx q[3];
rz(2.7240172) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

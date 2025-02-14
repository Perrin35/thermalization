OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3040721) q[0];
sx q[0];
rz(-2.4763595) q[0];
sx q[0];
rz(-1.2256149) q[0];
rz(2.4576814) q[1];
sx q[1];
rz(-0.39165762) q[1];
sx q[1];
rz(-2.439523) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64710275) q[0];
sx q[0];
rz(-2.4328874) q[0];
sx q[0];
rz(2.9604218) q[0];
rz(-pi) q[1];
rz(1.4058787) q[2];
sx q[2];
rz(-1.5497297) q[2];
sx q[2];
rz(1.5738459) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7309642) q[1];
sx q[1];
rz(-2.0515039) q[1];
sx q[1];
rz(2.5334873) q[1];
rz(-pi) q[2];
rz(0.87907531) q[3];
sx q[3];
rz(-1.4180776) q[3];
sx q[3];
rz(2.9756551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94032732) q[2];
sx q[2];
rz(-1.6028812) q[2];
sx q[2];
rz(0.24844696) q[2];
rz(-2.0072319) q[3];
sx q[3];
rz(-0.13392197) q[3];
sx q[3];
rz(-1.1321446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092875384) q[0];
sx q[0];
rz(-2.5889914) q[0];
sx q[0];
rz(-1.5930814) q[0];
rz(-0.50367194) q[1];
sx q[1];
rz(-1.9182938) q[1];
sx q[1];
rz(-0.37685397) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5706963) q[0];
sx q[0];
rz(-2.6691169) q[0];
sx q[0];
rz(0.04224472) q[0];
x q[1];
rz(-2.3289248) q[2];
sx q[2];
rz(-2.118131) q[2];
sx q[2];
rz(3.0603882) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8736564) q[1];
sx q[1];
rz(-2.5940478) q[1];
sx q[1];
rz(2.646562) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4941856) q[3];
sx q[3];
rz(-1.9025246) q[3];
sx q[3];
rz(-0.89523098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58439955) q[2];
sx q[2];
rz(-2.445745) q[2];
sx q[2];
rz(-2.2759571) q[2];
rz(1.4731167) q[3];
sx q[3];
rz(-0.98554635) q[3];
sx q[3];
rz(2.1540811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5675548) q[0];
sx q[0];
rz(-1.8542629) q[0];
sx q[0];
rz(-1.0512742) q[0];
rz(1.4569262) q[1];
sx q[1];
rz(-1.8345865) q[1];
sx q[1];
rz(-1.5164703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.595436) q[0];
sx q[0];
rz(-1.8519028) q[0];
sx q[0];
rz(-2.3670822) q[0];
rz(-pi) q[1];
rz(0.77645723) q[2];
sx q[2];
rz(-2.2946022) q[2];
sx q[2];
rz(-3.0645097) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1976002) q[1];
sx q[1];
rz(-1.1901344) q[1];
sx q[1];
rz(2.073111) q[1];
rz(1.5169237) q[3];
sx q[3];
rz(-1.5098715) q[3];
sx q[3];
rz(-2.1046154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7463344) q[2];
sx q[2];
rz(-1.0552152) q[2];
sx q[2];
rz(0.22221097) q[2];
rz(-2.639751) q[3];
sx q[3];
rz(-1.4072489) q[3];
sx q[3];
rz(-0.80880222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.89268452) q[0];
sx q[0];
rz(-0.48766708) q[0];
sx q[0];
rz(-1.1647613) q[0];
rz(-0.89272967) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(2.8597615) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.079275) q[0];
sx q[0];
rz(-1.4505511) q[0];
sx q[0];
rz(-2.6577302) q[0];
rz(-pi) q[1];
rz(0.64617363) q[2];
sx q[2];
rz(-2.3149688) q[2];
sx q[2];
rz(-1.3767565) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1077256) q[1];
sx q[1];
rz(-1.733779) q[1];
sx q[1];
rz(-0.070734338) q[1];
x q[2];
rz(0.025679703) q[3];
sx q[3];
rz(-1.1246944) q[3];
sx q[3];
rz(-3.1334973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.07936) q[2];
sx q[2];
rz(-2.5793109) q[2];
sx q[2];
rz(-2.6083561) q[2];
rz(1.0821292) q[3];
sx q[3];
rz(-2.5366668) q[3];
sx q[3];
rz(-2.3281039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.767652) q[0];
sx q[0];
rz(-0.39893183) q[0];
sx q[0];
rz(1.8178513) q[0];
rz(-2.3228877) q[1];
sx q[1];
rz(-0.74344126) q[1];
sx q[1];
rz(2.0735819) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0518347) q[0];
sx q[0];
rz(-1.4115507) q[0];
sx q[0];
rz(-1.490504) q[0];
rz(-0.6650659) q[2];
sx q[2];
rz(-2.1882727) q[2];
sx q[2];
rz(2.0047052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6883478) q[1];
sx q[1];
rz(-1.6736341) q[1];
sx q[1];
rz(2.4842841) q[1];
rz(1.8858484) q[3];
sx q[3];
rz(-0.49335155) q[3];
sx q[3];
rz(2.3137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1800804) q[2];
sx q[2];
rz(-1.1137806) q[2];
sx q[2];
rz(-0.86165825) q[2];
rz(-1.0263475) q[3];
sx q[3];
rz(-1.9828826) q[3];
sx q[3];
rz(-1.1177184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3468129) q[0];
sx q[0];
rz(-2.3220799) q[0];
sx q[0];
rz(-0.89282194) q[0];
rz(0.18094856) q[1];
sx q[1];
rz(-0.67829689) q[1];
sx q[1];
rz(-0.13042626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.923188) q[0];
sx q[0];
rz(-1.4124065) q[0];
sx q[0];
rz(3.0322101) q[0];
x q[1];
rz(0.60575859) q[2];
sx q[2];
rz(-1.1446627) q[2];
sx q[2];
rz(-2.5832165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0196258) q[1];
sx q[1];
rz(-1.6000012) q[1];
sx q[1];
rz(-1.3335511) q[1];
rz(-0.78108861) q[3];
sx q[3];
rz(-1.4295331) q[3];
sx q[3];
rz(-3.0716592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9170407) q[2];
sx q[2];
rz(-1.2265393) q[2];
sx q[2];
rz(-0.87635931) q[2];
rz(-1.8996436) q[3];
sx q[3];
rz(-0.94049898) q[3];
sx q[3];
rz(-2.131264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4001615) q[0];
sx q[0];
rz(-0.09859666) q[0];
sx q[0];
rz(-2.7778991) q[0];
rz(-0.74186507) q[1];
sx q[1];
rz(-1.0262841) q[1];
sx q[1];
rz(-0.40282869) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2296933) q[0];
sx q[0];
rz(-2.6586164) q[0];
sx q[0];
rz(-2.1080024) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3260457) q[2];
sx q[2];
rz(-2.6167653) q[2];
sx q[2];
rz(-1.8458402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5093358) q[1];
sx q[1];
rz(-1.4493353) q[1];
sx q[1];
rz(-3.0640825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36084036) q[3];
sx q[3];
rz(-2.4746102) q[3];
sx q[3];
rz(2.9426394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7824629) q[2];
sx q[2];
rz(-2.7733347) q[2];
sx q[2];
rz(1.5472319) q[2];
rz(-0.83526978) q[3];
sx q[3];
rz(-0.95649496) q[3];
sx q[3];
rz(-2.6529151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816724) q[0];
sx q[0];
rz(-1.8052354) q[0];
sx q[0];
rz(-0.51505995) q[0];
rz(-1.7022279) q[1];
sx q[1];
rz(-1.8481588) q[1];
sx q[1];
rz(2.7424367) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418509) q[0];
sx q[0];
rz(-3.0106595) q[0];
sx q[0];
rz(-3.127742) q[0];
rz(-pi) q[1];
rz(-0.84206669) q[2];
sx q[2];
rz(-1.7532323) q[2];
sx q[2];
rz(-1.4562171) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6281575) q[1];
sx q[1];
rz(-1.6353459) q[1];
sx q[1];
rz(1.3099531) q[1];
x q[2];
rz(-0.018980515) q[3];
sx q[3];
rz(-0.06566635) q[3];
sx q[3];
rz(1.1733023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1601552) q[2];
sx q[2];
rz(-1.258054) q[2];
sx q[2];
rz(2.0951648) q[2];
rz(0.6662755) q[3];
sx q[3];
rz(-0.25686887) q[3];
sx q[3];
rz(2.090914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94205725) q[0];
sx q[0];
rz(-3.0362447) q[0];
sx q[0];
rz(-2.5355205) q[0];
rz(-2.3145158) q[1];
sx q[1];
rz(-1.3304293) q[1];
sx q[1];
rz(-1.3333295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5904986) q[0];
sx q[0];
rz(-0.37305377) q[0];
sx q[0];
rz(-0.28753745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6624062) q[2];
sx q[2];
rz(-0.80305099) q[2];
sx q[2];
rz(-0.89862862) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5659396) q[1];
sx q[1];
rz(-2.2350532) q[1];
sx q[1];
rz(-0.14845005) q[1];
rz(-0.4584109) q[3];
sx q[3];
rz(-1.6473149) q[3];
sx q[3];
rz(-2.1297701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0102319) q[2];
sx q[2];
rz(-0.92326814) q[2];
sx q[2];
rz(-0.38491797) q[2];
rz(-2.5108003) q[3];
sx q[3];
rz(-1.2814949) q[3];
sx q[3];
rz(1.7990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25403062) q[0];
sx q[0];
rz(-1.7429202) q[0];
sx q[0];
rz(-1.4400462) q[0];
rz(-0.74132672) q[1];
sx q[1];
rz(-1.7404375) q[1];
sx q[1];
rz(-2.8828566) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7546623) q[0];
sx q[0];
rz(-1.566769) q[0];
sx q[0];
rz(-1.3248825) q[0];
x q[1];
rz(1.0503631) q[2];
sx q[2];
rz(-2.5602617) q[2];
sx q[2];
rz(-0.54301013) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96232167) q[1];
sx q[1];
rz(-1.2291161) q[1];
sx q[1];
rz(-0.82247295) q[1];
rz(-pi) q[2];
rz(-0.77787598) q[3];
sx q[3];
rz(-1.8612783) q[3];
sx q[3];
rz(-2.8992407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4197293) q[2];
sx q[2];
rz(-2.0294919) q[2];
sx q[2];
rz(-1.3191684) q[2];
rz(0.81926528) q[3];
sx q[3];
rz(-1.1300488) q[3];
sx q[3];
rz(-2.5949902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7405613) q[0];
sx q[0];
rz(-0.87420976) q[0];
sx q[0];
rz(2.5921205) q[0];
rz(1.4389379) q[1];
sx q[1];
rz(-0.66822744) q[1];
sx q[1];
rz(-2.6034036) q[1];
rz(3.13065) q[2];
sx q[2];
rz(-0.86075114) q[2];
sx q[2];
rz(1.1759947) q[2];
rz(-2.1323754) q[3];
sx q[3];
rz(-1.781096) q[3];
sx q[3];
rz(-2.5566035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

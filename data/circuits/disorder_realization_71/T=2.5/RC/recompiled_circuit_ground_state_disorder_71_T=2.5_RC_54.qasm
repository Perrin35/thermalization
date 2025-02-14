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
rz(3.8068258) q[0];
sx q[0];
rz(8.1991631) q[0];
rz(-0.6839112) q[1];
sx q[1];
rz(3.5332503) q[1];
sx q[1];
rz(11.864301) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618781) q[0];
sx q[0];
rz(-1.6883381) q[0];
sx q[0];
rz(-2.4410309) q[0];
rz(-1.735714) q[2];
sx q[2];
rz(-1.5497297) q[2];
sx q[2];
rz(1.5738459) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4106284) q[1];
sx q[1];
rz(-1.0900888) q[1];
sx q[1];
rz(-0.60810535) q[1];
rz(-pi) q[2];
rz(2.2625173) q[3];
sx q[3];
rz(-1.4180776) q[3];
sx q[3];
rz(-2.9756551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94032732) q[2];
sx q[2];
rz(-1.6028812) q[2];
sx q[2];
rz(-0.24844696) q[2];
rz(-1.1343608) q[3];
sx q[3];
rz(-0.13392197) q[3];
sx q[3];
rz(-2.0094481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0487173) q[0];
sx q[0];
rz(-2.5889914) q[0];
sx q[0];
rz(1.5930814) q[0];
rz(-0.50367194) q[1];
sx q[1];
rz(-1.2232989) q[1];
sx q[1];
rz(0.37685397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6183313) q[0];
sx q[0];
rz(-1.0987765) q[0];
sx q[0];
rz(-1.5492155) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84555618) q[2];
sx q[2];
rz(-2.2396002) q[2];
sx q[2];
rz(-1.9921274) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8736564) q[1];
sx q[1];
rz(-2.5940478) q[1];
sx q[1];
rz(-0.4950306) q[1];
rz(-0.33263388) q[3];
sx q[3];
rz(-1.4983699) q[3];
sx q[3];
rz(2.4410332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58439955) q[2];
sx q[2];
rz(-2.445745) q[2];
sx q[2];
rz(2.2759571) q[2];
rz(1.668476) q[3];
sx q[3];
rz(-0.98554635) q[3];
sx q[3];
rz(0.98751155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57403785) q[0];
sx q[0];
rz(-1.8542629) q[0];
sx q[0];
rz(-2.0903184) q[0];
rz(-1.4569262) q[1];
sx q[1];
rz(-1.3070062) q[1];
sx q[1];
rz(1.6251224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54615669) q[0];
sx q[0];
rz(-1.8519028) q[0];
sx q[0];
rz(-0.77451046) q[0];
rz(-pi) q[1];
rz(-0.90041884) q[2];
sx q[2];
rz(-2.134765) q[2];
sx q[2];
rz(-1.054686) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96788266) q[1];
sx q[1];
rz(-2.5213402) q[1];
sx q[1];
rz(2.2642676) q[1];
rz(-pi) q[2];
rz(1.5169237) q[3];
sx q[3];
rz(-1.6317211) q[3];
sx q[3];
rz(-1.0369773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7463344) q[2];
sx q[2];
rz(-1.0552152) q[2];
sx q[2];
rz(-2.9193817) q[2];
rz(2.639751) q[3];
sx q[3];
rz(-1.4072489) q[3];
sx q[3];
rz(-2.3327904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89268452) q[0];
sx q[0];
rz(-0.48766708) q[0];
sx q[0];
rz(1.9768313) q[0];
rz(-0.89272967) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(-0.28183118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062317693) q[0];
sx q[0];
rz(-1.4505511) q[0];
sx q[0];
rz(-2.6577302) q[0];
x q[1];
rz(0.64617363) q[2];
sx q[2];
rz(-0.82662383) q[2];
sx q[2];
rz(1.3767565) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4455681) q[1];
sx q[1];
rz(-0.17754517) q[1];
sx q[1];
rz(-1.976718) q[1];
rz(-pi) q[2];
rz(-0.025679703) q[3];
sx q[3];
rz(-1.1246944) q[3];
sx q[3];
rz(-0.0080953117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.07936) q[2];
sx q[2];
rz(-0.56228176) q[2];
sx q[2];
rz(0.53323659) q[2];
rz(2.0594635) q[3];
sx q[3];
rz(-2.5366668) q[3];
sx q[3];
rz(2.3281039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.767652) q[0];
sx q[0];
rz(-0.39893183) q[0];
sx q[0];
rz(1.8178513) q[0];
rz(-0.81870493) q[1];
sx q[1];
rz(-0.74344126) q[1];
sx q[1];
rz(-2.0735819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0897579) q[0];
sx q[0];
rz(-1.4115507) q[0];
sx q[0];
rz(-1.6510886) q[0];
rz(-pi) q[1];
rz(2.4765268) q[2];
sx q[2];
rz(-0.95331999) q[2];
sx q[2];
rz(-2.0047052) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9449759) q[1];
sx q[1];
rz(-2.2240337) q[1];
sx q[1];
rz(-1.4411627) q[1];
rz(-pi) q[2];
rz(0.16509861) q[3];
sx q[3];
rz(-2.0378651) q[3];
sx q[3];
rz(1.9593173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1800804) q[2];
sx q[2];
rz(-1.1137806) q[2];
sx q[2];
rz(-2.2799344) q[2];
rz(2.1152451) q[3];
sx q[3];
rz(-1.9828826) q[3];
sx q[3];
rz(-1.1177184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3468129) q[0];
sx q[0];
rz(-2.3220799) q[0];
sx q[0];
rz(-0.89282194) q[0];
rz(-0.18094856) q[1];
sx q[1];
rz(-2.4632958) q[1];
sx q[1];
rz(3.0111664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8065211) q[0];
sx q[0];
rz(-1.6788043) q[0];
sx q[0];
rz(1.4114702) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67306913) q[2];
sx q[2];
rz(-0.72491881) q[2];
sx q[2];
rz(0.47436213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45589009) q[1];
sx q[1];
rz(-1.3336542) q[1];
sx q[1];
rz(3.1115467) q[1];
rz(-pi) q[2];
rz(-2.942286) q[3];
sx q[3];
rz(-2.350507) q[3];
sx q[3];
rz(-1.6419322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22455198) q[2];
sx q[2];
rz(-1.9150534) q[2];
sx q[2];
rz(2.2652333) q[2];
rz(-1.8996436) q[3];
sx q[3];
rz(-2.2010937) q[3];
sx q[3];
rz(-1.0103286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7414311) q[0];
sx q[0];
rz(-0.09859666) q[0];
sx q[0];
rz(-0.36369351) q[0];
rz(-2.3997276) q[1];
sx q[1];
rz(-2.1153085) q[1];
sx q[1];
rz(-0.40282869) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8217709) q[0];
sx q[0];
rz(-1.9812225) q[0];
sx q[0];
rz(-0.26217802) q[0];
x q[1];
rz(-2.7637787) q[2];
sx q[2];
rz(-1.1973518) q[2];
sx q[2];
rz(-0.4682954) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6322569) q[1];
sx q[1];
rz(-1.4493353) q[1];
sx q[1];
rz(3.0640825) q[1];
rz(-2.5066911) q[3];
sx q[3];
rz(-1.7909808) q[3];
sx q[3];
rz(-1.4815154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3591298) q[2];
sx q[2];
rz(-2.7733347) q[2];
sx q[2];
rz(1.5472319) q[2];
rz(2.3063229) q[3];
sx q[3];
rz(-0.95649496) q[3];
sx q[3];
rz(0.48867759) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35992026) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(-0.51505995) q[0];
rz(-1.7022279) q[1];
sx q[1];
rz(-1.2934338) q[1];
sx q[1];
rz(-2.7424367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8137119) q[0];
sx q[0];
rz(-1.7017168) q[0];
sx q[0];
rz(1.5726202) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3005343) q[2];
sx q[2];
rz(-2.3944602) q[2];
sx q[2];
rz(0.085937339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6281575) q[1];
sx q[1];
rz(-1.6353459) q[1];
sx q[1];
rz(1.8316395) q[1];
rz(-pi) q[2];
rz(-0.065654556) q[3];
sx q[3];
rz(-1.5695509) q[3];
sx q[3];
rz(2.7630382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9814375) q[2];
sx q[2];
rz(-1.8835386) q[2];
sx q[2];
rz(-2.0951648) q[2];
rz(0.6662755) q[3];
sx q[3];
rz(-0.25686887) q[3];
sx q[3];
rz(2.090914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1995354) q[0];
sx q[0];
rz(-0.10534795) q[0];
sx q[0];
rz(-2.5355205) q[0];
rz(-0.82707682) q[1];
sx q[1];
rz(-1.8111633) q[1];
sx q[1];
rz(-1.3333295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2829958) q[0];
sx q[0];
rz(-1.2137611) q[0];
sx q[0];
rz(-1.4602565) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1252046) q[2];
sx q[2];
rz(-2.2632709) q[2];
sx q[2];
rz(-2.8852579) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.096790678) q[1];
sx q[1];
rz(-1.4540744) q[1];
sx q[1];
rz(2.2404284) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9700432) q[3];
sx q[3];
rz(-2.6772873) q[3];
sx q[3];
rz(0.40528471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0102319) q[2];
sx q[2];
rz(-0.92326814) q[2];
sx q[2];
rz(-2.7566747) q[2];
rz(2.5108003) q[3];
sx q[3];
rz(-1.8600978) q[3];
sx q[3];
rz(1.7990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25403062) q[0];
sx q[0];
rz(-1.3986724) q[0];
sx q[0];
rz(-1.7015464) q[0];
rz(-2.4002659) q[1];
sx q[1];
rz(-1.4011551) q[1];
sx q[1];
rz(0.25873605) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16782204) q[0];
sx q[0];
rz(-0.24594618) q[0];
sx q[0];
rz(-1.587338) q[0];
x q[1];
rz(-2.825794) q[2];
sx q[2];
rz(-1.074203) q[2];
sx q[2];
rz(1.1441355) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.179271) q[1];
sx q[1];
rz(-1.9124766) q[1];
sx q[1];
rz(2.3191197) q[1];
rz(-pi) q[2];
rz(2.7388938) q[3];
sx q[3];
rz(-0.81953632) q[3];
sx q[3];
rz(-2.0961026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72186333) q[2];
sx q[2];
rz(-1.1121007) q[2];
sx q[2];
rz(1.3191684) q[2];
rz(2.3223274) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(0.54660249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4010314) q[0];
sx q[0];
rz(-2.2673829) q[0];
sx q[0];
rz(-0.54947214) q[0];
rz(-1.4389379) q[1];
sx q[1];
rz(-2.4733652) q[1];
sx q[1];
rz(0.53818902) q[1];
rz(0.010942608) q[2];
sx q[2];
rz(-2.2808415) q[2];
sx q[2];
rz(-1.965598) q[2];
rz(1.9520252) q[3];
sx q[3];
rz(-0.59567957) q[3];
sx q[3];
rz(2.4761562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

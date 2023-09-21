OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(-0.95681325) q[0];
sx q[0];
rz(1.4332888) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8960436) q[0];
sx q[0];
rz(-0.25078068) q[0];
sx q[0];
rz(-1.8401237) q[0];
rz(1.1510017) q[2];
sx q[2];
rz(-1.3090054) q[2];
sx q[2];
rz(-2.7671438) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5073587) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(-1.5276315) q[1];
x q[2];
rz(-2.5979795) q[3];
sx q[3];
rz(-1.3136778) q[3];
sx q[3];
rz(1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-0.84709644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23628274) q[0];
sx q[0];
rz(-0.61203996) q[0];
sx q[0];
rz(2.406714) q[0];
x q[1];
rz(-1.2752091) q[2];
sx q[2];
rz(-0.81381961) q[2];
sx q[2];
rz(0.57397599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(0.54547711) q[1];
rz(-pi) q[2];
rz(-1.133425) q[3];
sx q[3];
rz(-1.8873894) q[3];
sx q[3];
rz(-0.40382995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29887154) q[0];
sx q[0];
rz(-1.5639595) q[0];
sx q[0];
rz(0.30928916) q[0];
x q[1];
rz(-2.7408319) q[2];
sx q[2];
rz(-0.71360613) q[2];
sx q[2];
rz(-1.00373) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1151162) q[1];
sx q[1];
rz(-0.83362245) q[1];
sx q[1];
rz(1.923418) q[1];
rz(0.35964386) q[3];
sx q[3];
rz(-0.76562866) q[3];
sx q[3];
rz(-1.5639203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30806914) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-0.72675881) q[2];
rz(-0.99749342) q[3];
sx q[3];
rz(-2.5163429) q[3];
sx q[3];
rz(-1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(2.3153268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59506455) q[0];
sx q[0];
rz(-2.1273158) q[0];
sx q[0];
rz(0.80842774) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(0.61566478) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1288209) q[1];
sx q[1];
rz(-1.9073309) q[1];
sx q[1];
rz(3.0613042) q[1];
rz(-pi) q[2];
rz(1.1816979) q[3];
sx q[3];
rz(-0.97550052) q[3];
sx q[3];
rz(-2.1967595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(-2.5240135) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78248258) q[0];
sx q[0];
rz(-1.5251274) q[0];
sx q[0];
rz(-3.111905) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64230201) q[2];
sx q[2];
rz(-1.2080492) q[2];
sx q[2];
rz(1.9833267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51296556) q[1];
sx q[1];
rz(-1.8497397) q[1];
sx q[1];
rz(0.97475027) q[1];
rz(0.17554749) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(0.73568425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(2.9079672) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-2.5036507) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(-0.18009137) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20443944) q[0];
sx q[0];
rz(-2.3817872) q[0];
sx q[0];
rz(-0.8888437) q[0];
rz(-pi) q[1];
rz(-1.300235) q[2];
sx q[2];
rz(-2.5205043) q[2];
sx q[2];
rz(-1.6294711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1485968) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(-2.6246043) q[1];
x q[2];
rz(-3.1230934) q[3];
sx q[3];
rz(-1.0780932) q[3];
sx q[3];
rz(2.8871356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-0.040239008) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38567625) q[0];
sx q[0];
rz(-1.3323116) q[0];
sx q[0];
rz(-3.0576865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9769386) q[2];
sx q[2];
rz(-1.5283661) q[2];
sx q[2];
rz(0.050780642) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-0.95655671) q[1];
sx q[1];
rz(2.9700301) q[1];
rz(-pi) q[2];
rz(-0.76544806) q[3];
sx q[3];
rz(-0.42740373) q[3];
sx q[3];
rz(-3.1020853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(-0.17383943) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.404495) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(-1.5054024) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36639402) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(-1.9154857) q[0];
x q[1];
rz(0.90754857) q[2];
sx q[2];
rz(-2.2758256) q[2];
sx q[2];
rz(-2.7570587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2124894) q[1];
sx q[1];
rz(-2.266532) q[1];
sx q[1];
rz(1.6131669) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96945854) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(0.45575842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4432916) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(2.495893) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.2649149) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1160693) q[0];
sx q[0];
rz(-2.2795296) q[0];
sx q[0];
rz(-2.619333) q[0];
rz(-pi) q[1];
rz(-2.0090946) q[2];
sx q[2];
rz(-1.0602789) q[2];
sx q[2];
rz(2.8709708) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1688655) q[1];
sx q[1];
rz(-1.7164413) q[1];
sx q[1];
rz(-1.5345854) q[1];
x q[2];
rz(1.3690788) q[3];
sx q[3];
rz(-1.1551876) q[3];
sx q[3];
rz(-0.49858529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-2.6169422) q[2];
rz(2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(-2.8961704) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.4642749) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9058162) q[0];
sx q[0];
rz(-0.24484867) q[0];
sx q[0];
rz(-2.66923) q[0];
x q[1];
rz(-0.22428959) q[2];
sx q[2];
rz(-0.99594342) q[2];
sx q[2];
rz(-1.6155417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.97800868) q[1];
sx q[1];
rz(-2.2469236) q[1];
sx q[1];
rz(2.7917557) q[1];
x q[2];
rz(0.085300307) q[3];
sx q[3];
rz(-1.7461667) q[3];
sx q[3];
rz(-2.8773957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2422553) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-2.0132813) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.8854234) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-2.0586661) q[2];
sx q[2];
rz(-2.2219873) q[2];
sx q[2];
rz(1.7359003) q[2];
rz(2.5916354) q[3];
sx q[3];
rz(-0.64823845) q[3];
sx q[3];
rz(1.3241495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

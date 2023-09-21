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
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(-1.9967611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.618453) q[0];
sx q[0];
rz(-1.3292399) q[0];
sx q[0];
rz(-0.068058204) q[0];
x q[1];
rz(2.152359) q[2];
sx q[2];
rz(-0.49058149) q[2];
sx q[2];
rz(1.7218334) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5073587) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(1.5276315) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54361312) q[3];
sx q[3];
rz(-1.3136778) q[3];
sx q[3];
rz(-1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.9887259) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(2.2944962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4439508) q[0];
sx q[0];
rz(-1.9662494) q[0];
sx q[0];
rz(0.4801428) q[0];
rz(-pi) q[1];
rz(-0.29909889) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(-0.15660827) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8250679) q[1];
sx q[1];
rz(-2.1148588) q[1];
sx q[1];
rz(-1.6506509) q[1];
rz(0.3470207) q[3];
sx q[3];
rz(-1.9850529) q[3];
sx q[3];
rz(-1.3115209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8427211) q[0];
sx q[0];
rz(-1.5639595) q[0];
sx q[0];
rz(2.8323035) q[0];
rz(-pi) q[1];
rz(2.7408319) q[2];
sx q[2];
rz(-2.4279865) q[2];
sx q[2];
rz(-1.00373) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5273683) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(0.36348344) q[1];
rz(-pi) q[2];
rz(-0.73266352) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(-2.8837567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(2.4148338) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(2.3153268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47047939) q[0];
sx q[0];
rz(-2.2320036) q[0];
sx q[0];
rz(-2.3040422) q[0];
x q[1];
rz(-1.516953) q[2];
sx q[2];
rz(-2.2701837) q[2];
sx q[2];
rz(0.92045036) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77376765) q[1];
sx q[1];
rz(-2.7959679) q[1];
sx q[1];
rz(1.3454382) q[1];
rz(-pi) q[2];
rz(0.51059254) q[3];
sx q[3];
rz(-2.4435352) q[3];
sx q[3];
rz(-1.576168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-0.91439247) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(0.61757913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78248258) q[0];
sx q[0];
rz(-1.5251274) q[0];
sx q[0];
rz(3.111905) q[0];
x q[1];
rz(0.64230201) q[2];
sx q[2];
rz(-1.2080492) q[2];
sx q[2];
rz(-1.158266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.87318201) q[1];
sx q[1];
rz(-1.0007443) q[1];
sx q[1];
rz(-0.33318712) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5579617) q[3];
sx q[3];
rz(-1.395263) q[3];
sx q[3];
rz(0.83735355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-0.63794199) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96930209) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(3.0517975) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(0.18009137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0952547) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(2.6020223) q[0];
rz(-pi) q[1];
rz(2.1744556) q[2];
sx q[2];
rz(-1.4146311) q[2];
sx q[2];
rz(-0.16317633) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1485968) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(-2.6246043) q[1];
rz(-pi) q[2];
rz(3.1230934) q[3];
sx q[3];
rz(-1.0780932) q[3];
sx q[3];
rz(0.25445709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(3.0684493) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-0.040239008) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0979472) q[0];
sx q[0];
rz(-2.8890434) q[0];
sx q[0];
rz(1.9027684) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9769386) q[2];
sx q[2];
rz(-1.5283661) q[2];
sx q[2];
rz(-0.050780642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8602627) q[1];
sx q[1];
rz(-1.7107692) q[1];
sx q[1];
rz(-0.94957385) q[1];
x q[2];
rz(1.2651029) q[3];
sx q[3];
rz(-1.8743268) q[3];
sx q[3];
rz(0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(-1.6361902) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36639402) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(-1.9154857) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90754857) q[2];
sx q[2];
rz(-0.86576701) q[2];
sx q[2];
rz(-2.7570587) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38547388) q[1];
sx q[1];
rz(-1.6033152) q[1];
sx q[1];
rz(0.69617747) q[1];
x q[2];
rz(2.1721341) q[3];
sx q[3];
rz(-1.8441895) q[3];
sx q[3];
rz(0.45575842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(-0.64569965) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18683534) q[0];
sx q[0];
rz(-1.9592013) q[0];
sx q[0];
rz(-2.3507621) q[0];
x q[1];
rz(0.64847704) q[2];
sx q[2];
rz(-0.65995526) q[2];
sx q[2];
rz(0.49396587) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1688655) q[1];
sx q[1];
rz(-1.4251514) q[1];
sx q[1];
rz(-1.5345854) q[1];
rz(1.7725138) q[3];
sx q[3];
rz(-1.1551876) q[3];
sx q[3];
rz(-2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(-0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
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
rz(-2.896744) q[0];
sx q[0];
rz(2.66923) q[0];
rz(-pi) q[1];
rz(2.1572838) q[2];
sx q[2];
rz(-1.7585635) q[2];
sx q[2];
rz(-0.078660065) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97800868) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-2.7917557) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7467935) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(1.0829265) q[2];
sx q[2];
rz(-2.2219873) q[2];
sx q[2];
rz(1.7359003) q[2];
rz(-1.9477378) q[3];
sx q[3];
rz(-2.1115163) q[3];
sx q[3];
rz(1.9797309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
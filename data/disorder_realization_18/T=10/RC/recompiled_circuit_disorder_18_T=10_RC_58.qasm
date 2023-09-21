OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(-1.4332888) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24554907) q[0];
sx q[0];
rz(-0.25078068) q[0];
sx q[0];
rz(-1.301469) q[0];
rz(-pi) q[1];
rz(0.28540622) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(-2.0602496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45935985) q[1];
sx q[1];
rz(-2.0233005) q[1];
sx q[1];
rz(-0.020999055) q[1];
rz(-pi) q[2];
rz(1.2727229) q[3];
sx q[3];
rz(-2.094659) q[3];
sx q[3];
rz(0.03447547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(1.2343181) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(-1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(2.2944962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69764187) q[0];
sx q[0];
rz(-1.1753433) q[0];
sx q[0];
rz(-2.6614499) q[0];
x q[1];
rz(-1.8663835) q[2];
sx q[2];
rz(-2.327773) q[2];
sx q[2];
rz(-2.5676167) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21287316) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(0.54547711) q[1];
rz(-pi) q[2];
rz(-0.3470207) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(-1.3115209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.117924) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(1.2891278) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(-1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-2.5701994) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8427211) q[0];
sx q[0];
rz(-1.5776331) q[0];
sx q[0];
rz(-0.30928916) q[0];
x q[1];
rz(-2.4685523) q[2];
sx q[2];
rz(-1.3125784) q[2];
sx q[2];
rz(0.25707993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5273683) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(-0.36348344) q[1];
x q[2];
rz(-2.4089291) q[3];
sx q[3];
rz(-1.817173) q[3];
sx q[3];
rz(0.25783595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30806914) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(0.72675881) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(2.1380651) q[0];
rz(-0.040680496) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(2.3153268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47047939) q[0];
sx q[0];
rz(-0.90958909) q[0];
sx q[0];
rz(0.83755042) q[0];
x q[1];
rz(-0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(-0.61566478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1288209) q[1];
sx q[1];
rz(-1.2342617) q[1];
sx q[1];
rz(3.0613042) q[1];
rz(-pi) q[2];
rz(-0.63185933) q[3];
sx q[3];
rz(-1.8903036) q[3];
sx q[3];
rz(2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(-0.91439247) q[2];
rz(0.10522035) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(-2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546346) q[0];
sx q[0];
rz(-1.5411396) q[0];
sx q[0];
rz(-1.6164854) q[0];
x q[1];
rz(-0.56474833) q[2];
sx q[2];
rz(-2.4167633) q[2];
sx q[2];
rz(0.85541475) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87318201) q[1];
sx q[1];
rz(-2.1408484) q[1];
sx q[1];
rz(-0.33318712) q[1];
x q[2];
rz(0.072237416) q[3];
sx q[3];
rz(-0.17599711) q[3];
sx q[3];
rz(2.3776059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.96930209) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(0.0897952) q[0];
rz(2.2604997) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0952547) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(-0.5395704) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8413576) q[2];
sx q[2];
rz(-0.62108835) q[2];
sx q[2];
rz(1.5121216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1485968) q[1];
sx q[1];
rz(-0.33278782) q[1];
sx q[1];
rz(0.51698835) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1230934) q[3];
sx q[3];
rz(-1.0780932) q[3];
sx q[3];
rz(2.8871356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2580516) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.9656666) q[2];
rz(-3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95773762) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(-1.4181597) q[0];
rz(1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0979472) q[0];
sx q[0];
rz(-2.8890434) q[0];
sx q[0];
rz(-1.9027684) q[0];
rz(-pi) q[1];
rz(1.1646541) q[2];
sx q[2];
rz(-1.5283661) q[2];
sx q[2];
rz(0.050780642) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0446339) q[1];
sx q[1];
rz(-0.63475906) q[1];
sx q[1];
rz(-1.8083014) q[1];
rz(-pi) q[2];
rz(-2.3761446) q[3];
sx q[3];
rz(-0.42740373) q[3];
sx q[3];
rz(-0.039507341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(1.7370976) q[0];
rz(1.2069758) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.5054024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0597499) q[0];
sx q[0];
rz(-1.8847701) q[0];
sx q[0];
rz(0.44072515) q[0];
x q[1];
rz(2.2340441) q[2];
sx q[2];
rz(-2.2758256) q[2];
sx q[2];
rz(2.7570587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7561188) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(-2.4454152) q[1];
x q[2];
rz(-2.1721341) q[3];
sx q[3];
rz(-1.8441895) q[3];
sx q[3];
rz(2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-0.62301821) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4432916) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(2.495893) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3919968) q[0];
sx q[0];
rz(-2.2889334) q[0];
sx q[0];
rz(-2.0977661) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5876797) q[2];
sx q[2];
rz(-1.950112) q[2];
sx q[2];
rz(-1.5253138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7281108) q[1];
sx q[1];
rz(-2.9915447) q[1];
sx q[1];
rz(-0.24197443) q[1];
rz(0.42616578) q[3];
sx q[3];
rz(-0.45939547) q[3];
sx q[3];
rz(-0.029749425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50679961) q[2];
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
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.6181347) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(-0.24542228) q[0];
rz(2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(-1.4642749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3905555) q[0];
sx q[0];
rz(-1.3532191) q[0];
sx q[0];
rz(1.4575973) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1572838) q[2];
sx q[2];
rz(-1.3830292) q[2];
sx q[2];
rz(0.078660065) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.163584) q[1];
sx q[1];
rz(-2.2469236) q[1];
sx q[1];
rz(-2.7917557) q[1];
rz(1.7467935) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(-1.2916816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-2.0132813) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.2561692) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(0.71184288) q[2];
sx q[2];
rz(-1.9528452) q[2];
sx q[2];
rz(-0.14609329) q[2];
rz(-1.1938548) q[3];
sx q[3];
rz(-1.0300763) q[3];
sx q[3];
rz(-1.1618617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

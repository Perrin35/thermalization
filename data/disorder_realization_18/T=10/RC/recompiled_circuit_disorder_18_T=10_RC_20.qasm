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
rz(5.3263721) q[0];
sx q[0];
rz(10.858067) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8960436) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(1.8401237) q[0];
rz(-1.990591) q[2];
sx q[2];
rz(-1.3090054) q[2];
sx q[2];
rz(-2.7671438) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0393385) q[1];
sx q[1];
rz(-1.551911) q[1];
sx q[1];
rz(-1.1182055) q[1];
rz(2.5979795) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(-1.4527814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(-1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-0.84709644) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4439508) q[0];
sx q[0];
rz(-1.9662494) q[0];
sx q[0];
rz(-2.6614499) q[0];
rz(-pi) q[1];
rz(0.77913021) q[2];
sx q[2];
rz(-1.3574294) q[2];
sx q[2];
rz(-1.9386292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6716799) q[1];
sx q[1];
rz(-2.5922853) q[1];
sx q[1];
rz(0.13110199) q[1];
rz(-pi) q[2];
rz(-2.0081677) q[3];
sx q[3];
rz(-1.2542033) q[3];
sx q[3];
rz(2.7377627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.4891362) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(-0.96673036) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2741094) q[0];
sx q[0];
rz(-1.2615146) q[0];
sx q[0];
rz(1.5779737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40076077) q[2];
sx q[2];
rz(-2.4279865) q[2];
sx q[2];
rz(2.1378627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0264764) q[1];
sx q[1];
rz(-0.83362245) q[1];
sx q[1];
rz(1.923418) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4089291) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(2.8837567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30806914) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(-0.72675881) q[2];
rz(2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
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
rz(-0.90958909) q[0];
sx q[0];
rz(2.3040422) q[0];
rz(-0.063886558) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(-2.137616) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5570045) q[1];
sx q[1];
rz(-1.6465721) q[1];
sx q[1];
rz(1.9083379) q[1];
x q[2];
rz(1.1816979) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(2.1967595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(-0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-0.61757913) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546346) q[0];
sx q[0];
rz(-1.600453) q[0];
sx q[0];
rz(1.6164854) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5768443) q[2];
sx q[2];
rz(-0.72482938) q[2];
sx q[2];
rz(0.85541475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4434112) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(1.0990259) q[1];
x q[2];
rz(-3.0693552) q[3];
sx q[3];
rz(-0.17599711) q[3];
sx q[3];
rz(2.3776059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(3.0517975) q[0];
rz(2.2604997) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(-2.9615013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0952547) q[0];
sx q[0];
rz(-1.006608) q[0];
sx q[0];
rz(-2.6020223) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.300235) q[2];
sx q[2];
rz(-0.62108835) q[2];
sx q[2];
rz(1.6294711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99299586) q[1];
sx q[1];
rz(-0.33278782) q[1];
sx q[1];
rz(0.51698835) q[1];
rz(1.0780219) q[3];
sx q[3];
rz(-1.5544976) q[3];
sx q[3];
rz(-1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(-1.1759261) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.4181597) q[0];
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
rz(1.9366074) q[0];
sx q[0];
rz(-1.6523223) q[0];
sx q[0];
rz(1.8100912) q[0];
x q[1];
rz(-1.6778498) q[2];
sx q[2];
rz(-2.7333626) q[2];
sx q[2];
rz(-1.5232616) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7525967) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(2.9700301) q[1];
x q[2];
rz(2.3761446) q[3];
sx q[3];
rz(-0.42740373) q[3];
sx q[3];
rz(-3.1020853) q[3];
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
rz(1.7447757) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.7370976) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(-1.6361902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7751986) q[0];
sx q[0];
rz(-1.1530071) q[0];
sx q[0];
rz(-1.9154857) q[0];
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
sx q[0];
rz(-pi/2) q[0];
rz(1.2124894) q[1];
sx q[1];
rz(-2.266532) q[1];
sx q[1];
rz(-1.6131669) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1721341) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(2.495893) q[0];
rz(-2.6422016) q[1];
sx q[1];
rz(-1.7130518) q[1];
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
rz(-1.1823913) q[0];
sx q[0];
rz(-0.79083058) q[0];
rz(2.0090946) q[2];
sx q[2];
rz(-1.0602789) q[2];
sx q[2];
rz(0.27062182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9727271) q[1];
sx q[1];
rz(-1.7164413) q[1];
sx q[1];
rz(-1.6070073) q[1];
rz(-2.7183652) q[3];
sx q[3];
rz(-1.755135) q[3];
sx q[3];
rz(-1.1545899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(-2.6169422) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.6181347) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(-2.8961704) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.4642749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3905555) q[0];
sx q[0];
rz(-1.3532191) q[0];
sx q[0];
rz(-1.6839954) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98430888) q[2];
sx q[2];
rz(-1.3830292) q[2];
sx q[2];
rz(-0.078660065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97800868) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-0.34983695) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0189832) q[3];
sx q[3];
rz(-0.19482329) q[3];
sx q[3];
rz(2.4217055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.8650919) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(2.0586661) q[2];
sx q[2];
rz(-0.91960533) q[2];
sx q[2];
rz(-1.4056924) q[2];
rz(2.5682156) q[3];
sx q[3];
rz(-1.2497414) q[3];
sx q[3];
rz(-2.9336815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.618453) q[0];
sx q[0];
rz(-1.3292399) q[0];
sx q[0];
rz(0.068058204) q[0];
rz(-pi) q[1];
rz(2.8561864) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(2.0602496) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5073587) q[1];
sx q[1];
rz(-2.6886352) q[1];
sx q[1];
rz(1.5276315) q[1];
rz(-0.47031109) q[3];
sx q[3];
rz(-2.5458126) q[3];
sx q[3];
rz(0.51629984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(1.2343181) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-2.3117075) q[0];
rz(-2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-0.84709644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9053099) q[0];
sx q[0];
rz(-2.5295527) q[0];
sx q[0];
rz(2.406714) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8424938) q[2];
sx q[2];
rz(-0.80183376) q[2];
sx q[2];
rz(0.15660827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3165247) q[1];
sx q[1];
rz(-2.1148588) q[1];
sx q[1];
rz(-1.6506509) q[1];
rz(-2.794572) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(1.3115209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-0.57139325) q[0];
rz(0.96673036) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29887154) q[0];
sx q[0];
rz(-1.5776331) q[0];
sx q[0];
rz(0.30928916) q[0];
rz(-pi) q[1];
rz(-0.40076077) q[2];
sx q[2];
rz(-0.71360613) q[2];
sx q[2];
rz(-2.1378627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1151162) q[1];
sx q[1];
rz(-0.83362245) q[1];
sx q[1];
rz(-1.2181746) q[1];
x q[2];
rz(-0.35964386) q[3];
sx q[3];
rz(-2.375964) q[3];
sx q[3];
rz(-1.5639203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(-0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.5184901) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
rz(-0.8262659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5465281) q[0];
sx q[0];
rz(-1.0142769) q[0];
sx q[0];
rz(-2.3331649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70010186) q[2];
sx q[2];
rz(-1.5296017) q[2];
sx q[2];
rz(0.61566478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0127718) q[1];
sx q[1];
rz(-1.2342617) q[1];
sx q[1];
rz(0.080288447) q[1];
rz(-pi) q[2];
rz(-1.1816979) q[3];
sx q[3];
rz(-0.97550052) q[3];
sx q[3];
rz(-0.94483313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(0.91439247) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(3.1359613) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3546346) q[0];
sx q[0];
rz(-1.600453) q[0];
sx q[0];
rz(-1.6164854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56474833) q[2];
sx q[2];
rz(-0.72482938) q[2];
sx q[2];
rz(-0.85541475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2684106) q[1];
sx q[1];
rz(-2.1408484) q[1];
sx q[1];
rz(2.8084055) q[1];
rz(-pi) q[2];
rz(-0.17554749) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(2.4059084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(-2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8343617) q[0];
sx q[0];
rz(-2.0198856) q[0];
sx q[0];
rz(0.93528549) q[0];
rz(-2.1744556) q[2];
sx q[2];
rz(-1.4146311) q[2];
sx q[2];
rz(0.16317633) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99299586) q[1];
sx q[1];
rz(-0.33278782) q[1];
sx q[1];
rz(-0.51698835) q[1];
rz(-pi) q[2];
rz(-0.018499231) q[3];
sx q[3];
rz(-2.0634994) q[3];
sx q[3];
rz(2.8871356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(-1.4181597) q[0];
rz(1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-0.040239008) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0979472) q[0];
sx q[0];
rz(-2.8890434) q[0];
sx q[0];
rz(-1.9027684) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1646541) q[2];
sx q[2];
rz(-1.6132266) q[2];
sx q[2];
rz(3.090812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8602627) q[1];
sx q[1];
rz(-1.4308235) q[1];
sx q[1];
rz(0.94957385) q[1];
rz(-pi) q[2];
rz(-0.76544806) q[3];
sx q[3];
rz(-0.42740373) q[3];
sx q[3];
rz(0.039507341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(1.404495) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(-1.5054024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907785) q[0];
sx q[0];
rz(-2.6065126) q[0];
sx q[0];
rz(0.65061609) q[0];
rz(-pi) q[1];
rz(0.82377388) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.4866231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9951524) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(-3.0909096) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0309615) q[3];
sx q[3];
rz(-2.4880829) q[3];
sx q[3];
rz(-1.489952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(0.15052477) q[2];
rz(1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4432916) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(-0.64569965) q[0];
rz(0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.8766778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0255233) q[0];
sx q[0];
rz(-2.2795296) q[0];
sx q[0];
rz(-0.52225964) q[0];
x q[1];
rz(2.4931156) q[2];
sx q[2];
rz(-0.65995526) q[2];
sx q[2];
rz(2.6476268) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39667323) q[1];
sx q[1];
rz(-1.5349689) q[1];
sx q[1];
rz(2.9958535) q[1];
rz(-pi) q[2];
rz(-2.7154269) q[3];
sx q[3];
rz(-0.45939547) q[3];
sx q[3];
rz(-0.029749425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-2.6169422) q[2];
rz(-0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6181347) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(-2.8961704) q[0];
rz(-2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.4642749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7510371) q[0];
sx q[0];
rz(-1.7883736) q[0];
sx q[0];
rz(1.6839954) q[0];
rz(-0.22428959) q[2];
sx q[2];
rz(-2.1456492) q[2];
sx q[2];
rz(1.6155417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.163584) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-2.7917557) q[1];
rz(-pi) q[2];
rz(-1.3947992) q[3];
sx q[3];
rz(-1.6547852) q[3];
sx q[3];
rz(1.2916816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8993373) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
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
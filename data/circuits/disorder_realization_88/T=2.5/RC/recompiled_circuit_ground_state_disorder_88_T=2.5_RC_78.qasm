OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2849046) q[0];
sx q[0];
rz(-1.8104799) q[0];
sx q[0];
rz(2.3321505) q[0];
rz(-2.5419905) q[1];
sx q[1];
rz(-1.124958) q[1];
sx q[1];
rz(-2.6198299) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7358462) q[0];
sx q[0];
rz(-3.0015115) q[0];
sx q[0];
rz(1.5790042) q[0];
rz(-1.8667562) q[2];
sx q[2];
rz(-2.1115536) q[2];
sx q[2];
rz(-2.7014033) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45568902) q[1];
sx q[1];
rz(-1.1215116) q[1];
sx q[1];
rz(-0.30251578) q[1];
rz(-pi) q[2];
rz(-0.14484804) q[3];
sx q[3];
rz(-2.600935) q[3];
sx q[3];
rz(2.7365505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.9894422) q[2];
sx q[2];
rz(-0.70789727) q[2];
sx q[2];
rz(1.36261) q[2];
rz(-1.6705492) q[3];
sx q[3];
rz(-1.5444642) q[3];
sx q[3];
rz(-2.7267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7489557) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(-1.0697399) q[0];
rz(-2.3906129) q[1];
sx q[1];
rz(-0.90194482) q[1];
sx q[1];
rz(2.353207) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8892944) q[0];
sx q[0];
rz(-1.3836593) q[0];
sx q[0];
rz(-0.62690134) q[0];
rz(2.4761202) q[2];
sx q[2];
rz(-1.5147527) q[2];
sx q[2];
rz(1.9201345) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.40882817) q[1];
sx q[1];
rz(-1.0548548) q[1];
sx q[1];
rz(0.22294362) q[1];
rz(-2.6308344) q[3];
sx q[3];
rz(-1.9572957) q[3];
sx q[3];
rz(-0.0051156839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1145733) q[2];
sx q[2];
rz(-0.35832778) q[2];
sx q[2];
rz(-2.2349854) q[2];
rz(1.00057) q[3];
sx q[3];
rz(-2.0228736) q[3];
sx q[3];
rz(-2.6226543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8949316) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(1.3038127) q[0];
rz(-1.1600102) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(-1.4132285) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50158635) q[0];
sx q[0];
rz(-1.4297856) q[0];
sx q[0];
rz(-1.4835351) q[0];
rz(-pi) q[1];
x q[1];
rz(2.505439) q[2];
sx q[2];
rz(-1.9463123) q[2];
sx q[2];
rz(-2.2454303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1294413) q[1];
sx q[1];
rz(-2.0291162) q[1];
sx q[1];
rz(-1.3740416) q[1];
x q[2];
rz(-3.0158442) q[3];
sx q[3];
rz(-2.2819321) q[3];
sx q[3];
rz(-1.8497194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7580737) q[2];
sx q[2];
rz(-1.5531837) q[2];
sx q[2];
rz(0.12913945) q[2];
rz(0.50390759) q[3];
sx q[3];
rz(-1.7241071) q[3];
sx q[3];
rz(0.57607877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1103766) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(0.84621286) q[0];
rz(0.57202488) q[1];
sx q[1];
rz(-1.5049479) q[1];
sx q[1];
rz(1.2342854) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6568778) q[0];
sx q[0];
rz(-1.6678255) q[0];
sx q[0];
rz(-2.808284) q[0];
rz(-pi) q[1];
rz(2.8568966) q[2];
sx q[2];
rz(-2.0139593) q[2];
sx q[2];
rz(2.5131132) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9920791) q[1];
sx q[1];
rz(-1.3617591) q[1];
sx q[1];
rz(-1.6531709) q[1];
rz(-pi) q[2];
rz(-2.0999072) q[3];
sx q[3];
rz(-0.84741454) q[3];
sx q[3];
rz(-0.46284562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7531551) q[2];
sx q[2];
rz(-2.2074102) q[2];
sx q[2];
rz(1.4446806) q[2];
rz(1.8709024) q[3];
sx q[3];
rz(-1.5710187) q[3];
sx q[3];
rz(1.9522033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5020849) q[0];
sx q[0];
rz(-1.6725699) q[0];
sx q[0];
rz(1.0953267) q[0];
rz(3.0785839) q[1];
sx q[1];
rz(-1.6687702) q[1];
sx q[1];
rz(-0.94620401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5216833) q[0];
sx q[0];
rz(-0.54542002) q[0];
sx q[0];
rz(1.6248466) q[0];
rz(-pi) q[1];
rz(0.49503742) q[2];
sx q[2];
rz(-2.4433141) q[2];
sx q[2];
rz(-2.3586522) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9621219) q[1];
sx q[1];
rz(-2.081916) q[1];
sx q[1];
rz(-1.5764135) q[1];
x q[2];
rz(-1.0466411) q[3];
sx q[3];
rz(-2.6970409) q[3];
sx q[3];
rz(1.9421645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63518628) q[2];
sx q[2];
rz(-2.0621641) q[2];
sx q[2];
rz(1.5080473) q[2];
rz(0.15277282) q[3];
sx q[3];
rz(-1.907405) q[3];
sx q[3];
rz(0.93301409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9043115) q[0];
sx q[0];
rz(-2.5560684) q[0];
sx q[0];
rz(-0.10865077) q[0];
rz(-0.12734224) q[1];
sx q[1];
rz(-1.2510977) q[1];
sx q[1];
rz(2.2551575) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3868635) q[0];
sx q[0];
rz(-0.51966705) q[0];
sx q[0];
rz(-3.1082252) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8511925) q[2];
sx q[2];
rz(-2.9290158) q[2];
sx q[2];
rz(0.55586284) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97798367) q[1];
sx q[1];
rz(-0.22630461) q[1];
sx q[1];
rz(-1.3127021) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16951321) q[3];
sx q[3];
rz(-1.4614786) q[3];
sx q[3];
rz(1.5284571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2966557) q[2];
sx q[2];
rz(-0.79534328) q[2];
sx q[2];
rz(-2.8080158) q[2];
rz(-2.2070456) q[3];
sx q[3];
rz(-0.75342527) q[3];
sx q[3];
rz(2.2540588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6880671) q[0];
sx q[0];
rz(-1.6418566) q[0];
sx q[0];
rz(2.8330579) q[0];
rz(-2.535179) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(2.6920614) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9251557) q[0];
sx q[0];
rz(-1.0569658) q[0];
sx q[0];
rz(-1.4975182) q[0];
x q[1];
rz(0.4662029) q[2];
sx q[2];
rz(-0.84377366) q[2];
sx q[2];
rz(-2.1920993) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9924589) q[1];
sx q[1];
rz(-0.95548742) q[1];
sx q[1];
rz(0.12507089) q[1];
rz(-pi) q[2];
rz(0.87089296) q[3];
sx q[3];
rz(-2.0051536) q[3];
sx q[3];
rz(0.30444452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0799847) q[2];
sx q[2];
rz(-1.4915165) q[2];
sx q[2];
rz(1.9645346) q[2];
rz(-2.4401149) q[3];
sx q[3];
rz(-2.8886815) q[3];
sx q[3];
rz(-2.285932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6845282) q[0];
sx q[0];
rz(-2.6236911) q[0];
sx q[0];
rz(-1.49217) q[0];
rz(1.0555142) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(0.42748705) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1558091) q[0];
sx q[0];
rz(-1.0571348) q[0];
sx q[0];
rz(-2.2751121) q[0];
rz(-pi) q[1];
rz(-0.34025365) q[2];
sx q[2];
rz(-0.45368735) q[2];
sx q[2];
rz(2.2668374) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44459773) q[1];
sx q[1];
rz(-2.8279404) q[1];
sx q[1];
rz(-0.072254953) q[1];
rz(-pi) q[2];
rz(1.8425138) q[3];
sx q[3];
rz(-1.875056) q[3];
sx q[3];
rz(0.4189241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9588354) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(-1.9367564) q[2];
rz(3.0086573) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.28250113) q[0];
sx q[0];
rz(-1.7318672) q[0];
sx q[0];
rz(-2.6891563) q[0];
rz(2.2826507) q[1];
sx q[1];
rz(-2.5622538) q[1];
sx q[1];
rz(1.6890242) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4389324) q[0];
sx q[0];
rz(-2.0468674) q[0];
sx q[0];
rz(-0.66113801) q[0];
rz(2.1253517) q[2];
sx q[2];
rz(-1.8565389) q[2];
sx q[2];
rz(0.72140297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0124546) q[1];
sx q[1];
rz(-0.5059537) q[1];
sx q[1];
rz(1.5652191) q[1];
x q[2];
rz(0.36859103) q[3];
sx q[3];
rz(-1.907371) q[3];
sx q[3];
rz(0.14813885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9026044) q[2];
sx q[2];
rz(-1.3024104) q[2];
sx q[2];
rz(-2.7461309) q[2];
rz(1.0836733) q[3];
sx q[3];
rz(-2.5189221) q[3];
sx q[3];
rz(-1.4130939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85635066) q[0];
sx q[0];
rz(-3.0265891) q[0];
sx q[0];
rz(-1.2913936) q[0];
rz(2.3639288) q[1];
sx q[1];
rz(-1.5592557) q[1];
sx q[1];
rz(1.8596328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0031105) q[0];
sx q[0];
rz(-1.3457451) q[0];
sx q[0];
rz(-0.22329325) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8516232) q[2];
sx q[2];
rz(-2.4735138) q[2];
sx q[2];
rz(-1.5886605) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3471372) q[1];
sx q[1];
rz(-1.0442808) q[1];
sx q[1];
rz(1.305425) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3814209) q[3];
sx q[3];
rz(-2.6145589) q[3];
sx q[3];
rz(2.1552212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97734863) q[2];
sx q[2];
rz(-2.7898596) q[2];
sx q[2];
rz(2.742761) q[2];
rz(2.2073958) q[3];
sx q[3];
rz(-0.7312921) q[3];
sx q[3];
rz(0.092223316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4451404) q[0];
sx q[0];
rz(-2.2284989) q[0];
sx q[0];
rz(-3.1365119) q[0];
rz(2.483881) q[1];
sx q[1];
rz(-1.7291768) q[1];
sx q[1];
rz(-1.9531858) q[1];
rz(-2.8270375) q[2];
sx q[2];
rz(-2.3125074) q[2];
sx q[2];
rz(-0.41994647) q[2];
rz(-1.6647958) q[3];
sx q[3];
rz(-2.5537659) q[3];
sx q[3];
rz(-1.4418494) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

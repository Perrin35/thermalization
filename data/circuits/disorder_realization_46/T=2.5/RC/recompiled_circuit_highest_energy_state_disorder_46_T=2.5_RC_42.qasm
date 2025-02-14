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
rz(0.47705874) q[0];
sx q[0];
rz(-0.30204371) q[0];
sx q[0];
rz(-1.2560691) q[0];
rz(-2.7657901) q[1];
sx q[1];
rz(-1.8355651) q[1];
sx q[1];
rz(-2.0130872) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9039584) q[0];
sx q[0];
rz(-0.65017831) q[0];
sx q[0];
rz(1.0976492) q[0];
x q[1];
rz(1.2641505) q[2];
sx q[2];
rz(-2.1696343) q[2];
sx q[2];
rz(1.6055589) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0839786) q[1];
sx q[1];
rz(-2.489016) q[1];
sx q[1];
rz(0.90820163) q[1];
x q[2];
rz(0.63367356) q[3];
sx q[3];
rz(-1.2094991) q[3];
sx q[3];
rz(0.47508815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.5350554) q[2];
sx q[2];
rz(-0.29416072) q[2];
sx q[2];
rz(-1.199022) q[2];
rz(2.8626275) q[3];
sx q[3];
rz(-1.6821034) q[3];
sx q[3];
rz(-0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8798384) q[0];
sx q[0];
rz(-0.45670515) q[0];
sx q[0];
rz(2.7754011) q[0];
rz(-1.8531894) q[1];
sx q[1];
rz(-2.2346456) q[1];
sx q[1];
rz(0.2624951) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27124559) q[0];
sx q[0];
rz(-0.34929212) q[0];
sx q[0];
rz(0.55089921) q[0];
x q[1];
rz(1.8196202) q[2];
sx q[2];
rz(-1.5142528) q[2];
sx q[2];
rz(2.4849934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4481363) q[1];
sx q[1];
rz(-0.46365204) q[1];
sx q[1];
rz(0.69000375) q[1];
x q[2];
rz(1.315824) q[3];
sx q[3];
rz(-1.5951968) q[3];
sx q[3];
rz(0.57282626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3890106) q[2];
sx q[2];
rz(-1.2685308) q[2];
sx q[2];
rz(2.9627723) q[2];
rz(-3.1255152) q[3];
sx q[3];
rz(-0.5355081) q[3];
sx q[3];
rz(0.9308365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7054684) q[0];
sx q[0];
rz(-0.12691623) q[0];
sx q[0];
rz(2.5674168) q[0];
rz(-2.9899959) q[1];
sx q[1];
rz(-2.6972289) q[1];
sx q[1];
rz(-1.9134329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3686217) q[0];
sx q[0];
rz(-2.6009692) q[0];
sx q[0];
rz(-1.2631046) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1636859) q[2];
sx q[2];
rz(-1.5139069) q[2];
sx q[2];
rz(-1.6757193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2252402) q[1];
sx q[1];
rz(-0.95969363) q[1];
sx q[1];
rz(3.0322187) q[1];
x q[2];
rz(-2.853573) q[3];
sx q[3];
rz(-2.0374808) q[3];
sx q[3];
rz(0.53190255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0000618) q[2];
sx q[2];
rz(-2.5059293) q[2];
sx q[2];
rz(-0.89111152) q[2];
rz(-2.7530503) q[3];
sx q[3];
rz(-1.6308234) q[3];
sx q[3];
rz(2.8019606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5163088) q[0];
sx q[0];
rz(-0.044965222) q[0];
sx q[0];
rz(0.48211023) q[0];
rz(-0.66857839) q[1];
sx q[1];
rz(-1.5336288) q[1];
sx q[1];
rz(0.63749981) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4607142) q[0];
sx q[0];
rz(-1.3377046) q[0];
sx q[0];
rz(0.97026578) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11580616) q[2];
sx q[2];
rz(-2.2409332) q[2];
sx q[2];
rz(1.1095123) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5341949) q[1];
sx q[1];
rz(-0.084102159) q[1];
sx q[1];
rz(0.41513319) q[1];
rz(-1.2145146) q[3];
sx q[3];
rz(-1.3351091) q[3];
sx q[3];
rz(1.5596149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.94347) q[2];
sx q[2];
rz(-2.1169457) q[2];
sx q[2];
rz(3.0812954) q[2];
rz(-0.30334011) q[3];
sx q[3];
rz(-0.077849418) q[3];
sx q[3];
rz(2.8764603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1497134) q[0];
sx q[0];
rz(-2.7880221) q[0];
sx q[0];
rz(-2.7283332) q[0];
rz(2.7464271) q[1];
sx q[1];
rz(-1.3576077) q[1];
sx q[1];
rz(-2.1392335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4226333) q[0];
sx q[0];
rz(-2.0074685) q[0];
sx q[0];
rz(-2.5366001) q[0];
rz(-pi) q[1];
rz(0.036040739) q[2];
sx q[2];
rz(-1.3610351) q[2];
sx q[2];
rz(-2.710611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52702466) q[1];
sx q[1];
rz(-0.90674149) q[1];
sx q[1];
rz(-0.27399691) q[1];
x q[2];
rz(1.8931562) q[3];
sx q[3];
rz(-1.3731706) q[3];
sx q[3];
rz(1.2569497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62785441) q[2];
sx q[2];
rz(-1.0971053) q[2];
sx q[2];
rz(-1.6045137) q[2];
rz(1.9419935) q[3];
sx q[3];
rz(-2.0204085) q[3];
sx q[3];
rz(-2.3698923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62726778) q[0];
sx q[0];
rz(-1.5031313) q[0];
sx q[0];
rz(2.7728873) q[0];
rz(0.74327028) q[1];
sx q[1];
rz(-2.5786046) q[1];
sx q[1];
rz(-0.35400131) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0403772) q[0];
sx q[0];
rz(-2.4987552) q[0];
sx q[0];
rz(-1.4835711) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2576302) q[2];
sx q[2];
rz(-2.0764362) q[2];
sx q[2];
rz(-2.2693116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94985574) q[1];
sx q[1];
rz(-0.60137784) q[1];
sx q[1];
rz(-1.8888425) q[1];
x q[2];
rz(2.5052222) q[3];
sx q[3];
rz(-1.2546347) q[3];
sx q[3];
rz(-0.13383257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.02504286) q[2];
sx q[2];
rz(-1.9321238) q[2];
sx q[2];
rz(0.0054736007) q[2];
rz(-0.51760751) q[3];
sx q[3];
rz(-2.7766683) q[3];
sx q[3];
rz(-3.1143809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7356877) q[0];
sx q[0];
rz(-0.55957496) q[0];
sx q[0];
rz(2.9285808) q[0];
rz(-1.4951911) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(-2.2921553) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6005858) q[0];
sx q[0];
rz(-1.6141506) q[0];
sx q[0];
rz(-0.61886532) q[0];
x q[1];
rz(-1.527571) q[2];
sx q[2];
rz(-1.5400949) q[2];
sx q[2];
rz(2.8376606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53687261) q[1];
sx q[1];
rz(-2.0402131) q[1];
sx q[1];
rz(1.0632656) q[1];
x q[2];
rz(-0.15504239) q[3];
sx q[3];
rz(-1.0477121) q[3];
sx q[3];
rz(-1.6301159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21462943) q[2];
sx q[2];
rz(-1.9081076) q[2];
sx q[2];
rz(-2.5978973) q[2];
rz(0.26116535) q[3];
sx q[3];
rz(-1.5800913) q[3];
sx q[3];
rz(-2.9847667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9549304) q[0];
sx q[0];
rz(-2.08827) q[0];
sx q[0];
rz(2.9539811) q[0];
rz(-1.1232417) q[1];
sx q[1];
rz(-2.3767327) q[1];
sx q[1];
rz(-0.16256464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.150972) q[0];
sx q[0];
rz(-2.2216068) q[0];
sx q[0];
rz(-1.1605754) q[0];
rz(-2.0948903) q[2];
sx q[2];
rz(-1.419029) q[2];
sx q[2];
rz(-0.87226112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9342207) q[1];
sx q[1];
rz(-1.6160864) q[1];
sx q[1];
rz(2.1846733) q[1];
rz(-pi) q[2];
rz(1.4054589) q[3];
sx q[3];
rz(-2.5389034) q[3];
sx q[3];
rz(1.7053982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8673914) q[2];
sx q[2];
rz(-2.826639) q[2];
sx q[2];
rz(-0.20142041) q[2];
rz(2.6650186) q[3];
sx q[3];
rz(-1.3786992) q[3];
sx q[3];
rz(2.1417638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9607214) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(2.5647105) q[0];
rz(-2.8374529) q[1];
sx q[1];
rz(-2.0174593) q[1];
sx q[1];
rz(-1.1041799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21872422) q[0];
sx q[0];
rz(-1.6593242) q[0];
sx q[0];
rz(1.9506111) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4494032) q[2];
sx q[2];
rz(-1.0052201) q[2];
sx q[2];
rz(0.34715101) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7930709) q[1];
sx q[1];
rz(-0.99731113) q[1];
sx q[1];
rz(-0.46646935) q[1];
rz(-pi) q[2];
x q[2];
rz(3.125413) q[3];
sx q[3];
rz(-1.9619298) q[3];
sx q[3];
rz(-1.1333381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96678174) q[2];
sx q[2];
rz(-0.2468144) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(-1.5934058) q[3];
sx q[3];
rz(-2.4631409) q[3];
sx q[3];
rz(2.9233176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9504647) q[0];
sx q[0];
rz(-0.8466962) q[0];
sx q[0];
rz(0.51093131) q[0];
rz(1.9966985) q[1];
sx q[1];
rz(-1.9547918) q[1];
sx q[1];
rz(-1.6330382) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.975113) q[0];
sx q[0];
rz(-0.60752019) q[0];
sx q[0];
rz(2.6954073) q[0];
rz(-pi) q[1];
rz(0.26692674) q[2];
sx q[2];
rz(-1.1014001) q[2];
sx q[2];
rz(-1.5933344) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.88860565) q[1];
sx q[1];
rz(-1.7481511) q[1];
sx q[1];
rz(-0.4528404) q[1];
x q[2];
rz(-1.3168961) q[3];
sx q[3];
rz(-1.2798556) q[3];
sx q[3];
rz(1.8559141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42023429) q[2];
sx q[2];
rz(-2.5098462) q[2];
sx q[2];
rz(-1.7101804) q[2];
rz(-3.1259649) q[3];
sx q[3];
rz(-1.2651919) q[3];
sx q[3];
rz(2.6354852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2331727) q[0];
sx q[0];
rz(-1.40042) q[0];
sx q[0];
rz(-0.9216876) q[0];
rz(-1.635101) q[1];
sx q[1];
rz(-1.2163305) q[1];
sx q[1];
rz(-0.46179927) q[1];
rz(-1.1882412) q[2];
sx q[2];
rz(-1.8391704) q[2];
sx q[2];
rz(3.0808629) q[2];
rz(-1.4931645) q[3];
sx q[3];
rz(-1.6059198) q[3];
sx q[3];
rz(-3.0618659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

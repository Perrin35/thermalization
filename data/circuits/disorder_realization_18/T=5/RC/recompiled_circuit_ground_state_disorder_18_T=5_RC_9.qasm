OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9590149) q[0];
sx q[0];
rz(-2.9380517) q[0];
sx q[0];
rz(1.2476873) q[0];
rz(1.3522476) q[1];
sx q[1];
rz(-0.69753733) q[1];
sx q[1];
rz(-0.69812671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1230948) q[0];
sx q[0];
rz(-2.4677489) q[0];
sx q[0];
rz(-0.68742623) q[0];
rz(-pi) q[1];
rz(2.6883672) q[2];
sx q[2];
rz(-1.5767323) q[2];
sx q[2];
rz(-0.56528948) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3580324) q[1];
sx q[1];
rz(-0.81746819) q[1];
sx q[1];
rz(1.1518351) q[1];
rz(-1.8258445) q[3];
sx q[3];
rz(-1.5849893) q[3];
sx q[3];
rz(3.1007781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84833604) q[2];
sx q[2];
rz(-1.9170599) q[2];
sx q[2];
rz(-2.2621034) q[2];
rz(2.4073811) q[3];
sx q[3];
rz(-1.1314393) q[3];
sx q[3];
rz(-0.82936275) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5010928) q[0];
sx q[0];
rz(-1.1859897) q[0];
sx q[0];
rz(-0.99047852) q[0];
rz(0.99539202) q[1];
sx q[1];
rz(-1.4442911) q[1];
sx q[1];
rz(-0.46517864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.746691) q[0];
sx q[0];
rz(-1.1550265) q[0];
sx q[0];
rz(-1.3711934) q[0];
x q[1];
rz(2.0700109) q[2];
sx q[2];
rz(-2.4518161) q[2];
sx q[2];
rz(0.76329714) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0439321) q[1];
sx q[1];
rz(-2.2273835) q[1];
sx q[1];
rz(0.93933479) q[1];
x q[2];
rz(0.2618913) q[3];
sx q[3];
rz(-1.7658111) q[3];
sx q[3];
rz(-3.0660925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5140932) q[2];
sx q[2];
rz(-1.2343312) q[2];
sx q[2];
rz(0.17847432) q[2];
rz(0.56525362) q[3];
sx q[3];
rz(-1.7491128) q[3];
sx q[3];
rz(-2.5843411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8239215) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(-0.70096651) q[0];
rz(-0.40755454) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(1.0754546) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0757119) q[0];
sx q[0];
rz(-1.9032818) q[0];
sx q[0];
rz(2.5677469) q[0];
rz(-1.8604061) q[2];
sx q[2];
rz(-0.67669213) q[2];
sx q[2];
rz(-1.2586359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6523931) q[1];
sx q[1];
rz(-1.4402188) q[1];
sx q[1];
rz(-1.9741535) q[1];
rz(-pi) q[2];
rz(1.386679) q[3];
sx q[3];
rz(-2.872481) q[3];
sx q[3];
rz(-1.6412013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.067523) q[2];
sx q[2];
rz(-1.3347551) q[2];
sx q[2];
rz(-1.734181) q[2];
rz(0.47809005) q[3];
sx q[3];
rz(-2.7985088) q[3];
sx q[3];
rz(0.34496719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57426977) q[0];
sx q[0];
rz(-1.747921) q[0];
sx q[0];
rz(-2.0950914) q[0];
rz(0.25431713) q[1];
sx q[1];
rz(-2.3260702) q[1];
sx q[1];
rz(0.3124803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82332084) q[0];
sx q[0];
rz(-1.7004564) q[0];
sx q[0];
rz(1.983485) q[0];
x q[1];
rz(-2.4813969) q[2];
sx q[2];
rz(-1.1822299) q[2];
sx q[2];
rz(1.6903433) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29280832) q[1];
sx q[1];
rz(-0.40221805) q[1];
sx q[1];
rz(2.4821698) q[1];
rz(-pi) q[2];
rz(1.0572079) q[3];
sx q[3];
rz(-1.2826589) q[3];
sx q[3];
rz(0.44618928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3089402) q[2];
sx q[2];
rz(-0.71844429) q[2];
sx q[2];
rz(0.18860513) q[2];
rz(0.23379937) q[3];
sx q[3];
rz(-0.89582396) q[3];
sx q[3];
rz(0.58376592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4466062) q[0];
sx q[0];
rz(-1.8291031) q[0];
sx q[0];
rz(0.0083228668) q[0];
rz(-2.7024929) q[1];
sx q[1];
rz(-1.3689901) q[1];
sx q[1];
rz(-1.2459374) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9326235) q[0];
sx q[0];
rz(-2.0216938) q[0];
sx q[0];
rz(1.3572925) q[0];
x q[1];
rz(-1.0143004) q[2];
sx q[2];
rz(-2.864553) q[2];
sx q[2];
rz(3.1173681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26995319) q[1];
sx q[1];
rz(-1.9325745) q[1];
sx q[1];
rz(2.8025318) q[1];
x q[2];
rz(-1.9381136) q[3];
sx q[3];
rz(-2.4532336) q[3];
sx q[3];
rz(1.7662774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40450725) q[2];
sx q[2];
rz(-3.0948907) q[2];
sx q[2];
rz(1.5076293) q[2];
rz(-2.5493933) q[3];
sx q[3];
rz(-1.3785572) q[3];
sx q[3];
rz(-0.73970214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1208039) q[0];
sx q[0];
rz(-1.4494267) q[0];
sx q[0];
rz(-2.8756496) q[0];
rz(-1.2159011) q[1];
sx q[1];
rz(-2.3561056) q[1];
sx q[1];
rz(1.1908092) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5191089) q[0];
sx q[0];
rz(-1.5734696) q[0];
sx q[0];
rz(0.83774211) q[0];
rz(-pi) q[1];
rz(3.1053084) q[2];
sx q[2];
rz(-2.1232833) q[2];
sx q[2];
rz(2.3123669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3336181) q[1];
sx q[1];
rz(-2.3876973) q[1];
sx q[1];
rz(0.011276007) q[1];
x q[2];
rz(1.5134638) q[3];
sx q[3];
rz(-1.2964037) q[3];
sx q[3];
rz(-2.969034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33209458) q[2];
sx q[2];
rz(-1.7569434) q[2];
sx q[2];
rz(-0.14796999) q[2];
rz(-2.68908) q[3];
sx q[3];
rz(-1.0871525) q[3];
sx q[3];
rz(0.40542671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0026523503) q[0];
sx q[0];
rz(-1.8754706) q[0];
sx q[0];
rz(-1.3462322) q[0];
rz(3.1345308) q[1];
sx q[1];
rz(-2.4738753) q[1];
sx q[1];
rz(-0.5074358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9503896) q[0];
sx q[0];
rz(-2.1090713) q[0];
sx q[0];
rz(-2.5274168) q[0];
rz(-pi) q[1];
rz(1.1603786) q[2];
sx q[2];
rz(-0.37330356) q[2];
sx q[2];
rz(-2.6609535) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9964357) q[1];
sx q[1];
rz(-1.5801589) q[1];
sx q[1];
rz(0.3939751) q[1];
rz(0.0094311992) q[3];
sx q[3];
rz(-1.8430897) q[3];
sx q[3];
rz(1.1444397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4355882) q[2];
sx q[2];
rz(-2.0861349) q[2];
sx q[2];
rz(2.9662507) q[2];
rz(0.38241479) q[3];
sx q[3];
rz(-1.9491111) q[3];
sx q[3];
rz(-0.46149883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14313702) q[0];
sx q[0];
rz(-1.0112421) q[0];
sx q[0];
rz(-1.829041) q[0];
rz(-3.0328499) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(-2.463602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10763845) q[0];
sx q[0];
rz(-2.0642515) q[0];
sx q[0];
rz(-3.1163868) q[0];
rz(2.3925875) q[2];
sx q[2];
rz(-1.555948) q[2];
sx q[2];
rz(-0.9204364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1428263) q[1];
sx q[1];
rz(-1.8590392) q[1];
sx q[1];
rz(1.0267797) q[1];
rz(-1.0394215) q[3];
sx q[3];
rz(-2.2153004) q[3];
sx q[3];
rz(2.927305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8678191) q[2];
sx q[2];
rz(-1.5742233) q[2];
sx q[2];
rz(-1.9847974) q[2];
rz(-0.50446883) q[3];
sx q[3];
rz(-1.0320832) q[3];
sx q[3];
rz(0.57197905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93541637) q[0];
sx q[0];
rz(-1.6125866) q[0];
sx q[0];
rz(-0.61844283) q[0];
rz(-0.84856021) q[1];
sx q[1];
rz(-0.51785523) q[1];
sx q[1];
rz(2.3458164) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.764641) q[0];
sx q[0];
rz(-1.5949773) q[0];
sx q[0];
rz(1.1005681) q[0];
x q[1];
rz(0.047330476) q[2];
sx q[2];
rz(-2.3162957) q[2];
sx q[2];
rz(-1.312544) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2256945) q[1];
sx q[1];
rz(-1.2978329) q[1];
sx q[1];
rz(1.8986439) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9913748) q[3];
sx q[3];
rz(-1.4338067) q[3];
sx q[3];
rz(1.610422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15647469) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(1.1953243) q[2];
rz(-0.34504238) q[3];
sx q[3];
rz(-2.2796567) q[3];
sx q[3];
rz(2.2752458) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4320375) q[0];
sx q[0];
rz(-2.9032752) q[0];
sx q[0];
rz(0.077614345) q[0];
rz(1.9084825) q[1];
sx q[1];
rz(-2.1014919) q[1];
sx q[1];
rz(2.3495823) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.663765) q[0];
sx q[0];
rz(-2.223513) q[0];
sx q[0];
rz(0.86413149) q[0];
rz(-1.7260688) q[2];
sx q[2];
rz(-1.4129935) q[2];
sx q[2];
rz(-1.4580245) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7866871) q[1];
sx q[1];
rz(-2.6514158) q[1];
sx q[1];
rz(-1.9016983) q[1];
rz(-pi) q[2];
x q[2];
rz(3.023572) q[3];
sx q[3];
rz(-2.3343918) q[3];
sx q[3];
rz(-2.0995875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21976694) q[2];
sx q[2];
rz(-0.91138387) q[2];
sx q[2];
rz(1.0731953) q[2];
rz(0.24751599) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(3.1246576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3794004) q[0];
sx q[0];
rz(-1.8397377) q[0];
sx q[0];
rz(2.4158438) q[0];
rz(2.2629867) q[1];
sx q[1];
rz(-0.85522063) q[1];
sx q[1];
rz(-1.9488889) q[1];
rz(1.0868308) q[2];
sx q[2];
rz(-1.7499583) q[2];
sx q[2];
rz(-1.3365895) q[2];
rz(2.2976919) q[3];
sx q[3];
rz(-1.5283199) q[3];
sx q[3];
rz(0.98315317) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

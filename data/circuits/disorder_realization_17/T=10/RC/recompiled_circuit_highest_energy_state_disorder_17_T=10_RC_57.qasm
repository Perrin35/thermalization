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
rz(1.849527) q[0];
sx q[0];
rz(5.4592291) q[0];
sx q[0];
rz(8.4973314) q[0];
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(-2.3271022) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8626635) q[0];
sx q[0];
rz(-1.9534847) q[0];
sx q[0];
rz(2.0576059) q[0];
rz(-pi) q[1];
rz(-2.1998946) q[2];
sx q[2];
rz(-0.88028958) q[2];
sx q[2];
rz(2.5494573) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42871745) q[1];
sx q[1];
rz(-0.71850417) q[1];
sx q[1];
rz(-2.3796622) q[1];
rz(-pi) q[2];
rz(0.82858927) q[3];
sx q[3];
rz(-2.0999206) q[3];
sx q[3];
rz(-0.55381004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77259511) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(0.43130809) q[2];
rz(-1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5876708) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(-1.824463) q[0];
rz(-2.092284) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(2.5228693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71402446) q[0];
sx q[0];
rz(-0.62494266) q[0];
sx q[0];
rz(-2.8929536) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87845408) q[2];
sx q[2];
rz(-1.2515504) q[2];
sx q[2];
rz(2.648166) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8730388) q[1];
sx q[1];
rz(-1.7015) q[1];
sx q[1];
rz(-1.9221582) q[1];
rz(1.9369164) q[3];
sx q[3];
rz(-2.0799326) q[3];
sx q[3];
rz(-0.81936554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72545663) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(-0.49187342) q[2];
rz(1.5791996) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(-2.5621342) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014515011) q[0];
sx q[0];
rz(-1.0579695) q[0];
sx q[0];
rz(2.8614817) q[0];
rz(3.0666871) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(1.3704971) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2477596) q[0];
sx q[0];
rz(-0.77447184) q[0];
sx q[0];
rz(1.2694556) q[0];
x q[1];
rz(-2.2417775) q[2];
sx q[2];
rz(-1.0428937) q[2];
sx q[2];
rz(-1.6646651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0947511) q[1];
sx q[1];
rz(-1.0221204) q[1];
sx q[1];
rz(1.8803174) q[1];
rz(0.13223572) q[3];
sx q[3];
rz(-1.5521129) q[3];
sx q[3];
rz(-1.4052237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8568153) q[2];
sx q[2];
rz(-1.6723526) q[2];
sx q[2];
rz(0.80186296) q[2];
rz(-0.91184584) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(-2.1626933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80482471) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(-2.5208852) q[0];
rz(-2.8630818) q[1];
sx q[1];
rz(-1.4418437) q[1];
sx q[1];
rz(-1.0034358) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28351706) q[0];
sx q[0];
rz(-1.2293613) q[0];
sx q[0];
rz(0.33591875) q[0];
x q[1];
rz(-1.150028) q[2];
sx q[2];
rz(-1.6347957) q[2];
sx q[2];
rz(2.6902932) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0520893) q[1];
sx q[1];
rz(-1.516253) q[1];
sx q[1];
rz(-2.8724345) q[1];
rz(-0.59856071) q[3];
sx q[3];
rz(-1.374647) q[3];
sx q[3];
rz(-0.3071827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3044546) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(3.0755074) q[2];
rz(-1.1262013) q[3];
sx q[3];
rz(-0.76032138) q[3];
sx q[3];
rz(2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5758301) q[0];
sx q[0];
rz(-3.059721) q[0];
sx q[0];
rz(-0.78731147) q[0];
rz(-0.85834223) q[1];
sx q[1];
rz(-0.97802496) q[1];
sx q[1];
rz(-2.0599005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72321945) q[0];
sx q[0];
rz(-2.0670509) q[0];
sx q[0];
rz(-2.4346274) q[0];
rz(1.101564) q[2];
sx q[2];
rz(-1.1936371) q[2];
sx q[2];
rz(-0.18845972) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.25095233) q[1];
sx q[1];
rz(-1.2861274) q[1];
sx q[1];
rz(-1.5642868) q[1];
rz(-pi) q[2];
rz(2.5852935) q[3];
sx q[3];
rz(-0.93141684) q[3];
sx q[3];
rz(-2.104708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16296884) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(-2.2990189) q[2];
rz(-2.8652371) q[3];
sx q[3];
rz(-2.1601951) q[3];
sx q[3];
rz(-1.3493376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2862947) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(1.1899765) q[0];
rz(-1.9189934) q[1];
sx q[1];
rz(-1.2346376) q[1];
sx q[1];
rz(-0.58293265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5258497) q[0];
sx q[0];
rz(-1.9186363) q[0];
sx q[0];
rz(2.6418532) q[0];
x q[1];
rz(-0.99722476) q[2];
sx q[2];
rz(-1.5590481) q[2];
sx q[2];
rz(-0.16596101) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8376298) q[1];
sx q[1];
rz(-2.8899341) q[1];
sx q[1];
rz(2.8117517) q[1];
x q[2];
rz(-0.56598466) q[3];
sx q[3];
rz(-2.4275587) q[3];
sx q[3];
rz(0.66431724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0588093) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(2.1018551) q[2];
rz(-2.5522363) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(-2.1350433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37650087) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(1.804922) q[0];
rz(-0.22557766) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(-2.371686) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25882803) q[0];
sx q[0];
rz(-2.6443091) q[0];
sx q[0];
rz(-2.8112341) q[0];
rz(2.5134447) q[2];
sx q[2];
rz(-2.2100115) q[2];
sx q[2];
rz(2.1232427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0437255) q[1];
sx q[1];
rz(-1.3732404) q[1];
sx q[1];
rz(1.6476589) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60524551) q[3];
sx q[3];
rz(-1.6351994) q[3];
sx q[3];
rz(1.4726255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.7485488) q[2];
sx q[2];
rz(2.1262271) q[2];
rz(2.9376302) q[3];
sx q[3];
rz(-1.5488307) q[3];
sx q[3];
rz(1.3913733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400951) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(-0.39328662) q[0];
rz(0.45010629) q[1];
sx q[1];
rz(-1.422912) q[1];
sx q[1];
rz(-0.030287655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6323551) q[0];
sx q[0];
rz(-1.8182753) q[0];
sx q[0];
rz(1.6028274) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.302839) q[2];
sx q[2];
rz(-0.7795802) q[2];
sx q[2];
rz(1.5991581) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62822484) q[1];
sx q[1];
rz(-1.0252044) q[1];
sx q[1];
rz(-1.5437158) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8775086) q[3];
sx q[3];
rz(-1.698016) q[3];
sx q[3];
rz(2.7252687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0344737) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(-2.0195473) q[2];
rz(-1.0551039) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(1.25157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17290641) q[0];
sx q[0];
rz(-2.1538669) q[0];
sx q[0];
rz(-2.3671142) q[0];
rz(0.6340181) q[1];
sx q[1];
rz(-1.0301544) q[1];
sx q[1];
rz(-1.0672306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241079) q[0];
sx q[0];
rz(-1.5536246) q[0];
sx q[0];
rz(-2.4033258) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7238317) q[2];
sx q[2];
rz(-0.31073292) q[2];
sx q[2];
rz(-2.8513214) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6065403) q[1];
sx q[1];
rz(-0.4036223) q[1];
sx q[1];
rz(-2.2500747) q[1];
rz(-pi) q[2];
rz(-2.0499633) q[3];
sx q[3];
rz(-1.7288343) q[3];
sx q[3];
rz(0.8534067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5344703) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(1.2657451) q[2];
rz(-0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-0.039073959) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48113111) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(1.2598502) q[0];
rz(1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(0.0025509603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1387864) q[0];
sx q[0];
rz(-1.4454137) q[0];
sx q[0];
rz(-0.12055293) q[0];
rz(1.506729) q[2];
sx q[2];
rz(-1.5612215) q[2];
sx q[2];
rz(-2.9970882) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3348744) q[1];
sx q[1];
rz(-0.93466991) q[1];
sx q[1];
rz(-1.0735372) q[1];
x q[2];
rz(-0.02648377) q[3];
sx q[3];
rz(-2.0051574) q[3];
sx q[3];
rz(-2.6514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.586414) q[2];
sx q[2];
rz(-1.5191583) q[2];
sx q[2];
rz(-2.4671538) q[2];
rz(2.0241578) q[3];
sx q[3];
rz(-1.0267461) q[3];
sx q[3];
rz(0.31291541) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085070327) q[0];
sx q[0];
rz(-1.1310348) q[0];
sx q[0];
rz(-0.52666589) q[0];
rz(0.3479192) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(2.6013026) q[2];
sx q[2];
rz(-2.7439596) q[2];
sx q[2];
rz(-0.31184218) q[2];
rz(1.4244637) q[3];
sx q[3];
rz(-1.5884593) q[3];
sx q[3];
rz(-2.9149169) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

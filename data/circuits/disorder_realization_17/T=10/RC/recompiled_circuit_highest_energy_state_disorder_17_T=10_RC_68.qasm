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
rz(-0.82395616) q[0];
sx q[0];
rz(-0.92744654) q[0];
rz(-4.1915512) q[1];
sx q[1];
rz(-0.97133049) q[1];
sx q[1];
rz(8.6102875) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2353219) q[0];
sx q[0];
rz(-2.5320279) q[0];
sx q[0];
rz(0.8602575) q[0];
rz(-pi) q[1];
rz(-0.94169803) q[2];
sx q[2];
rz(-0.88028958) q[2];
sx q[2];
rz(-2.5494573) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6223917) q[1];
sx q[1];
rz(-2.0425052) q[1];
sx q[1];
rz(2.5775108) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4708704) q[3];
sx q[3];
rz(-2.1938087) q[3];
sx q[3];
rz(0.58341208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77259511) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(0.43130809) q[2];
rz(1.4946233) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(-0.78540426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876708) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(-1.3171296) q[0];
rz(2.092284) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(-2.5228693) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241708) q[0];
sx q[0];
rz(-0.96786495) q[0];
sx q[0];
rz(-1.746491) q[0];
x q[1];
rz(-2.7359782) q[2];
sx q[2];
rz(-2.2219293) q[2];
sx q[2];
rz(1.8096015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2545027) q[1];
sx q[1];
rz(-1.9190333) q[1];
sx q[1];
rz(-3.00249) q[1];
rz(-pi) q[2];
rz(-1.9369164) q[3];
sx q[3];
rz(-2.0799326) q[3];
sx q[3];
rz(0.81936554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72545663) q[2];
sx q[2];
rz(-1.340103) q[2];
sx q[2];
rz(0.49187342) q[2];
rz(-1.5791996) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(-2.5621342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(0.28011093) q[0];
rz(-0.07490553) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(1.3704971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89383306) q[0];
sx q[0];
rz(-2.3671208) q[0];
sx q[0];
rz(1.2694556) q[0];
rz(0.63997322) q[2];
sx q[2];
rz(-1.0037862) q[2];
sx q[2];
rz(0.28653539) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.83082) q[1];
sx q[1];
rz(-1.3078863) q[1];
sx q[1];
rz(-2.5710158) q[1];
rz(-1.5519484) q[3];
sx q[3];
rz(-1.4385838) q[3];
sx q[3];
rz(2.978505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8568153) q[2];
sx q[2];
rz(-1.6723526) q[2];
sx q[2];
rz(2.3397297) q[2];
rz(-2.2297468) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(-0.97889939) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80482471) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(2.5208852) q[0];
rz(0.2785109) q[1];
sx q[1];
rz(-1.4418437) q[1];
sx q[1];
rz(2.1381569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0517563) q[0];
sx q[0];
rz(-0.47430719) q[0];
sx q[0];
rz(2.3186705) q[0];
rz(-pi) q[1];
rz(-1.9915646) q[2];
sx q[2];
rz(-1.6347957) q[2];
sx q[2];
rz(0.45129946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0520893) q[1];
sx q[1];
rz(-1.6253396) q[1];
sx q[1];
rz(-0.2691582) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3347606) q[3];
sx q[3];
rz(-2.1563362) q[3];
sx q[3];
rz(-1.3957617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8371381) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(1.1262013) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(-0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758301) q[0];
sx q[0];
rz(-3.059721) q[0];
sx q[0];
rz(-2.3542812) q[0];
rz(-2.2832504) q[1];
sx q[1];
rz(-0.97802496) q[1];
sx q[1];
rz(2.0599005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183732) q[0];
sx q[0];
rz(-1.0745418) q[0];
sx q[0];
rz(2.4346274) q[0];
x q[1];
rz(-2.7236347) q[2];
sx q[2];
rz(-2.0047028) q[2];
sx q[2];
rz(1.5746631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25095233) q[1];
sx q[1];
rz(-1.2861274) q[1];
sx q[1];
rz(1.5642868) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85162195) q[3];
sx q[3];
rz(-1.1332261) q[3];
sx q[3];
rz(2.9629666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9786238) q[2];
sx q[2];
rz(-1.9917515) q[2];
sx q[2];
rz(2.2990189) q[2];
rz(-2.8652371) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(1.3493376) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.855298) q[0];
sx q[0];
rz(-1.3110302) q[0];
sx q[0];
rz(-1.1899765) q[0];
rz(1.2225993) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(-2.55866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.603072) q[0];
sx q[0];
rz(-2.5412509) q[0];
sx q[0];
rz(2.4938514) q[0];
rz(-0.013986258) q[2];
sx q[2];
rz(-2.1443233) q[2];
sx q[2];
rz(-1.397246) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8376298) q[1];
sx q[1];
rz(-2.8899341) q[1];
sx q[1];
rz(-2.8117517) q[1];
x q[2];
rz(0.56598466) q[3];
sx q[3];
rz(-2.4275587) q[3];
sx q[3];
rz(-0.66431724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0588093) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(-1.0397376) q[2];
rz(2.5522363) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(2.1350433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650918) q[0];
sx q[0];
rz(-0.14987513) q[0];
sx q[0];
rz(-1.804922) q[0];
rz(0.22557766) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(-0.76990661) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63084376) q[0];
sx q[0];
rz(-1.1026369) q[0];
sx q[0];
rz(-1.3965142) q[0];
x q[1];
rz(-0.82775292) q[2];
sx q[2];
rz(-2.0619287) q[2];
sx q[2];
rz(-2.998005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.097867171) q[1];
sx q[1];
rz(-1.7683523) q[1];
sx q[1];
rz(-1.6476589) q[1];
rz(-0.60524551) q[3];
sx q[3];
rz(-1.5063933) q[3];
sx q[3];
rz(1.6689672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.7485488) q[2];
sx q[2];
rz(2.1262271) q[2];
rz(0.20396248) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(1.3913733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90149752) q[0];
sx q[0];
rz(-2.5101341) q[0];
sx q[0];
rz(-0.39328662) q[0];
rz(-2.6914864) q[1];
sx q[1];
rz(-1.422912) q[1];
sx q[1];
rz(3.111305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6323551) q[0];
sx q[0];
rz(-1.3233174) q[0];
sx q[0];
rz(-1.6028274) q[0];
rz(2.3322163) q[2];
sx q[2];
rz(-1.3835819) q[2];
sx q[2];
rz(-0.16448122) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68037625) q[1];
sx q[1];
rz(-2.595397) q[1];
sx q[1];
rz(-0.044574634) q[1];
rz(1.7025331) q[3];
sx q[3];
rz(-1.8326958) q[3];
sx q[3];
rz(-1.1887664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0344737) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(-2.0195473) q[2];
rz(-2.0864887) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17290641) q[0];
sx q[0];
rz(-0.98772573) q[0];
sx q[0];
rz(-2.3671142) q[0];
rz(-2.5075746) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(1.0672306) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7726583) q[0];
sx q[0];
rz(-0.83266363) q[0];
sx q[0];
rz(-1.5475818) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.417761) q[2];
sx q[2];
rz(-0.31073292) q[2];
sx q[2];
rz(0.29027127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18548705) q[1];
sx q[1];
rz(-1.2602578) q[1];
sx q[1];
rz(-0.26212543) q[1];
x q[2];
rz(1.9036071) q[3];
sx q[3];
rz(-0.5026256) q[3];
sx q[3];
rz(-0.42325936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5344703) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(1.8758476) q[2];
rz(-0.040146116) q[3];
sx q[3];
rz(-2.7669192) q[3];
sx q[3];
rz(0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6604615) q[0];
sx q[0];
rz(-0.88812319) q[0];
sx q[0];
rz(-1.2598502) q[0];
rz(1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(-3.1390417) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77230763) q[0];
sx q[0];
rz(-2.9678759) q[0];
sx q[0];
rz(-0.80887167) q[0];
rz(-1.719252) q[2];
sx q[2];
rz(-0.064777834) q[2];
sx q[2];
rz(-1.5744408) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3348744) q[1];
sx q[1];
rz(-2.2069227) q[1];
sx q[1];
rz(1.0735372) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1363013) q[3];
sx q[3];
rz(-1.5467724) q[3];
sx q[3];
rz(2.0498118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5551787) q[2];
sx q[2];
rz(-1.5191583) q[2];
sx q[2];
rz(-2.4671538) q[2];
rz(2.0241578) q[3];
sx q[3];
rz(-1.0267461) q[3];
sx q[3];
rz(-2.8286772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085070327) q[0];
sx q[0];
rz(-1.1310348) q[0];
sx q[0];
rz(-0.52666589) q[0];
rz(2.7936735) q[1];
sx q[1];
rz(-1.2384474) q[1];
sx q[1];
rz(1.7927982) q[1];
rz(-1.7835708) q[2];
sx q[2];
rz(-1.2322896) q[2];
sx q[2];
rz(-0.88862669) q[2];
rz(3.1237389) q[3];
sx q[3];
rz(-1.4244867) q[3];
sx q[3];
rz(-1.3467237) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

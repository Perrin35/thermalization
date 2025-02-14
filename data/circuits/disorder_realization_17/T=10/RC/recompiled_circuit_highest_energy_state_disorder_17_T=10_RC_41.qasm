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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0448582) q[0];
sx q[0];
rz(-1.1219026) q[0];
sx q[0];
rz(2.7142224) q[0];
x q[1];
rz(0.94169803) q[2];
sx q[2];
rz(-2.2613031) q[2];
sx q[2];
rz(0.59213537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3315288) q[1];
sx q[1];
rz(-1.074407) q[1];
sx q[1];
rz(-1.0277102) q[1];
rz(-pi) q[2];
rz(0.67072224) q[3];
sx q[3];
rz(-2.1938087) q[3];
sx q[3];
rz(-2.5581806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77259511) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(-0.43130809) q[2];
rz(1.4946233) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876708) q[0];
sx q[0];
rz(-0.19667721) q[0];
sx q[0];
rz(1.3171296) q[0];
rz(1.0493086) q[1];
sx q[1];
rz(-2.6056555) q[1];
sx q[1];
rz(-2.5228693) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174219) q[0];
sx q[0];
rz(-0.96786495) q[0];
sx q[0];
rz(1.746491) q[0];
rz(-pi) q[1];
rz(-0.87845408) q[2];
sx q[2];
rz(-1.2515504) q[2];
sx q[2];
rz(-2.648166) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8730388) q[1];
sx q[1];
rz(-1.4400926) q[1];
sx q[1];
rz(-1.9221582) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6027565) q[3];
sx q[3];
rz(-1.2528786) q[3];
sx q[3];
rz(-0.93618083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72545663) q[2];
sx q[2];
rz(-1.340103) q[2];
sx q[2];
rz(2.6497192) q[2];
rz(-1.5791996) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(-2.5621342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(-0.28011093) q[0];
rz(-3.0666871) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(-1.3704971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89383306) q[0];
sx q[0];
rz(-0.77447184) q[0];
sx q[0];
rz(-1.2694556) q[0];
x q[1];
rz(2.3241256) q[2];
sx q[2];
rz(-0.82767476) q[2];
sx q[2];
rz(-2.4823014) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0468415) q[1];
sx q[1];
rz(-1.0221204) q[1];
sx q[1];
rz(1.8803174) q[1];
rz(-pi) q[2];
rz(3.0008121) q[3];
sx q[3];
rz(-0.13354145) q[3];
sx q[3];
rz(-0.30511607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28477731) q[2];
sx q[2];
rz(-1.6723526) q[2];
sx q[2];
rz(-2.3397297) q[2];
rz(0.91184584) q[3];
sx q[3];
rz(-0.82404476) q[3];
sx q[3];
rz(0.97889939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80482471) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(-2.5208852) q[0];
rz(2.8630818) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(2.1381569) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28351706) q[0];
sx q[0];
rz(-1.2293613) q[0];
sx q[0];
rz(-2.8056739) q[0];
rz(-1.150028) q[2];
sx q[2];
rz(-1.506797) q[2];
sx q[2];
rz(0.45129946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.675337) q[1];
sx q[1];
rz(-1.8395443) q[1];
sx q[1];
rz(1.6273725) q[1];
rz(-pi) q[2];
rz(2.8025556) q[3];
sx q[3];
rz(-0.62612306) q[3];
sx q[3];
rz(-2.1563185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3044546) q[2];
sx q[2];
rz(-2.0292323) q[2];
sx q[2];
rz(-0.066085286) q[2];
rz(-1.1262013) q[3];
sx q[3];
rz(-0.76032138) q[3];
sx q[3];
rz(2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
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
rz(-2.1635677) q[1];
sx q[1];
rz(2.0599005) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3560548) q[0];
sx q[0];
rz(-0.83844664) q[0];
sx q[0];
rz(-2.4466956) q[0];
rz(-pi) q[1];
rz(2.0400286) q[2];
sx q[2];
rz(-1.9479556) q[2];
sx q[2];
rz(2.9531329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8199205) q[1];
sx q[1];
rz(-1.5770438) q[1];
sx q[1];
rz(2.856918) q[1];
rz(2.1882964) q[3];
sx q[3];
rz(-2.3205608) q[3];
sx q[3];
rz(-1.8428857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(2.2990189) q[2];
rz(-0.27635559) q[3];
sx q[3];
rz(-2.1601951) q[3];
sx q[3];
rz(-1.792255) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2862947) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(-1.9516161) q[0];
rz(1.2225993) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(0.58293265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5385206) q[0];
sx q[0];
rz(-0.6003418) q[0];
sx q[0];
rz(-0.64774127) q[0];
x q[1];
rz(3.1276064) q[2];
sx q[2];
rz(-2.1443233) q[2];
sx q[2];
rz(-1.397246) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5869904) q[1];
sx q[1];
rz(-1.6515367) q[1];
sx q[1];
rz(2.9029773) q[1];
rz(-2.575608) q[3];
sx q[3];
rz(-0.71403394) q[3];
sx q[3];
rz(0.66431724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0827834) q[2];
sx q[2];
rz(-1.8913816) q[2];
sx q[2];
rz(-2.1018551) q[2];
rz(0.58935634) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(-2.1350433) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650918) q[0];
sx q[0];
rz(-0.14987513) q[0];
sx q[0];
rz(1.804922) q[0];
rz(-0.22557766) q[1];
sx q[1];
rz(-2.070319) q[1];
sx q[1];
rz(-0.76990661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5107489) q[0];
sx q[0];
rz(-2.0389557) q[0];
sx q[0];
rz(1.7450784) q[0];
x q[1];
rz(2.2397584) q[2];
sx q[2];
rz(-2.2775501) q[2];
sx q[2];
rz(1.2400986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47178956) q[1];
sx q[1];
rz(-0.21179971) q[1];
sx q[1];
rz(2.7752911) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6490593) q[3];
sx q[3];
rz(-2.1746082) q[3];
sx q[3];
rz(2.9989238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62468195) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90149752) q[0];
sx q[0];
rz(-2.5101341) q[0];
sx q[0];
rz(-0.39328662) q[0];
rz(0.45010629) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(-3.111305) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5092376) q[0];
sx q[0];
rz(-1.3233174) q[0];
sx q[0];
rz(1.6028274) q[0];
rz(-1.302839) q[2];
sx q[2];
rz(-0.7795802) q[2];
sx q[2];
rz(1.5991581) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68037625) q[1];
sx q[1];
rz(-0.54619563) q[1];
sx q[1];
rz(-0.044574634) q[1];
rz(-pi) q[2];
rz(-1.7025331) q[3];
sx q[3];
rz(-1.3088969) q[3];
sx q[3];
rz(-1.1887664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1071189) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(-2.0195473) q[2];
rz(-2.0864887) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(-1.25157) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17290641) q[0];
sx q[0];
rz(-2.1538669) q[0];
sx q[0];
rz(2.3671142) q[0];
rz(-0.6340181) q[1];
sx q[1];
rz(-1.0301544) q[1];
sx q[1];
rz(1.0672306) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2174847) q[0];
sx q[0];
rz(-1.5536246) q[0];
sx q[0];
rz(0.73826684) q[0];
x q[1];
rz(1.8781233) q[2];
sx q[2];
rz(-1.5241703) q[2];
sx q[2];
rz(1.7152552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6744819) q[1];
sx q[1];
rz(-1.8201105) q[1];
sx q[1];
rz(1.2500019) q[1];
rz(-pi) q[2];
rz(-1.2379856) q[3];
sx q[3];
rz(-2.6389671) q[3];
sx q[3];
rz(0.42325936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5344703) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(1.8758476) q[2];
rz(0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-3.1025187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604615) q[0];
sx q[0];
rz(-0.88812319) q[0];
sx q[0];
rz(-1.8817425) q[0];
rz(1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(-3.1390417) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1387864) q[0];
sx q[0];
rz(-1.4454137) q[0];
sx q[0];
rz(3.0210397) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6348636) q[2];
sx q[2];
rz(-1.5612215) q[2];
sx q[2];
rz(-2.9970882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5469909) q[1];
sx q[1];
rz(-2.3560215) q[1];
sx q[1];
rz(-0.57348302) q[1];
rz(-pi) q[2];
rz(2.0052913) q[3];
sx q[3];
rz(-1.5948203) q[3];
sx q[3];
rz(-1.0917808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(0.67443887) q[2];
rz(-1.1174348) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(-0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6013026) q[2];
sx q[2];
rz(-0.39763309) q[2];
sx q[2];
rz(2.8297505) q[2];
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

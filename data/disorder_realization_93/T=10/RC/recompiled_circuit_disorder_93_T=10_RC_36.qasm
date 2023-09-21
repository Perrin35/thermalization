OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(-1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7605654) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(2.7678124) q[0];
x q[1];
rz(-2.997424) q[2];
sx q[2];
rz(-1.2906133) q[2];
sx q[2];
rz(2.7902214) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2687298) q[1];
sx q[1];
rz(-1.9951207) q[1];
sx q[1];
rz(2.5904168) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63763036) q[3];
sx q[3];
rz(-2.5507567) q[3];
sx q[3];
rz(1.5096111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(-1.9447928) q[0];
rz(-2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.4555567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6360977) q[0];
sx q[0];
rz(-0.60244766) q[0];
sx q[0];
rz(1.8683744) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6529085) q[2];
sx q[2];
rz(-1.1973235) q[2];
sx q[2];
rz(0.79370802) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93745366) q[1];
sx q[1];
rz(-2.0935645) q[1];
sx q[1];
rz(0.015603113) q[1];
x q[2];
rz(-0.060828408) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(2.9591566) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(-2.8296208) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.172539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55484178) q[0];
sx q[0];
rz(-1.2625853) q[0];
sx q[0];
rz(-0.59535938) q[0];
x q[1];
rz(1.8981947) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(-0.74795216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2492003) q[1];
sx q[1];
rz(-1.3274267) q[1];
sx q[1];
rz(-2.1057486) q[1];
x q[2];
rz(1.3451194) q[3];
sx q[3];
rz(-2.884622) q[3];
sx q[3];
rz(-0.32303177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-0.91119901) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38055414) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(-3.0134841) q[0];
rz(0.076106636) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(2.6180843) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1062766) q[0];
sx q[0];
rz(-1.2020532) q[0];
sx q[0];
rz(-0.62555255) q[0];
rz(1.4857616) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(2.0563682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7586786) q[1];
sx q[1];
rz(-0.74712336) q[1];
sx q[1];
rz(1.9488571) q[1];
x q[2];
rz(-3.0047699) q[3];
sx q[3];
rz(-2.1577583) q[3];
sx q[3];
rz(-2.6329991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(0.28856746) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(-0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(0.99194828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6742453) q[0];
sx q[0];
rz(-1.464198) q[0];
sx q[0];
rz(-0.17515134) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8736585) q[2];
sx q[2];
rz(-1.5793243) q[2];
sx q[2];
rz(-1.0964583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3757513) q[1];
sx q[1];
rz(-1.6110794) q[1];
sx q[1];
rz(2.8326616) q[1];
rz(1.4335853) q[3];
sx q[3];
rz(-0.76223323) q[3];
sx q[3];
rz(1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.7927992) q[0];
rz(0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8558559) q[0];
sx q[0];
rz(-2.5897313) q[0];
sx q[0];
rz(-2.7049271) q[0];
rz(-pi) q[1];
rz(1.233333) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(0.5459107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52275601) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(1.7230117) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.607843) q[3];
sx q[3];
rz(-2.1042049) q[3];
sx q[3];
rz(0.48983869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(-0.56387222) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(1.58889) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(-2.3197876) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1497027) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(-0.71233149) q[0];
rz(-2.2981811) q[2];
sx q[2];
rz(-1.1864098) q[2];
sx q[2];
rz(2.9316528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6268839) q[1];
sx q[1];
rz(-0.85299546) q[1];
sx q[1];
rz(1.2674598) q[1];
x q[2];
rz(-2.5542198) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(-1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(2.896893) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(-1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.8364505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1730723) q[0];
sx q[0];
rz(-2.1497796) q[0];
sx q[0];
rz(1.9645683) q[0];
rz(-pi) q[1];
rz(1.201198) q[2];
sx q[2];
rz(-2.2410789) q[2];
sx q[2];
rz(3.0073462) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95272428) q[1];
sx q[1];
rz(-0.75718588) q[1];
sx q[1];
rz(0.56307478) q[1];
rz(1.9100902) q[3];
sx q[3];
rz(-1.5929856) q[3];
sx q[3];
rz(2.0458178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.4902327) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-0.70294356) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2962869) q[0];
sx q[0];
rz(-1.9585113) q[0];
sx q[0];
rz(-1.182343) q[0];
rz(0.039673294) q[2];
sx q[2];
rz(-2.0797605) q[2];
sx q[2];
rz(-0.85658011) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34300464) q[1];
sx q[1];
rz(-1.4061905) q[1];
sx q[1];
rz(1.027147) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1861107) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(-2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(2.83589) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(0.46863619) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271284) q[0];
sx q[0];
rz(-0.44136029) q[0];
sx q[0];
rz(-2.6370185) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2655728) q[2];
sx q[2];
rz(-0.48601905) q[2];
sx q[2];
rz(-1.3181869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.66099) q[1];
sx q[1];
rz(-0.27462474) q[1];
sx q[1];
rz(-1.8250699) q[1];
rz(-pi) q[2];
rz(2.7077984) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(-2.1123561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6580711) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(-1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464012) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(-0.62190965) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-2.2767699) q[2];
sx q[2];
rz(-1.0100126) q[2];
sx q[2];
rz(1.0457912) q[2];
rz(-1.6735531) q[3];
sx q[3];
rz(-2.4874874) q[3];
sx q[3];
rz(-2.7174674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
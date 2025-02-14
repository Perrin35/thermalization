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
rz(0.33557284) q[0];
sx q[0];
rz(4.8290401) q[0];
sx q[0];
rz(7.3168559) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(2.1318336) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5752707) q[0];
sx q[0];
rz(-0.22804865) q[0];
sx q[0];
rz(0.98342498) q[0];
rz(-1.5255341) q[2];
sx q[2];
rz(-0.082556225) q[2];
sx q[2];
rz(-0.98708594) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.941135) q[1];
sx q[1];
rz(-1.4054646) q[1];
sx q[1];
rz(1.6881866) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5663773) q[3];
sx q[3];
rz(-1.7130245) q[3];
sx q[3];
rz(2.9517236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7734163) q[2];
sx q[2];
rz(-0.877031) q[2];
sx q[2];
rz(2.3167493) q[2];
rz(-0.33251897) q[3];
sx q[3];
rz(-2.2907292) q[3];
sx q[3];
rz(0.48377812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9124209) q[0];
sx q[0];
rz(-3.0569172) q[0];
sx q[0];
rz(2.605865) q[0];
rz(-2.3748705) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(-1.7492693) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9014382) q[0];
sx q[0];
rz(-1.9124228) q[0];
sx q[0];
rz(1.0761976) q[0];
x q[1];
rz(2.8004592) q[2];
sx q[2];
rz(-1.5798414) q[2];
sx q[2];
rz(-1.1760528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.630188) q[1];
sx q[1];
rz(-1.5492445) q[1];
sx q[1];
rz(2.3261916) q[1];
rz(-pi) q[2];
rz(1.4658607) q[3];
sx q[3];
rz(-1.4426878) q[3];
sx q[3];
rz(-2.6554633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9048189) q[2];
sx q[2];
rz(-0.83319131) q[2];
sx q[2];
rz(2.0429677) q[2];
rz(0.83313471) q[3];
sx q[3];
rz(-1.5980709) q[3];
sx q[3];
rz(1.0540849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5611834) q[0];
sx q[0];
rz(-1.9920749) q[0];
sx q[0];
rz(-2.4533601) q[0];
rz(-1.8904103) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(-1.2122663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040846856) q[0];
sx q[0];
rz(-1.0991352) q[0];
sx q[0];
rz(0.77420401) q[0];
rz(-pi) q[1];
x q[1];
rz(2.278944) q[2];
sx q[2];
rz(-0.35637636) q[2];
sx q[2];
rz(1.732812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1656629) q[1];
sx q[1];
rz(-1.5651795) q[1];
sx q[1];
rz(2.4826239) q[1];
rz(-1.6501649) q[3];
sx q[3];
rz(-1.575633) q[3];
sx q[3];
rz(-1.2629379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7324098) q[2];
sx q[2];
rz(-0.99486351) q[2];
sx q[2];
rz(0.44962064) q[2];
rz(-0.85363394) q[3];
sx q[3];
rz(-0.64695224) q[3];
sx q[3];
rz(-0.70931119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39321008) q[0];
sx q[0];
rz(-1.0164096) q[0];
sx q[0];
rz(0.40312314) q[0];
rz(2.368811) q[1];
sx q[1];
rz(-1.2330331) q[1];
sx q[1];
rz(0.82537878) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6306303) q[0];
sx q[0];
rz(-1.2494506) q[0];
sx q[0];
rz(-2.3031631) q[0];
rz(-pi) q[1];
rz(0.72751491) q[2];
sx q[2];
rz(-0.036002654) q[2];
sx q[2];
rz(1.9787479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1267438) q[1];
sx q[1];
rz(-0.9909317) q[1];
sx q[1];
rz(-0.090437263) q[1];
rz(-pi) q[2];
rz(0.89572923) q[3];
sx q[3];
rz(-0.78215996) q[3];
sx q[3];
rz(1.6767927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0394502) q[2];
sx q[2];
rz(-1.2692229) q[2];
sx q[2];
rz(2.6329182) q[2];
rz(1.2447119) q[3];
sx q[3];
rz(-2.0679943) q[3];
sx q[3];
rz(1.3332453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4077069) q[0];
sx q[0];
rz(-1.5579959) q[0];
sx q[0];
rz(-0.35633126) q[0];
rz(1.7286667) q[1];
sx q[1];
rz(-1.8316385) q[1];
sx q[1];
rz(-1.2948571) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6083487) q[0];
sx q[0];
rz(-0.49194592) q[0];
sx q[0];
rz(-2.66376) q[0];
x q[1];
rz(0.21759502) q[2];
sx q[2];
rz(-0.67309531) q[2];
sx q[2];
rz(0.714528) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7562189) q[1];
sx q[1];
rz(-1.1750147) q[1];
sx q[1];
rz(1.3211005) q[1];
x q[2];
rz(0.44900198) q[3];
sx q[3];
rz(-1.0706361) q[3];
sx q[3];
rz(-1.7412108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20092189) q[2];
sx q[2];
rz(-1.385043) q[2];
sx q[2];
rz(0.6400288) q[2];
rz(0.43404964) q[3];
sx q[3];
rz(-2.1899352) q[3];
sx q[3];
rz(1.1205589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5956748) q[0];
sx q[0];
rz(-2.7037103) q[0];
sx q[0];
rz(-2.3906999) q[0];
rz(-2.3789876) q[1];
sx q[1];
rz(-1.7060988) q[1];
sx q[1];
rz(-2.6079752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5845244) q[0];
sx q[0];
rz(-1.2831339) q[0];
sx q[0];
rz(0.30867048) q[0];
rz(-pi) q[1];
rz(-2.0178823) q[2];
sx q[2];
rz(-2.3229282) q[2];
sx q[2];
rz(2.1731203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68249629) q[1];
sx q[1];
rz(-2.4630416) q[1];
sx q[1];
rz(1.7595923) q[1];
x q[2];
rz(-2.153965) q[3];
sx q[3];
rz(-1.0915874) q[3];
sx q[3];
rz(-2.2892009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.077347191) q[2];
sx q[2];
rz(-2.7392445) q[2];
sx q[2];
rz(0.8858436) q[2];
rz(1.8387851) q[3];
sx q[3];
rz(-1.1833444) q[3];
sx q[3];
rz(2.6634789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9293905) q[0];
sx q[0];
rz(-2.2181856) q[0];
sx q[0];
rz(-0.59342629) q[0];
rz(-1.5133096) q[1];
sx q[1];
rz(-1.1791041) q[1];
sx q[1];
rz(-0.57366192) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79686576) q[0];
sx q[0];
rz(-2.0484951) q[0];
sx q[0];
rz(1.1487238) q[0];
rz(2.5638177) q[2];
sx q[2];
rz(-1.7794357) q[2];
sx q[2];
rz(-1.8633757) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5031858) q[1];
sx q[1];
rz(-1.3246623) q[1];
sx q[1];
rz(2.7315774) q[1];
rz(3.0629458) q[3];
sx q[3];
rz(-1.5203195) q[3];
sx q[3];
rz(-0.61928643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4297428) q[2];
sx q[2];
rz(-1.2028799) q[2];
sx q[2];
rz(1.708606) q[2];
rz(1.1687219) q[3];
sx q[3];
rz(-0.43787268) q[3];
sx q[3];
rz(1.6247113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2699921) q[0];
sx q[0];
rz(-2.2347436) q[0];
sx q[0];
rz(-0.20371833) q[0];
rz(-1.8261955) q[1];
sx q[1];
rz(-1.8351277) q[1];
sx q[1];
rz(1.809583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341605) q[0];
sx q[0];
rz(-2.1825509) q[0];
sx q[0];
rz(2.3513398) q[0];
rz(-pi) q[1];
rz(-0.56925781) q[2];
sx q[2];
rz(-2.1569886) q[2];
sx q[2];
rz(1.4365412) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97196619) q[1];
sx q[1];
rz(-2.1887767) q[1];
sx q[1];
rz(-2.0566259) q[1];
rz(-pi) q[2];
rz(1.733794) q[3];
sx q[3];
rz(-1.4428815) q[3];
sx q[3];
rz(2.7892512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77067644) q[2];
sx q[2];
rz(-0.081826536) q[2];
sx q[2];
rz(1.7601684) q[2];
rz(2.5904739) q[3];
sx q[3];
rz(-1.8128606) q[3];
sx q[3];
rz(2.3962925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119174) q[0];
sx q[0];
rz(-1.6918007) q[0];
sx q[0];
rz(-1.2592738) q[0];
rz(-1.3445541) q[1];
sx q[1];
rz(-2.5090736) q[1];
sx q[1];
rz(1.9154027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7525402) q[0];
sx q[0];
rz(-1.2864094) q[0];
sx q[0];
rz(3.000598) q[0];
x q[1];
rz(0.8441458) q[2];
sx q[2];
rz(-0.54960362) q[2];
sx q[2];
rz(-1.5908102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75424947) q[1];
sx q[1];
rz(-0.74973901) q[1];
sx q[1];
rz(-0.096258817) q[1];
rz(-2.1192083) q[3];
sx q[3];
rz(-0.58956857) q[3];
sx q[3];
rz(-2.0062674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.58174497) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(1.5398514) q[2];
rz(-1.3567443) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(-1.9073568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27720472) q[0];
sx q[0];
rz(-1.5978403) q[0];
sx q[0];
rz(-0.10449115) q[0];
rz(-1.744572) q[1];
sx q[1];
rz(-1.3518159) q[1];
sx q[1];
rz(1.5207965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6817271) q[0];
sx q[0];
rz(-2.2978373) q[0];
sx q[0];
rz(-1.8521502) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7303329) q[2];
sx q[2];
rz(-1.2345833) q[2];
sx q[2];
rz(1.5000864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0238259) q[1];
sx q[1];
rz(-1.858288) q[1];
sx q[1];
rz(-1.6790798) q[1];
rz(0.057985882) q[3];
sx q[3];
rz(-0.99000217) q[3];
sx q[3];
rz(-2.557366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70667679) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(-1.8527156) q[2];
rz(-1.0051109) q[3];
sx q[3];
rz(-2.0821327) q[3];
sx q[3];
rz(3.03481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2224251) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(-0.87986058) q[1];
sx q[1];
rz(-0.51645551) q[1];
sx q[1];
rz(-2.5194306) q[1];
rz(0.11779412) q[2];
sx q[2];
rz(-0.5024903) q[2];
sx q[2];
rz(2.6033664) q[2];
rz(1.6493116) q[3];
sx q[3];
rz(-2.169486) q[3];
sx q[3];
rz(-1.538856) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

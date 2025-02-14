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
rz(4.5728788) q[1];
sx q[1];
rz(3.6990777) q[1];
sx q[1];
rz(7.2929444) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57085867) q[0];
sx q[0];
rz(-1.6964127) q[0];
sx q[0];
rz(1.3799589) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4883243) q[2];
sx q[2];
rz(-1.5745275) q[2];
sx q[2];
rz(0.53860215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20045763) q[1];
sx q[1];
rz(-1.7361281) q[1];
sx q[1];
rz(-1.6881866) q[1];
rz(-pi) q[2];
x q[2];
rz(0.030850284) q[3];
sx q[3];
rz(-2.9992963) q[3];
sx q[3];
rz(-0.22103413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7734163) q[2];
sx q[2];
rz(-0.877031) q[2];
sx q[2];
rz(2.3167493) q[2];
rz(2.8090737) q[3];
sx q[3];
rz(-0.85086346) q[3];
sx q[3];
rz(2.6578145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9124209) q[0];
sx q[0];
rz(-0.084675463) q[0];
sx q[0];
rz(2.605865) q[0];
rz(2.3748705) q[1];
sx q[1];
rz(-1.352939) q[1];
sx q[1];
rz(1.3923233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24015445) q[0];
sx q[0];
rz(-1.9124228) q[0];
sx q[0];
rz(1.0761976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.027030305) q[2];
sx q[2];
rz(-0.34124869) q[2];
sx q[2];
rz(2.7213767) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.036511985) q[1];
sx q[1];
rz(-2.3859508) q[1];
sx q[1];
rz(1.6022268) q[1];
x q[2];
rz(3.0127834) q[3];
sx q[3];
rz(-1.6748689) q[3];
sx q[3];
rz(-1.0712119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9048189) q[2];
sx q[2];
rz(-0.83319131) q[2];
sx q[2];
rz(-1.0986249) q[2];
rz(-2.3084579) q[3];
sx q[3];
rz(-1.5980709) q[3];
sx q[3];
rz(-2.0875077) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58040923) q[0];
sx q[0];
rz(-1.1495178) q[0];
sx q[0];
rz(-2.4533601) q[0];
rz(1.8904103) q[1];
sx q[1];
rz(-0.65443188) q[1];
sx q[1];
rz(1.9293264) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1935342) q[0];
sx q[0];
rz(-0.89846629) q[0];
sx q[0];
rz(0.95114077) q[0];
x q[1];
rz(-2.9040292) q[2];
sx q[2];
rz(-1.8389987) q[2];
sx q[2];
rz(0.99246565) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9759298) q[1];
sx q[1];
rz(-1.5764131) q[1];
sx q[1];
rz(-2.4826239) q[1];
rz(-pi) q[2];
rz(-1.5098677) q[3];
sx q[3];
rz(-3.0620771) q[3];
sx q[3];
rz(2.7729976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4091829) q[2];
sx q[2];
rz(-2.1467291) q[2];
sx q[2];
rz(-2.691972) q[2];
rz(0.85363394) q[3];
sx q[3];
rz(-0.64695224) q[3];
sx q[3];
rz(0.70931119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.39321008) q[0];
sx q[0];
rz(-2.125183) q[0];
sx q[0];
rz(-2.7384695) q[0];
rz(2.368811) q[1];
sx q[1];
rz(-1.9085596) q[1];
sx q[1];
rz(-0.82537878) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8050418) q[0];
sx q[0];
rz(-2.2580349) q[0];
sx q[0];
rz(-2.7206793) q[0];
x q[1];
rz(3.1146997) q[2];
sx q[2];
rz(-1.5947358) q[2];
sx q[2];
rz(-0.31924143) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6352915) q[1];
sx q[1];
rz(-1.4951733) q[1];
sx q[1];
rz(2.1525394) q[1];
rz(2.585936) q[3];
sx q[3];
rz(-0.98816493) q[3];
sx q[3];
rz(0.61911303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0394502) q[2];
sx q[2];
rz(-1.2692229) q[2];
sx q[2];
rz(0.5086745) q[2];
rz(-1.2447119) q[3];
sx q[3];
rz(-1.0735984) q[3];
sx q[3];
rz(1.3332453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73388571) q[0];
sx q[0];
rz(-1.5579959) q[0];
sx q[0];
rz(-0.35633126) q[0];
rz(1.7286667) q[1];
sx q[1];
rz(-1.3099542) q[1];
sx q[1];
rz(1.2948571) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.533244) q[0];
sx q[0];
rz(-0.49194592) q[0];
sx q[0];
rz(2.66376) q[0];
x q[1];
rz(-2.4800972) q[2];
sx q[2];
rz(-1.7057888) q[2];
sx q[2];
rz(-2.1141426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2834245) q[1];
sx q[1];
rz(-1.8008261) q[1];
sx q[1];
rz(-0.40711222) q[1];
rz(2.6925907) q[3];
sx q[3];
rz(-2.0709566) q[3];
sx q[3];
rz(1.4003818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9406708) q[2];
sx q[2];
rz(-1.385043) q[2];
sx q[2];
rz(-2.5015639) q[2];
rz(0.43404964) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(2.0210338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5956748) q[0];
sx q[0];
rz(-2.7037103) q[0];
sx q[0];
rz(0.75089279) q[0];
rz(2.3789876) q[1];
sx q[1];
rz(-1.7060988) q[1];
sx q[1];
rz(2.6079752) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92350875) q[0];
sx q[0];
rz(-1.2752011) q[0];
sx q[0];
rz(1.8718998) q[0];
x q[1];
rz(2.0178823) q[2];
sx q[2];
rz(-0.81866449) q[2];
sx q[2];
rz(2.1731203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2184132) q[1];
sx q[1];
rz(-0.90649062) q[1];
sx q[1];
rz(-2.9914145) q[1];
x q[2];
rz(-0.814416) q[3];
sx q[3];
rz(-0.73668639) q[3];
sx q[3];
rz(-1.8131922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0642455) q[2];
sx q[2];
rz(-2.7392445) q[2];
sx q[2];
rz(2.2557491) q[2];
rz(-1.8387851) q[3];
sx q[3];
rz(-1.1833444) q[3];
sx q[3];
rz(0.47811374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2122022) q[0];
sx q[0];
rz(-0.92340702) q[0];
sx q[0];
rz(-2.5481664) q[0];
rz(1.5133096) q[1];
sx q[1];
rz(-1.9624886) q[1];
sx q[1];
rz(-0.57366192) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3447269) q[0];
sx q[0];
rz(-1.0930976) q[0];
sx q[0];
rz(1.1487238) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3232347) q[2];
sx q[2];
rz(-1.0070966) q[2];
sx q[2];
rz(-2.7147788) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5031858) q[1];
sx q[1];
rz(-1.8169303) q[1];
sx q[1];
rz(-2.7315774) q[1];
x q[2];
rz(-3.0629458) q[3];
sx q[3];
rz(-1.5203195) q[3];
sx q[3];
rz(0.61928643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4297428) q[2];
sx q[2];
rz(-1.2028799) q[2];
sx q[2];
rz(-1.708606) q[2];
rz(-1.1687219) q[3];
sx q[3];
rz(-2.70372) q[3];
sx q[3];
rz(-1.5168813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2699921) q[0];
sx q[0];
rz(-0.90684909) q[0];
sx q[0];
rz(-0.20371833) q[0];
rz(1.8261955) q[1];
sx q[1];
rz(-1.306465) q[1];
sx q[1];
rz(1.809583) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341605) q[0];
sx q[0];
rz(-0.95904175) q[0];
sx q[0];
rz(0.79025288) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90317995) q[2];
sx q[2];
rz(-1.1051264) q[2];
sx q[2];
rz(-2.6671034) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97196619) q[1];
sx q[1];
rz(-2.1887767) q[1];
sx q[1];
rz(-2.0566259) q[1];
x q[2];
rz(3.0119786) q[3];
sx q[3];
rz(-1.7324505) q[3];
sx q[3];
rz(1.9441138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3709162) q[2];
sx q[2];
rz(-0.081826536) q[2];
sx q[2];
rz(-1.7601684) q[2];
rz(-2.5904739) q[3];
sx q[3];
rz(-1.8128606) q[3];
sx q[3];
rz(0.74530017) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119174) q[0];
sx q[0];
rz(-1.4497919) q[0];
sx q[0];
rz(1.2592738) q[0];
rz(1.3445541) q[1];
sx q[1];
rz(-0.63251907) q[1];
sx q[1];
rz(-1.22619) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38905242) q[0];
sx q[0];
rz(-1.2864094) q[0];
sx q[0];
rz(-3.000598) q[0];
x q[1];
rz(-0.38649747) q[2];
sx q[2];
rz(-1.9718687) q[2];
sx q[2];
rz(2.3569466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75424947) q[1];
sx q[1];
rz(-0.74973901) q[1];
sx q[1];
rz(-0.096258817) q[1];
rz(-pi) q[2];
rz(-0.33554797) q[3];
sx q[3];
rz(-2.0651544) q[3];
sx q[3];
rz(1.7691607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5598477) q[2];
sx q[2];
rz(-2.2793844) q[2];
sx q[2];
rz(-1.5398514) q[2];
rz(1.7848484) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(1.2342359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8643879) q[0];
sx q[0];
rz(-1.5437523) q[0];
sx q[0];
rz(3.0371015) q[0];
rz(1.3970207) q[1];
sx q[1];
rz(-1.3518159) q[1];
sx q[1];
rz(-1.6207961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45986555) q[0];
sx q[0];
rz(-0.84375536) q[0];
sx q[0];
rz(1.8521502) q[0];
x q[1];
rz(1.7303329) q[2];
sx q[2];
rz(-1.9070093) q[2];
sx q[2];
rz(-1.6415063) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4222153) q[1];
sx q[1];
rz(-1.6746192) q[1];
sx q[1];
rz(0.28909282) q[1];
rz(2.1523624) q[3];
sx q[3];
rz(-1.6192659) q[3];
sx q[3];
rz(-1.0184108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4349159) q[2];
sx q[2];
rz(-2.2997663) q[2];
sx q[2];
rz(-1.288877) q[2];
rz(2.1364818) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(0.10678261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91916753) q[0];
sx q[0];
rz(-1.5594302) q[0];
sx q[0];
rz(2.9710309) q[0];
rz(0.87986058) q[1];
sx q[1];
rz(-2.6251371) q[1];
sx q[1];
rz(0.62216204) q[1];
rz(0.11779412) q[2];
sx q[2];
rz(-0.5024903) q[2];
sx q[2];
rz(2.6033664) q[2];
rz(-0.60012695) q[3];
sx q[3];
rz(-1.6356346) q[3];
sx q[3];
rz(0.076250565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

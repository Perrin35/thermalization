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
rz(-1.4541452) q[0];
sx q[0];
rz(-2.1079221) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(2.1318336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.570734) q[0];
sx q[0];
rz(-1.44518) q[0];
sx q[0];
rz(-1.7616338) q[0];
rz(3.1378488) q[2];
sx q[2];
rz(-1.4883248) q[2];
sx q[2];
rz(-1.0325026) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.941135) q[1];
sx q[1];
rz(-1.7361281) q[1];
sx q[1];
rz(1.453406) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14222957) q[3];
sx q[3];
rz(-1.5664219) q[3];
sx q[3];
rz(-1.7612918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7734163) q[2];
sx q[2];
rz(-0.877031) q[2];
sx q[2];
rz(2.3167493) q[2];
rz(-0.33251897) q[3];
sx q[3];
rz(-2.2907292) q[3];
sx q[3];
rz(-2.6578145) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2291718) q[0];
sx q[0];
rz(-0.084675463) q[0];
sx q[0];
rz(-2.605865) q[0];
rz(2.3748705) q[1];
sx q[1];
rz(-1.352939) q[1];
sx q[1];
rz(-1.7492693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.509393) q[0];
sx q[0];
rz(-2.0344816) q[0];
sx q[0];
rz(2.7576588) q[0];
rz(-pi) q[1];
rz(-0.34113348) q[2];
sx q[2];
rz(-1.5798414) q[2];
sx q[2];
rz(-1.1760528) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5114047) q[1];
sx q[1];
rz(-1.5492445) q[1];
sx q[1];
rz(0.81540108) q[1];
rz(1.6757319) q[3];
sx q[3];
rz(-1.6989048) q[3];
sx q[3];
rz(-2.6554633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9048189) q[2];
sx q[2];
rz(-2.3084013) q[2];
sx q[2];
rz(-2.0429677) q[2];
rz(-2.3084579) q[3];
sx q[3];
rz(-1.5435217) q[3];
sx q[3];
rz(-1.0540849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5611834) q[0];
sx q[0];
rz(-1.9920749) q[0];
sx q[0];
rz(2.4533601) q[0];
rz(-1.8904103) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(1.9293264) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935342) q[0];
sx q[0];
rz(-2.2431264) q[0];
sx q[0];
rz(0.95114077) q[0];
rz(-pi) q[1];
rz(1.8463676) q[2];
sx q[2];
rz(-1.7997121) q[2];
sx q[2];
rz(-2.6273397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1656629) q[1];
sx q[1];
rz(-1.5651795) q[1];
sx q[1];
rz(-0.65896874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.631725) q[3];
sx q[3];
rz(-0.079515545) q[3];
sx q[3];
rz(0.36859504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7324098) q[2];
sx q[2];
rz(-0.99486351) q[2];
sx q[2];
rz(0.44962064) q[2];
rz(-0.85363394) q[3];
sx q[3];
rz(-2.4946404) q[3];
sx q[3];
rz(0.70931119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7483826) q[0];
sx q[0];
rz(-1.0164096) q[0];
sx q[0];
rz(-2.7384695) q[0];
rz(2.368811) q[1];
sx q[1];
rz(-1.2330331) q[1];
sx q[1];
rz(0.82537878) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3365509) q[0];
sx q[0];
rz(-0.8835578) q[0];
sx q[0];
rz(-2.7206793) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1146997) q[2];
sx q[2];
rz(-1.5468569) q[2];
sx q[2];
rz(-2.8223512) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5063012) q[1];
sx q[1];
rz(-1.6464194) q[1];
sx q[1];
rz(-2.1525394) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2458634) q[3];
sx q[3];
rz(-2.3594327) q[3];
sx q[3];
rz(-1.6767927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1021425) q[2];
sx q[2];
rz(-1.8723698) q[2];
sx q[2];
rz(2.6329182) q[2];
rz(1.8968808) q[3];
sx q[3];
rz(-1.0735984) q[3];
sx q[3];
rz(-1.8083474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4077069) q[0];
sx q[0];
rz(-1.5835967) q[0];
sx q[0];
rz(-2.7852614) q[0];
rz(1.412926) q[1];
sx q[1];
rz(-1.3099542) q[1];
sx q[1];
rz(-1.2948571) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0771478) q[0];
sx q[0];
rz(-2.003621) q[0];
sx q[0];
rz(1.3291784) q[0];
x q[1];
rz(2.9239976) q[2];
sx q[2];
rz(-2.4684973) q[2];
sx q[2];
rz(-2.4270647) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7562189) q[1];
sx q[1];
rz(-1.1750147) q[1];
sx q[1];
rz(1.3211005) q[1];
rz(1.0255085) q[3];
sx q[3];
rz(-1.9615615) q[3];
sx q[3];
rz(-2.7440967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.20092189) q[2];
sx q[2];
rz(-1.385043) q[2];
sx q[2];
rz(-2.5015639) q[2];
rz(-0.43404964) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(1.1205589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5459179) q[0];
sx q[0];
rz(-2.7037103) q[0];
sx q[0];
rz(2.3906999) q[0];
rz(2.3789876) q[1];
sx q[1];
rz(-1.4354939) q[1];
sx q[1];
rz(-2.6079752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92350875) q[0];
sx q[0];
rz(-1.2752011) q[0];
sx q[0];
rz(-1.2696928) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7087174) q[2];
sx q[2];
rz(-2.2895669) q[2];
sx q[2];
rz(1.5804039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68249629) q[1];
sx q[1];
rz(-2.4630416) q[1];
sx q[1];
rz(1.7595923) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5848006) q[3];
sx q[3];
rz(-1.0602615) q[3];
sx q[3];
rz(0.42312121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.077347191) q[2];
sx q[2];
rz(-2.7392445) q[2];
sx q[2];
rz(-2.2557491) q[2];
rz(-1.8387851) q[3];
sx q[3];
rz(-1.1833444) q[3];
sx q[3];
rz(0.47811374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9293905) q[0];
sx q[0];
rz(-0.92340702) q[0];
sx q[0];
rz(-2.5481664) q[0];
rz(-1.5133096) q[1];
sx q[1];
rz(-1.9624886) q[1];
sx q[1];
rz(-2.5679307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57034111) q[0];
sx q[0];
rz(-1.9431337) q[0];
sx q[0];
rz(-2.6254197) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5638177) q[2];
sx q[2];
rz(-1.7794357) q[2];
sx q[2];
rz(-1.278217) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5031858) q[1];
sx q[1];
rz(-1.3246623) q[1];
sx q[1];
rz(0.41001525) q[1];
rz(1.6214293) q[3];
sx q[3];
rz(-1.4922499) q[3];
sx q[3];
rz(0.95548624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4297428) q[2];
sx q[2];
rz(-1.2028799) q[2];
sx q[2];
rz(-1.708606) q[2];
rz(-1.9728707) q[3];
sx q[3];
rz(-0.43787268) q[3];
sx q[3];
rz(1.6247113) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8716005) q[0];
sx q[0];
rz(-2.2347436) q[0];
sx q[0];
rz(-0.20371833) q[0];
rz(1.8261955) q[1];
sx q[1];
rz(-1.306465) q[1];
sx q[1];
rz(1.809583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7527134) q[0];
sx q[0];
rz(-2.1915771) q[0];
sx q[0];
rz(2.3546773) q[0];
x q[1];
rz(2.2384127) q[2];
sx q[2];
rz(-2.0364663) q[2];
sx q[2];
rz(2.6671034) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1696265) q[1];
sx q[1];
rz(-2.1887767) q[1];
sx q[1];
rz(2.0566259) q[1];
rz(-pi) q[2];
x q[2];
rz(1.733794) q[3];
sx q[3];
rz(-1.6987112) q[3];
sx q[3];
rz(-2.7892512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77067644) q[2];
sx q[2];
rz(-3.0597661) q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7119174) q[0];
sx q[0];
rz(-1.4497919) q[0];
sx q[0];
rz(-1.2592738) q[0];
rz(1.3445541) q[1];
sx q[1];
rz(-2.5090736) q[1];
sx q[1];
rz(-1.9154027) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9200468) q[0];
sx q[0];
rz(-1.4355005) q[0];
sx q[0];
rz(-1.8578803) q[0];
x q[1];
rz(0.38649747) q[2];
sx q[2];
rz(-1.9718687) q[2];
sx q[2];
rz(-2.3569466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88542467) q[1];
sx q[1];
rz(-2.3162335) q[1];
sx q[1];
rz(1.481545) q[1];
x q[2];
rz(-0.33554797) q[3];
sx q[3];
rz(-2.0651544) q[3];
sx q[3];
rz(-1.372432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5598477) q[2];
sx q[2];
rz(-2.2793844) q[2];
sx q[2];
rz(1.6017412) q[2];
rz(-1.3567443) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(1.2342359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(-0.27720472) q[0];
sx q[0];
rz(-1.5978403) q[0];
sx q[0];
rz(3.0371015) q[0];
rz(1.3970207) q[1];
sx q[1];
rz(-1.7897768) q[1];
sx q[1];
rz(1.6207961) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6817271) q[0];
sx q[0];
rz(-2.2978373) q[0];
sx q[0];
rz(1.8521502) q[0];
rz(-pi) q[1];
rz(0.42664032) q[2];
sx q[2];
rz(-0.37084118) q[2];
sx q[2];
rz(-1.1877354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0238259) q[1];
sx q[1];
rz(-1.858288) q[1];
sx q[1];
rz(1.4625129) q[1];
rz(-pi) q[2];
rz(-2.1523624) q[3];
sx q[3];
rz(-1.6192659) q[3];
sx q[3];
rz(-2.1231819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4349159) q[2];
sx q[2];
rz(-2.2997663) q[2];
sx q[2];
rz(-1.8527156) q[2];
rz(-2.1364818) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(-0.10678261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2224251) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(2.2617321) q[1];
sx q[1];
rz(-0.51645551) q[1];
sx q[1];
rz(-2.5194306) q[1];
rz(0.49956074) q[2];
sx q[2];
rz(-1.5141664) q[2];
sx q[2];
rz(0.92922633) q[2];
rz(0.11446791) q[3];
sx q[3];
rz(-0.60319067) q[3];
sx q[3];
rz(-1.4001455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

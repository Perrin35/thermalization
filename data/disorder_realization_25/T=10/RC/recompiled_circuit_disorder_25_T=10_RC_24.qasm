OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(-2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(-0.087892428) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0996617) q[0];
sx q[0];
rz(-0.66177216) q[0];
sx q[0];
rz(0.76627888) q[0];
rz(1.0432265) q[2];
sx q[2];
rz(-1.1967778) q[2];
sx q[2];
rz(-2.5803215) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80884113) q[1];
sx q[1];
rz(-2.07507) q[1];
sx q[1];
rz(-1.5682975) q[1];
rz(1.5486693) q[3];
sx q[3];
rz(-2.8828354) q[3];
sx q[3];
rz(-2.4165137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(-1.5287483) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(2.5879481) q[0];
rz(-1.2373295) q[1];
sx q[1];
rz(-1.7747223) q[1];
sx q[1];
rz(-1.233261) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20576661) q[0];
sx q[0];
rz(-0.92139771) q[0];
sx q[0];
rz(-2.4539024) q[0];
x q[1];
rz(0.1594752) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(-2.8603539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.065029) q[1];
sx q[1];
rz(-0.20670465) q[1];
sx q[1];
rz(1.6480584) q[1];
rz(1.3319098) q[3];
sx q[3];
rz(-0.91730648) q[3];
sx q[3];
rz(-1.4373923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0040940293) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(3.1393576) q[2];
rz(-0.8301174) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(-1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8925791) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(0.85154831) q[0];
rz(2.3705204) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(0.071203701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.203513) q[0];
sx q[0];
rz(-0.31103125) q[0];
sx q[0];
rz(0.87840338) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0729012) q[2];
sx q[2];
rz(-0.50902589) q[2];
sx q[2];
rz(-2.6976762) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.066312) q[1];
sx q[1];
rz(-1.7048786) q[1];
sx q[1];
rz(0.65384298) q[1];
rz(-pi) q[2];
rz(-2.2749388) q[3];
sx q[3];
rz(-2.1979077) q[3];
sx q[3];
rz(1.8303527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(-2.6181347) q[2];
rz(-0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.4098542) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(-0.82558924) q[0];
rz(1.1666974) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(1.694214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71890812) q[0];
sx q[0];
rz(-1.8005162) q[0];
sx q[0];
rz(-0.16750383) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7491313) q[2];
sx q[2];
rz(-1.5618088) q[2];
sx q[2];
rz(1.0967147) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96781603) q[1];
sx q[1];
rz(-2.2602091) q[1];
sx q[1];
rz(2.9348228) q[1];
x q[2];
rz(1.3284239) q[3];
sx q[3];
rz(-1.216785) q[3];
sx q[3];
rz(-2.6995475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(-2.9966667) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(-0.01468006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.1458364) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(0.68403912) q[0];
rz(-1.1902635) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.2129983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0096668) q[0];
sx q[0];
rz(-1.3103232) q[0];
sx q[0];
rz(2.7564604) q[0];
rz(-pi) q[1];
rz(2.2239387) q[2];
sx q[2];
rz(-1.0461763) q[2];
sx q[2];
rz(-1.7078924) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7291521) q[1];
sx q[1];
rz(-1.3153331) q[1];
sx q[1];
rz(-1.4490119) q[1];
rz(-pi) q[2];
rz(2.9641987) q[3];
sx q[3];
rz(-1.3131485) q[3];
sx q[3];
rz(0.28211668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(2.5115013) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1155788) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(-2.8552326) q[0];
rz(0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-0.58475959) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6385348) q[0];
sx q[0];
rz(-0.98031822) q[0];
sx q[0];
rz(-1.9496586) q[0];
rz(-pi) q[1];
rz(-2.5271687) q[2];
sx q[2];
rz(-1.5922058) q[2];
sx q[2];
rz(-1.9911839) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44536351) q[1];
sx q[1];
rz(-0.74659691) q[1];
sx q[1];
rz(-1.1397051) q[1];
rz(-pi) q[2];
rz(0.16634511) q[3];
sx q[3];
rz(-1.6206985) q[3];
sx q[3];
rz(2.9970616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9770603) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(-2.1785054) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(3.0528089) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(2.5928296) q[0];
rz(-3.0896297) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(-0.60639492) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43019766) q[0];
sx q[0];
rz(-2.093962) q[0];
sx q[0];
rz(-0.27336143) q[0];
rz(-2.7658471) q[2];
sx q[2];
rz(-2.1990015) q[2];
sx q[2];
rz(-1.9919765) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.28601521) q[1];
sx q[1];
rz(-2.3708323) q[1];
sx q[1];
rz(1.7325749) q[1];
x q[2];
rz(0.35550907) q[3];
sx q[3];
rz(-1.9572557) q[3];
sx q[3];
rz(0.10458065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8023119) q[0];
sx q[0];
rz(-0.34582129) q[0];
sx q[0];
rz(-0.75396496) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75922155) q[0];
sx q[0];
rz(-1.641944) q[0];
sx q[0];
rz(1.1573777) q[0];
rz(1.9367427) q[2];
sx q[2];
rz(-1.9908675) q[2];
sx q[2];
rz(-0.36047381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5880809) q[1];
sx q[1];
rz(-1.7953824) q[1];
sx q[1];
rz(2.5976546) q[1];
x q[2];
rz(-0.065683059) q[3];
sx q[3];
rz(-2.5816397) q[3];
sx q[3];
rz(-2.5064859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.078538744) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(-2.7770384) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(-0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446328) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(-1.0546168) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(1.1484336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371985) q[0];
sx q[0];
rz(-3.0657401) q[0];
sx q[0];
rz(1.8914468) q[0];
x q[1];
rz(0.36275136) q[2];
sx q[2];
rz(-2.0667549) q[2];
sx q[2];
rz(2.9438058) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0587522) q[1];
sx q[1];
rz(-2.0420923) q[1];
sx q[1];
rz(2.6575762) q[1];
rz(-1.5356482) q[3];
sx q[3];
rz(-1.9225549) q[3];
sx q[3];
rz(2.7621321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(0.77825528) q[2];
rz(0.19566472) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8675999) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(0.7397488) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(-0.69828066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2525576) q[0];
sx q[0];
rz(-1.6425942) q[0];
sx q[0];
rz(2.8385212) q[0];
rz(-pi) q[1];
rz(-1.2487222) q[2];
sx q[2];
rz(-3.0634355) q[2];
sx q[2];
rz(0.31101481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5729534) q[1];
sx q[1];
rz(-1.1146953) q[1];
sx q[1];
rz(0.24104636) q[1];
x q[2];
rz(-0.24153696) q[3];
sx q[3];
rz(-1.0165443) q[3];
sx q[3];
rz(-2.9709771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(1.9742112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2587851) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(0.48568934) q[2];
sx q[2];
rz(-0.66076856) q[2];
sx q[2];
rz(3.0271157) q[2];
rz(0.15016951) q[3];
sx q[3];
rz(-2.2554845) q[3];
sx q[3];
rz(3.1156202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(6.531125) q[0];
sx q[0];
rz(8.6046435) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(-0.087892428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9633006) q[0];
sx q[0];
rz(-1.1305729) q[0];
sx q[0];
rz(-0.51142366) q[0];
rz(-pi) q[1];
rz(2.0983661) q[2];
sx q[2];
rz(-1.1967778) q[2];
sx q[2];
rz(2.5803215) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81401285) q[1];
sx q[1];
rz(-2.6373133) q[1];
sx q[1];
rz(-3.1370647) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8294931) q[3];
sx q[3];
rz(-1.5764578) q[3];
sx q[3];
rz(-0.86710801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(-1.5287483) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.9083317) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20576661) q[0];
sx q[0];
rz(-0.92139771) q[0];
sx q[0];
rz(-0.68769023) q[0];
x q[1];
rz(-2.4670062) q[2];
sx q[2];
rz(-1.470675) q[2];
sx q[2];
rz(-1.9763725) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0765637) q[1];
sx q[1];
rz(-2.934888) q[1];
sx q[1];
rz(1.4935342) q[1];
rz(-2.4741461) q[3];
sx q[3];
rz(-1.7597886) q[3];
sx q[3];
rz(2.8611956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0040940293) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(0.0022350524) q[2];
rz(-2.3114752) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-0.85154831) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(-3.070389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6547346) q[0];
sx q[0];
rz(-1.808597) q[0];
sx q[0];
rz(0.20240692) q[0];
x q[1];
rz(-0.26239563) q[2];
sx q[2];
rz(-2.0121644) q[2];
sx q[2];
rz(1.0052094) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.60656602) q[1];
sx q[1];
rz(-0.92381322) q[1];
sx q[1];
rz(1.7391298) q[1];
rz(-2.2749388) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(1.31124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(-2.8097025) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7317384) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(-0.82558924) q[0];
rz(-1.9748953) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(-1.4473787) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3281876) q[0];
sx q[0];
rz(-1.4077328) q[0];
sx q[0];
rz(1.3379315) q[0];
rz(-0.0091323098) q[2];
sx q[2];
rz(-1.7491241) q[2];
sx q[2];
rz(-0.4724617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73562685) q[1];
sx q[1];
rz(-1.7298797) q[1];
sx q[1];
rz(-2.2707978) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36376603) q[3];
sx q[3];
rz(-1.7978661) q[3];
sx q[3];
rz(-1.0432537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1458364) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(2.4575535) q[0];
rz(1.1902635) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(1.9285944) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1367462) q[0];
sx q[0];
rz(-0.46127013) q[0];
sx q[0];
rz(-2.524551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.331779) q[2];
sx q[2];
rz(-0.81293101) q[2];
sx q[2];
rz(-0.44250689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2781093) q[1];
sx q[1];
rz(-0.28243318) q[1];
sx q[1];
rz(-0.43538283) q[1];
x q[2];
rz(-0.17739399) q[3];
sx q[3];
rz(-1.8284441) q[3];
sx q[3];
rz(2.859476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1155788) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(-0.28636006) q[0];
rz(-0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-2.5568331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6385348) q[0];
sx q[0];
rz(-0.98031822) q[0];
sx q[0];
rz(-1.9496586) q[0];
x q[1];
rz(0.61442394) q[2];
sx q[2];
rz(-1.5493869) q[2];
sx q[2];
rz(1.9911839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6905745) q[1];
sx q[1];
rz(-1.8585397) q[1];
sx q[1];
rz(2.2698127) q[1];
x q[2];
rz(0.16634511) q[3];
sx q[3];
rz(-1.5208941) q[3];
sx q[3];
rz(-2.9970616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.16453234) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(-0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(2.5928296) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(2.5351977) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8618146) q[0];
sx q[0];
rz(-1.3347515) q[0];
sx q[0];
rz(1.0311014) q[0];
x q[1];
rz(2.7658471) q[2];
sx q[2];
rz(-0.94259113) q[2];
sx q[2];
rz(1.1496161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.079271) q[1];
sx q[1];
rz(-2.3289526) q[1];
sx q[1];
rz(-2.9864242) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1612732) q[3];
sx q[3];
rz(-1.8990574) q[3];
sx q[3];
rz(-1.3271774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(2.6331804) q[2];
rz(-2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.8023119) q[0];
sx q[0];
rz(-0.34582129) q[0];
sx q[0];
rz(-0.75396496) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-2.9139013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4906625) q[0];
sx q[0];
rz(-0.41914808) q[0];
sx q[0];
rz(1.3952257) q[0];
rz(-pi) q[1];
rz(1.9367427) q[2];
sx q[2];
rz(-1.1507251) q[2];
sx q[2];
rz(-2.7811188) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55351174) q[1];
sx q[1];
rz(-1.7953824) q[1];
sx q[1];
rz(-0.54393804) q[1];
x q[2];
rz(-1.6119192) q[3];
sx q[3];
rz(-1.0121945) q[3];
sx q[3];
rz(0.55762824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(-0.83855808) q[2];
rz(-2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(-2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7446328) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(-0.80378419) q[0];
rz(2.0869758) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(1.1484336) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.046612) q[0];
sx q[0];
rz(-1.5946832) q[0];
sx q[0];
rz(-1.4987962) q[0];
rz(-pi) q[1];
rz(1.0461651) q[2];
sx q[2];
rz(-1.2534007) q[2];
sx q[2];
rz(1.5898926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.941274) q[1];
sx q[1];
rz(-0.66220821) q[1];
sx q[1];
rz(0.83076417) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35195841) q[3];
sx q[3];
rz(-1.5378012) q[3];
sx q[3];
rz(1.9381423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(-0.25319779) q[0];
rz(0.7397488) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(-2.443312) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2343276) q[0];
sx q[0];
rz(-2.8303879) q[0];
sx q[0];
rz(2.9051203) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.024784879) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(-0.63400808) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0638949) q[1];
sx q[1];
rz(-2.6297036) q[1];
sx q[1];
rz(-2.0236532) q[1];
x q[2];
rz(1.9395589) q[3];
sx q[3];
rz(-0.59951111) q[3];
sx q[3];
rz(0.60839073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11761052) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(-1.9742112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.2587851) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(-0.029126833) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(-2.6559033) q[2];
sx q[2];
rz(-0.66076856) q[2];
sx q[2];
rz(3.0271157) q[2];
rz(-2.2610353) q[3];
sx q[3];
rz(-1.6869443) q[3];
sx q[3];
rz(-1.501367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
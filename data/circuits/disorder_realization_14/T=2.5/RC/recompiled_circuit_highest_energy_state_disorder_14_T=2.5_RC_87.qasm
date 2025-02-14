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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(-2.8066714) q[0];
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(-0.60751539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4273194) q[0];
sx q[0];
rz(-1.8022984) q[0];
sx q[0];
rz(2.0718859) q[0];
x q[1];
rz(1.363475) q[2];
sx q[2];
rz(-2.6229515) q[2];
sx q[2];
rz(0.27388369) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3594397) q[1];
sx q[1];
rz(-0.50737587) q[1];
sx q[1];
rz(-2.8088267) q[1];
rz(-pi) q[2];
x q[2];
rz(2.257454) q[3];
sx q[3];
rz(-0.55301412) q[3];
sx q[3];
rz(-0.64698863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5241549) q[2];
sx q[2];
rz(-1.9257156) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(-0.080862008) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(1.9614356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2422159) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(-0.85897613) q[0];
rz(-1.2902749) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(-1.4069517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5840184) q[0];
sx q[0];
rz(-1.0913335) q[0];
sx q[0];
rz(-2.8842501) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.67703) q[2];
sx q[2];
rz(-0.33393814) q[2];
sx q[2];
rz(-0.33423697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3000112) q[1];
sx q[1];
rz(-1.661013) q[1];
sx q[1];
rz(-2.0848227) q[1];
x q[2];
rz(0.57699109) q[3];
sx q[3];
rz(-1.2170047) q[3];
sx q[3];
rz(0.64567034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0824288) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(1.3272237) q[2];
rz(-1.0446769) q[3];
sx q[3];
rz(-0.71377126) q[3];
sx q[3];
rz(-0.41675848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5663719) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(2.7161993) q[0];
rz(1.7644024) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(1.4345217) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1578428) q[0];
sx q[0];
rz(-2.2089777) q[0];
sx q[0];
rz(1.7179836) q[0];
rz(0.038952053) q[2];
sx q[2];
rz(-0.70868451) q[2];
sx q[2];
rz(1.612029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9454398) q[1];
sx q[1];
rz(-1.9827794) q[1];
sx q[1];
rz(2.4293071) q[1];
rz(-2.8355153) q[3];
sx q[3];
rz(-0.97460213) q[3];
sx q[3];
rz(-1.5706737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.44125685) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(1.8505992) q[2];
rz(-1.7839606) q[3];
sx q[3];
rz(-1.5057526) q[3];
sx q[3];
rz(-1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212379) q[0];
sx q[0];
rz(-1.0076032) q[0];
sx q[0];
rz(2.3714016) q[0];
rz(0.99984804) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(-1.8005449) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5504549) q[0];
sx q[0];
rz(-1.7831752) q[0];
sx q[0];
rz(0.55340931) q[0];
rz(-pi) q[1];
rz(1.0275339) q[2];
sx q[2];
rz(-2.3781666) q[2];
sx q[2];
rz(1.1929407) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8902258) q[1];
sx q[1];
rz(-0.40491762) q[1];
sx q[1];
rz(-0.7271073) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.099319066) q[3];
sx q[3];
rz(-1.2867974) q[3];
sx q[3];
rz(-1.8530396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3406713) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(3.0411804) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(1.8628619) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95555821) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(2.9146063) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(-1.7452128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5433301) q[0];
sx q[0];
rz(-2.6285183) q[0];
sx q[0];
rz(2.3725879) q[0];
rz(-pi) q[1];
rz(-0.41047217) q[2];
sx q[2];
rz(-1.7232401) q[2];
sx q[2];
rz(1.4727915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5310881) q[1];
sx q[1];
rz(-0.63734326) q[1];
sx q[1];
rz(0.12048851) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8245411) q[3];
sx q[3];
rz(-1.5932114) q[3];
sx q[3];
rz(-0.21460545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.22623006) q[2];
sx q[2];
rz(-1.1504983) q[2];
sx q[2];
rz(0.076233141) q[2];
rz(-1.7539615) q[3];
sx q[3];
rz(-0.97532719) q[3];
sx q[3];
rz(-1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6607894) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(-1.3355108) q[0];
rz(-2.238359) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(2.9551771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53509287) q[0];
sx q[0];
rz(-2.4267174) q[0];
sx q[0];
rz(0.8755279) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3720296) q[2];
sx q[2];
rz(-1.0093401) q[2];
sx q[2];
rz(-1.9622864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.820436) q[1];
sx q[1];
rz(-1.8119436) q[1];
sx q[1];
rz(0.92280573) q[1];
rz(-2.0669492) q[3];
sx q[3];
rz(-1.1285787) q[3];
sx q[3];
rz(-2.7520455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4621801) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(-1.3235486) q[2];
rz(1.9994252) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(-2.2511258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6776176) q[0];
sx q[0];
rz(-2.9587726) q[0];
sx q[0];
rz(-0.63419813) q[0];
rz(1.0448666) q[1];
sx q[1];
rz(-2.2506782) q[1];
sx q[1];
rz(-0.96010906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7564108) q[0];
sx q[0];
rz(-0.39968458) q[0];
sx q[0];
rz(1.5002285) q[0];
rz(-pi) q[1];
rz(-0.025770806) q[2];
sx q[2];
rz(-0.30116815) q[2];
sx q[2];
rz(2.0370551) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3297616) q[1];
sx q[1];
rz(-2.0961746) q[1];
sx q[1];
rz(-3.1282022) q[1];
rz(-0.76310254) q[3];
sx q[3];
rz(-1.0013442) q[3];
sx q[3];
rz(2.4214937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94379696) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(1.5459527) q[2];
rz(-1.4173077) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(-1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6087795) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(2.1667495) q[0];
rz(3.0335562) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.9727762) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9038531) q[0];
sx q[0];
rz(-2.0444336) q[0];
sx q[0];
rz(-2.6698378) q[0];
rz(1.5708099) q[2];
sx q[2];
rz(-0.80894366) q[2];
sx q[2];
rz(0.9521614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4296234) q[1];
sx q[1];
rz(-1.6584466) q[1];
sx q[1];
rz(-1.3369249) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5904558) q[3];
sx q[3];
rz(-1.0791313) q[3];
sx q[3];
rz(-1.8976854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10890266) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(2.4857793) q[2];
rz(3.0374895) q[3];
sx q[3];
rz(-1.7963573) q[3];
sx q[3];
rz(-2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26466894) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(1.4917829) q[0];
rz(-2.1022294) q[1];
sx q[1];
rz(-0.84016687) q[1];
sx q[1];
rz(-0.46844354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.428498) q[0];
sx q[0];
rz(-1.5665496) q[0];
sx q[0];
rz(-2.9539786) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8554815) q[2];
sx q[2];
rz(-1.9980901) q[2];
sx q[2];
rz(0.22186771) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2704308) q[1];
sx q[1];
rz(-2.3498145) q[1];
sx q[1];
rz(-1.7986078) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8599817) q[3];
sx q[3];
rz(-1.3002031) q[3];
sx q[3];
rz(0.56265807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.09482065) q[2];
sx q[2];
rz(-1.5609317) q[2];
sx q[2];
rz(-2.2136733) q[2];
rz(1.3377442) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(-1.4891589) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0874262) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(1.077865) q[0];
rz(2.7742591) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(1.8574538) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36791641) q[0];
sx q[0];
rz(-1.319029) q[0];
sx q[0];
rz(-2.0579703) q[0];
rz(-pi) q[1];
rz(-0.11576368) q[2];
sx q[2];
rz(-0.9204922) q[2];
sx q[2];
rz(-2.8011326) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2312647) q[1];
sx q[1];
rz(-2.961425) q[1];
sx q[1];
rz(2.641201) q[1];
x q[2];
rz(-1.7824836) q[3];
sx q[3];
rz(-0.84722391) q[3];
sx q[3];
rz(2.7016957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0014570634) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(2.7122811) q[2];
rz(0.15520994) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975288) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(2.2776729) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(-2.7511394) q[2];
sx q[2];
rz(-2.5627372) q[2];
sx q[2];
rz(2.550203) q[2];
rz(-0.94384296) q[3];
sx q[3];
rz(-1.2921492) q[3];
sx q[3];
rz(2.9155801) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

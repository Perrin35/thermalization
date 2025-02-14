OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72252005) q[0];
sx q[0];
rz(-0.68113405) q[0];
sx q[0];
rz(-1.0406159) q[0];
rz(2.430727) q[1];
sx q[1];
rz(3.791888) q[1];
sx q[1];
rz(11.314582) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2501564) q[0];
sx q[0];
rz(-0.52816872) q[0];
sx q[0];
rz(0.78869606) q[0];
x q[1];
rz(0.55177839) q[2];
sx q[2];
rz(-1.1349196) q[2];
sx q[2];
rz(-2.9575155) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33357269) q[1];
sx q[1];
rz(-0.90706149) q[1];
sx q[1];
rz(1.5395626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3146602) q[3];
sx q[3];
rz(-2.1029538) q[3];
sx q[3];
rz(1.102603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9616482) q[2];
sx q[2];
rz(-2.0521995) q[2];
sx q[2];
rz(-1.3570448) q[2];
rz(2.0876136) q[3];
sx q[3];
rz(-1.5318003) q[3];
sx q[3];
rz(-1.0668782) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8992952) q[0];
sx q[0];
rz(-0.59277788) q[0];
sx q[0];
rz(1.5574667) q[0];
rz(-1.8862995) q[1];
sx q[1];
rz(-0.86304945) q[1];
sx q[1];
rz(3.0994298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6145139) q[0];
sx q[0];
rz(-1.9111484) q[0];
sx q[0];
rz(0.69410364) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73017786) q[2];
sx q[2];
rz(-0.71766254) q[2];
sx q[2];
rz(2.4210986) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0706341) q[1];
sx q[1];
rz(-1.2786622) q[1];
sx q[1];
rz(2.9729384) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9412219) q[3];
sx q[3];
rz(-1.0215825) q[3];
sx q[3];
rz(2.7436243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92370799) q[2];
sx q[2];
rz(-1.8929241) q[2];
sx q[2];
rz(-2.5200747) q[2];
rz(-0.38226852) q[3];
sx q[3];
rz(-1.9999802) q[3];
sx q[3];
rz(-0.83466616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650448) q[0];
sx q[0];
rz(-1.6355729) q[0];
sx q[0];
rz(-1.3733093) q[0];
rz(0.37257591) q[1];
sx q[1];
rz(-0.5245477) q[1];
sx q[1];
rz(0.12038825) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7733369) q[0];
sx q[0];
rz(-1.0465413) q[0];
sx q[0];
rz(2.4144717) q[0];
rz(-pi) q[1];
rz(-3.1012717) q[2];
sx q[2];
rz(-1.2915269) q[2];
sx q[2];
rz(0.41353961) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1550202) q[1];
sx q[1];
rz(-0.84671686) q[1];
sx q[1];
rz(1.9968205) q[1];
rz(-pi) q[2];
rz(1.0387159) q[3];
sx q[3];
rz(-1.3339343) q[3];
sx q[3];
rz(-3.0719047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7239712) q[2];
sx q[2];
rz(-1.8310903) q[2];
sx q[2];
rz(-1.854678) q[2];
rz(-0.78479615) q[3];
sx q[3];
rz(-1.1104106) q[3];
sx q[3];
rz(2.3634461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6374636) q[0];
sx q[0];
rz(-2.6478719) q[0];
sx q[0];
rz(0.24566393) q[0];
rz(3.1371112) q[1];
sx q[1];
rz(-2.5878398) q[1];
sx q[1];
rz(0.063668879) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62718645) q[0];
sx q[0];
rz(-2.2133491) q[0];
sx q[0];
rz(-1.8802934) q[0];
rz(-pi) q[1];
rz(2.8174899) q[2];
sx q[2];
rz(-2.0162922) q[2];
sx q[2];
rz(-2.6112219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2331197) q[1];
sx q[1];
rz(-1.7571124) q[1];
sx q[1];
rz(-0.87636565) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7400209) q[3];
sx q[3];
rz(-2.662559) q[3];
sx q[3];
rz(0.019960545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8513546) q[2];
sx q[2];
rz(-1.7730224) q[2];
sx q[2];
rz(1.2151388) q[2];
rz(2.3859207) q[3];
sx q[3];
rz(-2.1561421) q[3];
sx q[3];
rz(-0.22172609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3168685) q[0];
sx q[0];
rz(-1.2687954) q[0];
sx q[0];
rz(-2.6487937) q[0];
rz(2.1700676) q[1];
sx q[1];
rz(-1.8172455) q[1];
sx q[1];
rz(2.9528101) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072666) q[0];
sx q[0];
rz(-1.2959238) q[0];
sx q[0];
rz(-1.3960394) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16021474) q[2];
sx q[2];
rz(-0.014774887) q[2];
sx q[2];
rz(-0.91195869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35313126) q[1];
sx q[1];
rz(-1.4945238) q[1];
sx q[1];
rz(2.0792873) q[1];
x q[2];
rz(0.64429342) q[3];
sx q[3];
rz(-1.3683934) q[3];
sx q[3];
rz(2.6930893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2762974) q[2];
sx q[2];
rz(-1.5027639) q[2];
sx q[2];
rz(3.0740645) q[2];
rz(1.4589795) q[3];
sx q[3];
rz(-0.71211165) q[3];
sx q[3];
rz(-2.0290831) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825298) q[0];
sx q[0];
rz(-1.9628061) q[0];
sx q[0];
rz(-1.4558526) q[0];
rz(2.9071232) q[1];
sx q[1];
rz(-0.5916943) q[1];
sx q[1];
rz(-0.04105982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9735255) q[0];
sx q[0];
rz(-2.4396908) q[0];
sx q[0];
rz(-3.0977416) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1424871) q[2];
sx q[2];
rz(-1.4465783) q[2];
sx q[2];
rz(1.6060843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.252255) q[1];
sx q[1];
rz(-0.66257157) q[1];
sx q[1];
rz(-0.8189447) q[1];
x q[2];
rz(-2.9931295) q[3];
sx q[3];
rz(-0.84772666) q[3];
sx q[3];
rz(-0.88048191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91512758) q[2];
sx q[2];
rz(-0.76639405) q[2];
sx q[2];
rz(-0.29092947) q[2];
rz(-2.9229524) q[3];
sx q[3];
rz(-2.4509957) q[3];
sx q[3];
rz(-2.2712928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9640294) q[0];
sx q[0];
rz(-2.18631) q[0];
sx q[0];
rz(-0.56354228) q[0];
rz(-0.82849416) q[1];
sx q[1];
rz(-0.68873134) q[1];
sx q[1];
rz(-0.73961893) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1391382) q[0];
sx q[0];
rz(-1.0268624) q[0];
sx q[0];
rz(-0.55940658) q[0];
rz(-pi) q[1];
rz(1.0534442) q[2];
sx q[2];
rz(-1.421325) q[2];
sx q[2];
rz(0.24565133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7259638) q[1];
sx q[1];
rz(-2.0705372) q[1];
sx q[1];
rz(2.8377297) q[1];
rz(-0.62862092) q[3];
sx q[3];
rz(-1.4328332) q[3];
sx q[3];
rz(-1.8688841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9699817) q[2];
sx q[2];
rz(-0.7526528) q[2];
sx q[2];
rz(2.8438711) q[2];
rz(3.1155078) q[3];
sx q[3];
rz(-1.1931158) q[3];
sx q[3];
rz(0.74371663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.96935) q[0];
sx q[0];
rz(-1.5765215) q[0];
sx q[0];
rz(-1.4132389) q[0];
rz(2.6353432) q[1];
sx q[1];
rz(-2.1126316) q[1];
sx q[1];
rz(-1.1594353) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0191553) q[0];
sx q[0];
rz(-1.8742689) q[0];
sx q[0];
rz(-1.6772683) q[0];
x q[1];
rz(-1.6165103) q[2];
sx q[2];
rz(-2.1518143) q[2];
sx q[2];
rz(-2.8948262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9196846) q[1];
sx q[1];
rz(-1.3338998) q[1];
sx q[1];
rz(-2.6568031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3881286) q[3];
sx q[3];
rz(-1.4770203) q[3];
sx q[3];
rz(-0.16760961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.96782288) q[2];
sx q[2];
rz(-0.92542595) q[2];
sx q[2];
rz(-1.7769163) q[2];
rz(-0.81554282) q[3];
sx q[3];
rz(-1.2929076) q[3];
sx q[3];
rz(1.7447785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5705465) q[0];
sx q[0];
rz(-1.0304352) q[0];
sx q[0];
rz(-1.9739157) q[0];
rz(0.22444621) q[1];
sx q[1];
rz(-0.75643221) q[1];
sx q[1];
rz(0.44463739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73564703) q[0];
sx q[0];
rz(-1.5155013) q[0];
sx q[0];
rz(-1.6223458) q[0];
x q[1];
rz(2.2811057) q[2];
sx q[2];
rz(-2.103946) q[2];
sx q[2];
rz(1.2459966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.81697631) q[1];
sx q[1];
rz(-1.908708) q[1];
sx q[1];
rz(-2.9959034) q[1];
rz(2.6921451) q[3];
sx q[3];
rz(-1.4420813) q[3];
sx q[3];
rz(1.4851324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.726752) q[2];
sx q[2];
rz(-1.8291992) q[2];
sx q[2];
rz(1.4078338) q[2];
rz(-1.6057181) q[3];
sx q[3];
rz(-1.9460257) q[3];
sx q[3];
rz(-2.2685952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6294412) q[0];
sx q[0];
rz(-0.81708556) q[0];
sx q[0];
rz(-1.0907115) q[0];
rz(1.0461944) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(2.2549021) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156977) q[0];
sx q[0];
rz(-0.52277589) q[0];
sx q[0];
rz(0.77976601) q[0];
x q[1];
rz(-2.1616814) q[2];
sx q[2];
rz(-1.8136029) q[2];
sx q[2];
rz(1.3473912) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7035148) q[1];
sx q[1];
rz(-1.5575927) q[1];
sx q[1];
rz(-2.7101906) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38409036) q[3];
sx q[3];
rz(-0.94058296) q[3];
sx q[3];
rz(1.5767136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1297168) q[2];
sx q[2];
rz(-2.1265714) q[2];
sx q[2];
rz(-1.8941619) q[2];
rz(-0.85793197) q[3];
sx q[3];
rz(-1.6777638) q[3];
sx q[3];
rz(1.5990545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7958551) q[0];
sx q[0];
rz(-1.416774) q[0];
sx q[0];
rz(-0.10373535) q[0];
rz(-0.70032447) q[1];
sx q[1];
rz(-1.2979957) q[1];
sx q[1];
rz(-1.8630984) q[1];
rz(-0.55194912) q[2];
sx q[2];
rz(-2.911534) q[2];
sx q[2];
rz(-0.96155675) q[2];
rz(-2.0250799) q[3];
sx q[3];
rz(-1.2496201) q[3];
sx q[3];
rz(2.7270185) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

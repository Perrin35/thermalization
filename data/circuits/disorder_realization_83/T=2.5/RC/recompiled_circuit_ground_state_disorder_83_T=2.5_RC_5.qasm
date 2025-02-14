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
rz(2.4604586) q[0];
sx q[0];
rz(13.606986) q[0];
rz(-0.71086565) q[1];
sx q[1];
rz(-0.65029538) q[1];
sx q[1];
rz(-1.8898036) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8914362) q[0];
sx q[0];
rz(-0.52816872) q[0];
sx q[0];
rz(-2.3528966) q[0];
x q[1];
rz(2.4151685) q[2];
sx q[2];
rz(-2.4527306) q[2];
sx q[2];
rz(2.3560694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7573479) q[1];
sx q[1];
rz(-2.4772344) q[1];
sx q[1];
rz(0.039907736) q[1];
rz(0.67485087) q[3];
sx q[3];
rz(-0.94776692) q[3];
sx q[3];
rz(-2.2365776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9616482) q[2];
sx q[2];
rz(-1.0893931) q[2];
sx q[2];
rz(-1.7845478) q[2];
rz(1.0539791) q[3];
sx q[3];
rz(-1.5318003) q[3];
sx q[3];
rz(-2.0747144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2422975) q[0];
sx q[0];
rz(-2.5488148) q[0];
sx q[0];
rz(-1.584126) q[0];
rz(-1.8862995) q[1];
sx q[1];
rz(-2.2785432) q[1];
sx q[1];
rz(-3.0994298) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164207) q[0];
sx q[0];
rz(-0.76043425) q[0];
sx q[0];
rz(0.50559931) q[0];
rz(-1.043528) q[2];
sx q[2];
rz(-1.0587436) q[2];
sx q[2];
rz(-1.5498424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0709585) q[1];
sx q[1];
rz(-1.8629304) q[1];
sx q[1];
rz(-0.16865428) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1290496) q[3];
sx q[3];
rz(-1.4002082) q[3];
sx q[3];
rz(-2.0743897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2178847) q[2];
sx q[2];
rz(-1.2486685) q[2];
sx q[2];
rz(-2.5200747) q[2];
rz(-0.38226852) q[3];
sx q[3];
rz(-1.9999802) q[3];
sx q[3];
rz(2.3069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650448) q[0];
sx q[0];
rz(-1.6355729) q[0];
sx q[0];
rz(1.3733093) q[0];
rz(0.37257591) q[1];
sx q[1];
rz(-0.5245477) q[1];
sx q[1];
rz(-3.0212044) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30992246) q[0];
sx q[0];
rz(-0.86747456) q[0];
sx q[0];
rz(-2.4256718) q[0];
x q[1];
rz(3.1012717) q[2];
sx q[2];
rz(-1.8500657) q[2];
sx q[2];
rz(0.41353961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2653343) q[1];
sx q[1];
rz(-1.2560532) q[1];
sx q[1];
rz(-2.3708484) q[1];
rz(-pi) q[2];
rz(1.1266842) q[3];
sx q[3];
rz(-0.57775195) q[3];
sx q[3];
rz(1.2611977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7239712) q[2];
sx q[2];
rz(-1.8310903) q[2];
sx q[2];
rz(1.2869147) q[2];
rz(-0.78479615) q[3];
sx q[3];
rz(-1.1104106) q[3];
sx q[3];
rz(-0.77814656) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50412905) q[0];
sx q[0];
rz(-2.6478719) q[0];
sx q[0];
rz(2.8959287) q[0];
rz(3.1371112) q[1];
sx q[1];
rz(-0.5537529) q[1];
sx q[1];
rz(3.0779238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.132936) q[0];
sx q[0];
rz(-1.8171165) q[0];
sx q[0];
rz(-2.4755347) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0374299) q[2];
sx q[2];
rz(-1.2793102) q[2];
sx q[2];
rz(-1.9574036) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88118756) q[1];
sx q[1];
rz(-0.71496012) q[1];
sx q[1];
rz(1.8572538) q[1];
x q[2];
rz(-1.3705092) q[3];
sx q[3];
rz(-1.1326579) q[3];
sx q[3];
rz(-2.6753257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29023805) q[2];
sx q[2];
rz(-1.7730224) q[2];
sx q[2];
rz(1.9264539) q[2];
rz(-0.75567192) q[3];
sx q[3];
rz(-2.1561421) q[3];
sx q[3];
rz(2.9198666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3168685) q[0];
sx q[0];
rz(-1.8727973) q[0];
sx q[0];
rz(-0.49279898) q[0];
rz(0.97152501) q[1];
sx q[1];
rz(-1.8172455) q[1];
sx q[1];
rz(0.18878254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072666) q[0];
sx q[0];
rz(-1.8456689) q[0];
sx q[0];
rz(-1.7455533) q[0];
x q[1];
rz(2.9813779) q[2];
sx q[2];
rz(-0.014774887) q[2];
sx q[2];
rz(-2.229634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7880612) q[1];
sx q[1];
rz(-2.6279098) q[1];
sx q[1];
rz(-1.7264926) q[1];
rz(-pi) q[2];
rz(-0.64429342) q[3];
sx q[3];
rz(-1.7731993) q[3];
sx q[3];
rz(-0.44850335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8652953) q[2];
sx q[2];
rz(-1.5027639) q[2];
sx q[2];
rz(-0.067528188) q[2];
rz(1.4589795) q[3];
sx q[3];
rz(-2.429481) q[3];
sx q[3];
rz(2.0290831) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825298) q[0];
sx q[0];
rz(-1.9628061) q[0];
sx q[0];
rz(1.68574) q[0];
rz(-2.9071232) q[1];
sx q[1];
rz(-2.5498984) q[1];
sx q[1];
rz(-0.04105982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11066786) q[0];
sx q[0];
rz(-0.86970702) q[0];
sx q[0];
rz(1.6078455) q[0];
x q[1];
rz(-2.9941999) q[2];
sx q[2];
rz(-1.0040548) q[2];
sx q[2];
rz(-3.0973377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0190553) q[1];
sx q[1];
rz(-1.1047939) q[1];
sx q[1];
rz(-2.6519777) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14846319) q[3];
sx q[3];
rz(-0.84772666) q[3];
sx q[3];
rz(-0.88048191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91512758) q[2];
sx q[2];
rz(-2.3751986) q[2];
sx q[2];
rz(-0.29092947) q[2];
rz(-2.9229524) q[3];
sx q[3];
rz(-2.4509957) q[3];
sx q[3];
rz(0.87029988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17756322) q[0];
sx q[0];
rz(-2.18631) q[0];
sx q[0];
rz(-0.56354228) q[0];
rz(0.82849416) q[1];
sx q[1];
rz(-2.4528613) q[1];
sx q[1];
rz(2.4019737) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2639573) q[0];
sx q[0];
rz(-2.3821914) q[0];
sx q[0];
rz(-2.2910221) q[0];
rz(-pi) q[1];
rz(-1.0534442) q[2];
sx q[2];
rz(-1.421325) q[2];
sx q[2];
rz(-0.24565133) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1465211) q[1];
sx q[1];
rz(-2.5634829) q[1];
sx q[1];
rz(1.0694631) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23187238) q[3];
sx q[3];
rz(-0.64157569) q[3];
sx q[3];
rz(-2.6565463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.171611) q[2];
sx q[2];
rz(-0.7526528) q[2];
sx q[2];
rz(0.29772154) q[2];
rz(0.026084829) q[3];
sx q[3];
rz(-1.9484768) q[3];
sx q[3];
rz(-2.397876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.96935) q[0];
sx q[0];
rz(-1.5650711) q[0];
sx q[0];
rz(-1.4132389) q[0];
rz(2.6353432) q[1];
sx q[1];
rz(-1.0289611) q[1];
sx q[1];
rz(-1.9821573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9205639) q[0];
sx q[0];
rz(-0.32106303) q[0];
sx q[0];
rz(-2.8144224) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.069483453) q[2];
sx q[2];
rz(-0.58260703) q[2];
sx q[2];
rz(0.16361388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9156933) q[1];
sx q[1];
rz(-2.0409313) q[1];
sx q[1];
rz(-1.8371831) q[1];
rz(-2.3881286) q[3];
sx q[3];
rz(-1.6645724) q[3];
sx q[3];
rz(-2.973983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96782288) q[2];
sx q[2];
rz(-2.2161667) q[2];
sx q[2];
rz(1.7769163) q[2];
rz(-2.3260498) q[3];
sx q[3];
rz(-1.2929076) q[3];
sx q[3];
rz(-1.7447785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5705465) q[0];
sx q[0];
rz(-2.1111574) q[0];
sx q[0];
rz(1.1676769) q[0];
rz(2.9171464) q[1];
sx q[1];
rz(-0.75643221) q[1];
sx q[1];
rz(2.6969553) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015405795) q[0];
sx q[0];
rz(-0.07557902) q[0];
sx q[0];
rz(-0.74962693) q[0];
x q[1];
rz(-2.4801587) q[2];
sx q[2];
rz(-0.97451658) q[2];
sx q[2];
rz(-3.0542947) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9079355) q[1];
sx q[1];
rz(-0.3668712) q[1];
sx q[1];
rz(-1.9625825) q[1];
rz(-pi) q[2];
rz(-0.44944758) q[3];
sx q[3];
rz(-1.6995113) q[3];
sx q[3];
rz(1.6564603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4148407) q[2];
sx q[2];
rz(-1.8291992) q[2];
sx q[2];
rz(-1.4078338) q[2];
rz(-1.5358745) q[3];
sx q[3];
rz(-1.1955669) q[3];
sx q[3];
rz(0.87299743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51215148) q[0];
sx q[0];
rz(-0.81708556) q[0];
sx q[0];
rz(1.0907115) q[0];
rz(2.0953983) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(0.88669056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156977) q[0];
sx q[0];
rz(-2.6188168) q[0];
sx q[0];
rz(-0.77976601) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8517285) q[2];
sx q[2];
rz(-0.99946195) q[2];
sx q[2];
rz(0.38331616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1613933) q[1];
sx q[1];
rz(-2.7100013) q[1];
sx q[1];
rz(3.1100247) q[1];
x q[2];
rz(-2.0453457) q[3];
sx q[3];
rz(-0.72418776) q[3];
sx q[3];
rz(0.97557035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1297168) q[2];
sx q[2];
rz(-2.1265714) q[2];
sx q[2];
rz(-1.2474308) q[2];
rz(-0.85793197) q[3];
sx q[3];
rz(-1.6777638) q[3];
sx q[3];
rz(-1.5425382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7958551) q[0];
sx q[0];
rz(-1.416774) q[0];
sx q[0];
rz(-0.10373535) q[0];
rz(0.70032447) q[1];
sx q[1];
rz(-1.843597) q[1];
sx q[1];
rz(1.2784943) q[1];
rz(-2.9447484) q[2];
sx q[2];
rz(-1.4509401) q[2];
sx q[2];
rz(0.069139253) q[2];
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

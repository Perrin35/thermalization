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
rz(2.1009768) q[0];
rz(2.430727) q[1];
sx q[1];
rz(3.791888) q[1];
sx q[1];
rz(11.314582) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2501564) q[0];
sx q[0];
rz(-2.6134239) q[0];
sx q[0];
rz(2.3528966) q[0];
x q[1];
rz(-2.5898143) q[2];
sx q[2];
rz(-1.1349196) q[2];
sx q[2];
rz(-2.9575155) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.80802) q[1];
sx q[1];
rz(-0.90706149) q[1];
sx q[1];
rz(-1.6020301) q[1];
x q[2];
rz(-0.8269324) q[3];
sx q[3];
rz(-2.1029538) q[3];
sx q[3];
rz(2.0389897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17994443) q[2];
sx q[2];
rz(-2.0521995) q[2];
sx q[2];
rz(-1.7845478) q[2];
rz(1.0539791) q[3];
sx q[3];
rz(-1.5318003) q[3];
sx q[3];
rz(1.0668782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2422975) q[0];
sx q[0];
rz(-2.5488148) q[0];
sx q[0];
rz(1.5574667) q[0];
rz(1.2552931) q[1];
sx q[1];
rz(-2.2785432) q[1];
sx q[1];
rz(-3.0994298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6145139) q[0];
sx q[0];
rz(-1.9111484) q[0];
sx q[0];
rz(2.447489) q[0];
rz(-0.57664906) q[2];
sx q[2];
rz(-1.1167142) q[2];
sx q[2];
rz(-0.25694914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54916064) q[1];
sx q[1];
rz(-1.7322408) q[1];
sx q[1];
rz(1.8668951) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8852127) q[3];
sx q[3];
rz(-0.58107983) q[3];
sx q[3];
rz(-2.3725841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2178847) q[2];
sx q[2];
rz(-1.8929241) q[2];
sx q[2];
rz(-2.5200747) q[2];
rz(0.38226852) q[3];
sx q[3];
rz(-1.9999802) q[3];
sx q[3];
rz(-2.3069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37654787) q[0];
sx q[0];
rz(-1.5060197) q[0];
sx q[0];
rz(-1.7682834) q[0];
rz(-2.7690167) q[1];
sx q[1];
rz(-0.5245477) q[1];
sx q[1];
rz(0.12038825) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7733369) q[0];
sx q[0];
rz(-1.0465413) q[0];
sx q[0];
rz(-2.4144717) q[0];
x q[1];
rz(-1.2913114) q[2];
sx q[2];
rz(-1.6095543) q[2];
sx q[2];
rz(-1.1683769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38594151) q[1];
sx q[1];
rz(-2.3214968) q[1];
sx q[1];
rz(2.7044317) q[1];
rz(-pi) q[2];
rz(-1.1266842) q[3];
sx q[3];
rz(-0.57775195) q[3];
sx q[3];
rz(1.8803949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7239712) q[2];
sx q[2];
rz(-1.8310903) q[2];
sx q[2];
rz(-1.2869147) q[2];
rz(-2.3567965) q[3];
sx q[3];
rz(-2.0311821) q[3];
sx q[3];
rz(-0.77814656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50412905) q[0];
sx q[0];
rz(-2.6478719) q[0];
sx q[0];
rz(-0.24566393) q[0];
rz(0.0044814666) q[1];
sx q[1];
rz(-0.5537529) q[1];
sx q[1];
rz(0.063668879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5144062) q[0];
sx q[0];
rz(-0.92824358) q[0];
sx q[0];
rz(1.2612993) q[0];
rz(-pi) q[1];
rz(-0.32410279) q[2];
sx q[2];
rz(-1.1253005) q[2];
sx q[2];
rz(2.6112219) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88118756) q[1];
sx q[1];
rz(-2.4266325) q[1];
sx q[1];
rz(1.2843389) q[1];
x q[2];
rz(0.40157179) q[3];
sx q[3];
rz(-0.47903362) q[3];
sx q[3];
rz(-0.019960545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8513546) q[2];
sx q[2];
rz(-1.3685702) q[2];
sx q[2];
rz(-1.2151388) q[2];
rz(0.75567192) q[3];
sx q[3];
rz(-0.98545051) q[3];
sx q[3];
rz(2.9198666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82472411) q[0];
sx q[0];
rz(-1.2687954) q[0];
sx q[0];
rz(2.6487937) q[0];
rz(-0.97152501) q[1];
sx q[1];
rz(-1.3243472) q[1];
sx q[1];
rz(0.18878254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072666) q[0];
sx q[0];
rz(-1.8456689) q[0];
sx q[0];
rz(1.3960394) q[0];
rz(2.9813779) q[2];
sx q[2];
rz(-0.014774887) q[2];
sx q[2];
rz(-2.229634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7884614) q[1];
sx q[1];
rz(-1.6470688) q[1];
sx q[1];
rz(2.0792873) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64429342) q[3];
sx q[3];
rz(-1.3683934) q[3];
sx q[3];
rz(-0.44850335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2762974) q[2];
sx q[2];
rz(-1.6388288) q[2];
sx q[2];
rz(-0.067528188) q[2];
rz(-1.6826132) q[3];
sx q[3];
rz(-0.71211165) q[3];
sx q[3];
rz(1.1125096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7825298) q[0];
sx q[0];
rz(-1.9628061) q[0];
sx q[0];
rz(-1.68574) q[0];
rz(-2.9071232) q[1];
sx q[1];
rz(-0.5916943) q[1];
sx q[1];
rz(-3.1005328) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9735255) q[0];
sx q[0];
rz(-2.4396908) q[0];
sx q[0];
rz(-3.0977416) q[0];
x q[1];
rz(0.99910555) q[2];
sx q[2];
rz(-1.4465783) q[2];
sx q[2];
rz(1.5355084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8893376) q[1];
sx q[1];
rz(-0.66257157) q[1];
sx q[1];
rz(2.322648) q[1];
rz(-2.2993574) q[3];
sx q[3];
rz(-1.6819309) q[3];
sx q[3];
rz(-0.59167093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91512758) q[2];
sx q[2];
rz(-0.76639405) q[2];
sx q[2];
rz(0.29092947) q[2];
rz(2.9229524) q[3];
sx q[3];
rz(-0.69059697) q[3];
sx q[3];
rz(-2.2712928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9640294) q[0];
sx q[0];
rz(-0.95528269) q[0];
sx q[0];
rz(-2.5780504) q[0];
rz(-0.82849416) q[1];
sx q[1];
rz(-2.4528613) q[1];
sx q[1];
rz(-2.4019737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8816895) q[0];
sx q[0];
rz(-1.0994403) q[0];
sx q[0];
rz(2.1905633) q[0];
x q[1];
rz(-2.0881485) q[2];
sx q[2];
rz(-1.421325) q[2];
sx q[2];
rz(0.24565133) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4156289) q[1];
sx q[1];
rz(-2.0705372) q[1];
sx q[1];
rz(-0.30386297) q[1];
rz(-pi) q[2];
rz(-1.4007934) q[3];
sx q[3];
rz(-2.1925266) q[3];
sx q[3];
rz(0.19843693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9699817) q[2];
sx q[2];
rz(-2.3889399) q[2];
sx q[2];
rz(-0.29772154) q[2];
rz(-3.1155078) q[3];
sx q[3];
rz(-1.9484768) q[3];
sx q[3];
rz(-2.397876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.96935) q[0];
sx q[0];
rz(-1.5650711) q[0];
sx q[0];
rz(1.7283537) q[0];
rz(-2.6353432) q[1];
sx q[1];
rz(-2.1126316) q[1];
sx q[1];
rz(1.1594353) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9205639) q[0];
sx q[0];
rz(-2.8205296) q[0];
sx q[0];
rz(2.8144224) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5250823) q[2];
sx q[2];
rz(-0.9897784) q[2];
sx q[2];
rz(2.8948262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7680507) q[1];
sx q[1];
rz(-2.6061879) q[1];
sx q[1];
rz(2.6635976) q[1];
rz(1.4425464) q[3];
sx q[3];
rz(-0.82144605) q[3];
sx q[3];
rz(1.3155703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1737698) q[2];
sx q[2];
rz(-0.92542595) q[2];
sx q[2];
rz(-1.7769163) q[2];
rz(2.3260498) q[3];
sx q[3];
rz(-1.848685) q[3];
sx q[3];
rz(1.3968141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5705465) q[0];
sx q[0];
rz(-2.1111574) q[0];
sx q[0];
rz(1.1676769) q[0];
rz(-0.22444621) q[1];
sx q[1];
rz(-0.75643221) q[1];
sx q[1];
rz(-0.44463739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015405795) q[0];
sx q[0];
rz(-0.07557902) q[0];
sx q[0];
rz(2.3919657) q[0];
x q[1];
rz(-2.3063956) q[2];
sx q[2];
rz(-2.2822501) q[2];
sx q[2];
rz(-2.2830613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81697631) q[1];
sx q[1];
rz(-1.2328846) q[1];
sx q[1];
rz(2.9959034) q[1];
rz(-pi) q[2];
rz(0.28953449) q[3];
sx q[3];
rz(-2.6752895) q[3];
sx q[3];
rz(-0.34573761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.726752) q[2];
sx q[2];
rz(-1.3123935) q[2];
sx q[2];
rz(-1.7337588) q[2];
rz(1.5358745) q[3];
sx q[3];
rz(-1.9460257) q[3];
sx q[3];
rz(-2.2685952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6294412) q[0];
sx q[0];
rz(-2.3245071) q[0];
sx q[0];
rz(1.0907115) q[0];
rz(-2.0953983) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(2.2549021) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7162298) q[0];
sx q[0];
rz(-1.2078478) q[0];
sx q[0];
rz(1.9557529) q[0];
rz(-pi) q[1];
rz(-0.97991122) q[2];
sx q[2];
rz(-1.8136029) q[2];
sx q[2];
rz(1.7942015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.43807784) q[1];
sx q[1];
rz(-1.5575927) q[1];
sx q[1];
rz(-2.7101906) q[1];
x q[2];
rz(-0.38409036) q[3];
sx q[3];
rz(-0.94058296) q[3];
sx q[3];
rz(-1.5767136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0118759) q[2];
sx q[2];
rz(-1.0150212) q[2];
sx q[2];
rz(-1.2474308) q[2];
rz(0.85793197) q[3];
sx q[3];
rz(-1.6777638) q[3];
sx q[3];
rz(1.5425382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7958551) q[0];
sx q[0];
rz(-1.7248187) q[0];
sx q[0];
rz(3.0378573) q[0];
rz(-2.4412682) q[1];
sx q[1];
rz(-1.843597) q[1];
sx q[1];
rz(1.2784943) q[1];
rz(2.5896435) q[2];
sx q[2];
rz(-2.911534) q[2];
sx q[2];
rz(-0.96155675) q[2];
rz(-0.35459749) q[3];
sx q[3];
rz(-2.000256) q[3];
sx q[3];
rz(1.3091814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

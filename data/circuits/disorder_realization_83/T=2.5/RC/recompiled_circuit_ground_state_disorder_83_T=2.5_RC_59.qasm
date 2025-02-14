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
rz(-2.7530818) q[0];
sx q[0];
rz(-1.2076993) q[0];
sx q[0];
rz(-1.9632505) q[0];
rz(-2.5898143) q[2];
sx q[2];
rz(-2.006673) q[2];
sx q[2];
rz(2.9575155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.80802) q[1];
sx q[1];
rz(-2.2345312) q[1];
sx q[1];
rz(-1.5395626) q[1];
x q[2];
rz(0.67485087) q[3];
sx q[3];
rz(-2.1938257) q[3];
sx q[3];
rz(-0.90501508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.17994443) q[2];
sx q[2];
rz(-1.0893931) q[2];
sx q[2];
rz(-1.3570448) q[2];
rz(2.0876136) q[3];
sx q[3];
rz(-1.6097924) q[3];
sx q[3];
rz(-2.0747144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.8992952) q[0];
sx q[0];
rz(-0.59277788) q[0];
sx q[0];
rz(-1.584126) q[0];
rz(-1.2552931) q[1];
sx q[1];
rz(-2.2785432) q[1];
sx q[1];
rz(-0.042162808) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3688646) q[0];
sx q[0];
rz(-2.2180551) q[0];
sx q[0];
rz(-2.002541) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0980646) q[2];
sx q[2];
rz(-1.0587436) q[2];
sx q[2];
rz(1.5498424) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.592432) q[1];
sx q[1];
rz(-1.4093519) q[1];
sx q[1];
rz(-1.2746975) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0125431) q[3];
sx q[3];
rz(-1.7413845) q[3];
sx q[3];
rz(2.0743897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2178847) q[2];
sx q[2];
rz(-1.8929241) q[2];
sx q[2];
rz(-0.62151796) q[2];
rz(0.38226852) q[3];
sx q[3];
rz(-1.9999802) q[3];
sx q[3];
rz(0.83466616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650448) q[0];
sx q[0];
rz(-1.6355729) q[0];
sx q[0];
rz(1.7682834) q[0];
rz(-2.7690167) q[1];
sx q[1];
rz(-0.5245477) q[1];
sx q[1];
rz(-3.0212044) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7733369) q[0];
sx q[0];
rz(-2.0950514) q[0];
sx q[0];
rz(0.7271209) q[0];
rz(-pi) q[1];
rz(-1.7104501) q[2];
sx q[2];
rz(-2.8595028) q[2];
sx q[2];
rz(2.8733746) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87625831) q[1];
sx q[1];
rz(-1.2560532) q[1];
sx q[1];
rz(0.77074428) q[1];
rz(-2.0149085) q[3];
sx q[3];
rz(-0.57775195) q[3];
sx q[3];
rz(1.2611977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7239712) q[2];
sx q[2];
rz(-1.8310903) q[2];
sx q[2];
rz(1.854678) q[2];
rz(2.3567965) q[3];
sx q[3];
rz(-1.1104106) q[3];
sx q[3];
rz(2.3634461) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50412905) q[0];
sx q[0];
rz(-2.6478719) q[0];
sx q[0];
rz(0.24566393) q[0];
rz(-3.1371112) q[1];
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
rz(2.0086567) q[0];
sx q[0];
rz(-1.3244761) q[0];
sx q[0];
rz(2.4755347) q[0];
x q[1];
rz(1.1041627) q[2];
sx q[2];
rz(-1.8622824) q[2];
sx q[2];
rz(-1.9574036) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2331197) q[1];
sx q[1];
rz(-1.7571124) q[1];
sx q[1];
rz(2.265227) q[1];
x q[2];
rz(-0.44594619) q[3];
sx q[3];
rz(-1.7519439) q[3];
sx q[3];
rz(-1.9511478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8513546) q[2];
sx q[2];
rz(-1.7730224) q[2];
sx q[2];
rz(1.2151388) q[2];
rz(0.75567192) q[3];
sx q[3];
rz(-0.98545051) q[3];
sx q[3];
rz(-0.22172609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82472411) q[0];
sx q[0];
rz(-1.2687954) q[0];
sx q[0];
rz(-2.6487937) q[0];
rz(0.97152501) q[1];
sx q[1];
rz(-1.3243472) q[1];
sx q[1];
rz(-0.18878254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59072666) q[0];
sx q[0];
rz(-1.2959238) q[0];
sx q[0];
rz(-1.3960394) q[0];
rz(-pi) q[1];
rz(-0.16021474) q[2];
sx q[2];
rz(-0.014774887) q[2];
sx q[2];
rz(0.91195869) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1752191) q[1];
sx q[1];
rz(-2.0776675) q[1];
sx q[1];
rz(0.08726774) q[1];
rz(2.8123662) q[3];
sx q[3];
rz(-2.4705973) q[3];
sx q[3];
rz(-1.7578357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8652953) q[2];
sx q[2];
rz(-1.5027639) q[2];
sx q[2];
rz(-0.067528188) q[2];
rz(1.4589795) q[3];
sx q[3];
rz(-0.71211165) q[3];
sx q[3];
rz(-2.0290831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825298) q[0];
sx q[0];
rz(-1.9628061) q[0];
sx q[0];
rz(-1.4558526) q[0];
rz(-2.9071232) q[1];
sx q[1];
rz(-0.5916943) q[1];
sx q[1];
rz(-3.1005328) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9735255) q[0];
sx q[0];
rz(-0.70190185) q[0];
sx q[0];
rz(0.043851093) q[0];
x q[1];
rz(-2.1424871) q[2];
sx q[2];
rz(-1.6950144) q[2];
sx q[2];
rz(1.6060843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1225374) q[1];
sx q[1];
rz(-2.0367987) q[1];
sx q[1];
rz(2.6519777) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14846319) q[3];
sx q[3];
rz(-2.293866) q[3];
sx q[3];
rz(-2.2611107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91512758) q[2];
sx q[2];
rz(-0.76639405) q[2];
sx q[2];
rz(-2.8506632) q[2];
rz(-0.21864024) q[3];
sx q[3];
rz(-2.4509957) q[3];
sx q[3];
rz(2.2712928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9640294) q[0];
sx q[0];
rz(-2.18631) q[0];
sx q[0];
rz(0.56354228) q[0];
rz(0.82849416) q[1];
sx q[1];
rz(-2.4528613) q[1];
sx q[1];
rz(2.4019737) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1391382) q[0];
sx q[0];
rz(-1.0268624) q[0];
sx q[0];
rz(-2.5821861) q[0];
rz(-pi) q[1];
x q[1];
rz(2.970026) q[2];
sx q[2];
rz(-2.0818149) q[2];
sx q[2];
rz(-1.2406065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7259638) q[1];
sx q[1];
rz(-2.0705372) q[1];
sx q[1];
rz(-2.8377297) q[1];
x q[2];
rz(0.23187238) q[3];
sx q[3];
rz(-0.64157569) q[3];
sx q[3];
rz(-2.6565463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9699817) q[2];
sx q[2];
rz(-2.3889399) q[2];
sx q[2];
rz(0.29772154) q[2];
rz(0.026084829) q[3];
sx q[3];
rz(-1.1931158) q[3];
sx q[3];
rz(2.397876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.96935) q[0];
sx q[0];
rz(-1.5650711) q[0];
sx q[0];
rz(1.4132389) q[0];
rz(-0.50624943) q[1];
sx q[1];
rz(-2.1126316) q[1];
sx q[1];
rz(1.9821573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12243733) q[0];
sx q[0];
rz(-1.8742689) q[0];
sx q[0];
rz(-1.4643244) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5250823) q[2];
sx q[2];
rz(-2.1518143) q[2];
sx q[2];
rz(-0.24676649) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2258994) q[1];
sx q[1];
rz(-2.0409313) q[1];
sx q[1];
rz(-1.8371831) q[1];
x q[2];
rz(1.4425464) q[3];
sx q[3];
rz(-2.3201466) q[3];
sx q[3];
rz(-1.3155703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1737698) q[2];
sx q[2];
rz(-0.92542595) q[2];
sx q[2];
rz(-1.7769163) q[2];
rz(-2.3260498) q[3];
sx q[3];
rz(-1.2929076) q[3];
sx q[3];
rz(1.3968141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-2.3851604) q[1];
sx q[1];
rz(-0.44463739) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035919) q[0];
sx q[0];
rz(-1.5193257) q[0];
sx q[0];
rz(-0.055368467) q[0];
rz(-2.3063956) q[2];
sx q[2];
rz(-2.2822501) q[2];
sx q[2];
rz(0.85853133) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4363775) q[1];
sx q[1];
rz(-1.4333998) q[1];
sx q[1];
rz(-1.912053) q[1];
x q[2];
rz(2.6921451) q[3];
sx q[3];
rz(-1.6995113) q[3];
sx q[3];
rz(1.6564603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4148407) q[2];
sx q[2];
rz(-1.3123935) q[2];
sx q[2];
rz(1.4078338) q[2];
rz(1.6057181) q[3];
sx q[3];
rz(-1.1955669) q[3];
sx q[3];
rz(-2.2685952) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6294412) q[0];
sx q[0];
rz(-2.3245071) q[0];
sx q[0];
rz(2.0508811) q[0];
rz(-1.0461944) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(-2.2549021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7156977) q[0];
sx q[0];
rz(-2.6188168) q[0];
sx q[0];
rz(-2.3618266) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28986411) q[2];
sx q[2];
rz(-2.1421307) q[2];
sx q[2];
rz(0.38331616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0149518) q[1];
sx q[1];
rz(-2.0021582) q[1];
sx q[1];
rz(-1.5562612) q[1];
rz(1.0962469) q[3];
sx q[3];
rz(-2.4174049) q[3];
sx q[3];
rz(2.1660223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0118759) q[2];
sx q[2];
rz(-2.1265714) q[2];
sx q[2];
rz(1.8941619) q[2];
rz(0.85793197) q[3];
sx q[3];
rz(-1.4638289) q[3];
sx q[3];
rz(-1.5425382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7958551) q[0];
sx q[0];
rz(-1.7248187) q[0];
sx q[0];
rz(3.0378573) q[0];
rz(-0.70032447) q[1];
sx q[1];
rz(-1.2979957) q[1];
sx q[1];
rz(-1.8630984) q[1];
rz(2.5896435) q[2];
sx q[2];
rz(-2.911534) q[2];
sx q[2];
rz(-0.96155675) q[2];
rz(2.7869952) q[3];
sx q[3];
rz(-2.000256) q[3];
sx q[3];
rz(1.3091814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

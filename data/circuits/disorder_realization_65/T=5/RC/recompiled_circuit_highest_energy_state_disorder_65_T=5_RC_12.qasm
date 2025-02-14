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
rz(-1.3574358) q[0];
sx q[0];
rz(-2.5166002) q[0];
sx q[0];
rz(-1.5224737) q[0];
rz(1.1433262) q[1];
sx q[1];
rz(3.7781236) q[1];
sx q[1];
rz(14.349024) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56464897) q[0];
sx q[0];
rz(-1.2822106) q[0];
sx q[0];
rz(2.1493069) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9404561) q[2];
sx q[2];
rz(-1.8471175) q[2];
sx q[2];
rz(-1.2075961) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3431817) q[1];
sx q[1];
rz(-2.792387) q[1];
sx q[1];
rz(-0.83425547) q[1];
rz(-0.10349689) q[3];
sx q[3];
rz(-2.0465133) q[3];
sx q[3];
rz(-0.58799839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6978567) q[2];
sx q[2];
rz(-2.3190658) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(-0.13992986) q[3];
sx q[3];
rz(-1.5014476) q[3];
sx q[3];
rz(-0.073277624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78854617) q[0];
sx q[0];
rz(-1.5209501) q[0];
sx q[0];
rz(-2.7781558) q[0];
rz(2.4711171) q[1];
sx q[1];
rz(-1.8164219) q[1];
sx q[1];
rz(2.061981) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99936679) q[0];
sx q[0];
rz(-1.3329263) q[0];
sx q[0];
rz(-2.1145942) q[0];
x q[1];
rz(-2.1035467) q[2];
sx q[2];
rz(-1.4584625) q[2];
sx q[2];
rz(-2.7137604) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9448589) q[1];
sx q[1];
rz(-2.8962039) q[1];
sx q[1];
rz(-1.3887482) q[1];
x q[2];
rz(-0.3877181) q[3];
sx q[3];
rz(-0.61655515) q[3];
sx q[3];
rz(-1.8247557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4662027) q[2];
sx q[2];
rz(-0.73489302) q[2];
sx q[2];
rz(-2.7030763) q[2];
rz(-1.0112666) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(1.6605759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.6587875) q[0];
sx q[0];
rz(-2.6130982) q[0];
sx q[0];
rz(-0.82861376) q[0];
rz(1.2476745) q[1];
sx q[1];
rz(-1.9745461) q[1];
sx q[1];
rz(-0.25965986) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4858711) q[0];
sx q[0];
rz(-1.8410793) q[0];
sx q[0];
rz(2.2673658) q[0];
x q[1];
rz(2.7607418) q[2];
sx q[2];
rz(-1.9162188) q[2];
sx q[2];
rz(-0.38174111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34777113) q[1];
sx q[1];
rz(-2.035537) q[1];
sx q[1];
rz(2.3427367) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91783701) q[3];
sx q[3];
rz(-1.3363095) q[3];
sx q[3];
rz(3.0378436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.20716771) q[2];
sx q[2];
rz(-1.577689) q[2];
sx q[2];
rz(-0.820532) q[2];
rz(-2.3151243) q[3];
sx q[3];
rz(-2.1491437) q[3];
sx q[3];
rz(-1.3124527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25927037) q[0];
sx q[0];
rz(-2.3311908) q[0];
sx q[0];
rz(-2.209254) q[0];
rz(-1.8238292) q[1];
sx q[1];
rz(-1.980314) q[1];
sx q[1];
rz(0.12277776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2185635) q[0];
sx q[0];
rz(-1.6549131) q[0];
sx q[0];
rz(0.53972678) q[0];
rz(-2.8343646) q[2];
sx q[2];
rz(-2.3005552) q[2];
sx q[2];
rz(-0.73452362) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7386094) q[1];
sx q[1];
rz(-1.5697641) q[1];
sx q[1];
rz(2.3139075) q[1];
rz(-pi) q[2];
rz(-1.8639542) q[3];
sx q[3];
rz(-2.1159626) q[3];
sx q[3];
rz(2.4923862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5968898) q[2];
sx q[2];
rz(-2.2201316) q[2];
sx q[2];
rz(-0.60453647) q[2];
rz(0.70945159) q[3];
sx q[3];
rz(-0.76990288) q[3];
sx q[3];
rz(0.82367045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9619047) q[0];
sx q[0];
rz(-3.0178495) q[0];
sx q[0];
rz(-1.3667579) q[0];
rz(-2.1429515) q[1];
sx q[1];
rz(-0.81592453) q[1];
sx q[1];
rz(-2.6008115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75434573) q[0];
sx q[0];
rz(-1.3636813) q[0];
sx q[0];
rz(-1.7133173) q[0];
rz(-2.3983058) q[2];
sx q[2];
rz(-1.0096482) q[2];
sx q[2];
rz(2.4744792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97637002) q[1];
sx q[1];
rz(-1.762855) q[1];
sx q[1];
rz(2.8403175) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28022061) q[3];
sx q[3];
rz(-0.58455211) q[3];
sx q[3];
rz(1.9909797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46405408) q[2];
sx q[2];
rz(-2.3276734) q[2];
sx q[2];
rz(0.82826725) q[2];
rz(1.6895435) q[3];
sx q[3];
rz(-1.4938846) q[3];
sx q[3];
rz(-2.233708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.117711) q[0];
sx q[0];
rz(-1.9458867) q[0];
sx q[0];
rz(-1.9697795) q[0];
rz(-0.68012971) q[1];
sx q[1];
rz(-1.9155733) q[1];
sx q[1];
rz(-0.5562869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85324641) q[0];
sx q[0];
rz(-2.6814406) q[0];
sx q[0];
rz(-0.96921885) q[0];
rz(-3.1388144) q[2];
sx q[2];
rz(-1.4242715) q[2];
sx q[2];
rz(0.79039449) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43832512) q[1];
sx q[1];
rz(-2.5653953) q[1];
sx q[1];
rz(-3.0238815) q[1];
rz(2.5285491) q[3];
sx q[3];
rz(-2.4073232) q[3];
sx q[3];
rz(-3.120852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3351626) q[2];
sx q[2];
rz(-1.8589636) q[2];
sx q[2];
rz(-0.20509091) q[2];
rz(-0.80275503) q[3];
sx q[3];
rz(-2.0358678) q[3];
sx q[3];
rz(-2.5870489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71451181) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(-0.22739534) q[0];
rz(0.79044509) q[1];
sx q[1];
rz(-0.77722725) q[1];
sx q[1];
rz(2.3048185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1632501) q[0];
sx q[0];
rz(-2.0866835) q[0];
sx q[0];
rz(1.31869) q[0];
rz(-pi) q[1];
rz(-2.5535832) q[2];
sx q[2];
rz(-1.8711897) q[2];
sx q[2];
rz(-0.38640768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10936508) q[1];
sx q[1];
rz(-1.1676452) q[1];
sx q[1];
rz(-1.9126066) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9430591) q[3];
sx q[3];
rz(-2.296631) q[3];
sx q[3];
rz(-1.5133295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5922015) q[2];
sx q[2];
rz(-1.4791919) q[2];
sx q[2];
rz(-0.57360348) q[2];
rz(-0.32522374) q[3];
sx q[3];
rz(-2.2461788) q[3];
sx q[3];
rz(2.7371244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0322872) q[0];
sx q[0];
rz(-2.9635297) q[0];
sx q[0];
rz(2.4626515) q[0];
rz(0.13488723) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(1.4237684) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7622691) q[0];
sx q[0];
rz(-2.0998619) q[0];
sx q[0];
rz(-0.53405116) q[0];
rz(-pi) q[1];
rz(0.94332327) q[2];
sx q[2];
rz(-2.3296142) q[2];
sx q[2];
rz(1.2710149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1011185) q[1];
sx q[1];
rz(-2.3777739) q[1];
sx q[1];
rz(-1.5274672) q[1];
rz(-1.0342977) q[3];
sx q[3];
rz(-0.52742672) q[3];
sx q[3];
rz(-2.3303383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5474825) q[2];
sx q[2];
rz(-0.91071931) q[2];
sx q[2];
rz(2.9098848) q[2];
rz(-2.1314651) q[3];
sx q[3];
rz(-1.2332656) q[3];
sx q[3];
rz(0.89099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6758839) q[0];
sx q[0];
rz(-0.75513419) q[0];
sx q[0];
rz(2.6278507) q[0];
rz(-1.4328009) q[1];
sx q[1];
rz(-1.865973) q[1];
sx q[1];
rz(1.6339711) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20376539) q[0];
sx q[0];
rz(-1.7047593) q[0];
sx q[0];
rz(-1.521827) q[0];
x q[1];
rz(-2.5812923) q[2];
sx q[2];
rz(-2.8447084) q[2];
sx q[2];
rz(-0.96499824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0158314) q[1];
sx q[1];
rz(-1.4471701) q[1];
sx q[1];
rz(-3.0112253) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1166508) q[3];
sx q[3];
rz(-1.8224622) q[3];
sx q[3];
rz(-0.97398678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12019176) q[2];
sx q[2];
rz(-2.2970565) q[2];
sx q[2];
rz(0.25944844) q[2];
rz(-1.9868959) q[3];
sx q[3];
rz(-0.76609937) q[3];
sx q[3];
rz(-1.2999138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65086377) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(1.4890626) q[0];
rz(2.5015855) q[1];
sx q[1];
rz(-1.0970486) q[1];
sx q[1];
rz(0.99001137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31371597) q[0];
sx q[0];
rz(-1.4826688) q[0];
sx q[0];
rz(-0.14493305) q[0];
x q[1];
rz(-0.76837627) q[2];
sx q[2];
rz(-1.9208761) q[2];
sx q[2];
rz(-2.1742252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1460625) q[1];
sx q[1];
rz(-1.7319253) q[1];
sx q[1];
rz(1.1481768) q[1];
x q[2];
rz(0.02496533) q[3];
sx q[3];
rz(-1.348266) q[3];
sx q[3];
rz(2.3937267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9330357) q[2];
sx q[2];
rz(-1.2988337) q[2];
sx q[2];
rz(0.2253069) q[2];
rz(-0.92132583) q[3];
sx q[3];
rz(-0.39042979) q[3];
sx q[3];
rz(-0.050617378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0113572) q[0];
sx q[0];
rz(-0.36778944) q[0];
sx q[0];
rz(0.854048) q[0];
rz(-1.8232952) q[1];
sx q[1];
rz(-1.5250991) q[1];
sx q[1];
rz(-0.93223882) q[1];
rz(1.9226045) q[2];
sx q[2];
rz(-0.27588898) q[2];
sx q[2];
rz(0.5994748) q[2];
rz(-1.2750219) q[3];
sx q[3];
rz(-0.56850962) q[3];
sx q[3];
rz(-3.1103984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

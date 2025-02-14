OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0140822) q[0];
sx q[0];
rz(-0.2427225) q[0];
sx q[0];
rz(2.2665562) q[0];
rz(-1.8889282) q[1];
sx q[1];
rz(-1.6369605) q[1];
sx q[1];
rz(2.5445282) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93347862) q[0];
sx q[0];
rz(-2.7467318) q[0];
sx q[0];
rz(-1.3482679) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4389285) q[2];
sx q[2];
rz(-1.4070373) q[2];
sx q[2];
rz(0.57715774) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7272612) q[1];
sx q[1];
rz(-0.67318577) q[1];
sx q[1];
rz(-0.080869599) q[1];
x q[2];
rz(0.050757082) q[3];
sx q[3];
rz(-1.4043706) q[3];
sx q[3];
rz(0.084648477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.732932) q[2];
sx q[2];
rz(-0.084736846) q[2];
sx q[2];
rz(-1.6049467) q[2];
rz(1.599865) q[3];
sx q[3];
rz(-2.7432975) q[3];
sx q[3];
rz(-2.1275684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9059471) q[0];
sx q[0];
rz(-2.5962317) q[0];
sx q[0];
rz(1.9769309) q[0];
rz(-2.2230478) q[1];
sx q[1];
rz(-0.48857498) q[1];
sx q[1];
rz(2.4426544) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0889269) q[0];
sx q[0];
rz(-0.97481004) q[0];
sx q[0];
rz(1.2794897) q[0];
rz(0.87911682) q[2];
sx q[2];
rz(-1.7335636) q[2];
sx q[2];
rz(0.04479822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1876549) q[1];
sx q[1];
rz(-1.6587757) q[1];
sx q[1];
rz(1.6768558) q[1];
x q[2];
rz(1.5941991) q[3];
sx q[3];
rz(-1.7197945) q[3];
sx q[3];
rz(3.0210173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32236448) q[2];
sx q[2];
rz(-1.4388204) q[2];
sx q[2];
rz(1.9361852) q[2];
rz(2.9820005) q[3];
sx q[3];
rz(-1.7528088) q[3];
sx q[3];
rz(3.0157183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50652248) q[0];
sx q[0];
rz(-0.078190088) q[0];
sx q[0];
rz(-3.0248094) q[0];
rz(1.1306688) q[1];
sx q[1];
rz(-3.0394381) q[1];
sx q[1];
rz(-3.0932313) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.48232) q[0];
sx q[0];
rz(-3.0249615) q[0];
sx q[0];
rz(-0.92088033) q[0];
rz(-pi) q[1];
rz(3.02028) q[2];
sx q[2];
rz(-1.5229791) q[2];
sx q[2];
rz(-1.3791093) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7680184) q[1];
sx q[1];
rz(-1.2559203) q[1];
sx q[1];
rz(-0.86656226) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4428407) q[3];
sx q[3];
rz(-1.5317917) q[3];
sx q[3];
rz(-1.4641264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2241263) q[2];
sx q[2];
rz(-1.5175061) q[2];
sx q[2];
rz(-0.65829128) q[2];
rz(1.0607464) q[3];
sx q[3];
rz(-2.3297533) q[3];
sx q[3];
rz(0.15061024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0105522) q[0];
sx q[0];
rz(-0.083715938) q[0];
sx q[0];
rz(0.87378275) q[0];
rz(-2.1281758) q[1];
sx q[1];
rz(-0.0031009379) q[1];
sx q[1];
rz(-0.82469624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31451685) q[0];
sx q[0];
rz(-1.6028048) q[0];
sx q[0];
rz(1.5289983) q[0];
x q[1];
rz(1.0295139) q[2];
sx q[2];
rz(-0.80179399) q[2];
sx q[2];
rz(3.1355592) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5910278) q[1];
sx q[1];
rz(-0.30375215) q[1];
sx q[1];
rz(-1.135056) q[1];
rz(-pi) q[2];
rz(-0.6225067) q[3];
sx q[3];
rz(-1.1460163) q[3];
sx q[3];
rz(0.82889601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.970845) q[2];
sx q[2];
rz(-2.8837995) q[2];
sx q[2];
rz(0.26138678) q[2];
rz(1.9445253) q[3];
sx q[3];
rz(-1.7944929) q[3];
sx q[3];
rz(2.4813467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18545397) q[0];
sx q[0];
rz(-3.0281797) q[0];
sx q[0];
rz(1.0979106) q[0];
rz(-0.55770355) q[1];
sx q[1];
rz(-3.1126366) q[1];
sx q[1];
rz(1.2835693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4056617) q[0];
sx q[0];
rz(-0.064583555) q[0];
sx q[0];
rz(0.92900999) q[0];
x q[1];
rz(0.26153841) q[2];
sx q[2];
rz(-1.4976302) q[2];
sx q[2];
rz(1.9436702) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9030212) q[1];
sx q[1];
rz(-2.0687177) q[1];
sx q[1];
rz(1.7037039) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5398272) q[3];
sx q[3];
rz(-0.93671173) q[3];
sx q[3];
rz(-0.82179797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4881318) q[2];
sx q[2];
rz(-0.9119432) q[2];
sx q[2];
rz(0.27257356) q[2];
rz(-1.163698) q[3];
sx q[3];
rz(-1.3397763) q[3];
sx q[3];
rz(0.94289601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1777451) q[0];
sx q[0];
rz(-3.0680532) q[0];
sx q[0];
rz(-2.9195926) q[0];
rz(-1.47413) q[1];
sx q[1];
rz(-0.024802955) q[1];
sx q[1];
rz(-2.7892392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2896293) q[0];
sx q[0];
rz(-1.5709988) q[0];
sx q[0];
rz(1.5727497) q[0];
rz(-pi) q[1];
rz(0.9325005) q[2];
sx q[2];
rz(-0.83915448) q[2];
sx q[2];
rz(2.7256903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2655991) q[1];
sx q[1];
rz(-0.76294661) q[1];
sx q[1];
rz(-0.89150064) q[1];
rz(1.1389987) q[3];
sx q[3];
rz(-0.96240265) q[3];
sx q[3];
rz(2.268057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12677975) q[2];
sx q[2];
rz(-1.1819906) q[2];
sx q[2];
rz(-2.354055) q[2];
rz(-3.0102503) q[3];
sx q[3];
rz(-2.1613439) q[3];
sx q[3];
rz(-2.378715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545559) q[0];
sx q[0];
rz(-0.035722345) q[0];
sx q[0];
rz(0.63294739) q[0];
rz(-2.5542906) q[1];
sx q[1];
rz(-0.018922806) q[1];
sx q[1];
rz(2.6515554) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49289277) q[0];
sx q[0];
rz(-2.9577575) q[0];
sx q[0];
rz(2.0931758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7964726) q[2];
sx q[2];
rz(-1.39087) q[2];
sx q[2];
rz(-2.3210482) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8146805) q[1];
sx q[1];
rz(-2.1974954) q[1];
sx q[1];
rz(1.9912849) q[1];
x q[2];
rz(2.1774749) q[3];
sx q[3];
rz(-2.5593966) q[3];
sx q[3];
rz(-2.4011998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5309165) q[2];
sx q[2];
rz(-2.5848415) q[2];
sx q[2];
rz(2.4497633) q[2];
rz(1.0308712) q[3];
sx q[3];
rz(-2.2014047) q[3];
sx q[3];
rz(-2.809439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0523025) q[0];
sx q[0];
rz(-0.060929935) q[0];
sx q[0];
rz(2.6887509) q[0];
rz(-0.37372681) q[1];
sx q[1];
rz(-0.12416298) q[1];
sx q[1];
rz(-2.9492818) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7953291) q[0];
sx q[0];
rz(-1.6843616) q[0];
sx q[0];
rz(1.5862443) q[0];
x q[1];
rz(1.1636803) q[2];
sx q[2];
rz(-1.0958899) q[2];
sx q[2];
rz(2.6559336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2365742) q[1];
sx q[1];
rz(-0.55779167) q[1];
sx q[1];
rz(1.94005) q[1];
x q[2];
rz(1.2471036) q[3];
sx q[3];
rz(-1.7028042) q[3];
sx q[3];
rz(-2.4292601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.33334357) q[2];
sx q[2];
rz(-3.0875751) q[2];
sx q[2];
rz(-0.56600189) q[2];
rz(2.443215) q[3];
sx q[3];
rz(-1.799182) q[3];
sx q[3];
rz(-0.061605569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1762539) q[0];
sx q[0];
rz(-2.7955671) q[0];
sx q[0];
rz(1.6842496) q[0];
rz(0.95218843) q[1];
sx q[1];
rz(-0.0038650611) q[1];
sx q[1];
rz(1.7239408) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3838634) q[0];
sx q[0];
rz(-1.8762021) q[0];
sx q[0];
rz(-0.1975493) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.023073828) q[2];
sx q[2];
rz(-1.0468113) q[2];
sx q[2];
rz(0.11701458) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7240282) q[1];
sx q[1];
rz(-1.7463435) q[1];
sx q[1];
rz(-1.6035622) q[1];
x q[2];
rz(-1.7030956) q[3];
sx q[3];
rz(-1.5217363) q[3];
sx q[3];
rz(0.2459615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5119322) q[2];
sx q[2];
rz(-2.8544482) q[2];
sx q[2];
rz(0.49893898) q[2];
rz(0.75368369) q[3];
sx q[3];
rz(-2.6028809) q[3];
sx q[3];
rz(0.57873571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8653487) q[0];
sx q[0];
rz(-2.8911599) q[0];
sx q[0];
rz(0.28755406) q[0];
rz(0.7175256) q[1];
sx q[1];
rz(-0.10401195) q[1];
sx q[1];
rz(-1.2984553) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13253015) q[0];
sx q[0];
rz(-2.1611737) q[0];
sx q[0];
rz(-0.68254614) q[0];
rz(-2.1537142) q[2];
sx q[2];
rz(-2.0955502) q[2];
sx q[2];
rz(-0.74452213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76261452) q[1];
sx q[1];
rz(-1.839338) q[1];
sx q[1];
rz(1.9767799) q[1];
x q[2];
rz(1.6631546) q[3];
sx q[3];
rz(-1.256681) q[3];
sx q[3];
rz(0.11387728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31460497) q[2];
sx q[2];
rz(-1.8525476) q[2];
sx q[2];
rz(-2.2513466) q[2];
rz(-3.1173949) q[3];
sx q[3];
rz(-0.20166339) q[3];
sx q[3];
rz(-2.3307687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.191962) q[0];
sx q[0];
rz(-2.5856615) q[0];
sx q[0];
rz(-0.91188201) q[0];
rz(0.98056071) q[1];
sx q[1];
rz(-2.801827) q[1];
sx q[1];
rz(-2.7580072) q[1];
rz(-1.7928852) q[2];
sx q[2];
rz(-2.8206456) q[2];
sx q[2];
rz(-1.7802466) q[2];
rz(-1.8611363) q[3];
sx q[3];
rz(-1.3383337) q[3];
sx q[3];
rz(0.70575502) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

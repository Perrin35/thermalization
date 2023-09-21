OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(4.8882422) q[0];
sx q[0];
rz(12.56765) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2389195) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(2.8047049) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84471976) q[2];
sx q[2];
rz(-1.8978999) q[2];
sx q[2];
rz(2.6127882) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85841132) q[1];
sx q[1];
rz(-1.5309257) q[1];
sx q[1];
rz(1.4285425) q[1];
rz(2.6585456) q[3];
sx q[3];
rz(-2.8222198) q[3];
sx q[3];
rz(-1.3787624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.8700245) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.7659448) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028582024) q[0];
sx q[0];
rz(-0.41886371) q[0];
sx q[0];
rz(0.21582614) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8091082) q[2];
sx q[2];
rz(-1.6147699) q[2];
sx q[2];
rz(0.51648742) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.980629) q[1];
sx q[1];
rz(-1.7501608) q[1];
sx q[1];
rz(0.37100002) q[1];
rz(-pi) q[2];
rz(0.35956412) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-2.0294702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-2.6518872) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.906357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2429758) q[0];
sx q[0];
rz(-1.1968062) q[0];
sx q[0];
rz(0.6181194) q[0];
rz(1.3182993) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(2.8351438) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0864799) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(-0.26992814) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0838305) q[3];
sx q[3];
rz(-0.999513) q[3];
sx q[3];
rz(0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-0.80379379) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9712517) q[0];
sx q[0];
rz(-0.60702885) q[0];
sx q[0];
rz(1.5081362) q[0];
rz(1.7159963) q[2];
sx q[2];
rz(-1.928066) q[2];
sx q[2];
rz(1.6574588) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.90384342) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(-0.018149908) q[1];
rz(0.5991163) q[3];
sx q[3];
rz(-2.0261507) q[3];
sx q[3];
rz(1.8653387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-3.0692549) q[2];
rz(0.37483254) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(-1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.9015076) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70532521) q[0];
sx q[0];
rz(-1.3316532) q[0];
sx q[0];
rz(1.5843841) q[0];
rz(-pi) q[1];
rz(1.8334332) q[2];
sx q[2];
rz(-1.8034168) q[2];
sx q[2];
rz(-0.55839415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48092914) q[1];
sx q[1];
rz(-0.96748057) q[1];
sx q[1];
rz(-1.6515345) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9661071) q[3];
sx q[3];
rz(-1.9786414) q[3];
sx q[3];
rz(0.84433489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.038625) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.4098343) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96419656) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(-3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.6606768) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4303362) q[0];
sx q[0];
rz(-0.38861409) q[0];
sx q[0];
rz(1.3769763) q[0];
rz(-pi) q[1];
rz(-3.1408177) q[2];
sx q[2];
rz(-2.0159855) q[2];
sx q[2];
rz(-1.0926525) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8112091) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(-1.1793544) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6310286) q[3];
sx q[3];
rz(-0.67959736) q[3];
sx q[3];
rz(-0.70590245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.4690171) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.8168824) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4949188) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(0.17636756) q[0];
x q[1];
rz(-0.84682805) q[2];
sx q[2];
rz(-2.4412051) q[2];
sx q[2];
rz(2.2221153) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.99609375) q[1];
sx q[1];
rz(-0.56963527) q[1];
sx q[1];
rz(0.74665229) q[1];
x q[2];
rz(2.2277101) q[3];
sx q[3];
rz(-1.8808639) q[3];
sx q[3];
rz(3.0403746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(-2.1929072) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.2896279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063110654) q[0];
sx q[0];
rz(-2.1266298) q[0];
sx q[0];
rz(-2.1102935) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81823924) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(2.408839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21266567) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(3.0287663) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4054221) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6219014) q[0];
sx q[0];
rz(-1.1150517) q[0];
sx q[0];
rz(-2.4569608) q[0];
rz(-1.7868144) q[2];
sx q[2];
rz(-1.8945433) q[2];
sx q[2];
rz(0.87460364) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64687318) q[1];
sx q[1];
rz(-0.44534007) q[1];
sx q[1];
rz(-2.7237707) q[1];
rz(-pi) q[2];
rz(-1.4275527) q[3];
sx q[3];
rz(-2.4510265) q[3];
sx q[3];
rz(0.43977246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-2.4003417) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(-2.5575976) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-2.6328971) q[0];
rz(-3.0264061) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(2.4597816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1334575) q[0];
sx q[0];
rz(-1.4032161) q[0];
sx q[0];
rz(2.125678) q[0];
x q[1];
rz(0.68306142) q[2];
sx q[2];
rz(-0.94711727) q[2];
sx q[2];
rz(1.9649486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9224285) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(-1.0132755) q[1];
x q[2];
rz(1.0107332) q[3];
sx q[3];
rz(-1.1162973) q[3];
sx q[3];
rz(-1.0487995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(-0.60733168) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(1.0495937) q[2];
sx q[2];
rz(-0.97432077) q[2];
sx q[2];
rz(-1.5885098) q[2];
rz(2.6736005) q[3];
sx q[3];
rz(-2.7157126) q[3];
sx q[3];
rz(-1.4117905) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

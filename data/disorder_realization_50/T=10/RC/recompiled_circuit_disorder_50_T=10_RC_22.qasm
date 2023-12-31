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
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89844184) q[0];
sx q[0];
rz(-1.4282707) q[0];
sx q[0];
rz(1.5205163) q[0];
rz(-pi) q[1];
rz(1.0983659) q[2];
sx q[2];
rz(-0.78394267) q[2];
sx q[2];
rz(0.69477615) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1577658) q[1];
sx q[1];
rz(-0.14769927) q[1];
sx q[1];
rz(1.2965135) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8567021) q[3];
sx q[3];
rz(-1.7171515) q[3];
sx q[3];
rz(-2.871606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4202704) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(1.7659448) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1130106) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(2.9257665) q[0];
x q[1];
rz(-1.3324845) q[2];
sx q[2];
rz(-1.6147699) q[2];
sx q[2];
rz(-0.51648742) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8010506) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(-1.3786475) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1249173) q[3];
sx q[3];
rz(-1.2614054) q[3];
sx q[3];
rz(0.644899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(-1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(2.9009853) q[0];
rz(0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.906357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2429758) q[0];
sx q[0];
rz(-1.1968062) q[0];
sx q[0];
rz(-2.5234733) q[0];
rz(-pi) q[1];
rz(-0.58358242) q[2];
sx q[2];
rz(-0.43791134) q[2];
sx q[2];
rz(2.8180518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55360824) q[1];
sx q[1];
rz(-1.8318892) q[1];
sx q[1];
rz(-1.3081461) q[1];
x q[2];
rz(1.0838305) q[3];
sx q[3];
rz(-0.999513) q[3];
sx q[3];
rz(-0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075994611) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(-0.94961387) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-0.11985699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45194295) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(2.1769051) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4255964) q[2];
sx q[2];
rz(-1.2135266) q[2];
sx q[2];
rz(1.6574588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90384342) q[1];
sx q[1];
rz(-1.3386054) q[1];
sx q[1];
rz(-3.1234427) q[1];
rz(-0.5991163) q[3];
sx q[3];
rz(-2.0261507) q[3];
sx q[3];
rz(-1.8653387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(3.0692549) q[2];
rz(2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(-1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(1.240085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70532521) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(-1.5843841) q[0];
x q[1];
rz(-2.901022) q[2];
sx q[2];
rz(-1.3153937) q[2];
sx q[2];
rz(-2.0672928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1357437) q[1];
sx q[1];
rz(-1.6372576) q[1];
sx q[1];
rz(-0.60484109) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1754856) q[3];
sx q[3];
rz(-1.1629512) q[3];
sx q[3];
rz(-2.2972578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.96419656) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(0.05274996) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(-1.6606768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67985095) q[0];
sx q[0];
rz(-1.6438419) q[0];
sx q[0];
rz(-1.1887656) q[0];
x q[1];
rz(-1.125607) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(-2.6631151) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33038352) q[1];
sx q[1];
rz(-2.5685709) q[1];
sx q[1];
rz(-1.9622383) q[1];
rz(-pi) q[2];
rz(-3.092993) q[3];
sx q[3];
rz(-2.2489293) q[3];
sx q[3];
rz(0.62852678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(-0.56224242) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19787191) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.8168824) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727407) q[0];
sx q[0];
rz(-1.4434837) q[0];
sx q[0];
rz(-1.5479814) q[0];
x q[1];
rz(2.2947646) q[2];
sx q[2];
rz(-2.4412051) q[2];
sx q[2];
rz(-0.91947739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.087254698) q[1];
sx q[1];
rz(-1.9458276) q[1];
sx q[1];
rz(0.43941984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0539829) q[3];
sx q[3];
rz(-0.71648589) q[3];
sx q[3];
rz(-2.0487752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5614732) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396486) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(2.6810714) q[0];
rz(-0.1000239) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.2896279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78532082) q[0];
sx q[0];
rz(-2.3873781) q[0];
sx q[0];
rz(-2.450599) q[0];
rz(-0.71596594) q[2];
sx q[2];
rz(-0.95574524) q[2];
sx q[2];
rz(-1.3032608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21266567) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(-0.11282632) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21059489) q[3];
sx q[3];
rz(-0.67376332) q[3];
sx q[3];
rz(-1.3099561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(-2.0588493) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(1.224068) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6219014) q[0];
sx q[0];
rz(-1.1150517) q[0];
sx q[0];
rz(0.68463188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8106868) q[2];
sx q[2];
rz(-1.3661642) q[2];
sx q[2];
rz(2.3757039) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0375367) q[1];
sx q[1];
rz(-1.9754585) q[1];
sx q[1];
rz(1.3794823) q[1];
rz(-pi) q[2];
rz(0.88528021) q[3];
sx q[3];
rz(-1.4797398) q[3];
sx q[3];
rz(-1.0202927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.042645) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(-2.4597816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54036056) q[0];
sx q[0];
rz(-1.0245748) q[0];
sx q[0];
rz(-0.19646125) q[0];
x q[1];
rz(-2.290906) q[2];
sx q[2];
rz(-0.88973532) q[2];
sx q[2];
rz(-2.1249287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9224285) q[1];
sx q[1];
rz(-2.104496) q[1];
sx q[1];
rz(-1.0132755) q[1];
x q[2];
rz(2.1308594) q[3];
sx q[3];
rz(-2.0252953) q[3];
sx q[3];
rz(2.0927932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(2.8005023) q[2];
rz(-1.0836481) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(2.9705689) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9733799) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(0.60733168) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-1.0495937) q[2];
sx q[2];
rz(-2.1672719) q[2];
sx q[2];
rz(1.5530829) q[2];
rz(-0.38469436) q[3];
sx q[3];
rz(-1.7582498) q[3];
sx q[3];
rz(2.869217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

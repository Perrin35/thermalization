OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(3.0043998) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(-0.52991968) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13745452) q[0];
sx q[0];
rz(-1.19085) q[0];
sx q[0];
rz(0.11560346) q[0];
x q[1];
rz(-2.2719703) q[2];
sx q[2];
rz(-2.6343971) q[2];
sx q[2];
rz(-1.6091572) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73218988) q[1];
sx q[1];
rz(-2.4362872) q[1];
sx q[1];
rz(2.3975055) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3210117) q[3];
sx q[3];
rz(-0.82818177) q[3];
sx q[3];
rz(-3.0299377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(1.5256933) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(-3.120378) q[0];
rz(-1.9477828) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-2.3056727) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056136925) q[0];
sx q[0];
rz(-1.5888927) q[0];
sx q[0];
rz(-1.4319112) q[0];
rz(-pi) q[1];
rz(-0.95401986) q[2];
sx q[2];
rz(-2.4980133) q[2];
sx q[2];
rz(-1.2611024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6902496) q[1];
sx q[1];
rz(-2.2022044) q[1];
sx q[1];
rz(0.33957014) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4017176) q[3];
sx q[3];
rz(-1.8679108) q[3];
sx q[3];
rz(1.5996931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2772284) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(1.345984) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(-2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-2.0879478) q[0];
rz(-1.2288278) q[1];
sx q[1];
rz(-1.5412953) q[1];
sx q[1];
rz(-2.704481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8374098) q[0];
sx q[0];
rz(-1.8619616) q[0];
sx q[0];
rz(3.0512179) q[0];
rz(-pi) q[1];
rz(2.0929298) q[2];
sx q[2];
rz(-2.6792567) q[2];
sx q[2];
rz(-2.1585652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8846109) q[1];
sx q[1];
rz(-1.7092488) q[1];
sx q[1];
rz(0.42582663) q[1];
rz(-0.06023076) q[3];
sx q[3];
rz(-2.0647991) q[3];
sx q[3];
rz(1.6935108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.019471021) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(-2.1195228) q[2];
rz(1.9034889) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6195174) q[0];
sx q[0];
rz(-1.2457122) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(0.13521067) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-0.19128004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0606196) q[0];
sx q[0];
rz(-1.4220211) q[0];
sx q[0];
rz(3.0540375) q[0];
x q[1];
rz(-0.27387597) q[2];
sx q[2];
rz(-2.2857776) q[2];
sx q[2];
rz(0.82211923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31186715) q[1];
sx q[1];
rz(-2.3448181) q[1];
sx q[1];
rz(1.4273248) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53381613) q[3];
sx q[3];
rz(-1.8752408) q[3];
sx q[3];
rz(2.387405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68025756) q[2];
sx q[2];
rz(-0.98538435) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(-2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8028832) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(2.5849735) q[0];
rz(0.11511766) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-0.97250485) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8079677) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(-1.6385965) q[0];
x q[1];
rz(0.82053484) q[2];
sx q[2];
rz(-1.8168601) q[2];
sx q[2];
rz(-2.8503502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.433686) q[1];
sx q[1];
rz(-1.3771332) q[1];
sx q[1];
rz(-0.14240264) q[1];
rz(-1.6226107) q[3];
sx q[3];
rz(-0.6165781) q[3];
sx q[3];
rz(-0.50293621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3187023) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.5931607) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(-2.1877066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71972972) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(-1.7156037) q[0];
rz(2.0772207) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-0.37429601) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4842589) q[0];
sx q[0];
rz(-1.1614292) q[0];
sx q[0];
rz(-1.9166458) q[0];
rz(-0.49403814) q[2];
sx q[2];
rz(-1.3400153) q[2];
sx q[2];
rz(1.5649232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.018505521) q[1];
sx q[1];
rz(-1.0235041) q[1];
sx q[1];
rz(1.2612543) q[1];
rz(-2.4738594) q[3];
sx q[3];
rz(-2.1453834) q[3];
sx q[3];
rz(-0.4414861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2465683) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(-0.95823112) q[2];
rz(2.9124177) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-0.62098256) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7763057) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(-2.2348485) q[0];
rz(2.0523741) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(1.3100756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7999254) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(-3.0068586) q[0];
rz(-pi) q[1];
rz(-1.9016978) q[2];
sx q[2];
rz(-2.2563997) q[2];
sx q[2];
rz(0.14724018) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3771364) q[1];
sx q[1];
rz(-2.7546282) q[1];
sx q[1];
rz(-2.6167234) q[1];
x q[2];
rz(-2.5861916) q[3];
sx q[3];
rz(-1.5474833) q[3];
sx q[3];
rz(-0.59613746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92581302) q[2];
sx q[2];
rz(-2.6999707) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-0.057597615) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9782372) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(2.9220007) q[0];
rz(2.638468) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(0.84987744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4655315) q[0];
sx q[0];
rz(-1.7277576) q[0];
sx q[0];
rz(1.7992875) q[0];
rz(-pi) q[1];
rz(-0.88768994) q[2];
sx q[2];
rz(-1.6437093) q[2];
sx q[2];
rz(-2.2788252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0975768) q[1];
sx q[1];
rz(-2.3313287) q[1];
sx q[1];
rz(-3.0939328) q[1];
rz(-pi) q[2];
rz(2.0073118) q[3];
sx q[3];
rz(-0.28537649) q[3];
sx q[3];
rz(-0.34261045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5510817) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(1.760651) q[2];
rz(0.75602174) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(-2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-1.6324156) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(-0.40503043) q[0];
rz(2.6889154) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(1.8639494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.02361) q[0];
sx q[0];
rz(-1.8007468) q[0];
sx q[0];
rz(-0.50595565) q[0];
rz(-pi) q[1];
rz(0.69127609) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(2.861475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1956049) q[1];
sx q[1];
rz(-1.1254278) q[1];
sx q[1];
rz(0.040954879) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9926662) q[3];
sx q[3];
rz(-0.78242362) q[3];
sx q[3];
rz(-2.9187834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3383125) q[2];
sx q[2];
rz(-0.40643224) q[2];
sx q[2];
rz(-2.7837616) q[2];
rz(1.7221649) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(1.0740124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(-0.11225587) q[0];
rz(2.2414801) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(2.9311438) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68073273) q[0];
sx q[0];
rz(-2.0740777) q[0];
sx q[0];
rz(1.7778648) q[0];
rz(1.1591572) q[2];
sx q[2];
rz(-1.7214081) q[2];
sx q[2];
rz(0.48356907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3810972) q[1];
sx q[1];
rz(-1.3638745) q[1];
sx q[1];
rz(2.8269672) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8080766) q[3];
sx q[3];
rz(-1.8144326) q[3];
sx q[3];
rz(2.7436649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9615053) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(1.5562742) q[2];
rz(-1.2735584) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5719941) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(-2.3251484) q[1];
sx q[1];
rz(-1.888231) q[1];
sx q[1];
rz(2.9838557) q[1];
rz(1.2522092) q[2];
sx q[2];
rz(-2.6503485) q[2];
sx q[2];
rz(2.0624401) q[2];
rz(2.3002426) q[3];
sx q[3];
rz(-1.8302866) q[3];
sx q[3];
rz(2.4091099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
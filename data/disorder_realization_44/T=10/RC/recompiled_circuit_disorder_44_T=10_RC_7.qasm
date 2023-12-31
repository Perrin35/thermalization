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
rz(-0.13719288) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(2.611673) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7512902) q[0];
sx q[0];
rz(-1.6781224) q[0];
sx q[0];
rz(-1.1885378) q[0];
rz(-pi) q[1];
rz(-0.86962236) q[2];
sx q[2];
rz(-2.6343971) q[2];
sx q[2];
rz(-1.5324355) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5296386) q[1];
sx q[1];
rz(-1.0736335) q[1];
sx q[1];
rz(1.0477209) q[1];
rz(1.3210117) q[3];
sx q[3];
rz(-0.82818177) q[3];
sx q[3];
rz(3.0299377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7933554) q[0];
sx q[0];
rz(-1.1084778) q[0];
sx q[0];
rz(0.021214699) q[0];
rz(-1.9477828) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-2.3056727) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5171889) q[0];
sx q[0];
rz(-1.7096585) q[0];
sx q[0];
rz(-0.018272321) q[0];
rz(-0.95401986) q[2];
sx q[2];
rz(-2.4980133) q[2];
sx q[2];
rz(-1.2611024) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.086120124) q[1];
sx q[1];
rz(-1.2985843) q[1];
sx q[1];
rz(2.230456) q[1];
x q[2];
rz(1.1777982) q[3];
sx q[3];
rz(-2.2713695) q[3];
sx q[3];
rz(0.29004471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8643643) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.7956087) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(-0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-2.0879478) q[0];
rz(1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(0.4371117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3041829) q[0];
sx q[0];
rz(-1.279631) q[0];
sx q[0];
rz(0.090374723) q[0];
rz(-pi) q[1];
rz(-2.0929298) q[2];
sx q[2];
rz(-0.46233593) q[2];
sx q[2];
rz(0.98302746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5324085) q[1];
sx q[1];
rz(-2.6951365) q[1];
sx q[1];
rz(0.32534728) q[1];
x q[2];
rz(-1.4594853) q[3];
sx q[3];
rz(-2.6442332) q[3];
sx q[3];
rz(-1.5745844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1221216) q[2];
sx q[2];
rz(-2.3601668) q[2];
sx q[2];
rz(2.1195228) q[2];
rz(-1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(-0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.6195174) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(-3.006382) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(-2.9503126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080973074) q[0];
sx q[0];
rz(-1.7195716) q[0];
sx q[0];
rz(-0.087555126) q[0];
rz(-pi) q[1];
rz(-2.8677167) q[2];
sx q[2];
rz(-0.855815) q[2];
sx q[2];
rz(0.82211923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8297255) q[1];
sx q[1];
rz(-2.3448181) q[1];
sx q[1];
rz(1.4273248) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9207657) q[3];
sx q[3];
rz(-2.0776437) q[3];
sx q[3];
rz(2.1496273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(1.0106687) q[2];
rz(2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8028832) q[0];
sx q[0];
rz(-2.8864679) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(0.97250485) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2302549) q[0];
sx q[0];
rz(-1.5033493) q[0];
sx q[0];
rz(3.0393242) q[0];
rz(-pi) q[1];
rz(0.33072492) q[2];
sx q[2];
rz(-0.84825584) q[2];
sx q[2];
rz(1.5028138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16469615) q[1];
sx q[1];
rz(-1.4310734) q[1];
sx q[1];
rz(1.3752027) q[1];
rz(-3.1049018) q[3];
sx q[3];
rz(-2.1864236) q[3];
sx q[3];
rz(-2.7021367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82289034) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(1.5931607) q[2];
rz(-1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.4218629) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(0.37429601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3930646) q[0];
sx q[0];
rz(-2.6120798) q[0];
sx q[0];
rz(2.4783496) q[0];
x q[1];
rz(-2.6475545) q[2];
sx q[2];
rz(-1.8015773) q[2];
sx q[2];
rz(1.5649232) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6090138) q[1];
sx q[1];
rz(-2.5206869) q[1];
sx q[1];
rz(-0.46355526) q[1];
rz(-pi) q[2];
rz(-0.80769844) q[3];
sx q[3];
rz(-2.2904615) q[3];
sx q[3];
rz(2.6157275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2465683) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(0.95823112) q[2];
rz(0.22917497) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36528698) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(2.2348485) q[0];
rz(-2.0523741) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.3100756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10446564) q[0];
sx q[0];
rz(-1.5195527) q[0];
sx q[0];
rz(2.7535704) q[0];
rz(-pi) q[1];
rz(-2.7635049) q[2];
sx q[2];
rz(-2.3921161) q[2];
sx q[2];
rz(-0.34989244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9359365) q[1];
sx q[1];
rz(-1.9034791) q[1];
sx q[1];
rz(1.7722305) q[1];
rz(-pi) q[2];
rz(0.55540107) q[3];
sx q[3];
rz(-1.5474833) q[3];
sx q[3];
rz(2.5454552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92581302) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(2.8619134) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(3.083995) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16335547) q[0];
sx q[0];
rz(-1.6219448) q[0];
sx q[0];
rz(2.9220007) q[0];
rz(0.50312463) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69707623) q[0];
sx q[0];
rz(-2.8651617) q[0];
sx q[0];
rz(2.1806549) q[0];
rz(-1.6860028) q[2];
sx q[2];
rz(-2.4552279) q[2];
sx q[2];
rz(-2.3442868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0975768) q[1];
sx q[1];
rz(-0.8102639) q[1];
sx q[1];
rz(-3.0939328) q[1];
rz(3.0181846) q[3];
sx q[3];
rz(-1.3128237) q[3];
sx q[3];
rz(2.3464399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59051096) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(1.3809416) q[2];
rz(-2.3855709) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(-0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5091771) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(-2.7365622) q[0];
rz(-2.6889154) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(-1.8639494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.02361) q[0];
sx q[0];
rz(-1.3408459) q[0];
sx q[0];
rz(-2.635637) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56565506) q[2];
sx q[2];
rz(-2.3662162) q[2];
sx q[2];
rz(-0.8650118) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1956049) q[1];
sx q[1];
rz(-2.0161649) q[1];
sx q[1];
rz(-3.1006378) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3074179) q[3];
sx q[3];
rz(-1.8636384) q[3];
sx q[3];
rz(-1.0398231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8032802) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(0.35783106) q[2];
rz(1.7221649) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(1.0740124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91838592) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(0.11225587) q[0];
rz(0.90011251) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(0.21044883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68073273) q[0];
sx q[0];
rz(-1.0675149) q[0];
sx q[0];
rz(1.3637278) q[0];
rz(-pi) q[1];
rz(-1.9333282) q[2];
sx q[2];
rz(-2.7047485) q[2];
sx q[2];
rz(2.3853962) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3810972) q[1];
sx q[1];
rz(-1.7777182) q[1];
sx q[1];
rz(-0.31462545) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8912192) q[3];
sx q[3];
rz(-1.3406521) q[3];
sx q[3];
rz(-1.2311414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18008733) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(1.5853184) q[2];
rz(1.2735584) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5719941) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(1.2522092) q[2];
sx q[2];
rz(-2.6503485) q[2];
sx q[2];
rz(2.0624401) q[2];
rz(0.84135009) q[3];
sx q[3];
rz(-1.311306) q[3];
sx q[3];
rz(-0.73248274) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-2.8357808) q[0];
sx q[0];
rz(-0.6337136) q[0];
sx q[0];
rz(-0.60854882) q[0];
rz(2.4701056) q[1];
sx q[1];
rz(-0.66317135) q[1];
sx q[1];
rz(-1.4594249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08069399) q[0];
sx q[0];
rz(-2.0804724) q[0];
sx q[0];
rz(0.26943785) q[0];
rz(-pi) q[1];
rz(-2.0107277) q[2];
sx q[2];
rz(-2.6509883) q[2];
sx q[2];
rz(3.0705796) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3069923) q[1];
sx q[1];
rz(-1.1728334) q[1];
sx q[1];
rz(-2.2521521) q[1];
x q[2];
rz(2.3816649) q[3];
sx q[3];
rz(-1.6960521) q[3];
sx q[3];
rz(1.4533991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7451611) q[2];
sx q[2];
rz(-0.67239422) q[2];
sx q[2];
rz(2.3154837) q[2];
rz(-0.36469001) q[3];
sx q[3];
rz(-1.6148022) q[3];
sx q[3];
rz(0.29432347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18225886) q[0];
sx q[0];
rz(-1.8987645) q[0];
sx q[0];
rz(0.44573319) q[0];
rz(-2.5511197) q[1];
sx q[1];
rz(-1.0346552) q[1];
sx q[1];
rz(0.42764923) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70156258) q[0];
sx q[0];
rz(-2.0272581) q[0];
sx q[0];
rz(3.0892526) q[0];
rz(-0.98154624) q[2];
sx q[2];
rz(-0.84651154) q[2];
sx q[2];
rz(1.2638448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6972097) q[1];
sx q[1];
rz(-2.1661609) q[1];
sx q[1];
rz(0.31636379) q[1];
rz(-pi) q[2];
rz(0.79645653) q[3];
sx q[3];
rz(-1.5743619) q[3];
sx q[3];
rz(-2.2623747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4599956) q[2];
sx q[2];
rz(-0.57969105) q[2];
sx q[2];
rz(-2.106529) q[2];
rz(-1.5486108) q[3];
sx q[3];
rz(-2.9620453) q[3];
sx q[3];
rz(-2.2523994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11397938) q[0];
sx q[0];
rz(-0.005682156) q[0];
sx q[0];
rz(0.11153829) q[0];
rz(1.5882209) q[1];
sx q[1];
rz(-2.7370743) q[1];
sx q[1];
rz(2.9389971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3039379) q[0];
sx q[0];
rz(-1.9662153) q[0];
sx q[0];
rz(-2.9962792) q[0];
x q[1];
rz(-2.3417046) q[2];
sx q[2];
rz(-2.7383907) q[2];
sx q[2];
rz(-3.0416807) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2592273) q[1];
sx q[1];
rz(-0.86556065) q[1];
sx q[1];
rz(2.3036876) q[1];
x q[2];
rz(2.5509994) q[3];
sx q[3];
rz(-1.9265129) q[3];
sx q[3];
rz(-1.6999262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69278875) q[2];
sx q[2];
rz(-2.0710129) q[2];
sx q[2];
rz(-2.1602574) q[2];
rz(-2.8547309) q[3];
sx q[3];
rz(-0.57932866) q[3];
sx q[3];
rz(-2.9441492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2364748) q[0];
sx q[0];
rz(-0.050967) q[0];
sx q[0];
rz(-0.68848759) q[0];
rz(-2.9972637) q[1];
sx q[1];
rz(-1.1878443) q[1];
sx q[1];
rz(0.7483288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9084602) q[0];
sx q[0];
rz(-1.5262506) q[0];
sx q[0];
rz(-2.4219803) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2690721) q[2];
sx q[2];
rz(-1.8278484) q[2];
sx q[2];
rz(-1.014705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6627745) q[1];
sx q[1];
rz(-1.8689026) q[1];
sx q[1];
rz(-1.7653905) q[1];
x q[2];
rz(2.2714628) q[3];
sx q[3];
rz(-1.3324454) q[3];
sx q[3];
rz(2.7173661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0904842) q[2];
sx q[2];
rz(-1.313831) q[2];
sx q[2];
rz(-1.8947821) q[2];
rz(0.11940739) q[3];
sx q[3];
rz(-1.9236919) q[3];
sx q[3];
rz(0.30521211) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3731821) q[0];
sx q[0];
rz(-0.38653448) q[0];
sx q[0];
rz(0.61491948) q[0];
rz(0.83456314) q[1];
sx q[1];
rz(-1.7465697) q[1];
sx q[1];
rz(-2.3051197) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30395384) q[0];
sx q[0];
rz(-1.5364354) q[0];
sx q[0];
rz(-0.0092577309) q[0];
x q[1];
rz(2.4725288) q[2];
sx q[2];
rz(-2.0537801) q[2];
sx q[2];
rz(1.3945182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53764105) q[1];
sx q[1];
rz(-1.0102235) q[1];
sx q[1];
rz(-3.0106697) q[1];
rz(2.3616124) q[3];
sx q[3];
rz(-2.2298919) q[3];
sx q[3];
rz(-2.7798228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5323083) q[2];
sx q[2];
rz(-2.534635) q[2];
sx q[2];
rz(-0.92030805) q[2];
rz(-2.3406384) q[3];
sx q[3];
rz(-2.6194173) q[3];
sx q[3];
rz(-0.33867684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9698708) q[0];
sx q[0];
rz(-2.0750177) q[0];
sx q[0];
rz(-2.875476) q[0];
rz(-1.4316106) q[1];
sx q[1];
rz(-1.1191198) q[1];
sx q[1];
rz(0.59763479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1475567) q[0];
sx q[0];
rz(-1.7606252) q[0];
sx q[0];
rz(-0.61845431) q[0];
rz(-pi) q[1];
rz(0.061068717) q[2];
sx q[2];
rz(-2.960223) q[2];
sx q[2];
rz(2.2812115) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2955851) q[1];
sx q[1];
rz(-2.4425828) q[1];
sx q[1];
rz(3/(11*pi)) q[1];
x q[2];
rz(-1.131874) q[3];
sx q[3];
rz(-0.35177975) q[3];
sx q[3];
rz(-1.1949415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0133682) q[2];
sx q[2];
rz(-2.8743663) q[2];
sx q[2];
rz(0.5989778) q[2];
rz(-0.61601764) q[3];
sx q[3];
rz(-0.44121656) q[3];
sx q[3];
rz(-1.0429355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1004341) q[0];
sx q[0];
rz(-2.2293595) q[0];
sx q[0];
rz(2.8225733) q[0];
rz(-0.11095412) q[1];
sx q[1];
rz(-1.7496505) q[1];
sx q[1];
rz(0.14161938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.148478) q[0];
sx q[0];
rz(-1.6671868) q[0];
sx q[0];
rz(-2.4575811) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1371731) q[2];
sx q[2];
rz(-1.0194155) q[2];
sx q[2];
rz(-1.0221635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1846022) q[1];
sx q[1];
rz(-0.58035589) q[1];
sx q[1];
rz(0.34837153) q[1];
rz(-pi) q[2];
rz(-1.3124489) q[3];
sx q[3];
rz(-0.53505361) q[3];
sx q[3];
rz(0.38219562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.55613279) q[2];
sx q[2];
rz(-2.1161049) q[2];
sx q[2];
rz(0.11036135) q[2];
rz(-2.6335671) q[3];
sx q[3];
rz(-3.09943) q[3];
sx q[3];
rz(1.149811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13462774) q[0];
sx q[0];
rz(-2.4899794) q[0];
sx q[0];
rz(1.2192669) q[0];
rz(0.02267516) q[1];
sx q[1];
rz(-2.2494648) q[1];
sx q[1];
rz(-2.5882744) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39529453) q[0];
sx q[0];
rz(-0.91522258) q[0];
sx q[0];
rz(3.0297999) q[0];
rz(-pi) q[1];
rz(0.99626096) q[2];
sx q[2];
rz(-1.2424412) q[2];
sx q[2];
rz(-0.36267126) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7003498) q[1];
sx q[1];
rz(-1.9158331) q[1];
sx q[1];
rz(-2.5923968) q[1];
rz(-pi) q[2];
rz(0.88388367) q[3];
sx q[3];
rz(-2.7629768) q[3];
sx q[3];
rz(0.0061638262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40552178) q[2];
sx q[2];
rz(-2.5260479) q[2];
sx q[2];
rz(-1.3496529) q[2];
rz(-3.0810629) q[3];
sx q[3];
rz(-0.90977257) q[3];
sx q[3];
rz(-0.6282261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26199207) q[0];
sx q[0];
rz(-2.1475726) q[0];
sx q[0];
rz(0.0066198786) q[0];
rz(2.5026542) q[1];
sx q[1];
rz(-0.21955755) q[1];
sx q[1];
rz(0.42983291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2868067) q[0];
sx q[0];
rz(-1.9298975) q[0];
sx q[0];
rz(-0.39639985) q[0];
rz(1.0581994) q[2];
sx q[2];
rz(-0.48658961) q[2];
sx q[2];
rz(-1.8712107) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4749516) q[1];
sx q[1];
rz(-0.78814298) q[1];
sx q[1];
rz(1.5402543) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58542975) q[3];
sx q[3];
rz(-1.6282035) q[3];
sx q[3];
rz(-1.2186183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85475737) q[2];
sx q[2];
rz(-1.1755875) q[2];
sx q[2];
rz(-0.01469928) q[2];
rz(1.1356575) q[3];
sx q[3];
rz(-2.0292) q[3];
sx q[3];
rz(-1.3057825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056674615) q[0];
sx q[0];
rz(-0.25923964) q[0];
sx q[0];
rz(-1.8490476) q[0];
rz(0.1560642) q[1];
sx q[1];
rz(-2.3239457) q[1];
sx q[1];
rz(-2.3110716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95076319) q[0];
sx q[0];
rz(-0.84759334) q[0];
sx q[0];
rz(3.0087185) q[0];
rz(-pi) q[1];
rz(2.8409675) q[2];
sx q[2];
rz(-0.81659277) q[2];
sx q[2];
rz(-2.7378094) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76354181) q[1];
sx q[1];
rz(-1.4908067) q[1];
sx q[1];
rz(-1.3536118) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51955207) q[3];
sx q[3];
rz(-0.41659912) q[3];
sx q[3];
rz(2.8263457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0326651) q[2];
sx q[2];
rz(-0.69643164) q[2];
sx q[2];
rz(-2.6984974) q[2];
rz(0.14919925) q[3];
sx q[3];
rz(-0.34330338) q[3];
sx q[3];
rz(-0.33299115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.14015848) q[0];
sx q[0];
rz(-0.59068155) q[0];
sx q[0];
rz(-0.97472192) q[0];
rz(0.99826605) q[1];
sx q[1];
rz(-1.9010192) q[1];
sx q[1];
rz(-0.98767282) q[1];
rz(2.7985991) q[2];
sx q[2];
rz(-1.7609114) q[2];
sx q[2];
rz(-1.5821725) q[2];
rz(-1.042749) q[3];
sx q[3];
rz(-0.73560148) q[3];
sx q[3];
rz(-0.92675496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

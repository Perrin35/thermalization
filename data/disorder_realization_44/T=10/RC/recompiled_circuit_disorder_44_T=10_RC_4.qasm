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
rz(2.611673) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13745452) q[0];
sx q[0];
rz(-1.19085) q[0];
sx q[0];
rz(0.11560346) q[0];
x q[1];
rz(2.797384) q[2];
sx q[2];
rz(-1.1905626) q[2];
sx q[2];
rz(-2.3772079) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73218988) q[1];
sx q[1];
rz(-0.70530546) q[1];
sx q[1];
rz(0.7440872) q[1];
rz(1.8205809) q[3];
sx q[3];
rz(-2.3134109) q[3];
sx q[3];
rz(-0.111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4522176) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.5161139) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(-1.6158993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(0.021214699) q[0];
rz(-1.9477828) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(-0.83591998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6244038) q[0];
sx q[0];
rz(-1.4319341) q[0];
sx q[0];
rz(3.1233203) q[0];
rz(-pi) q[1];
rz(-2.732227) q[2];
sx q[2];
rz(-1.0592807) q[2];
sx q[2];
rz(-0.53586938) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6902496) q[1];
sx q[1];
rz(-0.93938821) q[1];
sx q[1];
rz(0.33957014) q[1];
rz(-pi) q[2];
rz(2.4017176) q[3];
sx q[3];
rz(-1.2736819) q[3];
sx q[3];
rz(1.5996931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2772284) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(1.345984) q[2];
rz(-0.35955444) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(2.0879478) q[0];
rz(-1.2288278) q[1];
sx q[1];
rz(-1.6002974) q[1];
sx q[1];
rz(-0.4371117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489704) q[0];
sx q[0];
rz(-1.4842352) q[0];
sx q[0];
rz(1.2785046) q[0];
x q[1];
rz(1.0486629) q[2];
sx q[2];
rz(-0.46233593) q[2];
sx q[2];
rz(0.98302746) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6091842) q[1];
sx q[1];
rz(-2.6951365) q[1];
sx q[1];
rz(0.32534728) q[1];
rz(2.0655572) q[3];
sx q[3];
rz(-1.5177739) q[3];
sx q[3];
rz(3.0474636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1221216) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(-2.1195228) q[2];
rz(-1.9034889) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6195174) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(3.006382) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(2.9503126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5028338) q[0];
sx q[0];
rz(-1.4842108) q[0];
sx q[0];
rz(1.7201352) q[0];
rz(-0.27387597) q[2];
sx q[2];
rz(-2.2857776) q[2];
sx q[2];
rz(0.82211923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0335238) q[1];
sx q[1];
rz(-2.3571157) q[1];
sx q[1];
rz(2.9963521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55320923) q[3];
sx q[3];
rz(-2.534453) q[3];
sx q[3];
rz(2.7943484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4613351) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(-2.130924) q[2];
rz(2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(2.9060569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(0.97250485) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8079677) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(1.6385965) q[0];
rz(1.9237256) q[2];
sx q[2];
rz(-0.7820411) q[2];
sx q[2];
rz(2.1176586) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9768965) q[1];
sx q[1];
rz(-1.7105192) q[1];
sx q[1];
rz(1.3752027) q[1];
x q[2];
rz(-1.5189819) q[3];
sx q[3];
rz(-2.5250146) q[3];
sx q[3];
rz(2.6386564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82289034) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.548432) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(1.4259889) q[0];
rz(-2.0772207) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-2.7672966) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74852809) q[0];
sx q[0];
rz(-0.52951282) q[0];
sx q[0];
rz(2.4783496) q[0];
rz(-pi) q[1];
rz(2.6475545) q[2];
sx q[2];
rz(-1.3400153) q[2];
sx q[2];
rz(-1.5766694) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53257886) q[1];
sx q[1];
rz(-2.5206869) q[1];
sx q[1];
rz(0.46355526) q[1];
rz(-2.3338942) q[3];
sx q[3];
rz(-2.2904615) q[3];
sx q[3];
rz(0.52586517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2465683) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(2.1833615) q[2];
rz(-0.22917497) q[3];
sx q[3];
rz(-1.4567679) q[3];
sx q[3];
rz(-0.62098256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36528698) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(-2.2348485) q[0];
rz(2.0523741) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(1.3100756) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6543286) q[0];
sx q[0];
rz(-1.1833106) q[0];
sx q[0];
rz(1.6261473) q[0];
x q[1];
rz(0.37808772) q[2];
sx q[2];
rz(-2.3921161) q[2];
sx q[2];
rz(-0.34989244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8430427) q[1];
sx q[1];
rz(-1.7610465) q[1];
sx q[1];
rz(-0.33903867) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5982315) q[3];
sx q[3];
rz(-2.1260288) q[3];
sx q[3];
rz(-0.96019402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2157796) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(1.6581992) q[2];
rz(-2.8619134) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-3.083995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16335547) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(2.9220007) q[0];
rz(0.50312463) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726639) q[0];
sx q[0];
rz(-1.3451631) q[0];
sx q[0];
rz(0.16107852) q[0];
rz(-pi) q[1];
rz(-0.093896534) q[2];
sx q[2];
rz(-0.88985032) q[2];
sx q[2];
rz(0.64881334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.97493193) q[1];
sx q[1];
rz(-0.76172511) q[1];
sx q[1];
rz(1.5207661) q[1];
x q[2];
rz(-1.8306584) q[3];
sx q[3];
rz(-1.4514918) q[3];
sx q[3];
rz(-0.80727778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5510817) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(-1.760651) q[2];
rz(2.3855709) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.6324156) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(-2.7365622) q[0];
rz(-0.45267725) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(1.2776432) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67280806) q[0];
sx q[0];
rz(-2.0622258) q[0];
sx q[0];
rz(-1.3093033) q[0];
rz(-pi) q[1];
rz(-2.5759376) q[2];
sx q[2];
rz(-2.3662162) q[2];
sx q[2];
rz(0.8650118) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9459878) q[1];
sx q[1];
rz(-2.0161649) q[1];
sx q[1];
rz(-3.1006378) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7550335) q[3];
sx q[3];
rz(-2.2694526) q[3];
sx q[3];
rz(-0.78702918) q[3];
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
rz(1.4194277) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(1.0740124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91838592) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(3.0293368) q[0];
rz(-0.90011251) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(-2.9311438) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99104098) q[0];
sx q[0];
rz(-1.3897087) q[0];
sx q[0];
rz(2.6291356) q[0];
rz(-pi) q[1];
rz(-1.9333282) q[2];
sx q[2];
rz(-2.7047485) q[2];
sx q[2];
rz(2.3853962) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3810972) q[1];
sx q[1];
rz(-1.3638745) q[1];
sx q[1];
rz(0.31462545) q[1];
x q[2];
rz(-0.25037346) q[3];
sx q[3];
rz(-1.8009406) q[3];
sx q[3];
rz(1.2311414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.18008733) q[2];
sx q[2];
rz(-0.6663565) q[2];
sx q[2];
rz(1.5853184) q[2];
rz(1.8680343) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(2.5719941) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(-0.81644425) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(1.2522092) q[2];
sx q[2];
rz(-2.6503485) q[2];
sx q[2];
rz(2.0624401) q[2];
rz(0.3420842) q[3];
sx q[3];
rz(-0.87089201) q[3];
sx q[3];
rz(1.0637829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
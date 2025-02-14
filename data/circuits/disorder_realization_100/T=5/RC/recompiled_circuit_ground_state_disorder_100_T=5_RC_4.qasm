OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72803175) q[0];
sx q[0];
rz(2.1821238) q[0];
sx q[0];
rz(9.9766599) q[0];
rz(-0.38504398) q[1];
sx q[1];
rz(4.9209891) q[1];
sx q[1];
rz(9.8434386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6543579) q[0];
sx q[0];
rz(-2.4208768) q[0];
sx q[0];
rz(-2.0134175) q[0];
rz(-0.44960449) q[2];
sx q[2];
rz(-1.4105964) q[2];
sx q[2];
rz(2.2267603) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.234698) q[1];
sx q[1];
rz(-1.3677246) q[1];
sx q[1];
rz(-1.6613735) q[1];
rz(-0.017114279) q[3];
sx q[3];
rz(-1.601222) q[3];
sx q[3];
rz(-0.78998128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5358413) q[2];
sx q[2];
rz(-1.8391106) q[2];
sx q[2];
rz(-0.38899404) q[2];
rz(-0.89748663) q[3];
sx q[3];
rz(-2.5806081) q[3];
sx q[3];
rz(-1.8628023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2314583) q[0];
sx q[0];
rz(-0.95160216) q[0];
sx q[0];
rz(-2.8125473) q[0];
rz(1.344205) q[1];
sx q[1];
rz(-2.3842594) q[1];
sx q[1];
rz(-3.0175041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0253925) q[0];
sx q[0];
rz(-1.3611462) q[0];
sx q[0];
rz(-1.8127182) q[0];
rz(-pi) q[1];
rz(0.28606881) q[2];
sx q[2];
rz(-1.7033615) q[2];
sx q[2];
rz(0.90434597) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6581003) q[1];
sx q[1];
rz(-2.0112733) q[1];
sx q[1];
rz(0.3962724) q[1];
rz(-pi) q[2];
rz(-1.8120519) q[3];
sx q[3];
rz(-0.71863758) q[3];
sx q[3];
rz(-2.1332316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8590392) q[2];
sx q[2];
rz(-2.6586847) q[2];
sx q[2];
rz(2.5340951) q[2];
rz(-0.88827682) q[3];
sx q[3];
rz(-1.6162623) q[3];
sx q[3];
rz(1.4265149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4597976) q[0];
sx q[0];
rz(-1.2920222) q[0];
sx q[0];
rz(0.23455308) q[0];
rz(-1.0598496) q[1];
sx q[1];
rz(-1.9748961) q[1];
sx q[1];
rz(0.92811981) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7674196) q[0];
sx q[0];
rz(-0.95967996) q[0];
sx q[0];
rz(2.0815297) q[0];
rz(-0.71100997) q[2];
sx q[2];
rz(-0.59525049) q[2];
sx q[2];
rz(-0.75314116) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4473572) q[1];
sx q[1];
rz(-0.35590812) q[1];
sx q[1];
rz(-2.7463972) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1909799) q[3];
sx q[3];
rz(-0.39576021) q[3];
sx q[3];
rz(-2.2896965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5782535) q[2];
sx q[2];
rz(-1.4207062) q[2];
sx q[2];
rz(1.2467747) q[2];
rz(2.9747544) q[3];
sx q[3];
rz(-0.92883795) q[3];
sx q[3];
rz(1.6125352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4406776) q[0];
sx q[0];
rz(-3.0755141) q[0];
sx q[0];
rz(0.75827688) q[0];
rz(1.5658763) q[1];
sx q[1];
rz(-0.58454746) q[1];
sx q[1];
rz(2.27104) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1681002) q[0];
sx q[0];
rz(-2.3273558) q[0];
sx q[0];
rz(-0.6310985) q[0];
rz(2.5640268) q[2];
sx q[2];
rz(-2.019276) q[2];
sx q[2];
rz(2.4048793) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.16337559) q[1];
sx q[1];
rz(-1.4184457) q[1];
sx q[1];
rz(1.1635029) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6279041) q[3];
sx q[3];
rz(-1.6293173) q[3];
sx q[3];
rz(2.1128775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42252758) q[2];
sx q[2];
rz(-2.1336522) q[2];
sx q[2];
rz(-0.75971216) q[2];
rz(-0.64940137) q[3];
sx q[3];
rz(-1.5348744) q[3];
sx q[3];
rz(1.5788797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11288697) q[0];
sx q[0];
rz(-0.6539456) q[0];
sx q[0];
rz(1.9586067) q[0];
rz(-0.68823367) q[1];
sx q[1];
rz(-1.8249244) q[1];
sx q[1];
rz(-2.7722955) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1210021) q[0];
sx q[0];
rz(-0.90743104) q[0];
sx q[0];
rz(-1.6055907) q[0];
rz(0.30226548) q[2];
sx q[2];
rz(-2.4068953) q[2];
sx q[2];
rz(-2.1718028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5253882) q[1];
sx q[1];
rz(-1.3862351) q[1];
sx q[1];
rz(0.86854684) q[1];
rz(-0.33325382) q[3];
sx q[3];
rz(-1.2566393) q[3];
sx q[3];
rz(-0.21018782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94283048) q[2];
sx q[2];
rz(-1.3254712) q[2];
sx q[2];
rz(-2.14373) q[2];
rz(-0.61740795) q[3];
sx q[3];
rz(-2.182775) q[3];
sx q[3];
rz(-0.82120419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038641039) q[0];
sx q[0];
rz(-1.8674253) q[0];
sx q[0];
rz(1.8983023) q[0];
rz(2.359911) q[1];
sx q[1];
rz(-1.2739173) q[1];
sx q[1];
rz(2.8807358) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2170625) q[0];
sx q[0];
rz(-0.53981656) q[0];
sx q[0];
rz(-2.7227719) q[0];
rz(-pi) q[1];
rz(0.37520295) q[2];
sx q[2];
rz(-1.1625566) q[2];
sx q[2];
rz(-0.21438504) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8159224) q[1];
sx q[1];
rz(-1.4192033) q[1];
sx q[1];
rz(1.1071015) q[1];
rz(-pi) q[2];
rz(2.2024037) q[3];
sx q[3];
rz(-1.7429461) q[3];
sx q[3];
rz(0.84014713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6692052) q[2];
sx q[2];
rz(-0.86296764) q[2];
sx q[2];
rz(-0.48259398) q[2];
rz(0.11416301) q[3];
sx q[3];
rz(-1.0897021) q[3];
sx q[3];
rz(2.193006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77940762) q[0];
sx q[0];
rz(-3.0857093) q[0];
sx q[0];
rz(-2.8756397) q[0];
rz(-2.7159122) q[1];
sx q[1];
rz(-1.2454147) q[1];
sx q[1];
rz(-0.46806213) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6702995) q[0];
sx q[0];
rz(-0.54139304) q[0];
sx q[0];
rz(1.584021) q[0];
rz(-1.9114248) q[2];
sx q[2];
rz(-1.4832895) q[2];
sx q[2];
rz(-0.82078314) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99331964) q[1];
sx q[1];
rz(-2.5778505) q[1];
sx q[1];
rz(-2.0658653) q[1];
rz(1.0506094) q[3];
sx q[3];
rz(-1.597763) q[3];
sx q[3];
rz(0.20668465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88428664) q[2];
sx q[2];
rz(-1.2904737) q[2];
sx q[2];
rz(2.5878944) q[2];
rz(2.2181559) q[3];
sx q[3];
rz(-0.90673509) q[3];
sx q[3];
rz(-3.0807909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6786574) q[0];
sx q[0];
rz(-0.23660062) q[0];
sx q[0];
rz(-0.67805725) q[0];
rz(2.1642115) q[1];
sx q[1];
rz(-1.9934318) q[1];
sx q[1];
rz(1.4378907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554094) q[0];
sx q[0];
rz(-2.0269505) q[0];
sx q[0];
rz(0.21502226) q[0];
rz(-2.5776093) q[2];
sx q[2];
rz(-1.5684557) q[2];
sx q[2];
rz(2.7067647) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6594636) q[1];
sx q[1];
rz(-1.8244484) q[1];
sx q[1];
rz(-1.8096022) q[1];
x q[2];
rz(2.1071042) q[3];
sx q[3];
rz(-0.89867979) q[3];
sx q[3];
rz(-1.514707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4108654) q[2];
sx q[2];
rz(-2.5812456) q[2];
sx q[2];
rz(0.66696683) q[2];
rz(-3.0242331) q[3];
sx q[3];
rz(-1.8662235) q[3];
sx q[3];
rz(1.7608775) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19408195) q[0];
sx q[0];
rz(-1.5503333) q[0];
sx q[0];
rz(0.34341735) q[0];
rz(-1.6471479) q[1];
sx q[1];
rz(-2.1518555) q[1];
sx q[1];
rz(-0.48430482) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0750879) q[0];
sx q[0];
rz(-2.4874788) q[0];
sx q[0];
rz(-0.44371407) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3745851) q[2];
sx q[2];
rz(-2.5756774) q[2];
sx q[2];
rz(-0.2919251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6086413) q[1];
sx q[1];
rz(-1.7923454) q[1];
sx q[1];
rz(-0.59355841) q[1];
x q[2];
rz(-0.24933322) q[3];
sx q[3];
rz(-1.7847381) q[3];
sx q[3];
rz(0.46505022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9453498) q[2];
sx q[2];
rz(-2.0413155) q[2];
sx q[2];
rz(-0.84570447) q[2];
rz(1.9218933) q[3];
sx q[3];
rz(-1.889735) q[3];
sx q[3];
rz(-0.20488258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6824816) q[0];
sx q[0];
rz(-0.24523188) q[0];
sx q[0];
rz(-2.5526175) q[0];
rz(-2.466195) q[1];
sx q[1];
rz(-0.94869906) q[1];
sx q[1];
rz(1.4640456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9927514) q[0];
sx q[0];
rz(-2.0939079) q[0];
sx q[0];
rz(0.4776202) q[0];
x q[1];
rz(-2.3039742) q[2];
sx q[2];
rz(-1.5116201) q[2];
sx q[2];
rz(-2.0967332) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88661609) q[1];
sx q[1];
rz(-0.58838298) q[1];
sx q[1];
rz(-0.083440668) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0328224) q[3];
sx q[3];
rz(-2.2115876) q[3];
sx q[3];
rz(-2.6907211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7318763) q[2];
sx q[2];
rz(-1.7475374) q[2];
sx q[2];
rz(2.0173006) q[2];
rz(-0.8479979) q[3];
sx q[3];
rz(-2.336899) q[3];
sx q[3];
rz(3.1112352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59415862) q[0];
sx q[0];
rz(-2.0068598) q[0];
sx q[0];
rz(1.6225847) q[0];
rz(2.2183954) q[1];
sx q[1];
rz(-2.1122439) q[1];
sx q[1];
rz(-1.5069638) q[1];
rz(-1.3849003) q[2];
sx q[2];
rz(-1.8757314) q[2];
sx q[2];
rz(1.3396946) q[2];
rz(0.93529978) q[3];
sx q[3];
rz(-2.7705396) q[3];
sx q[3];
rz(-2.0075575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

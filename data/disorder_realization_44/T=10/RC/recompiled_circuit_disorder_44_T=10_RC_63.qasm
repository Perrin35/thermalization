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
rz(6.7232806) q[0];
sx q[0];
rz(6.4203782) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(-0.52991968) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7007028) q[0];
sx q[0];
rz(-2.7452677) q[0];
sx q[0];
rz(1.2896145) q[0];
rz(-pi) q[1];
rz(2.797384) q[2];
sx q[2];
rz(-1.1905626) q[2];
sx q[2];
rz(0.76438475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.73218988) q[1];
sx q[1];
rz(-2.4362872) q[1];
sx q[1];
rz(0.7440872) q[1];
rz(0.75833851) q[3];
sx q[3];
rz(-1.3876649) q[3];
sx q[3];
rz(-1.5116215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4522176) q[2];
sx q[2];
rz(-1.8415425) q[2];
sx q[2];
rz(-2.8049862) q[2];
rz(1.6254788) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.34823725) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(-3.120378) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(-0.83591998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056136925) q[0];
sx q[0];
rz(-1.5526999) q[0];
sx q[0];
rz(-1.4319112) q[0];
x q[1];
rz(-2.1875728) q[2];
sx q[2];
rz(-0.6435794) q[2];
sx q[2];
rz(-1.2611024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6902496) q[1];
sx q[1];
rz(-0.93938821) q[1];
sx q[1];
rz(-0.33957014) q[1];
x q[2];
rz(-2.7153035) q[3];
sx q[3];
rz(-0.78668919) q[3];
sx q[3];
rz(0.28144893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8643643) q[2];
sx q[2];
rz(-1.1479062) q[2];
sx q[2];
rz(-1.7956087) q[2];
rz(-2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(2.0879478) q[0];
rz(-1.2288278) q[1];
sx q[1];
rz(-1.5412953) q[1];
sx q[1];
rz(0.4371117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1431883) q[0];
sx q[0];
rz(-0.30448738) q[0];
sx q[0];
rz(1.8633153) q[0];
rz(-0.24361165) q[2];
sx q[2];
rz(-1.1738452) q[2];
sx q[2];
rz(1.5872019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25698173) q[1];
sx q[1];
rz(-1.4323438) q[1];
sx q[1];
rz(-2.715766) q[1];
x q[2];
rz(1.4594853) q[3];
sx q[3];
rz(-2.6442332) q[3];
sx q[3];
rz(-1.5670083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.019471021) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(-1.9034889) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(-2.7220272) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(-0.98130256) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-2.0573261) q[1];
sx q[1];
rz(2.9503126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6879038) q[0];
sx q[0];
rz(-2.9691302) q[0];
sx q[0];
rz(1.0426636) q[0];
rz(-pi) q[1];
rz(1.26881) q[2];
sx q[2];
rz(-2.3846845) q[2];
sx q[2];
rz(-0.41727558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10806882) q[1];
sx q[1];
rz(-2.3571157) q[1];
sx q[1];
rz(-2.9963521) q[1];
rz(-pi) q[2];
rz(-2.5883834) q[3];
sx q[3];
rz(-2.534453) q[3];
sx q[3];
rz(0.34724423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(-1.0106687) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(-0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-2.8028832) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(-0.97250485) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8079677) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(-1.6385965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33072492) q[2];
sx q[2];
rz(-2.2933368) q[2];
sx q[2];
rz(-1.6387788) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9768965) q[1];
sx q[1];
rz(-1.4310734) q[1];
sx q[1];
rz(-1.7663899) q[1];
rz(1.6226107) q[3];
sx q[3];
rz(-2.5250146) q[3];
sx q[3];
rz(2.6386564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3187023) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.548432) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-0.32309353) q[3];
sx q[3];
rz(-0.95388609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4218629) q[0];
sx q[0];
rz(-1.2635764) q[0];
sx q[0];
rz(-1.4259889) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-1.0168889) q[1];
sx q[1];
rz(-2.7672966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22898856) q[0];
sx q[0];
rz(-1.2545663) q[0];
sx q[0];
rz(0.4321179) q[0];
x q[1];
rz(1.3099953) q[2];
sx q[2];
rz(-1.0909832) q[2];
sx q[2];
rz(3.0248883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7172076) q[1];
sx q[1];
rz(-1.3076412) q[1];
sx q[1];
rz(-0.56916635) q[1];
x q[2];
rz(-0.80769844) q[3];
sx q[3];
rz(-2.2904615) q[3];
sx q[3];
rz(-0.52586517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2465683) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(2.1833615) q[2];
rz(-2.9124177) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(2.5206101) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36528698) q[0];
sx q[0];
rz(-1.9488652) q[0];
sx q[0];
rz(-0.90674415) q[0];
rz(1.0892185) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(-1.8315171) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3416672) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(-0.13473405) q[0];
x q[1];
rz(-0.37808772) q[2];
sx q[2];
rz(-0.74947658) q[2];
sx q[2];
rz(-0.34989244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8430427) q[1];
sx q[1];
rz(-1.7610465) q[1];
sx q[1];
rz(2.802554) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.044192627) q[3];
sx q[3];
rz(-2.5857539) q[3];
sx q[3];
rz(1.0122055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2157796) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(-1.6581992) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(0.057597615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.2917152) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6760611) q[0];
sx q[0];
rz(-1.413835) q[0];
sx q[0];
rz(-1.7992875) q[0];
x q[1];
rz(1.6860028) q[2];
sx q[2];
rz(-0.68636471) q[2];
sx q[2];
rz(-2.3442868) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1666607) q[1];
sx q[1];
rz(-0.76172511) q[1];
sx q[1];
rz(-1.5207661) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8306584) q[3];
sx q[3];
rz(-1.4514918) q[3];
sx q[3];
rz(-0.80727778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5510817) q[2];
sx q[2];
rz(-1.8293646) q[2];
sx q[2];
rz(1.760651) q[2];
rz(-0.75602174) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(-0.40503043) q[0];
rz(-0.45267725) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(-1.8639494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4687846) q[0];
sx q[0];
rz(-1.0793669) q[0];
sx q[0];
rz(-1.3093033) q[0];
x q[1];
rz(-2.5759376) q[2];
sx q[2];
rz(-0.7753765) q[2];
sx q[2];
rz(-0.8650118) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35754044) q[1];
sx q[1];
rz(-1.5338384) q[1];
sx q[1];
rz(-1.1251015) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3074179) q[3];
sx q[3];
rz(-1.2779543) q[3];
sx q[3];
rz(-1.0398231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3383125) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(-2.7837616) q[2];
rz(-1.4194277) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(2.0675802) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(3.0293368) q[0];
rz(-2.2414801) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(0.21044883) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8716547) q[0];
sx q[0];
rz(-2.6007814) q[0];
sx q[0];
rz(-2.7842115) q[0];
rz(-2.9774882) q[2];
sx q[2];
rz(-1.1640942) q[2];
sx q[2];
rz(1.1526398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3887742) q[1];
sx q[1];
rz(-0.37467271) q[1];
sx q[1];
rz(2.545536) q[1];
rz(-pi) q[2];
rz(2.3841303) q[3];
sx q[3];
rz(-2.8031581) q[3];
sx q[3];
rz(0.38871845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18008733) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(1.5562742) q[2];
rz(-1.2735584) q[3];
sx q[3];
rz(-2.5189416) q[3];
sx q[3];
rz(0.20726985) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56959854) q[0];
sx q[0];
rz(-2.270569) q[0];
sx q[0];
rz(1.7763174) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(1.1007166) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(2.7995085) q[3];
sx q[3];
rz(-2.2707006) q[3];
sx q[3];
rz(-2.0778098) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
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
rz(0.91520619) q[0];
sx q[0];
rz(2.6829166) q[0];
sx q[0];
rz(12.248326) q[0];
rz(1.3660499) q[1];
sx q[1];
rz(-2.4429758) q[1];
sx q[1];
rz(0.06279343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88616497) q[0];
sx q[0];
rz(-1.0647827) q[0];
sx q[0];
rz(0.54765986) q[0];
x q[1];
rz(-2.7032825) q[2];
sx q[2];
rz(-2.2026014) q[2];
sx q[2];
rz(0.28011986) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5972728) q[1];
sx q[1];
rz(-1.8493127) q[1];
sx q[1];
rz(2.0803948) q[1];
rz(-pi) q[2];
rz(2.1476521) q[3];
sx q[3];
rz(-1.2341043) q[3];
sx q[3];
rz(-1.5693046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.723145) q[2];
sx q[2];
rz(-1.3157088) q[2];
sx q[2];
rz(-2.1212228) q[2];
rz(0.72888914) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(1.6093904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51434022) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(0.038272055) q[0];
rz(-3.017784) q[1];
sx q[1];
rz(-0.47272155) q[1];
sx q[1];
rz(-1.5708057) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8077652) q[0];
sx q[0];
rz(-0.021919576) q[0];
sx q[0];
rz(2.0448565) q[0];
rz(-1.6699617) q[2];
sx q[2];
rz(-0.80201521) q[2];
sx q[2];
rz(2.2834509) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60197406) q[1];
sx q[1];
rz(-1.5718196) q[1];
sx q[1];
rz(-0.00042331528) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6917211) q[3];
sx q[3];
rz(-2.0837415) q[3];
sx q[3];
rz(0.67739048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52246419) q[2];
sx q[2];
rz(-1.3294486) q[2];
sx q[2];
rz(-0.5788571) q[2];
rz(-1.045643) q[3];
sx q[3];
rz(-0.78543109) q[3];
sx q[3];
rz(-0.75700179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-2.4724562) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(0.45397595) q[0];
rz(0.16481915) q[1];
sx q[1];
rz(-0.77607981) q[1];
sx q[1];
rz(-1.4705315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1810581) q[0];
sx q[0];
rz(-1.9583324) q[0];
sx q[0];
rz(-2.9563532) q[0];
rz(2.0963832) q[2];
sx q[2];
rz(-1.8586577) q[2];
sx q[2];
rz(0.20341104) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3537703) q[1];
sx q[1];
rz(-2.2332237) q[1];
sx q[1];
rz(-1.9595998) q[1];
rz(-2.4999714) q[3];
sx q[3];
rz(-0.76867662) q[3];
sx q[3];
rz(-2.8856483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37956023) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(1.3564159) q[2];
rz(-0.49946579) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(-1.0511901) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22719638) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(2.0830578) q[0];
rz(-1.1610441) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(3.0243691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0441705) q[0];
sx q[0];
rz(-1.6099842) q[0];
sx q[0];
rz(0.54255658) q[0];
rz(-pi) q[1];
rz(2.2160596) q[2];
sx q[2];
rz(-1.4129854) q[2];
sx q[2];
rz(0.35343808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30018729) q[1];
sx q[1];
rz(-0.8148199) q[1];
sx q[1];
rz(-2.7542265) q[1];
rz(-pi) q[2];
rz(-1.1612558) q[3];
sx q[3];
rz(-0.97417384) q[3];
sx q[3];
rz(2.2724702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78493541) q[2];
sx q[2];
rz(-1.6931809) q[2];
sx q[2];
rz(1.7735205) q[2];
rz(-1.8562227) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(-2.4135597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0990937) q[0];
sx q[0];
rz(-0.92785257) q[0];
sx q[0];
rz(1.4121144) q[0];
rz(1.0026576) q[1];
sx q[1];
rz(-2.899677) q[1];
sx q[1];
rz(1.5589421) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2501734) q[0];
sx q[0];
rz(-1.0153664) q[0];
sx q[0];
rz(-2.069418) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71835235) q[2];
sx q[2];
rz(-1.6334772) q[2];
sx q[2];
rz(-2.2557392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0274899) q[1];
sx q[1];
rz(-1.5943502) q[1];
sx q[1];
rz(-1.0204169) q[1];
rz(-pi) q[2];
rz(-1.2491426) q[3];
sx q[3];
rz(-1.6603004) q[3];
sx q[3];
rz(-2.9200875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58854181) q[2];
sx q[2];
rz(-1.43247) q[2];
sx q[2];
rz(0.39349619) q[2];
rz(-2.1816175) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(0.967832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680962) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(0.18336503) q[0];
rz(2.6496437) q[1];
sx q[1];
rz(-0.79194561) q[1];
sx q[1];
rz(-1.2194182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6693319) q[0];
sx q[0];
rz(-0.12024256) q[0];
sx q[0];
rz(-0.73445102) q[0];
x q[1];
rz(-1.9398795) q[2];
sx q[2];
rz(-0.08412349) q[2];
sx q[2];
rz(0.28757986) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9522493) q[1];
sx q[1];
rz(-0.5533411) q[1];
sx q[1];
rz(2.4207741) q[1];
rz(1.4063356) q[3];
sx q[3];
rz(-2.6643848) q[3];
sx q[3];
rz(-2.8754614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(-0.08237002) q[2];
rz(-0.92665893) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(1.5026708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53967404) q[0];
sx q[0];
rz(-0.32783666) q[0];
sx q[0];
rz(0.71640054) q[0];
rz(2.7377103) q[1];
sx q[1];
rz(-1.5682181) q[1];
sx q[1];
rz(1.4366038) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8536896) q[0];
sx q[0];
rz(-2.7349964) q[0];
sx q[0];
rz(-1.2000183) q[0];
x q[1];
rz(1.4384957) q[2];
sx q[2];
rz(-2.0952333) q[2];
sx q[2];
rz(-0.99614267) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75175367) q[1];
sx q[1];
rz(-1.1734278) q[1];
sx q[1];
rz(-1.6786871) q[1];
x q[2];
rz(-0.8088423) q[3];
sx q[3];
rz(-1.0029965) q[3];
sx q[3];
rz(-2.0140935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3639823) q[2];
sx q[2];
rz(-0.97475514) q[2];
sx q[2];
rz(-3.0762365) q[2];
rz(-2.2743547) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.48677483) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(-0.73100334) q[0];
rz(2.3664318) q[1];
sx q[1];
rz(-2.8072) q[1];
sx q[1];
rz(2.7269272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1461244) q[0];
sx q[0];
rz(-1.0421317) q[0];
sx q[0];
rz(2.29236) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.783466) q[2];
sx q[2];
rz(-2.0913107) q[2];
sx q[2];
rz(-1.6047603) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2036671) q[1];
sx q[1];
rz(-1.3199184) q[1];
sx q[1];
rz(2.832042) q[1];
rz(-pi) q[2];
rz(-0.13402588) q[3];
sx q[3];
rz(-1.7334337) q[3];
sx q[3];
rz(-2.7664879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23184648) q[2];
sx q[2];
rz(-2.6555588) q[2];
sx q[2];
rz(-3.0492142) q[2];
rz(0.77629027) q[3];
sx q[3];
rz(-1.679136) q[3];
sx q[3];
rz(2.9379454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6542776) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(3.1238632) q[0];
rz(-2.0954633) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(-1.0606934) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6700376) q[0];
sx q[0];
rz(-3.0714543) q[0];
sx q[0];
rz(-0.20151968) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1589375) q[2];
sx q[2];
rz(-1.3545631) q[2];
sx q[2];
rz(2.4311291) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2613147) q[1];
sx q[1];
rz(-0.78431097) q[1];
sx q[1];
rz(-1.3256636) q[1];
rz(-pi) q[2];
rz(-0.10194998) q[3];
sx q[3];
rz(-0.52604874) q[3];
sx q[3];
rz(2.7395484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1549418) q[2];
sx q[2];
rz(-2.3255746) q[2];
sx q[2];
rz(0.6443392) q[2];
rz(0.84195697) q[3];
sx q[3];
rz(-1.5717477) q[3];
sx q[3];
rz(-0.90698609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3802721) q[0];
sx q[0];
rz(-0.43571061) q[0];
sx q[0];
rz(-2.6897588) q[0];
rz(1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(1.5712646) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1399978) q[0];
sx q[0];
rz(-2.9712935) q[0];
sx q[0];
rz(3.1055121) q[0];
rz(3.0250835) q[2];
sx q[2];
rz(-2.0429869) q[2];
sx q[2];
rz(-1.0509059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.328969) q[1];
sx q[1];
rz(-1.4419893) q[1];
sx q[1];
rz(-1.179943) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4191946) q[3];
sx q[3];
rz(-1.5766451) q[3];
sx q[3];
rz(-1.4129782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5181804) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(1.5709467) q[2];
rz(-0.10562854) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-2.6251729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.327772) q[0];
sx q[0];
rz(-1.0186503) q[0];
sx q[0];
rz(-3.0042197) q[0];
rz(0.2188006) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(1.995718) q[2];
sx q[2];
rz(-2.6617202) q[2];
sx q[2];
rz(-0.94750994) q[2];
rz(-1.9365334) q[3];
sx q[3];
rz(-2.1954721) q[3];
sx q[3];
rz(1.1385067) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-2.0353844) q[0];
sx q[0];
rz(2.455403) q[0];
sx q[0];
rz(10.578293) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(0.44841132) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8657239) q[0];
sx q[0];
rz(-0.65058904) q[0];
sx q[0];
rz(3.0537729) q[0];
x q[1];
rz(1.469606) q[2];
sx q[2];
rz(-1.3162344) q[2];
sx q[2];
rz(2.2806666) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37473512) q[1];
sx q[1];
rz(-2.285706) q[1];
sx q[1];
rz(2.292407) q[1];
rz(3.0541522) q[3];
sx q[3];
rz(-2.6387847) q[3];
sx q[3];
rz(-1.8016165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6260234) q[2];
sx q[2];
rz(-1.7531351) q[2];
sx q[2];
rz(-0.62421978) q[2];
rz(-2.7624687) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(2.8930801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666331) q[0];
sx q[0];
rz(-2.3269854) q[0];
sx q[0];
rz(2.1061184) q[0];
rz(1.8677208) q[1];
sx q[1];
rz(-2.404411) q[1];
sx q[1];
rz(2.6587291) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57377061) q[0];
sx q[0];
rz(-1.616524) q[0];
sx q[0];
rz(3.0450495) q[0];
rz(2.5245978) q[2];
sx q[2];
rz(-2.9386387) q[2];
sx q[2];
rz(-0.49398068) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9804364) q[1];
sx q[1];
rz(-2.709678) q[1];
sx q[1];
rz(-1.5932139) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8206536) q[3];
sx q[3];
rz(-1.648099) q[3];
sx q[3];
rz(-2.2805205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.024293385) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(-1.8710322) q[2];
rz(0.67000669) q[3];
sx q[3];
rz(-2.5195401) q[3];
sx q[3];
rz(-2.2445934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(-2.4025412) q[0];
rz(2.1801379) q[1];
sx q[1];
rz(-0.32197222) q[1];
sx q[1];
rz(2.3846073) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77853816) q[0];
sx q[0];
rz(-0.86150384) q[0];
sx q[0];
rz(-2.7947133) q[0];
x q[1];
rz(-0.44481014) q[2];
sx q[2];
rz(-0.74193566) q[2];
sx q[2];
rz(-0.43659376) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5730505) q[1];
sx q[1];
rz(-1.7049763) q[1];
sx q[1];
rz(-2.5106984) q[1];
rz(1.3540165) q[3];
sx q[3];
rz(-2.8269973) q[3];
sx q[3];
rz(-1.6198688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6387393) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(2.4364831) q[2];
rz(-2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(-2.7307935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0099156378) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(1.9158844) q[0];
rz(1.907584) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(3.0753678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4915337) q[0];
sx q[0];
rz(-2.8797132) q[0];
sx q[0];
rz(1.9240379) q[0];
x q[1];
rz(0.29860626) q[2];
sx q[2];
rz(-0.82393194) q[2];
sx q[2];
rz(0.9048942) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.855763) q[1];
sx q[1];
rz(-1.0914088) q[1];
sx q[1];
rz(-0.3605607) q[1];
rz(-pi) q[2];
rz(2.8230571) q[3];
sx q[3];
rz(-2.7580166) q[3];
sx q[3];
rz(2.2936402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0951198) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(-0.44060102) q[2];
rz(2.1353841) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(-2.3073176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42136583) q[0];
sx q[0];
rz(-2.328673) q[0];
sx q[0];
rz(1.6356069) q[0];
rz(-1.5274564) q[1];
sx q[1];
rz(-1.2200049) q[1];
sx q[1];
rz(-1.6857326) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57350006) q[0];
sx q[0];
rz(-1.4254036) q[0];
sx q[0];
rz(0.83515867) q[0];
rz(0.22730374) q[2];
sx q[2];
rz(-1.6844201) q[2];
sx q[2];
rz(-2.794968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8090022) q[1];
sx q[1];
rz(-2.2223516) q[1];
sx q[1];
rz(1.5963735) q[1];
rz(-pi) q[2];
rz(1.4167352) q[3];
sx q[3];
rz(-0.96821456) q[3];
sx q[3];
rz(-1.3101206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27318925) q[2];
sx q[2];
rz(-1.5634147) q[2];
sx q[2];
rz(-0.91147649) q[2];
rz(-0.8738001) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(3.1349283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8580496) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(-0.13033303) q[0];
rz(-2.6538972) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(2.3977051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7672507) q[0];
sx q[0];
rz(-0.09627405) q[0];
sx q[0];
rz(-1.0262579) q[0];
rz(-1.3792453) q[2];
sx q[2];
rz(-0.63427502) q[2];
sx q[2];
rz(0.77855643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2762017) q[1];
sx q[1];
rz(-1.8542093) q[1];
sx q[1];
rz(1.0083501) q[1];
rz(-0.33244953) q[3];
sx q[3];
rz(-0.5088734) q[3];
sx q[3];
rz(-0.32989855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0222212) q[2];
sx q[2];
rz(-2.2890942) q[2];
sx q[2];
rz(-0.96088299) q[2];
rz(-0.39572257) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(0.1161639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845067) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(-2.0945666) q[0];
rz(2.537435) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(-0.39271694) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34907863) q[0];
sx q[0];
rz(-0.70988467) q[0];
sx q[0];
rz(-1.2997846) q[0];
rz(3.0839594) q[2];
sx q[2];
rz(-2.1890292) q[2];
sx q[2];
rz(-1.0317486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3369357) q[1];
sx q[1];
rz(-1.5439022) q[1];
sx q[1];
rz(2.8662445) q[1];
x q[2];
rz(-1.3052528) q[3];
sx q[3];
rz(-1.2528386) q[3];
sx q[3];
rz(-1.4972655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.087223209) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(-1.8325904) q[2];
rz(-1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8488309) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(-0.30717474) q[0];
rz(1.3245026) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(-2.2408392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4221638) q[0];
sx q[0];
rz(-1.8329282) q[0];
sx q[0];
rz(1.5533226) q[0];
rz(-2.6866954) q[2];
sx q[2];
rz(-2.1730246) q[2];
sx q[2];
rz(-2.086161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5430205) q[1];
sx q[1];
rz(-1.2391587) q[1];
sx q[1];
rz(-2.806862) q[1];
rz(-pi) q[2];
rz(-2.9489904) q[3];
sx q[3];
rz(-2.5012272) q[3];
sx q[3];
rz(-0.34512025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8073392) q[2];
sx q[2];
rz(-3.0901577) q[2];
sx q[2];
rz(0.61484289) q[2];
rz(-1.0963415) q[3];
sx q[3];
rz(-2.5494826) q[3];
sx q[3];
rz(-2.443327) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39356247) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(0.50931859) q[0];
rz(-2.9604984) q[1];
sx q[1];
rz(-1.4053922) q[1];
sx q[1];
rz(-2.1720355) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90617243) q[0];
sx q[0];
rz(-2.3114822) q[0];
sx q[0];
rz(1.3614015) q[0];
rz(-pi) q[1];
rz(-0.050886919) q[2];
sx q[2];
rz(-1.0641838) q[2];
sx q[2];
rz(1.2801054) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6020376) q[1];
sx q[1];
rz(-1.6786715) q[1];
sx q[1];
rz(-2.8007568) q[1];
rz(2.2990555) q[3];
sx q[3];
rz(-1.304783) q[3];
sx q[3];
rz(-2.0913948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8934882) q[2];
sx q[2];
rz(-2.9073145) q[2];
sx q[2];
rz(-1.0815557) q[2];
rz(0.23165101) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.76081) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(-2.494452) q[0];
rz(0.45516792) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(2.3796577) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9058911) q[0];
sx q[0];
rz(-0.76607031) q[0];
sx q[0];
rz(-0.42278843) q[0];
rz(2.3222515) q[2];
sx q[2];
rz(-2.0692673) q[2];
sx q[2];
rz(2.8070531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.476452) q[1];
sx q[1];
rz(-2.2287205) q[1];
sx q[1];
rz(1.5275147) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1338828) q[3];
sx q[3];
rz(-1.1033909) q[3];
sx q[3];
rz(-1.7860408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.053293856) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(-2.2376412) q[2];
rz(1.5229185) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.224613) q[0];
sx q[0];
rz(-1.1876748) q[0];
sx q[0];
rz(-1.351958) q[0];
rz(-3.0147973) q[1];
sx q[1];
rz(-2.3190111) q[1];
sx q[1];
rz(0.93228985) q[1];
rz(2.8895072) q[2];
sx q[2];
rz(-1.132195) q[2];
sx q[2];
rz(1.8762527) q[2];
rz(-1.5679172) q[3];
sx q[3];
rz(-2.0007912) q[3];
sx q[3];
rz(3.1040647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

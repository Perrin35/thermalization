OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(-2.2025684) q[0];
sx q[0];
rz(0.0052069081) q[0];
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(-1.9519238) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0931041) q[0];
sx q[0];
rz(-2.5079873) q[0];
sx q[0];
rz(-2.5392169) q[0];
x q[1];
rz(-0.66733811) q[2];
sx q[2];
rz(-2.9138406) q[2];
sx q[2];
rz(1.8337133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3559239) q[1];
sx q[1];
rz(-2.7969116) q[1];
sx q[1];
rz(1.1204526) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7351904) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(0.5973967) q[2];
rz(-1.3655837) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.8610154) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(3.0283668) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9202168) q[0];
sx q[0];
rz(-0.42718233) q[0];
sx q[0];
rz(1.4772052) q[0];
rz(-pi) q[1];
rz(1.5905154) q[2];
sx q[2];
rz(-1.4975417) q[2];
sx q[2];
rz(0.8578701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3791703) q[1];
sx q[1];
rz(-1.1602853) q[1];
sx q[1];
rz(1.1344086) q[1];
rz(0.23785915) q[3];
sx q[3];
rz(-0.97994643) q[3];
sx q[3];
rz(-1.8464586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-2.9193027) q[2];
rz(2.9120581) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(-1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-3.0139794) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1193082) q[0];
sx q[0];
rz(-2.2420609) q[0];
sx q[0];
rz(0.22247252) q[0];
rz(-pi) q[1];
rz(-0.80673809) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(-2.9850609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0763921) q[1];
sx q[1];
rz(-0.25307357) q[1];
sx q[1];
rz(-1.2364926) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8887799) q[3];
sx q[3];
rz(-0.95889927) q[3];
sx q[3];
rz(1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3601274) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(0.310251) q[2];
rz(-0.34960738) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67868245) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(2.983685) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(0.37240949) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87293738) q[0];
sx q[0];
rz(-0.71632179) q[0];
sx q[0];
rz(1.3852081) q[0];
x q[1];
rz(2.1521058) q[2];
sx q[2];
rz(-1.2876236) q[2];
sx q[2];
rz(-0.95209939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1540893) q[1];
sx q[1];
rz(-1.1945063) q[1];
sx q[1];
rz(-0.16102468) q[1];
rz(-pi) q[2];
rz(2.1129205) q[3];
sx q[3];
rz(-0.91915932) q[3];
sx q[3];
rz(1.2937014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(0.27238971) q[2];
rz(0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(-0.4549543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(2.0102665) q[0];
rz(-0.5979901) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-2.4647443) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880496) q[0];
sx q[0];
rz(-1.4918259) q[0];
sx q[0];
rz(-2.700564) q[0];
x q[1];
rz(-2.4104775) q[2];
sx q[2];
rz(-1.7321246) q[2];
sx q[2];
rz(2.4246755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21351335) q[1];
sx q[1];
rz(-2.885474) q[1];
sx q[1];
rz(0.10180267) q[1];
rz(-pi) q[2];
x q[2];
rz(0.058651794) q[3];
sx q[3];
rz(-0.39415112) q[3];
sx q[3];
rz(-0.7002206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(0.14906135) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(1.4591249) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-0.93313342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695078) q[0];
sx q[0];
rz(-2.5818995) q[0];
sx q[0];
rz(-2.1493388) q[0];
rz(-pi) q[1];
rz(0.60284166) q[2];
sx q[2];
rz(-1.8800929) q[2];
sx q[2];
rz(-1.3476635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3999403) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(-1.9386577) q[1];
rz(2.9643781) q[3];
sx q[3];
rz(-2.2209077) q[3];
sx q[3];
rz(0.23055102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3559945) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(-0.28665001) q[2];
rz(0.24946985) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0934802) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(-1.4666784) q[0];
rz(-1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(-1.0010304) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769343) q[0];
sx q[0];
rz(-1.5334198) q[0];
sx q[0];
rz(1.2901165) q[0];
rz(1.940735) q[2];
sx q[2];
rz(-1.8897893) q[2];
sx q[2];
rz(-1.621304) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66623464) q[1];
sx q[1];
rz(-1.1821786) q[1];
sx q[1];
rz(2.9734441) q[1];
rz(-pi) q[2];
rz(-1.3592968) q[3];
sx q[3];
rz(-2.2106417) q[3];
sx q[3];
rz(1.1246455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(-3.0380847) q[2];
rz(0.59182709) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(-0.6119588) q[0];
rz(-1.7991964) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-0.38988316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108111) q[0];
sx q[0];
rz(-2.1446052) q[0];
sx q[0];
rz(1.1382889) q[0];
rz(-pi) q[1];
rz(1.1578324) q[2];
sx q[2];
rz(-2.0418842) q[2];
sx q[2];
rz(-0.72380356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0583565) q[1];
sx q[1];
rz(-2.5857158) q[1];
sx q[1];
rz(-0.91169375) q[1];
rz(-pi) q[2];
x q[2];
rz(2.249986) q[3];
sx q[3];
rz(-1.3636175) q[3];
sx q[3];
rz(-2.2694015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3020246) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(0.71425444) q[2];
rz(0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488688) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(2.3642448) q[0];
rz(-2.2946987) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(2.8651967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4682639) q[0];
sx q[0];
rz(-1.4110761) q[0];
sx q[0];
rz(2.7295652) q[0];
rz(-pi) q[1];
rz(1.7644291) q[2];
sx q[2];
rz(-1.7948705) q[2];
sx q[2];
rz(-2.5335238) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8917577) q[1];
sx q[1];
rz(-1.3502305) q[1];
sx q[1];
rz(2.9794429) q[1];
rz(-pi) q[2];
rz(2.5329481) q[3];
sx q[3];
rz(-1.3666271) q[3];
sx q[3];
rz(2.5592786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(-3.1048807) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(-1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-2.7888443) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477662) q[0];
sx q[0];
rz(-2.2802135) q[0];
sx q[0];
rz(-1.312027) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2376386) q[2];
sx q[2];
rz(-0.91234708) q[2];
sx q[2];
rz(1.2308987) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5926338) q[1];
sx q[1];
rz(-2.0083798) q[1];
sx q[1];
rz(-1.2195107) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.685931) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(1.957422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-0.34118787) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108903) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(-2.4717992) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(0.40217051) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(-2.1445198) q[3];
sx q[3];
rz(-2.1763747) q[3];
sx q[3];
rz(2.5381266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

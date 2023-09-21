OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(-0.0052069081) q[0];
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(1.1896689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7548646) q[0];
sx q[0];
rz(-1.0611738) q[0];
sx q[0];
rz(1.1763563) q[0];
rz(-2.9615133) q[2];
sx q[2];
rz(-1.4305978) q[2];
sx q[2];
rz(2.7498498) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7835044) q[1];
sx q[1];
rz(-1.4231829) q[1];
sx q[1];
rz(1.258177) q[1];
x q[2];
rz(-0.40640229) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(1.3960081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71620119) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(2.544196) q[2];
rz(1.3655837) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(-1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4269203) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(-3.0283668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069488) q[0];
sx q[0];
rz(-1.5320677) q[0];
sx q[0];
rz(-1.1452655) q[0];
rz(-pi) q[1];
rz(-1.5510773) q[2];
sx q[2];
rz(-1.4975417) q[2];
sx q[2];
rz(-2.2837226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51616878) q[1];
sx q[1];
rz(-2.5516769) q[1];
sx q[1];
rz(-2.3708458) q[1];
x q[2];
rz(2.1749288) q[3];
sx q[3];
rz(-1.3738487) q[3];
sx q[3];
rz(0.40991022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.018628) q[2];
sx q[2];
rz(-0.48015067) q[2];
sx q[2];
rz(2.9193027) q[2];
rz(-2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46626058) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(-2.0942028) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(-3.0139794) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4532967) q[0];
sx q[0];
rz(-1.3971546) q[0];
sx q[0];
rz(0.88734532) q[0];
x q[1];
rz(-2.3348546) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(-0.15653175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.322791) q[1];
sx q[1];
rz(-1.6530419) q[1];
sx q[1];
rz(-1.8104042) q[1];
rz(-pi) q[2];
rz(-0.25281275) q[3];
sx q[3];
rz(-2.1826934) q[3];
sx q[3];
rz(-1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(0.34960738) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-0.67868245) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(-2.7754916) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(0.37240949) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2686553) q[0];
sx q[0];
rz(-2.4252709) q[0];
sx q[0];
rz(-1.7563845) q[0];
rz(-pi) q[1];
rz(-2.0580975) q[2];
sx q[2];
rz(-0.63939017) q[2];
sx q[2];
rz(-1.0207748) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6652674) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(-1.1900351) q[1];
rz(-2.1129205) q[3];
sx q[3];
rz(-2.2224333) q[3];
sx q[3];
rz(1.2937014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7238414) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(0.27238971) q[2];
rz(-2.2327936) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(-2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(1.1313261) q[0];
rz(2.5436026) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15166053) q[0];
sx q[0];
rz(-0.44758546) q[0];
sx q[0];
rz(-0.18330343) q[0];
rz(0.23907451) q[2];
sx q[2];
rz(-0.74547807) q[2];
sx q[2];
rz(2.4649232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9280793) q[1];
sx q[1];
rz(-2.885474) q[1];
sx q[1];
rz(-3.03979) q[1];
x q[2];
rz(2.7480514) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(-2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1260219) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(2.3964264) q[1];
sx q[1];
rz(-2.7909653) q[1];
sx q[1];
rz(2.2084592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4695078) q[0];
sx q[0];
rz(-0.55969319) q[0];
sx q[0];
rz(-2.1493388) q[0];
rz(-pi) q[1];
rz(-1.2007347) q[2];
sx q[2];
rz(-1.0002631) q[2];
sx q[2];
rz(-2.7119315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16490368) q[1];
sx q[1];
rz(-2.7737244) q[1];
sx q[1];
rz(-1.5644172) q[1];
x q[2];
rz(1.3429787) q[3];
sx q[3];
rz(-0.67043793) q[3];
sx q[3];
rz(-0.05712856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3559945) q[2];
sx q[2];
rz(-2.9170673) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0934802) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(-1.4666784) q[0];
rz(1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(1.0010304) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46465835) q[0];
sx q[0];
rz(-1.5334198) q[0];
sx q[0];
rz(1.8514762) q[0];
rz(-pi) q[1];
rz(-1.2008576) q[2];
sx q[2];
rz(-1.2518034) q[2];
sx q[2];
rz(-1.5202886) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8965473) q[1];
sx q[1];
rz(-0.42173112) q[1];
sx q[1];
rz(-1.1827724) q[1];
rz(-pi) q[2];
rz(-2.8666897) q[3];
sx q[3];
rz(-0.66920815) q[3];
sx q[3];
rz(-1.6717403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(3.0380847) q[2];
rz(0.59182709) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(-2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(-0.6119588) q[0];
rz(-1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-2.7517095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0355426) q[0];
sx q[0];
rz(-0.70362008) q[0];
sx q[0];
rz(0.57530595) q[0];
x q[1];
rz(2.4742485) q[2];
sx q[2];
rz(-2.5255425) q[2];
sx q[2];
rz(-1.6499856) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4854359) q[1];
sx q[1];
rz(-1.1404783) q[1];
sx q[1];
rz(2.7780611) q[1];
rz(-pi) q[2];
rz(2.8777458) q[3];
sx q[3];
rz(-0.90875328) q[3];
sx q[3];
rz(2.2784233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8395681) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-0.71425444) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(-2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(0.27639595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24647507) q[0];
sx q[0];
rz(-0.44024375) q[0];
sx q[0];
rz(-2.759139) q[0];
rz(-0.70119621) q[2];
sx q[2];
rz(-2.8465135) q[2];
sx q[2];
rz(1.3311177) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8564057) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(-1.3473947) q[1];
x q[2];
rz(2.7941462) q[3];
sx q[3];
rz(-2.5037519) q[3];
sx q[3];
rz(-0.70536246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(3.1048807) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(2.7888443) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(1.030285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.83537) q[0];
sx q[0];
rz(-1.3754002) q[0];
sx q[0];
rz(-2.4154001) q[0];
rz(-pi) q[1];
rz(-2.7416517) q[2];
sx q[2];
rz(-0.72657864) q[2];
sx q[2];
rz(2.4253997) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5926338) q[1];
sx q[1];
rz(-2.0083798) q[1];
sx q[1];
rz(-1.9220819) q[1];
rz(1.9029721) q[3];
sx q[3];
rz(-2.1579451) q[3];
sx q[3];
rz(1.4004933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108903) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(2.4717992) q[1];
sx q[1];
rz(-2.5302946) q[1];
sx q[1];
rz(-1.5763462) q[1];
rz(1.7651991) q[2];
sx q[2];
rz(-1.9971078) q[2];
sx q[2];
rz(2.8304921) q[2];
rz(0.66486799) q[3];
sx q[3];
rz(-2.3330199) q[3];
sx q[3];
rz(-1.4521269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
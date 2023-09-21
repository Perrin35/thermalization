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
rz(3.1363857) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(-1.985328) q[1];
sx q[1];
rz(-1.1896689) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.386728) q[0];
sx q[0];
rz(-2.0804188) q[0];
sx q[0];
rz(1.1763563) q[0];
rz(0.66733811) q[2];
sx q[2];
rz(-0.22775209) q[2];
sx q[2];
rz(1.8337133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3559239) q[1];
sx q[1];
rz(-2.7969116) q[1];
sx q[1];
rz(-2.0211401) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7351904) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(-1.3960081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71620119) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(0.5973967) q[2];
rz(1.776009) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(-0.11322583) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118606) q[0];
sx q[0];
rz(-1.9959873) q[0];
sx q[0];
rz(3.0990764) q[0];
rz(-pi) q[1];
rz(1.5905154) q[2];
sx q[2];
rz(-1.644051) q[2];
sx q[2];
rz(-0.8578701) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7659521) q[1];
sx q[1];
rz(-1.1728219) q[1];
sx q[1];
rz(2.6938733) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1749288) q[3];
sx q[3];
rz(-1.767744) q[3];
sx q[3];
rz(-2.7316824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.018628) q[2];
sx q[2];
rz(-0.48015067) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(-2.9120581) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
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
rz(0.12761322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.815044) q[0];
sx q[0];
rz(-0.70171261) q[0];
sx q[0];
rz(-1.2998507) q[0];
x q[1];
rz(-0.87330841) q[2];
sx q[2];
rz(-2.1224988) q[2];
sx q[2];
rz(1.973195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0652005) q[1];
sx q[1];
rz(-2.8885191) q[1];
sx q[1];
rz(1.2364926) q[1];
rz(0.94362887) q[3];
sx q[3];
rz(-1.3645932) q[3];
sx q[3];
rz(-0.19087736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(-2.7919853) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(0.37240949) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62896699) q[0];
sx q[0];
rz(-0.86932875) q[0];
sx q[0];
rz(-0.15928282) q[0];
x q[1];
rz(0.33505586) q[2];
sx q[2];
rz(-2.1261566) q[2];
sx q[2];
rz(-2.7044538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47632521) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(1.9515576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59471547) q[3];
sx q[3];
rz(-2.3200431) q[3];
sx q[3];
rz(-0.51175129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(2.2327936) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(-2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(-2.0102665) q[0];
rz(-2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-2.4647443) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8880496) q[0];
sx q[0];
rz(-1.4918259) q[0];
sx q[0];
rz(2.700564) q[0];
rz(-pi) q[1];
rz(-2.9025181) q[2];
sx q[2];
rz(-2.3961146) q[2];
sx q[2];
rz(-2.4649232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9280793) q[1];
sx q[1];
rz(-0.25611862) q[1];
sx q[1];
rz(-3.03979) q[1];
x q[2];
rz(1.5951717) q[3];
sx q[3];
rz(-1.9642324) q[3];
sx q[3];
rz(-0.76373053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(-2.9925313) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(-1.0173652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(2.2084592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81290302) q[0];
sx q[0];
rz(-2.0314386) q[0];
sx q[0];
rz(0.33005379) q[0];
rz(-pi) q[1];
rz(2.538751) q[2];
sx q[2];
rz(-1.8800929) q[2];
sx q[2];
rz(1.3476635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9698525) q[1];
sx q[1];
rz(-1.2029359) q[1];
sx q[1];
rz(-0.0024585558) q[1];
rz(-pi) q[2];
rz(-0.17721456) q[3];
sx q[3];
rz(-0.92068499) q[3];
sx q[3];
rz(2.9110416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(-2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0934802) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(1.4666784) q[0];
rz(2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(1.0010304) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2350378) q[0];
sx q[0];
rz(-2.8585003) q[0];
sx q[0];
rz(-1.4366158) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2008576) q[2];
sx q[2];
rz(-1.2518034) q[2];
sx q[2];
rz(1.621304) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96879362) q[1];
sx q[1];
rz(-1.4152923) q[1];
sx q[1];
rz(1.1771727) q[1];
rz(1.3592968) q[3];
sx q[3];
rz(-0.930951) q[3];
sx q[3];
rz(-2.0169472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.137407) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(-3.0380847) q[2];
rz(-2.5497656) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-0.83759585) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(0.6119588) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-2.7517095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33078157) q[0];
sx q[0];
rz(-0.99698742) q[0];
sx q[0];
rz(-1.1382889) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4742485) q[2];
sx q[2];
rz(-0.61605011) q[2];
sx q[2];
rz(-1.4916071) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0583565) q[1];
sx q[1];
rz(-2.5857158) q[1];
sx q[1];
rz(-0.91169375) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8777458) q[3];
sx q[3];
rz(-2.2328394) q[3];
sx q[3];
rz(-0.86316934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8395681) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(2.4273382) q[2];
rz(0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-0.57202655) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927239) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(-2.2946987) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(0.27639595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8951176) q[0];
sx q[0];
rz(-2.7013489) q[0];
sx q[0];
rz(2.759139) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70119621) q[2];
sx q[2];
rz(-2.8465135) q[2];
sx q[2];
rz(1.3311177) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8564057) q[1];
sx q[1];
rz(-1.4126084) q[1];
sx q[1];
rz(1.7941979) q[1];
rz(1.3235839) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(-2.0127399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85236621) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(-3.0140871) q[2];
rz(0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7413095) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(2.1113077) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093826483) q[0];
sx q[0];
rz(-2.2802135) q[0];
sx q[0];
rz(-1.8295656) q[0];
rz(0.68600168) q[2];
sx q[2];
rz(-1.8324319) q[2];
sx q[2];
rz(-0.54856448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2620169) q[1];
sx q[1];
rz(-0.55393065) q[1];
sx q[1];
rz(-0.63417706) q[1];
rz(-0.61334765) q[3];
sx q[3];
rz(-1.2958591) q[3];
sx q[3];
rz(3.1230694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-0.34118787) q[2];
rz(0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
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
rz(0.13070233) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(-2.4717992) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(1.7651991) q[2];
sx q[2];
rz(-1.9971078) q[2];
sx q[2];
rz(2.8304921) q[2];
rz(-2.4521811) q[3];
sx q[3];
rz(-2.0333615) q[3];
sx q[3];
rz(0.61483308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
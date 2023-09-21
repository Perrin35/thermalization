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
rz(1.7967254) q[1];
sx q[1];
rz(4.2978573) q[1];
sx q[1];
rz(8.2351091) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016276377) q[0];
sx q[0];
rz(-1.2287041) q[0];
sx q[0];
rz(-2.5972511) q[0];
x q[1];
rz(2.9615133) q[2];
sx q[2];
rz(-1.4305978) q[2];
sx q[2];
rz(0.3917429) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7835044) q[1];
sx q[1];
rz(-1.4231829) q[1];
sx q[1];
rz(-1.258177) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7351904) q[3];
sx q[3];
rz(-2.1592525) q[3];
sx q[3];
rz(-1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(0.5973967) q[2];
rz(-1.776009) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(0.11322583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2213759) q[0];
sx q[0];
rz(-2.7144103) q[0];
sx q[0];
rz(1.4772052) q[0];
rz(-pi) q[1];
rz(2.8791061) q[2];
sx q[2];
rz(-0.075857698) q[2];
sx q[2];
rz(0.59466098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7624224) q[1];
sx q[1];
rz(-1.1602853) q[1];
sx q[1];
rz(-1.1344086) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9037335) q[3];
sx q[3];
rz(-0.97994643) q[3];
sx q[3];
rz(-1.295134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(-2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46626058) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(1.0473898) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(-3.0139794) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532967) q[0];
sx q[0];
rz(-1.3971546) q[0];
sx q[0];
rz(-0.88734532) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3348546) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(2.9850609) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.81880169) q[1];
sx q[1];
rz(-1.4885508) q[1];
sx q[1];
rz(1.8104042) q[1];
rz(2.8887799) q[3];
sx q[3];
rz(-2.1826934) q[3];
sx q[3];
rz(-1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3601274) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(0.310251) q[2];
rz(0.34960738) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(2.983685) q[0];
rz(2.7754916) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(2.7691832) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5126257) q[0];
sx q[0];
rz(-0.86932875) q[0];
sx q[0];
rz(-0.15928282) q[0];
rz(-pi) q[1];
rz(2.8065368) q[2];
sx q[2];
rz(-1.0154361) q[2];
sx q[2];
rz(-2.7044538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73788961) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(-1.1853192) q[1];
x q[2];
rz(2.4140671) q[3];
sx q[3];
rz(-1.993506) q[3];
sx q[3];
rz(-2.5142575) q[3];
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
rz(-0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(0.4549543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.882778) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(2.0102665) q[0];
rz(-0.5979901) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-0.67684832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15166053) q[0];
sx q[0];
rz(-0.44758546) q[0];
sx q[0];
rz(2.9582892) q[0];
rz(-pi) q[1];
rz(2.9025181) q[2];
sx q[2];
rz(-2.3961146) q[2];
sx q[2];
rz(2.4649232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21351335) q[1];
sx q[1];
rz(-2.885474) q[1];
sx q[1];
rz(3.03979) q[1];
rz(-2.7480514) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(-2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.9963025) q[2];
rz(0.14906135) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(-1.0173652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260219) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(2.2084592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60676735) q[0];
sx q[0];
rz(-1.8653499) q[0];
sx q[0];
rz(-2.0539001) q[0];
x q[1];
rz(-1.9408579) q[2];
sx q[2];
rz(-2.1413295) q[2];
sx q[2];
rz(0.42966118) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.16490368) q[1];
sx q[1];
rz(-0.36786825) q[1];
sx q[1];
rz(1.5771754) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91306367) q[3];
sx q[3];
rz(-1.7115895) q[3];
sx q[3];
rz(-1.6933683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3559945) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(-0.24946985) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(-0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0934802) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(-1.4666784) q[0];
rz(-2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(2.1405623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46465835) q[0];
sx q[0];
rz(-1.5334198) q[0];
sx q[0];
rz(-1.8514762) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2008576) q[2];
sx q[2];
rz(-1.8897893) q[2];
sx q[2];
rz(1.621304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.172799) q[1];
sx q[1];
rz(-1.7263004) q[1];
sx q[1];
rz(-1.96442) q[1];
rz(-1.7822958) q[3];
sx q[3];
rz(-0.930951) q[3];
sx q[3];
rz(1.1246455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.004185685) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(-0.10350791) q[2];
rz(-2.5497656) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(-0.38988316) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.10605) q[0];
sx q[0];
rz(-0.70362008) q[0];
sx q[0];
rz(-0.57530595) q[0];
x q[1];
rz(0.50750081) q[2];
sx q[2];
rz(-1.9365053) q[2];
sx q[2];
rz(-2.4909004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65615678) q[1];
sx q[1];
rz(-1.1404783) q[1];
sx q[1];
rz(-0.3635316) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(-0.3020246) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(-2.4273382) q[2];
rz(0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927239) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(0.77734787) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-2.8651967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8951176) q[0];
sx q[0];
rz(-2.7013489) q[0];
sx q[0];
rz(-2.759139) q[0];
rz(2.4403964) q[2];
sx q[2];
rz(-0.29507911) q[2];
sx q[2];
rz(-1.3311177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8918708) q[1];
sx q[1];
rz(-2.8686214) q[1];
sx q[1];
rz(0.94675468) q[1];
rz(-pi) q[2];
rz(-0.3474465) q[3];
sx q[3];
rz(-0.63784079) q[3];
sx q[3];
rz(0.70536246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(-1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093826483) q[0];
sx q[0];
rz(-2.2802135) q[0];
sx q[0];
rz(1.312027) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39994098) q[2];
sx q[2];
rz(-0.72657864) q[2];
sx q[2];
rz(2.4253997) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2620169) q[1];
sx q[1];
rz(-0.55393065) q[1];
sx q[1];
rz(-0.63417706) q[1];
rz(-pi) q[2];
rz(1.9029721) q[3];
sx q[3];
rz(-0.98364753) q[3];
sx q[3];
rz(1.7410994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-2.7080766) q[2];
sx q[2];
rz(-1.3939861) q[2];
sx q[2];
rz(1.3409333) q[2];
rz(2.4767247) q[3];
sx q[3];
rz(-0.80857279) q[3];
sx q[3];
rz(1.6894658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
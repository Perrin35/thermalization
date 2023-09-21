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
rz(-1.9519238) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0931041) q[0];
sx q[0];
rz(-0.63360533) q[0];
sx q[0];
rz(-2.5392169) q[0];
rz(0.66733811) q[2];
sx q[2];
rz(-2.9138406) q[2];
sx q[2];
rz(1.3078794) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2602106) q[1];
sx q[1];
rz(-1.8799026) q[1];
sx q[1];
rz(0.15501546) q[1];
rz(-pi) q[2];
rz(-2.1055929) q[3];
sx q[3];
rz(-2.4403604) q[3];
sx q[3];
rz(-1.086077) q[3];
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
rz(-0.5973967) q[2];
rz(1.3655837) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(1.8610154) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4269203) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(0.11322583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0229867) q[0];
sx q[0];
rz(-1.9959873) q[0];
sx q[0];
rz(-3.0990764) q[0];
rz(-pi) q[1];
rz(3.0683238) q[2];
sx q[2];
rz(-1.5904625) q[2];
sx q[2];
rz(0.71436963) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37564056) q[1];
sx q[1];
rz(-1.9687708) q[1];
sx q[1];
rz(2.6938733) q[1];
rz(-pi) q[2];
rz(-2.9037335) q[3];
sx q[3];
rz(-0.97994643) q[3];
sx q[3];
rz(1.295134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1229646) q[2];
sx q[2];
rz(-0.48015067) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46626058) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(-2.0942028) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(3.0139794) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022284431) q[0];
sx q[0];
rz(-2.2420609) q[0];
sx q[0];
rz(0.22247252) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87330841) q[2];
sx q[2];
rz(-1.0190939) q[2];
sx q[2];
rz(1.1683977) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0652005) q[1];
sx q[1];
rz(-2.8885191) q[1];
sx q[1];
rz(-1.9051001) q[1];
rz(-pi) q[2];
rz(1.913194) q[3];
sx q[3];
rz(-0.6558334) q[3];
sx q[3];
rz(2.0369903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(-0.34960738) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(-1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(2.7691832) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62896699) q[0];
sx q[0];
rz(-2.2722639) q[0];
sx q[0];
rz(-0.15928282) q[0];
rz(-2.1521058) q[2];
sx q[2];
rz(-1.8539691) q[2];
sx q[2];
rz(2.1894933) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47632521) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(1.9515576) q[1];
x q[2];
rz(2.1129205) q[3];
sx q[3];
rz(-2.2224333) q[3];
sx q[3];
rz(1.8478912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(-0.90879905) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(-0.4549543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(-1.1313261) q[0];
rz(-2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2535431) q[0];
sx q[0];
rz(-1.4918259) q[0];
sx q[0];
rz(0.44102863) q[0];
rz(-1.3555688) q[2];
sx q[2];
rz(-2.2903246) q[2];
sx q[2];
rz(2.1446251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0332898) q[1];
sx q[1];
rz(-1.8255594) q[1];
sx q[1];
rz(1.5441896) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5464209) q[3];
sx q[3];
rz(-1.9642324) q[3];
sx q[3];
rz(-2.3778621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.9963025) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(-0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(0.93313342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60676735) q[0];
sx q[0];
rz(-1.2762428) q[0];
sx q[0];
rz(1.0876925) q[0];
rz(-pi) q[1];
rz(0.51322333) q[2];
sx q[2];
rz(-0.66868082) q[2];
sx q[2];
rz(0.19323397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7416523) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(1.2029349) q[1];
rz(0.91306367) q[3];
sx q[3];
rz(-1.4300031) q[3];
sx q[3];
rz(-1.4482244) q[3];
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
rz(-2.8921228) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0934802) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(-1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(-1.0010304) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2350378) q[0];
sx q[0];
rz(-0.28309238) q[0];
sx q[0];
rz(1.4366158) q[0];
rz(-pi) q[1];
rz(0.83058968) q[2];
sx q[2];
rz(-2.6579654) q[2];
sx q[2];
rz(-2.4110576) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24504532) q[1];
sx q[1];
rz(-0.42173112) q[1];
sx q[1];
rz(-1.9588203) q[1];
rz(0.6506728) q[3];
sx q[3];
rz(-1.4015897) q[3];
sx q[3];
rz(-0.31864877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.004185685) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(3.0380847) q[2];
rz(-2.5497656) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(0.6119588) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(2.7517095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.10605) q[0];
sx q[0];
rz(-2.4379726) q[0];
sx q[0];
rz(2.5662867) q[0];
x q[1];
rz(0.6673442) q[2];
sx q[2];
rz(-0.61605011) q[2];
sx q[2];
rz(1.4916071) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0695614) q[1];
sx q[1];
rz(-1.2417292) q[1];
sx q[1];
rz(1.1142932) q[1];
rz(-pi) q[2];
rz(0.89160664) q[3];
sx q[3];
rz(-1.3636175) q[3];
sx q[3];
rz(-0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8395681) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(0.71425444) q[2];
rz(-0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927239) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(0.77734787) q[0];
rz(0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-0.27639595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24647507) q[0];
sx q[0];
rz(-2.7013489) q[0];
sx q[0];
rz(-0.38245364) q[0];
rz(-2.4403964) q[2];
sx q[2];
rz(-2.8465135) q[2];
sx q[2];
rz(-1.3311177) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8918708) q[1];
sx q[1];
rz(-0.2729713) q[1];
sx q[1];
rz(0.94675468) q[1];
rz(-pi) q[2];
rz(-0.60864457) q[3];
sx q[3];
rz(-1.7749655) q[3];
sx q[3];
rz(0.58231402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(-0.036711983) q[3];
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
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47980967) q[0];
sx q[0];
rz(-2.3942238) q[0];
sx q[0];
rz(-2.8519147) q[0];
rz(2.7416517) q[2];
sx q[2];
rz(-0.72657864) q[2];
sx q[2];
rz(0.71619294) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5489588) q[1];
sx q[1];
rz(-2.0083798) q[1];
sx q[1];
rz(-1.9220819) q[1];
x q[2];
rz(-2.528245) q[3];
sx q[3];
rz(-1.8457335) q[3];
sx q[3];
rz(3.1230694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(2.8004048) q[2];
rz(0.7406922) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(3.0497131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13070233) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(-2.7394221) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(-0.99707281) q[3];
sx q[3];
rz(-0.96521796) q[3];
sx q[3];
rz(-0.60346606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

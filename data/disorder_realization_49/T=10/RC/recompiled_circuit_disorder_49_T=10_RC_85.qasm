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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016276377) q[0];
sx q[0];
rz(-1.9128886) q[0];
sx q[0];
rz(0.54434158) q[0];
x q[1];
rz(1.7132684) q[2];
sx q[2];
rz(-1.3925029) q[2];
sx q[2];
rz(-1.1536191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2602106) q[1];
sx q[1];
rz(-1.26169) q[1];
sx q[1];
rz(0.15501546) q[1];
rz(-pi) q[2];
rz(-2.7351904) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(-1.3960081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(-0.5973967) q[2];
rz(1.776009) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(-0.33357099) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(-0.11322583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118606) q[0];
sx q[0];
rz(-1.1456053) q[0];
sx q[0];
rz(3.0990764) q[0];
rz(-pi) q[1];
rz(-3.0683238) q[2];
sx q[2];
rz(-1.5511302) q[2];
sx q[2];
rz(-2.427223) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7624224) q[1];
sx q[1];
rz(-1.9813073) q[1];
sx q[1];
rz(1.1344086) q[1];
rz(0.23785915) q[3];
sx q[3];
rz(-2.1616462) q[3];
sx q[3];
rz(1.8464586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.018628) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.46626058) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(-0.88071841) q[0];
rz(1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(3.0139794) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32654861) q[0];
sx q[0];
rz(-0.70171261) q[0];
sx q[0];
rz(-1.841742) q[0];
rz(0.87330841) q[2];
sx q[2];
rz(-2.1224988) q[2];
sx q[2];
rz(1.1683977) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.322791) q[1];
sx q[1];
rz(-1.4885508) q[1];
sx q[1];
rz(-1.8104042) q[1];
rz(-pi) q[2];
rz(1.2283986) q[3];
sx q[3];
rz(-0.6558334) q[3];
sx q[3];
rz(1.1046023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-2.8313417) q[2];
rz(0.34960738) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.7278956) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4629102) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(2.983685) q[0];
rz(2.7754916) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(2.7691832) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62896699) q[0];
sx q[0];
rz(-0.86932875) q[0];
sx q[0];
rz(-0.15928282) q[0];
rz(2.1521058) q[2];
sx q[2];
rz(-1.2876236) q[2];
sx q[2];
rz(2.1894933) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.73788961) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(1.9562734) q[1];
rz(2.1129205) q[3];
sx q[3];
rz(-0.91915932) q[3];
sx q[3];
rz(-1.8478912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(-0.27238971) q[2];
rz(0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.882778) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(2.0102665) q[0];
rz(2.5436026) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35447435) q[0];
sx q[0];
rz(-1.1312383) q[0];
sx q[0];
rz(-1.6580824) q[0];
rz(-pi) q[1];
rz(-1.7860239) q[2];
sx q[2];
rz(-2.2903246) q[2];
sx q[2];
rz(-2.1446251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21351335) q[1];
sx q[1];
rz(-0.25611862) q[1];
sx q[1];
rz(-3.03979) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5951717) q[3];
sx q[3];
rz(-1.1773603) q[3];
sx q[3];
rz(2.3778621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(-2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(-2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-2.7909653) q[1];
sx q[1];
rz(-2.2084592) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6720848) q[0];
sx q[0];
rz(-0.55969319) q[0];
sx q[0];
rz(2.1493388) q[0];
rz(-2.538751) q[2];
sx q[2];
rz(-1.8800929) q[2];
sx q[2];
rz(-1.3476635) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7416523) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(1.9386577) q[1];
rz(0.91306367) q[3];
sx q[3];
rz(-1.4300031) q[3];
sx q[3];
rz(1.6933683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(0.24946985) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(-2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(-2.1215227) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(-2.1405623) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046227) q[0];
sx q[0];
rz(-1.2903178) q[0];
sx q[0];
rz(3.1026955) q[0];
x q[1];
rz(1.940735) q[2];
sx q[2];
rz(-1.8897893) q[2];
sx q[2];
rz(-1.621304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.172799) q[1];
sx q[1];
rz(-1.7263004) q[1];
sx q[1];
rz(-1.96442) q[1];
x q[2];
rz(-1.3592968) q[3];
sx q[3];
rz(-2.2106417) q[3];
sx q[3];
rz(1.1246455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(-0.10350791) q[2];
rz(-2.5497656) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(-2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.9806769) q[1];
sx q[1];
rz(-2.7517095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.10605) q[0];
sx q[0];
rz(-2.4379726) q[0];
sx q[0];
rz(-0.57530595) q[0];
rz(-2.6340918) q[2];
sx q[2];
rz(-1.2050873) q[2];
sx q[2];
rz(2.4909004) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0695614) q[1];
sx q[1];
rz(-1.2417292) q[1];
sx q[1];
rz(-2.0272994) q[1];
rz(1.893703) q[3];
sx q[3];
rz(-0.70526037) q[3];
sx q[3];
rz(-0.44912072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-0.71425444) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927239) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(-2.3642448) q[0];
rz(2.2946987) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(-0.27639595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24647507) q[0];
sx q[0];
rz(-2.7013489) q[0];
sx q[0];
rz(-0.38245364) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3771636) q[2];
sx q[2];
rz(-1.3467222) q[2];
sx q[2];
rz(0.60806882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2497219) q[1];
sx q[1];
rz(-2.8686214) q[1];
sx q[1];
rz(-0.94675468) q[1];
x q[2];
rz(1.3235839) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(-2.0127399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2892264) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(-3.1048807) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(-1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(1.030285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477662) q[0];
sx q[0];
rz(-0.86137912) q[0];
sx q[0];
rz(1.312027) q[0];
x q[1];
rz(1.9039541) q[2];
sx q[2];
rz(-0.91234708) q[2];
sx q[2];
rz(1.2308987) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0093617) q[1];
sx q[1];
rz(-1.8877601) q[1];
sx q[1];
rz(-2.6793381) q[1];
x q[2];
rz(0.45566166) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(1.957422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(0.34118787) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-0.091879524) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13070233) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(2.4717992) q[1];
sx q[1];
rz(-2.5302946) q[1];
sx q[1];
rz(-1.5763462) q[1];
rz(0.40217051) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(0.99707281) q[3];
sx q[3];
rz(-2.1763747) q[3];
sx q[3];
rz(2.5381266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

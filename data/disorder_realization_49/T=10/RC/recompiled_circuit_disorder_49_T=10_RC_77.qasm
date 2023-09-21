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
rz(1.1896689) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016276377) q[0];
sx q[0];
rz(-1.9128886) q[0];
sx q[0];
rz(-2.5972511) q[0];
x q[1];
rz(2.4742545) q[2];
sx q[2];
rz(-2.9138406) q[2];
sx q[2];
rz(1.8337133) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3559239) q[1];
sx q[1];
rz(-0.34468109) q[1];
sx q[1];
rz(2.0211401) q[1];
x q[2];
rz(-0.40640229) q[3];
sx q[3];
rz(-2.1592525) q[3];
sx q[3];
rz(1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.71620119) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(2.544196) q[2];
rz(-1.3655837) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(1.8610154) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(3.0283668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2213759) q[0];
sx q[0];
rz(-2.7144103) q[0];
sx q[0];
rz(-1.4772052) q[0];
x q[1];
rz(1.5510773) q[2];
sx q[2];
rz(-1.4975417) q[2];
sx q[2];
rz(-0.8578701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7659521) q[1];
sx q[1];
rz(-1.9687708) q[1];
sx q[1];
rz(-0.44771938) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9085957) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(0.88463569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(-1.0473898) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(-0.12761322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32654861) q[0];
sx q[0];
rz(-2.43988) q[0];
sx q[0];
rz(-1.841742) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3348546) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(-0.15653175) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81880169) q[1];
sx q[1];
rz(-1.4885508) q[1];
sx q[1];
rz(-1.8104042) q[1];
rz(-pi) q[2];
rz(0.25281275) q[3];
sx q[3];
rz(-2.1826934) q[3];
sx q[3];
rz(1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(0.310251) q[2];
rz(2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(-1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67868245) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(0.15790766) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(-0.37240949) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62896699) q[0];
sx q[0];
rz(-0.86932875) q[0];
sx q[0];
rz(-2.9823098) q[0];
x q[1];
rz(-1.0834951) q[2];
sx q[2];
rz(-2.5022025) q[2];
sx q[2];
rz(-1.0207748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73788961) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(-1.1853192) q[1];
rz(2.4140671) q[3];
sx q[3];
rz(-1.993506) q[3];
sx q[3];
rz(-2.5142575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7238414) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(-0.27238971) q[2];
rz(-0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(0.4549543) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(-1.1313261) q[0];
rz(0.5979901) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-2.4647443) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2535431) q[0];
sx q[0];
rz(-1.4918259) q[0];
sx q[0];
rz(-2.700564) q[0];
rz(-pi) q[1];
rz(-1.7860239) q[2];
sx q[2];
rz(-0.85126801) q[2];
sx q[2];
rz(-0.99696751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.685806) q[1];
sx q[1];
rz(-1.5965441) q[1];
sx q[1];
rz(0.2548494) q[1];
rz(0.39354126) q[3];
sx q[3];
rz(-1.593309) q[3];
sx q[3];
rz(-0.81641203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0040434917) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(-0.14906135) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(1.4591249) q[0];
rz(-2.3964264) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-0.93313342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6720848) q[0];
sx q[0];
rz(-0.55969319) q[0];
sx q[0];
rz(-0.99225386) q[0];
x q[1];
rz(-0.60284166) q[2];
sx q[2];
rz(-1.2614998) q[2];
sx q[2];
rz(1.7939292) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7416523) q[1];
sx q[1];
rz(-1.5730904) q[1];
sx q[1];
rz(-1.9386577) q[1];
rz(-pi) q[2];
x q[2];
rz(2.228529) q[3];
sx q[3];
rz(-1.7115895) q[3];
sx q[3];
rz(-1.4482244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(-2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-0.048112415) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.4666784) q[0];
rz(-2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(2.1405623) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6769343) q[0];
sx q[0];
rz(-1.5334198) q[0];
sx q[0];
rz(-1.2901165) q[0];
rz(-pi) q[1];
rz(1.2008576) q[2];
sx q[2];
rz(-1.8897893) q[2];
sx q[2];
rz(-1.5202886) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.475358) q[1];
sx q[1];
rz(-1.1821786) q[1];
sx q[1];
rz(-2.9734441) q[1];
x q[2];
rz(-0.6506728) q[3];
sx q[3];
rz(-1.7400029) q[3];
sx q[3];
rz(-0.31864877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.004185685) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(0.10350791) q[2];
rz(-2.5497656) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-0.83759585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(-2.5296339) q[0];
rz(-1.3423963) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(-0.38988316) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99446699) q[0];
sx q[0];
rz(-1.9305221) q[0];
sx q[0];
rz(2.5229182) q[0];
rz(-pi) q[1];
rz(-1.1578324) q[2];
sx q[2];
rz(-2.0418842) q[2];
sx q[2];
rz(-2.4177891) q[2];
rz(-pi) q[3];
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
rz(-0.26384683) q[3];
sx q[3];
rz(-2.2328394) q[3];
sx q[3];
rz(-2.2784233) q[3];
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
rz(2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927239) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(2.8651967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696635) q[0];
sx q[0];
rz(-1.164325) q[0];
sx q[0];
rz(-1.3967692) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9133965) q[2];
sx q[2];
rz(-1.3820634) q[2];
sx q[2];
rz(-0.91918321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8917577) q[1];
sx q[1];
rz(-1.3502305) q[1];
sx q[1];
rz(0.16214976) q[1];
rz(2.7941462) q[3];
sx q[3];
rz(-2.5037519) q[3];
sx q[3];
rz(2.4362302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(-0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-0.35274831) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(-1.030285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47980967) q[0];
sx q[0];
rz(-0.74736887) q[0];
sx q[0];
rz(-0.28967793) q[0];
rz(-pi) q[1];
rz(-0.68600168) q[2];
sx q[2];
rz(-1.8324319) q[2];
sx q[2];
rz(0.54856448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2620169) q[1];
sx q[1];
rz(-0.55393065) q[1];
sx q[1];
rz(0.63417706) q[1];
rz(1.9029721) q[3];
sx q[3];
rz(-0.98364753) q[3];
sx q[3];
rz(-1.4004933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-0.34118787) q[2];
rz(-2.4009005) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(-0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13070233) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(-2.7080766) q[2];
sx q[2];
rz(-1.3939861) q[2];
sx q[2];
rz(1.3409333) q[2];
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
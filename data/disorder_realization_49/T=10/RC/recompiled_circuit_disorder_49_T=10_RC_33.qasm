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
rz(4.0806169) q[0];
sx q[0];
rz(9.4299849) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(4.2978573) q[1];
sx q[1];
rz(8.2351091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016276377) q[0];
sx q[0];
rz(-1.9128886) q[0];
sx q[0];
rz(-0.54434158) q[0];
rz(-pi) q[1];
x q[1];
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
x q[0];
rz(-0.78566879) q[1];
sx q[1];
rz(-0.34468109) q[1];
sx q[1];
rz(2.0211401) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0359997) q[3];
sx q[3];
rz(-0.70123226) q[3];
sx q[3];
rz(-2.0555156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
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
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(1.0936273) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(0.11322583) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0229867) q[0];
sx q[0];
rz(-1.1456053) q[0];
sx q[0];
rz(-0.042516275) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8791061) q[2];
sx q[2];
rz(-0.075857698) q[2];
sx q[2];
rz(2.5469317) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7624224) q[1];
sx q[1];
rz(-1.9813073) q[1];
sx q[1];
rz(1.1344086) q[1];
rz(-pi) q[2];
rz(1.232997) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(2.256957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.018628) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(-2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(3.0632609) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46626058) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(-1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(0.12761322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1193082) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(-0.22247252) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67655501) q[2];
sx q[2];
rz(-0.99202079) q[2];
sx q[2];
rz(-2.3253331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0652005) q[1];
sx q[1];
rz(-2.8885191) q[1];
sx q[1];
rz(-1.9051001) q[1];
rz(-1.913194) q[3];
sx q[3];
rz(-2.4857593) q[3];
sx q[3];
rz(-1.1046023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3601274) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(2.8313417) q[2];
rz(2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67868245) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(2.7691832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62896699) q[0];
sx q[0];
rz(-0.86932875) q[0];
sx q[0];
rz(0.15928282) q[0];
rz(2.8065368) q[2];
sx q[2];
rz(-1.0154361) q[2];
sx q[2];
rz(-2.7044538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73788961) q[1];
sx q[1];
rz(-2.7338114) q[1];
sx q[1];
rz(1.1853192) q[1];
rz(1.0286721) q[3];
sx q[3];
rz(-2.2224333) q[3];
sx q[3];
rz(-1.8478912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(-2.2327936) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(-2.0102665) q[0];
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
rz(1.2535431) q[0];
sx q[0];
rz(-1.6497668) q[0];
sx q[0];
rz(-0.44102863) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7860239) q[2];
sx q[2];
rz(-2.2903246) q[2];
sx q[2];
rz(-2.1446251) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4557867) q[1];
sx q[1];
rz(-1.5450486) q[1];
sx q[1];
rz(0.2548494) q[1];
x q[2];
rz(0.39354126) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(0.81641203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.9963025) q[2];
rz(-2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(-2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260219) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-2.7909653) q[1];
sx q[1];
rz(-2.2084592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81290302) q[0];
sx q[0];
rz(-2.0314386) q[0];
sx q[0];
rz(-2.8115389) q[0];
rz(-pi) q[1];
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
rz(-pi) q[0];
x q[0];
rz(-2.9698525) q[1];
sx q[1];
rz(-1.9386567) q[1];
sx q[1];
rz(-0.0024585558) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3429787) q[3];
sx q[3];
rz(-2.4711547) q[3];
sx q[3];
rz(-3.0844641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(-0.28665001) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.048112415) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.4666784) q[0];
rz(2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(-2.1405623) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046227) q[0];
sx q[0];
rz(-1.2903178) q[0];
sx q[0];
rz(-3.1026955) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8011488) q[2];
sx q[2];
rz(-1.2203487) q[2];
sx q[2];
rz(0.0705138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24504532) q[1];
sx q[1];
rz(-0.42173112) q[1];
sx q[1];
rz(1.9588203) q[1];
rz(-pi) q[2];
rz(1.7822958) q[3];
sx q[3];
rz(-2.2106417) q[3];
sx q[3];
rz(1.1246455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.004185685) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(-0.10350791) q[2];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(0.6119588) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(-0.38988316) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1471257) q[0];
sx q[0];
rz(-1.2110706) q[0];
sx q[0];
rz(0.61867449) q[0];
x q[1];
rz(-0.50750081) q[2];
sx q[2];
rz(-1.9365053) q[2];
sx q[2];
rz(-0.65069228) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0695614) q[1];
sx q[1];
rz(-1.8998635) q[1];
sx q[1];
rz(2.0272994) q[1];
x q[2];
rz(-0.26384683) q[3];
sx q[3];
rz(-0.90875328) q[3];
sx q[3];
rz(2.2784233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
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
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(2.2946987) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(-0.27639595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24647507) q[0];
sx q[0];
rz(-2.7013489) q[0];
sx q[0];
rz(-2.759139) q[0];
rz(1.3771636) q[2];
sx q[2];
rz(-1.3467222) q[2];
sx q[2];
rz(-2.5335238) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8918708) q[1];
sx q[1];
rz(-2.8686214) q[1];
sx q[1];
rz(0.94675468) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3474465) q[3];
sx q[3];
rz(-2.5037519) q[3];
sx q[3];
rz(0.70536246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(0.12750553) q[2];
rz(-3.1048807) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-2.7888443) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(1.030285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.661783) q[0];
sx q[0];
rz(-2.3942238) q[0];
sx q[0];
rz(0.28967793) q[0];
x q[1];
rz(0.68600168) q[2];
sx q[2];
rz(-1.3091607) q[2];
sx q[2];
rz(-2.5930282) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87957571) q[1];
sx q[1];
rz(-0.55393065) q[1];
sx q[1];
rz(-2.5074156) q[1];
rz(2.528245) q[3];
sx q[3];
rz(-1.8457335) q[3];
sx q[3];
rz(0.018523286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(2.4009005) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(3.0497131) q[3];
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
x q[1];
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
rz(2.4521811) q[3];
sx q[3];
rz(-1.1082311) q[3];
sx q[3];
rz(-2.5267596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

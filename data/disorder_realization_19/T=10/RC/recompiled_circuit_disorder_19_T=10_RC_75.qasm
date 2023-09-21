OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(4.5933525) q[0];
sx q[0];
rz(10.070355) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1323224) q[0];
sx q[0];
rz(-1.2148804) q[0];
sx q[0];
rz(0.15412553) q[0];
x q[1];
rz(-0.069734863) q[2];
sx q[2];
rz(-0.94143922) q[2];
sx q[2];
rz(-0.33915181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59939304) q[1];
sx q[1];
rz(-1.7855872) q[1];
sx q[1];
rz(-2.3989124) q[1];
rz(-pi) q[2];
rz(2.6333991) q[3];
sx q[3];
rz(-0.71532202) q[3];
sx q[3];
rz(0.040166044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92007414) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(-2.3085964) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-0.54498589) q[0];
rz(0.90822059) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-2.3166336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075515) q[0];
sx q[0];
rz(-1.3836765) q[0];
sx q[0];
rz(-1.9812816) q[0];
x q[1];
rz(0.40290515) q[2];
sx q[2];
rz(-1.5803792) q[2];
sx q[2];
rz(-2.8275037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7488885) q[1];
sx q[1];
rz(-0.57447937) q[1];
sx q[1];
rz(-0.21484612) q[1];
x q[2];
rz(-0.19552688) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(0.60002458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(0.17641243) q[2];
rz(-0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(-3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-0.46974716) q[0];
rz(-1.5247955) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(3.1030531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8124213) q[0];
sx q[0];
rz(-0.0077795452) q[0];
sx q[0];
rz(-1.6532941) q[0];
rz(-0.58263393) q[2];
sx q[2];
rz(-0.27432549) q[2];
sx q[2];
rz(0.81253101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5650428) q[1];
sx q[1];
rz(-1.0817263) q[1];
sx q[1];
rz(0.069688571) q[1];
x q[2];
rz(0.97088082) q[3];
sx q[3];
rz(-2.0754793) q[3];
sx q[3];
rz(-1.2952627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.58275756) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(0.91252404) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(0.34969774) q[0];
rz(1.8967459) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-2.9464338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0450268) q[0];
sx q[0];
rz(-1.4782527) q[0];
sx q[0];
rz(1.6916313) q[0];
x q[1];
rz(-0.48093502) q[2];
sx q[2];
rz(-1.0785042) q[2];
sx q[2];
rz(-0.42602691) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8735503) q[1];
sx q[1];
rz(-0.51635427) q[1];
sx q[1];
rz(2.0622846) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5545198) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-2.7635014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17079167) q[0];
sx q[0];
rz(-1.0266725) q[0];
sx q[0];
rz(-2.8118954) q[0];
x q[1];
rz(-0.80412229) q[2];
sx q[2];
rz(-0.99493775) q[2];
sx q[2];
rz(2.5428307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.34199076) q[1];
sx q[1];
rz(-1.9365891) q[1];
sx q[1];
rz(-2.0288543) q[1];
rz(2.7630745) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(-2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(0.88796973) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75509214) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(-2.7639672) q[0];
rz(2.8185484) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(-0.91517085) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5765502) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(-2.8898426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2890131) q[2];
sx q[2];
rz(-0.70194178) q[2];
sx q[2];
rz(-0.41350565) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8296216) q[1];
sx q[1];
rz(-1.8835856) q[1];
sx q[1];
rz(-1.642986) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6842501) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(-1.3464348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(-0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(-1.908318) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.4253915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6293684) q[0];
sx q[0];
rz(-0.47591305) q[0];
sx q[0];
rz(-2.7251564) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5064429) q[2];
sx q[2];
rz(-1.1255985) q[2];
sx q[2];
rz(1.839947) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2196301) q[1];
sx q[1];
rz(-1.1929409) q[1];
sx q[1];
rz(-0.65803836) q[1];
rz(-pi) q[2];
rz(1.5464877) q[3];
sx q[3];
rz(-2.6377502) q[3];
sx q[3];
rz(2.3712096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(-1.210775) q[2];
rz(0.0018421729) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84412557) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(-0.67529768) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(-2.3892367) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90826666) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(2.6106735) q[0];
rz(-pi) q[1];
rz(-0.92213995) q[2];
sx q[2];
rz(-1.6862009) q[2];
sx q[2];
rz(0.61308544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4560495) q[1];
sx q[1];
rz(-1.6838264) q[1];
sx q[1];
rz(-3.0526524) q[1];
x q[2];
rz(-1.7767056) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(-2.0549783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-3.1414462) q[2];
rz(2.032062) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(2.1355656) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.3148274) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00025230322) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(-1.89639) q[0];
rz(-pi) q[1];
rz(1.7039653) q[2];
sx q[2];
rz(-1.9317992) q[2];
sx q[2];
rz(1.5683057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4970376) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(2.6130996) q[1];
x q[2];
rz(0.0074578961) q[3];
sx q[3];
rz(-2.3141626) q[3];
sx q[3];
rz(-2.5936562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(-0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7219287) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(-2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(2.749696) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058404) q[0];
sx q[0];
rz(-1.579111) q[0];
sx q[0];
rz(1.0132829) q[0];
rz(-pi) q[1];
rz(-1.395412) q[2];
sx q[2];
rz(-2.4782964) q[2];
sx q[2];
rz(1.7593918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0185768) q[1];
sx q[1];
rz(-1.7565109) q[1];
sx q[1];
rz(-2.3260444) q[1];
x q[2];
rz(1.8896905) q[3];
sx q[3];
rz(-2.1074169) q[3];
sx q[3];
rz(-1.5434138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.60951704) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.8023087) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(-0.17386439) q[2];
sx q[2];
rz(-2.3542913) q[2];
sx q[2];
rz(-2.3402294) q[2];
rz(0.34655456) q[3];
sx q[3];
rz(-1.0481491) q[3];
sx q[3];
rz(-0.92878503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

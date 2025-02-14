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
rz(-2.8831557) q[0];
sx q[0];
rz(-2.8539477) q[0];
sx q[0];
rz(-2.0339461) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(0.62825656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4391651) q[0];
sx q[0];
rz(-1.747923) q[0];
sx q[0];
rz(-3.0217373) q[0];
rz(-1.9626748) q[2];
sx q[2];
rz(-2.2547743) q[2];
sx q[2];
rz(0.98525233) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3648277) q[1];
sx q[1];
rz(-1.3421632) q[1];
sx q[1];
rz(-1.6079812) q[1];
rz(-pi) q[2];
rz(-1.0367582) q[3];
sx q[3];
rz(-1.9566571) q[3];
sx q[3];
rz(-0.88632583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2884752) q[2];
sx q[2];
rz(-1.4318117) q[2];
sx q[2];
rz(0.34118578) q[2];
rz(0.8902542) q[3];
sx q[3];
rz(-1.8743926) q[3];
sx q[3];
rz(2.6064579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.5219236) q[0];
sx q[0];
rz(-2.4578019) q[0];
sx q[0];
rz(0.71440119) q[0];
rz(1.5996541) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(0.094706789) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6833362) q[0];
sx q[0];
rz(-2.8074675) q[0];
sx q[0];
rz(-1.279356) q[0];
rz(0.33659597) q[2];
sx q[2];
rz(-2.1456652) q[2];
sx q[2];
rz(-0.47913715) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3871555) q[1];
sx q[1];
rz(-1.913694) q[1];
sx q[1];
rz(1.8217525) q[1];
x q[2];
rz(2.8625032) q[3];
sx q[3];
rz(-0.52625873) q[3];
sx q[3];
rz(2.581439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16227214) q[2];
sx q[2];
rz(-1.2459735) q[2];
sx q[2];
rz(-2.6675513) q[2];
rz(-0.76215172) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(-0.4711841) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8629465) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(-2.2886724) q[0];
rz(-0.22311738) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(0.51582897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11027656) q[0];
sx q[0];
rz(-3.0678684) q[0];
sx q[0];
rz(0.81839968) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.91024) q[2];
sx q[2];
rz(-2.4759549) q[2];
sx q[2];
rz(-0.15979824) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78385508) q[1];
sx q[1];
rz(-2.2180495) q[1];
sx q[1];
rz(-0.81070645) q[1];
rz(-0.34308691) q[3];
sx q[3];
rz(-2.2394387) q[3];
sx q[3];
rz(2.9647788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3333266) q[2];
sx q[2];
rz(-2.1648679) q[2];
sx q[2];
rz(0.094956368) q[2];
rz(0.44148463) q[3];
sx q[3];
rz(-1.7850826) q[3];
sx q[3];
rz(-1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3917291) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(1.0561426) q[0];
rz(-0.9437584) q[1];
sx q[1];
rz(-2.8902003) q[1];
sx q[1];
rz(1.6897078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7565931) q[0];
sx q[0];
rz(-0.90314048) q[0];
sx q[0];
rz(0.7636015) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46197666) q[2];
sx q[2];
rz(-1.9236132) q[2];
sx q[2];
rz(1.8227673) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6744212) q[1];
sx q[1];
rz(-1.3302263) q[1];
sx q[1];
rz(1.6075385) q[1];
rz(-pi) q[2];
rz(1.5783674) q[3];
sx q[3];
rz(-1.4601344) q[3];
sx q[3];
rz(-0.31294926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(-0.70017868) q[2];
rz(2.2705196) q[3];
sx q[3];
rz(-1.0659404) q[3];
sx q[3];
rz(1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(-2.6569271) q[0];
rz(-2.2635745) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(2.7427618) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8190552) q[0];
sx q[0];
rz(-2.5761671) q[0];
sx q[0];
rz(-0.25382332) q[0];
rz(-0.70339922) q[2];
sx q[2];
rz(-1.3079155) q[2];
sx q[2];
rz(2.2180598) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5113885) q[1];
sx q[1];
rz(-0.92899073) q[1];
sx q[1];
rz(-2.6766269) q[1];
rz(-pi) q[2];
x q[2];
rz(0.871931) q[3];
sx q[3];
rz(-1.1355577) q[3];
sx q[3];
rz(2.1520681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6681119) q[2];
sx q[2];
rz(-1.6010189) q[2];
sx q[2];
rz(2.2637746) q[2];
rz(-0.38464883) q[3];
sx q[3];
rz(-0.84482241) q[3];
sx q[3];
rz(-1.1788684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8574852) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(-0.2704764) q[0];
rz(-2.8349304) q[1];
sx q[1];
rz(-2.6276734) q[1];
sx q[1];
rz(-2.564863) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3424042) q[0];
sx q[0];
rz(-1.6584089) q[0];
sx q[0];
rz(-1.7566998) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88097303) q[2];
sx q[2];
rz(-1.8554129) q[2];
sx q[2];
rz(1.6135482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5698147) q[1];
sx q[1];
rz(-1.0228436) q[1];
sx q[1];
rz(-0.46787365) q[1];
x q[2];
rz(3.1012975) q[3];
sx q[3];
rz(-2.8065805) q[3];
sx q[3];
rz(-0.96714902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7755255) q[2];
sx q[2];
rz(-2.0719353) q[2];
sx q[2];
rz(2.1712187) q[2];
rz(-0.35106418) q[3];
sx q[3];
rz(-2.909436) q[3];
sx q[3];
rz(2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.905726) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(2.6406777) q[0];
rz(0.72055912) q[1];
sx q[1];
rz(-0.84545285) q[1];
sx q[1];
rz(-0.93773425) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9792694) q[0];
sx q[0];
rz(-2.9340194) q[0];
sx q[0];
rz(-1.7718149) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9483614) q[2];
sx q[2];
rz(-2.015471) q[2];
sx q[2];
rz(1.3832165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5825469) q[1];
sx q[1];
rz(-1.2101047) q[1];
sx q[1];
rz(2.2466895) q[1];
rz(-pi) q[2];
rz(2.4547415) q[3];
sx q[3];
rz(-1.5631286) q[3];
sx q[3];
rz(0.075014278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0096036) q[2];
sx q[2];
rz(-0.14543532) q[2];
sx q[2];
rz(-0.84849882) q[2];
rz(-2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(-1.2257303) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7312412) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(-0.89286667) q[0];
rz(-3.0800173) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(-2.9796013) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238857) q[0];
sx q[0];
rz(-3.0382503) q[0];
sx q[0];
rz(1.7211821) q[0];
rz(1.4561557) q[2];
sx q[2];
rz(-0.87213665) q[2];
sx q[2];
rz(-1.6841152) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87016856) q[1];
sx q[1];
rz(-1.7080016) q[1];
sx q[1];
rz(0.31183621) q[1];
rz(3.0477338) q[3];
sx q[3];
rz(-1.4673309) q[3];
sx q[3];
rz(2.2458592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.091247678) q[2];
sx q[2];
rz(-2.7140736) q[2];
sx q[2];
rz(0.72489911) q[2];
rz(-0.27724087) q[3];
sx q[3];
rz(-2.1767949) q[3];
sx q[3];
rz(-0.086517081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8898833) q[0];
sx q[0];
rz(-2.4885663) q[0];
sx q[0];
rz(0.025803331) q[0];
rz(-1.7426527) q[1];
sx q[1];
rz(-1.6394697) q[1];
sx q[1];
rz(0.10765156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7079332) q[0];
sx q[0];
rz(-1.7042394) q[0];
sx q[0];
rz(-0.99876499) q[0];
rz(-pi) q[1];
rz(0.4716267) q[2];
sx q[2];
rz(-1.1749069) q[2];
sx q[2];
rz(-0.63360032) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.746826) q[1];
sx q[1];
rz(-1.0977543) q[1];
sx q[1];
rz(1.5302782) q[1];
rz(1.7259476) q[3];
sx q[3];
rz(-3.0265813) q[3];
sx q[3];
rz(-1.0105159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4699576) q[2];
sx q[2];
rz(-1.2592955) q[2];
sx q[2];
rz(-1.0939481) q[2];
rz(-0.57540244) q[3];
sx q[3];
rz(-2.3016774) q[3];
sx q[3];
rz(2.0247139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6082918) q[0];
sx q[0];
rz(-0.58654439) q[0];
sx q[0];
rz(-2.4285512) q[0];
rz(2.9167922) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(-1.0932659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521234) q[0];
sx q[0];
rz(-2.2575245) q[0];
sx q[0];
rz(-2.5144469) q[0];
rz(0.20168882) q[2];
sx q[2];
rz(-1.8400157) q[2];
sx q[2];
rz(1.9112916) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.374129) q[1];
sx q[1];
rz(-2.0258365) q[1];
sx q[1];
rz(1.6947075) q[1];
rz(-pi) q[2];
rz(1.8328463) q[3];
sx q[3];
rz(-2.7152677) q[3];
sx q[3];
rz(2.654512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0010058086) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(-0.33983964) q[2];
rz(3.0991683) q[3];
sx q[3];
rz(-1.2845311) q[3];
sx q[3];
rz(-0.62780082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18702678) q[0];
sx q[0];
rz(-1.5196336) q[0];
sx q[0];
rz(-1.2484311) q[0];
rz(2.3595702) q[1];
sx q[1];
rz(-0.93688688) q[1];
sx q[1];
rz(-0.72601906) q[1];
rz(2.4154631) q[2];
sx q[2];
rz(-0.48116044) q[2];
sx q[2];
rz(3.0143123) q[2];
rz(-1.9966765) q[3];
sx q[3];
rz(-1.799634) q[3];
sx q[3];
rz(-1.4992731) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

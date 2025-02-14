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
rz(-1.9189605) q[0];
sx q[0];
rz(-0.55813342) q[0];
sx q[0];
rz(-2.4827935) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(-0.69488156) q[1];
sx q[1];
rz(-3.053009) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026042875) q[0];
sx q[0];
rz(-2.6691873) q[0];
sx q[0];
rz(1.1360854) q[0];
rz(-pi) q[1];
rz(-2.1363356) q[2];
sx q[2];
rz(-1.9449678) q[2];
sx q[2];
rz(-2.1924874) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0483612) q[1];
sx q[1];
rz(-2.1856896) q[1];
sx q[1];
rz(-0.1860861) q[1];
rz(-2.5837059) q[3];
sx q[3];
rz(-0.87432623) q[3];
sx q[3];
rz(0.59244746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6170071) q[2];
sx q[2];
rz(-1.1716537) q[2];
sx q[2];
rz(2.6535772) q[2];
rz(-0.46767849) q[3];
sx q[3];
rz(-2.6165163) q[3];
sx q[3];
rz(-2.974143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2117598) q[0];
sx q[0];
rz(-2.3303895) q[0];
sx q[0];
rz(1.8875341) q[0];
rz(0.60011855) q[1];
sx q[1];
rz(-1.3328726) q[1];
sx q[1];
rz(-0.41807237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9771043) q[0];
sx q[0];
rz(-2.573624) q[0];
sx q[0];
rz(-0.33645679) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.977191) q[2];
sx q[2];
rz(-0.65509701) q[2];
sx q[2];
rz(0.10318211) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.969287) q[1];
sx q[1];
rz(-1.2551771) q[1];
sx q[1];
rz(0.068788485) q[1];
rz(-pi) q[2];
rz(0.63184072) q[3];
sx q[3];
rz(-2.2152293) q[3];
sx q[3];
rz(1.0570248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40689251) q[2];
sx q[2];
rz(-2.0105346) q[2];
sx q[2];
rz(-0.1184173) q[2];
rz(-2.7393869) q[3];
sx q[3];
rz(-1.6536568) q[3];
sx q[3];
rz(-1.8124628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44620946) q[0];
sx q[0];
rz(-2.8948247) q[0];
sx q[0];
rz(1.9870019) q[0];
rz(-2.0788976) q[1];
sx q[1];
rz(-2.1176391) q[1];
sx q[1];
rz(1.9564995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4615153) q[0];
sx q[0];
rz(-0.74225589) q[0];
sx q[0];
rz(1.1552325) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1923033) q[2];
sx q[2];
rz(-0.54735294) q[2];
sx q[2];
rz(-2.0277835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.628988) q[1];
sx q[1];
rz(-1.9459241) q[1];
sx q[1];
rz(-2.8400322) q[1];
rz(1.3122769) q[3];
sx q[3];
rz(-1.962477) q[3];
sx q[3];
rz(-1.2092115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9065173) q[2];
sx q[2];
rz(-0.86248988) q[2];
sx q[2];
rz(1.6974576) q[2];
rz(2.2200572) q[3];
sx q[3];
rz(-2.5369365) q[3];
sx q[3];
rz(-3.0864033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4398572) q[0];
sx q[0];
rz(-3.0210962) q[0];
sx q[0];
rz(3.0635656) q[0];
rz(1.1174508) q[1];
sx q[1];
rz(-1.8945339) q[1];
sx q[1];
rz(1.4720346) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32737752) q[0];
sx q[0];
rz(-2.3348885) q[0];
sx q[0];
rz(-0.14847688) q[0];
rz(-0.53773625) q[2];
sx q[2];
rz(-1.0555648) q[2];
sx q[2];
rz(0.73440427) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2591272) q[1];
sx q[1];
rz(-1.5992016) q[1];
sx q[1];
rz(-1.9309429) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4185798) q[3];
sx q[3];
rz(-1.8545056) q[3];
sx q[3];
rz(-0.32645389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83634079) q[2];
sx q[2];
rz(-0.64771104) q[2];
sx q[2];
rz(-2.9913537) q[2];
rz(2.5495106) q[3];
sx q[3];
rz(-1.838107) q[3];
sx q[3];
rz(-1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93064654) q[0];
sx q[0];
rz(-0.35363126) q[0];
sx q[0];
rz(-0.76552248) q[0];
rz(-2.3402479) q[1];
sx q[1];
rz(-1.067602) q[1];
sx q[1];
rz(1.1917005) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0203637) q[0];
sx q[0];
rz(-2.1018711) q[0];
sx q[0];
rz(-1.9241821) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0982047) q[2];
sx q[2];
rz(-2.049438) q[2];
sx q[2];
rz(-0.20028534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4420061) q[1];
sx q[1];
rz(-0.68574673) q[1];
sx q[1];
rz(-1.6722635) q[1];
rz(1.3016316) q[3];
sx q[3];
rz(-1.3084305) q[3];
sx q[3];
rz(3.0848665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2661065) q[2];
sx q[2];
rz(-2.2849639) q[2];
sx q[2];
rz(2.6066656) q[2];
rz(-2.5180425) q[3];
sx q[3];
rz(-1.1112735) q[3];
sx q[3];
rz(0.88304869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542338) q[0];
sx q[0];
rz(-0.41340241) q[0];
sx q[0];
rz(2.2075388) q[0];
rz(1.2454698) q[1];
sx q[1];
rz(-1.586069) q[1];
sx q[1];
rz(0.62560558) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0348822) q[0];
sx q[0];
rz(-2.3090274) q[0];
sx q[0];
rz(-2.9009079) q[0];
x q[1];
rz(-3.0482381) q[2];
sx q[2];
rz(-0.90269222) q[2];
sx q[2];
rz(-0.36076383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4491419) q[1];
sx q[1];
rz(-2.3896273) q[1];
sx q[1];
rz(-0.43159952) q[1];
x q[2];
rz(-2.2295093) q[3];
sx q[3];
rz(-2.8588223) q[3];
sx q[3];
rz(1.2706437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0921359) q[2];
sx q[2];
rz(-0.9149887) q[2];
sx q[2];
rz(-3.0165239) q[2];
rz(-0.72758979) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(-2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8759988) q[0];
sx q[0];
rz(-0.47530526) q[0];
sx q[0];
rz(-1.2566316) q[0];
rz(-1.1988634) q[1];
sx q[1];
rz(-1.7190944) q[1];
sx q[1];
rz(-1.8667603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0555818) q[0];
sx q[0];
rz(-2.0025829) q[0];
sx q[0];
rz(2.5837269) q[0];
rz(-2.1076249) q[2];
sx q[2];
rz(-1.2481999) q[2];
sx q[2];
rz(0.19746298) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15423552) q[1];
sx q[1];
rz(-1.2200095) q[1];
sx q[1];
rz(-2.039775) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2009462) q[3];
sx q[3];
rz(-0.33888926) q[3];
sx q[3];
rz(1.3047993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4510497) q[2];
sx q[2];
rz(-1.0131016) q[2];
sx q[2];
rz(-0.96735442) q[2];
rz(1.5585772) q[3];
sx q[3];
rz(-1.6396921) q[3];
sx q[3];
rz(-0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79307443) q[0];
sx q[0];
rz(-2.5418042) q[0];
sx q[0];
rz(0.10979688) q[0];
rz(-1.9684567) q[1];
sx q[1];
rz(-0.99183142) q[1];
sx q[1];
rz(2.0593624) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5062456) q[0];
sx q[0];
rz(-1.325334) q[0];
sx q[0];
rz(2.3574791) q[0];
rz(1.3668109) q[2];
sx q[2];
rz(-2.6480541) q[2];
sx q[2];
rz(0.94369027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6746862) q[1];
sx q[1];
rz(-1.162858) q[1];
sx q[1];
rz(2.3323184) q[1];
rz(-pi) q[2];
rz(-2.5030062) q[3];
sx q[3];
rz(-1.878384) q[3];
sx q[3];
rz(-1.547055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3829696) q[2];
sx q[2];
rz(-1.0926282) q[2];
sx q[2];
rz(0.36337241) q[2];
rz(0.97909561) q[3];
sx q[3];
rz(-2.4978814) q[3];
sx q[3];
rz(2.6399829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8363504) q[0];
sx q[0];
rz(-0.79638052) q[0];
sx q[0];
rz(-0.35261944) q[0];
rz(-1.3615707) q[1];
sx q[1];
rz(-1.1905328) q[1];
sx q[1];
rz(-0.1415267) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3985719) q[0];
sx q[0];
rz(-1.5367265) q[0];
sx q[0];
rz(0.16966804) q[0];
rz(-pi) q[1];
rz(0.63371079) q[2];
sx q[2];
rz(-2.4437208) q[2];
sx q[2];
rz(0.77264589) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5302116) q[1];
sx q[1];
rz(-2.5712643) q[1];
sx q[1];
rz(2.4228816) q[1];
rz(-pi) q[2];
rz(2.448672) q[3];
sx q[3];
rz(-2.3669405) q[3];
sx q[3];
rz(1.2415394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3129468) q[2];
sx q[2];
rz(-1.6673648) q[2];
sx q[2];
rz(1.0860156) q[2];
rz(-2.9588251) q[3];
sx q[3];
rz(-2.1153617) q[3];
sx q[3];
rz(-0.97203794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52520853) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(0.64424789) q[0];
rz(3.0746025) q[1];
sx q[1];
rz(-1.3863775) q[1];
sx q[1];
rz(0.11360528) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4098334) q[0];
sx q[0];
rz(-1.5749965) q[0];
sx q[0];
rz(1.3467711) q[0];
rz(0.082661672) q[2];
sx q[2];
rz(-1.4689406) q[2];
sx q[2];
rz(1.5422623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0225102) q[1];
sx q[1];
rz(-1.7412698) q[1];
sx q[1];
rz(-1.4524231) q[1];
x q[2];
rz(-0.04256256) q[3];
sx q[3];
rz(-1.213338) q[3];
sx q[3];
rz(-2.8379311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99864787) q[2];
sx q[2];
rz(-1.2544268) q[2];
sx q[2];
rz(-0.76510731) q[2];
rz(-2.9633925) q[3];
sx q[3];
rz(-1.5811812) q[3];
sx q[3];
rz(2.3044738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2160303) q[0];
sx q[0];
rz(-2.3860274) q[0];
sx q[0];
rz(-0.26269333) q[0];
rz(2.8021011) q[1];
sx q[1];
rz(-1.2295634) q[1];
sx q[1];
rz(-2.0801574) q[1];
rz(-0.65503623) q[2];
sx q[2];
rz(-1.8476386) q[2];
sx q[2];
rz(2.4551433) q[2];
rz(2.3332023) q[3];
sx q[3];
rz(-0.42338531) q[3];
sx q[3];
rz(0.12082621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

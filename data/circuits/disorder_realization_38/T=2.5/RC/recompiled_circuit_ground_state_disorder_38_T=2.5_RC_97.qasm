OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.030169686) q[0];
sx q[0];
rz(4.5017894) q[0];
sx q[0];
rz(5.9202249) q[0];
rz(1.9650004) q[1];
sx q[1];
rz(1.7791553) q[1];
sx q[1];
rz(11.122009) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43508726) q[0];
sx q[0];
rz(-1.4772863) q[0];
sx q[0];
rz(-2.963752) q[0];
x q[1];
rz(0.40048234) q[2];
sx q[2];
rz(-2.8304511) q[2];
sx q[2];
rz(-0.36382407) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0054659) q[1];
sx q[1];
rz(-0.44584063) q[1];
sx q[1];
rz(-1.3492829) q[1];
rz(-pi) q[2];
rz(3.1009272) q[3];
sx q[3];
rz(-0.96974428) q[3];
sx q[3];
rz(-1.9448626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1732257) q[2];
sx q[2];
rz(-2.3543365) q[2];
sx q[2];
rz(-3.0464029) q[2];
rz(1.7732874) q[3];
sx q[3];
rz(-0.94822001) q[3];
sx q[3];
rz(-0.29909721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66265166) q[0];
sx q[0];
rz(-2.3421685) q[0];
sx q[0];
rz(-0.69315243) q[0];
rz(-0.33723801) q[1];
sx q[1];
rz(-1.8353029) q[1];
sx q[1];
rz(2.3487924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900729) q[0];
sx q[0];
rz(-1.8197104) q[0];
sx q[0];
rz(-0.78686692) q[0];
x q[1];
rz(0.22979615) q[2];
sx q[2];
rz(-1.0655128) q[2];
sx q[2];
rz(-1.2838703) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0568367) q[1];
sx q[1];
rz(-0.82377071) q[1];
sx q[1];
rz(-0.39577248) q[1];
rz(1.0569237) q[3];
sx q[3];
rz(-2.3213291) q[3];
sx q[3];
rz(-2.6721374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1482131) q[2];
sx q[2];
rz(-0.17017636) q[2];
sx q[2];
rz(-1.4519838) q[2];
rz(0.65886894) q[3];
sx q[3];
rz(-1.4631319) q[3];
sx q[3];
rz(0.39670593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0186998) q[0];
sx q[0];
rz(-2.1227699) q[0];
sx q[0];
rz(-3.0320211) q[0];
rz(-2.295629) q[1];
sx q[1];
rz(-0.94596091) q[1];
sx q[1];
rz(2.3931961) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40846857) q[0];
sx q[0];
rz(-1.4811902) q[0];
sx q[0];
rz(1.0138847) q[0];
rz(-pi) q[1];
rz(0.72785925) q[2];
sx q[2];
rz(-0.72304356) q[2];
sx q[2];
rz(-2.7287908) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8019164) q[1];
sx q[1];
rz(-1.5775562) q[1];
sx q[1];
rz(-1.6338145) q[1];
x q[2];
rz(-2.9633065) q[3];
sx q[3];
rz(-2.1349838) q[3];
sx q[3];
rz(-2.215204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1154068) q[2];
sx q[2];
rz(-1.4226961) q[2];
sx q[2];
rz(2.4288948) q[2];
rz(1.6648434) q[3];
sx q[3];
rz(-1.289117) q[3];
sx q[3];
rz(-3.0570928) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1348212) q[0];
sx q[0];
rz(-2.0016142) q[0];
sx q[0];
rz(3.0546597) q[0];
rz(-2.8839819) q[1];
sx q[1];
rz(-2.0517709) q[1];
sx q[1];
rz(1.849252) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5851) q[0];
sx q[0];
rz(-0.79946639) q[0];
sx q[0];
rz(1.9071346) q[0];
rz(0.33165641) q[2];
sx q[2];
rz(-1.4185959) q[2];
sx q[2];
rz(-2.4401522) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.72491652) q[1];
sx q[1];
rz(-1.2236164) q[1];
sx q[1];
rz(-0.93840878) q[1];
rz(-pi) q[2];
rz(-0.35859006) q[3];
sx q[3];
rz(-1.4309644) q[3];
sx q[3];
rz(-0.38966013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8452235) q[2];
sx q[2];
rz(-1.1901647) q[2];
sx q[2];
rz(-2.2947218) q[2];
rz(1.752468) q[3];
sx q[3];
rz(-1.4877157) q[3];
sx q[3];
rz(-0.59051591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95319372) q[0];
sx q[0];
rz(-2.0162835) q[0];
sx q[0];
rz(-0.97287792) q[0];
rz(-2.1994622) q[1];
sx q[1];
rz(-2.3140488) q[1];
sx q[1];
rz(3.1408302) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98648345) q[0];
sx q[0];
rz(-2.0707978) q[0];
sx q[0];
rz(2.1613104) q[0];
x q[1];
rz(0.43652759) q[2];
sx q[2];
rz(-0.68219705) q[2];
sx q[2];
rz(1.6641219) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.768133) q[1];
sx q[1];
rz(-2.4985857) q[1];
sx q[1];
rz(-1.0932342) q[1];
x q[2];
rz(2.5464779) q[3];
sx q[3];
rz(-1.5059808) q[3];
sx q[3];
rz(2.7832372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7611277) q[2];
sx q[2];
rz(-1.4030115) q[2];
sx q[2];
rz(0.38499704) q[2];
rz(-0.72832406) q[3];
sx q[3];
rz(-0.60898048) q[3];
sx q[3];
rz(0.42496067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7321135) q[0];
sx q[0];
rz(-0.80687579) q[0];
sx q[0];
rz(-2.683486) q[0];
rz(1.4234281) q[1];
sx q[1];
rz(-0.79045311) q[1];
sx q[1];
rz(-0.093553392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6038743) q[0];
sx q[0];
rz(-0.97386003) q[0];
sx q[0];
rz(-0.56054795) q[0];
x q[1];
rz(-1.3504998) q[2];
sx q[2];
rz(-1.7335906) q[2];
sx q[2];
rz(-0.3045813) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97620353) q[1];
sx q[1];
rz(-1.3553344) q[1];
sx q[1];
rz(0.54241753) q[1];
rz(-pi) q[2];
rz(3.1074353) q[3];
sx q[3];
rz(-0.56006008) q[3];
sx q[3];
rz(0.36049662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2576922) q[2];
sx q[2];
rz(-0.42282405) q[2];
sx q[2];
rz(-0.60919961) q[2];
rz(0.60112634) q[3];
sx q[3];
rz(-1.6060035) q[3];
sx q[3];
rz(2.6570184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5179829) q[0];
sx q[0];
rz(-1.4731151) q[0];
sx q[0];
rz(-1.8286888) q[0];
rz(3.0598705) q[1];
sx q[1];
rz(-1.3879958) q[1];
sx q[1];
rz(1.0645701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1601105) q[0];
sx q[0];
rz(-0.89799268) q[0];
sx q[0];
rz(2.299911) q[0];
rz(-0.95439744) q[2];
sx q[2];
rz(-1.5570882) q[2];
sx q[2];
rz(-2.990266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2456817) q[1];
sx q[1];
rz(-0.57990042) q[1];
sx q[1];
rz(-1.4407451) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96780583) q[3];
sx q[3];
rz(-1.9120312) q[3];
sx q[3];
rz(-1.5016709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99819034) q[2];
sx q[2];
rz(-0.4021796) q[2];
sx q[2];
rz(-1.1581988) q[2];
rz(0.68425933) q[3];
sx q[3];
rz(-1.8161769) q[3];
sx q[3];
rz(-1.6392596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.174468) q[0];
sx q[0];
rz(-3.0376349) q[0];
sx q[0];
rz(-2.6167468) q[0];
rz(1.1809433) q[1];
sx q[1];
rz(-2.3044105) q[1];
sx q[1];
rz(-0.7276853) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0200266) q[0];
sx q[0];
rz(-1.7870197) q[0];
sx q[0];
rz(2.2452767) q[0];
x q[1];
rz(-2.4196376) q[2];
sx q[2];
rz(-1.180207) q[2];
sx q[2];
rz(0.73467322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3062417) q[1];
sx q[1];
rz(-1.1351651) q[1];
sx q[1];
rz(1.2795875) q[1];
rz(-pi) q[2];
rz(-0.10730174) q[3];
sx q[3];
rz(-2.2224217) q[3];
sx q[3];
rz(-0.95355031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9119447) q[2];
sx q[2];
rz(-1.0934528) q[2];
sx q[2];
rz(1.3463119) q[2];
rz(0.66156578) q[3];
sx q[3];
rz(-1.9965568) q[3];
sx q[3];
rz(-3.055011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0410556) q[0];
sx q[0];
rz(-2.4424398) q[0];
sx q[0];
rz(-2.2030785) q[0];
rz(1.8427294) q[1];
sx q[1];
rz(-2.2524736) q[1];
sx q[1];
rz(-0.13359698) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7093606) q[0];
sx q[0];
rz(-1.2424766) q[0];
sx q[0];
rz(-1.6698056) q[0];
rz(-pi) q[1];
rz(1.8739545) q[2];
sx q[2];
rz(-1.016262) q[2];
sx q[2];
rz(-1.6078439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2707227) q[1];
sx q[1];
rz(-0.36859757) q[1];
sx q[1];
rz(2.9245118) q[1];
rz(-pi) q[2];
rz(-1.9383921) q[3];
sx q[3];
rz(-0.9617241) q[3];
sx q[3];
rz(-1.6742791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9962697) q[2];
sx q[2];
rz(-1.0066373) q[2];
sx q[2];
rz(2.2553196) q[2];
rz(-2.8730872) q[3];
sx q[3];
rz(-2.0566745) q[3];
sx q[3];
rz(1.5119542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4292384) q[0];
sx q[0];
rz(-2.7463284) q[0];
sx q[0];
rz(-0.23671737) q[0];
rz(0.16353823) q[1];
sx q[1];
rz(-1.5312559) q[1];
sx q[1];
rz(-2.9530361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8343413) q[0];
sx q[0];
rz(-2.3545958) q[0];
sx q[0];
rz(-2.045162) q[0];
rz(0.47958572) q[2];
sx q[2];
rz(-1.34362) q[2];
sx q[2];
rz(-0.6152252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.982354) q[1];
sx q[1];
rz(-2.6272863) q[1];
sx q[1];
rz(1.6883172) q[1];
rz(-pi) q[2];
x q[2];
rz(1.439752) q[3];
sx q[3];
rz(-2.0072972) q[3];
sx q[3];
rz(0.97697645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4208372) q[2];
sx q[2];
rz(-1.5034224) q[2];
sx q[2];
rz(-0.76589626) q[2];
rz(3.1229373) q[3];
sx q[3];
rz(-1.3466287) q[3];
sx q[3];
rz(-0.53001058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6388549) q[0];
sx q[0];
rz(-1.9866332) q[0];
sx q[0];
rz(-1.6992983) q[0];
rz(2.8522708) q[1];
sx q[1];
rz(-1.0565636) q[1];
sx q[1];
rz(0.4531959) q[1];
rz(-0.53884698) q[2];
sx q[2];
rz(-1.4543903) q[2];
sx q[2];
rz(-1.7001154) q[2];
rz(-1.4342593) q[3];
sx q[3];
rz(-0.47341163) q[3];
sx q[3];
rz(0.51281131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

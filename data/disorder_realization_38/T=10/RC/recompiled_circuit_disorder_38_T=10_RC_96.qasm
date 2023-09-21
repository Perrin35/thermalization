OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(-2.1698706) q[1];
sx q[1];
rz(0.5501774) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8367856) q[0];
sx q[0];
rz(-1.3512304) q[0];
sx q[0];
rz(-3.1138793) q[0];
rz(2.7031541) q[2];
sx q[2];
rz(-1.0796094) q[2];
sx q[2];
rz(-3.0868798) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8093811) q[1];
sx q[1];
rz(-2.2978133) q[1];
sx q[1];
rz(1.0512645) q[1];
x q[2];
rz(-0.95991858) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(-1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41574079) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(2.501781) q[2];
rz(-0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(0.71331435) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.5126022) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0523895) q[0];
sx q[0];
rz(-0.69141483) q[0];
sx q[0];
rz(2.0877439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.011629148) q[2];
sx q[2];
rz(-2.8577869) q[2];
sx q[2];
rz(1.1139882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7035547) q[1];
sx q[1];
rz(-2.2433271) q[1];
sx q[1];
rz(-2.9707675) q[1];
rz(-pi) q[2];
rz(0.87725957) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(3.0960494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9397395) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(0.29552466) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(2.3957516) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94507664) q[0];
sx q[0];
rz(-1.3498422) q[0];
sx q[0];
rz(-3.0483079) q[0];
rz(-1.4065811) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(-0.28902136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5513788) q[1];
sx q[1];
rz(-0.56400245) q[1];
sx q[1];
rz(-2.6022634) q[1];
rz(0.45332076) q[3];
sx q[3];
rz(-1.7954149) q[3];
sx q[3];
rz(-2.2090705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(-0.68874613) q[2];
rz(0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(-3.0396089) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(0.27329683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647588) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(1.8434974) q[0];
rz(0.64812135) q[2];
sx q[2];
rz(-2.0378049) q[2];
sx q[2];
rz(-1.8304706) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3015775) q[1];
sx q[1];
rz(-1.3372278) q[1];
sx q[1];
rz(-2.9905) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9076505) q[3];
sx q[3];
rz(-2.4787239) q[3];
sx q[3];
rz(-0.58594221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79167241) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-2.1544429) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(2.2228918) q[0];
x q[1];
rz(3.1290595) q[2];
sx q[2];
rz(-2.4977411) q[2];
sx q[2];
rz(2.8963793) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3293216) q[1];
sx q[1];
rz(-1.2580039) q[1];
sx q[1];
rz(2.4815464) q[1];
x q[2];
rz(-3.0994814) q[3];
sx q[3];
rz(-1.1782421) q[3];
sx q[3];
rz(-2.7819355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-0.21128543) q[2];
rz(0.42090297) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(-0.063407272) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(1.2794367) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1267705) q[0];
sx q[0];
rz(-1.556067) q[0];
sx q[0];
rz(1.2889839) q[0];
rz(-0.40462599) q[2];
sx q[2];
rz(-2.4936516) q[2];
sx q[2];
rz(0.65587703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64843236) q[1];
sx q[1];
rz(-1.6006159) q[1];
sx q[1];
rz(-2.0926507) q[1];
rz(-pi) q[2];
rz(-1.389723) q[3];
sx q[3];
rz(-1.4317792) q[3];
sx q[3];
rz(2.3982323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-0.77077579) q[2];
rz(-1.6714913) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(-3.1380222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647033) q[0];
sx q[0];
rz(-1.2029552) q[0];
sx q[0];
rz(0.44643114) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3891719) q[2];
sx q[2];
rz(-1.8634897) q[2];
sx q[2];
rz(-0.21586403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7230941) q[1];
sx q[1];
rz(-1.8685891) q[1];
sx q[1];
rz(2.198092) q[1];
x q[2];
rz(-1.2083863) q[3];
sx q[3];
rz(-0.21167314) q[3];
sx q[3];
rz(0.025312245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(2.9428633) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(2.7238817) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51628095) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(-0.40813804) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(0.67869854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0146101) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(1.2922657) q[0];
rz(-pi) q[1];
rz(-0.46531123) q[2];
sx q[2];
rz(-0.88854549) q[2];
sx q[2];
rz(-2.0899783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.855367) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(1.0463868) q[1];
rz(-1.7292761) q[3];
sx q[3];
rz(-1.1500119) q[3];
sx q[3];
rz(2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.17710182) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(2.881799) q[0];
rz(0.70867509) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567112) q[0];
sx q[0];
rz(-1.3538133) q[0];
sx q[0];
rz(2.7816539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2277463) q[2];
sx q[2];
rz(-2.4545049) q[2];
sx q[2];
rz(-2.3234141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37104169) q[1];
sx q[1];
rz(-2.7275804) q[1];
sx q[1];
rz(-0.71055926) q[1];
x q[2];
rz(2.5832289) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(-2.5496897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32593411) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(-2.7549426) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-0.3607761) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26786131) q[0];
sx q[0];
rz(-2.972313) q[0];
sx q[0];
rz(1.4219567) q[0];
rz(1.0656409) q[2];
sx q[2];
rz(-2.897122) q[2];
sx q[2];
rz(0.62015647) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.39721397) q[1];
sx q[1];
rz(-2.386464) q[1];
sx q[1];
rz(-0.20882864) q[1];
x q[2];
rz(-1.7182699) q[3];
sx q[3];
rz(-1.7383988) q[3];
sx q[3];
rz(2.9180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(2.5676981) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5671134) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(1.3654937) q[2];
sx q[2];
rz(-1.6663972) q[2];
sx q[2];
rz(-1.321928) q[2];
rz(2.3902262) q[3];
sx q[3];
rz(-2.1204227) q[3];
sx q[3];
rz(-1.942996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
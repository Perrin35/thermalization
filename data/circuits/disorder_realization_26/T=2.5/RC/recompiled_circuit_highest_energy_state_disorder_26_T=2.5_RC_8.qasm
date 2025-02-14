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
rz(0.8910203) q[0];
sx q[0];
rz(-1.2863337) q[0];
sx q[0];
rz(-2.2556055) q[0];
rz(0.52357829) q[1];
sx q[1];
rz(-2.5819267) q[1];
sx q[1];
rz(-1.3203415) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65151764) q[0];
sx q[0];
rz(-2.1044255) q[0];
sx q[0];
rz(-0.20488157) q[0];
rz(-0.34446005) q[2];
sx q[2];
rz(-2.3340324) q[2];
sx q[2];
rz(-1.5986313) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7976171) q[1];
sx q[1];
rz(-1.5098176) q[1];
sx q[1];
rz(2.134255) q[1];
rz(-pi) q[2];
rz(-2.3137847) q[3];
sx q[3];
rz(-2.1850065) q[3];
sx q[3];
rz(0.18060623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3628799) q[2];
sx q[2];
rz(-1.820463) q[2];
sx q[2];
rz(2.2134181) q[2];
rz(-1.5859531) q[3];
sx q[3];
rz(-1.0471683) q[3];
sx q[3];
rz(-1.9716523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.25598079) q[0];
sx q[0];
rz(-0.33807999) q[0];
sx q[0];
rz(2.0804491) q[0];
rz(-1.2127009) q[1];
sx q[1];
rz(-2.3586528) q[1];
sx q[1];
rz(1.8276259) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43976682) q[0];
sx q[0];
rz(-1.9226388) q[0];
sx q[0];
rz(2.3934796) q[0];
rz(-1.018386) q[2];
sx q[2];
rz(-0.72478849) q[2];
sx q[2];
rz(1.2040491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3281341) q[1];
sx q[1];
rz(-0.5603921) q[1];
sx q[1];
rz(-1.8916393) q[1];
rz(-pi) q[2];
rz(-0.98269083) q[3];
sx q[3];
rz(-2.324571) q[3];
sx q[3];
rz(2.6213561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5709915) q[2];
sx q[2];
rz(-2.1676895) q[2];
sx q[2];
rz(-2.2368597) q[2];
rz(-1.5915126) q[3];
sx q[3];
rz(-1.855987) q[3];
sx q[3];
rz(-1.8486283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631113) q[0];
sx q[0];
rz(-0.08992973) q[0];
sx q[0];
rz(-0.91823804) q[0];
rz(-0.69084424) q[1];
sx q[1];
rz(-1.7450688) q[1];
sx q[1];
rz(2.3450559) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7762866) q[0];
sx q[0];
rz(-1.7788634) q[0];
sx q[0];
rz(-0.093642276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.042823533) q[2];
sx q[2];
rz(-1.6302675) q[2];
sx q[2];
rz(2.8958993) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84948375) q[1];
sx q[1];
rz(-0.70569084) q[1];
sx q[1];
rz(-1.73805) q[1];
rz(-pi) q[2];
rz(0.34557202) q[3];
sx q[3];
rz(-2.1949147) q[3];
sx q[3];
rz(-2.6081599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1306661) q[2];
sx q[2];
rz(-0.9504168) q[2];
sx q[2];
rz(-1.8939023) q[2];
rz(2.583336) q[3];
sx q[3];
rz(-1.4562675) q[3];
sx q[3];
rz(3.1202417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70304865) q[0];
sx q[0];
rz(-0.049228638) q[0];
sx q[0];
rz(-0.7274279) q[0];
rz(-2.9797331) q[1];
sx q[1];
rz(-1.1715803) q[1];
sx q[1];
rz(-1.1579827) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2271661) q[0];
sx q[0];
rz(-3.0241443) q[0];
sx q[0];
rz(-2.1690452) q[0];
rz(-pi) q[1];
rz(-3.1406265) q[2];
sx q[2];
rz(-1.7100934) q[2];
sx q[2];
rz(-1.7666221) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0772845) q[1];
sx q[1];
rz(-1.3798215) q[1];
sx q[1];
rz(3.0160267) q[1];
rz(-pi) q[2];
rz(1.4876266) q[3];
sx q[3];
rz(-1.0462049) q[3];
sx q[3];
rz(0.34357854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27266476) q[2];
sx q[2];
rz(-1.9395892) q[2];
sx q[2];
rz(1.391601) q[2];
rz(2.9113655) q[3];
sx q[3];
rz(-1.1743436) q[3];
sx q[3];
rz(-2.9496884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4153083) q[0];
sx q[0];
rz(-3.0400161) q[0];
sx q[0];
rz(-1.104785) q[0];
rz(0.11100189) q[1];
sx q[1];
rz(-2.0260725) q[1];
sx q[1];
rz(1.312779) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5376624) q[0];
sx q[0];
rz(-2.8712508) q[0];
sx q[0];
rz(-0.84690799) q[0];
rz(-1.6592624) q[2];
sx q[2];
rz(-2.1328204) q[2];
sx q[2];
rz(-2.1254181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.032636383) q[1];
sx q[1];
rz(-1.7708888) q[1];
sx q[1];
rz(1.1632827) q[1];
x q[2];
rz(-0.85185527) q[3];
sx q[3];
rz(-1.8924215) q[3];
sx q[3];
rz(2.3648398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.50292) q[2];
sx q[2];
rz(-2.0827677) q[2];
sx q[2];
rz(-1.0271094) q[2];
rz(0.98661679) q[3];
sx q[3];
rz(-2.8592181) q[3];
sx q[3];
rz(1.4237039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33243772) q[0];
sx q[0];
rz(-1.9092535) q[0];
sx q[0];
rz(-0.80068457) q[0];
rz(0.77146161) q[1];
sx q[1];
rz(-2.0636676) q[1];
sx q[1];
rz(1.7652184) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0906935) q[0];
sx q[0];
rz(-0.92779033) q[0];
sx q[0];
rz(2.0741295) q[0];
rz(3.1065953) q[2];
sx q[2];
rz(-2.0366324) q[2];
sx q[2];
rz(-2.310487) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6668762) q[1];
sx q[1];
rz(-2.2762269) q[1];
sx q[1];
rz(-1.5882467) q[1];
rz(-2.3153051) q[3];
sx q[3];
rz(-1.0184231) q[3];
sx q[3];
rz(-0.9780405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8419522) q[2];
sx q[2];
rz(-1.2066634) q[2];
sx q[2];
rz(3.1237777) q[2];
rz(2.8282015) q[3];
sx q[3];
rz(-2.119795) q[3];
sx q[3];
rz(0.71948403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7844836) q[0];
sx q[0];
rz(-1.5672368) q[0];
sx q[0];
rz(-0.16726476) q[0];
rz(0.74370614) q[1];
sx q[1];
rz(-2.2552762) q[1];
sx q[1];
rz(-0.13626616) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0038827) q[0];
sx q[0];
rz(-1.2943177) q[0];
sx q[0];
rz(2.2599758) q[0];
x q[1];
rz(0.17301128) q[2];
sx q[2];
rz(-1.271651) q[2];
sx q[2];
rz(2.0188221) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4140461) q[1];
sx q[1];
rz(-0.36484499) q[1];
sx q[1];
rz(2.9529497) q[1];
x q[2];
rz(-1.6266009) q[3];
sx q[3];
rz(-1.9161092) q[3];
sx q[3];
rz(-1.2962411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83884376) q[2];
sx q[2];
rz(-2.9911797) q[2];
sx q[2];
rz(2.0602843) q[2];
rz(0.14351621) q[3];
sx q[3];
rz(-1.84294) q[3];
sx q[3];
rz(-1.6941841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74893394) q[0];
sx q[0];
rz(-1.8086139) q[0];
sx q[0];
rz(0.64312154) q[0];
rz(0.92203036) q[1];
sx q[1];
rz(-2.8755201) q[1];
sx q[1];
rz(-1.6304852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.155678) q[0];
sx q[0];
rz(-1.2307457) q[0];
sx q[0];
rz(2.9513533) q[0];
rz(1.8505356) q[2];
sx q[2];
rz(-1.6806753) q[2];
sx q[2];
rz(0.85816511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0865759) q[1];
sx q[1];
rz(-0.19682238) q[1];
sx q[1];
rz(-1.5301401) q[1];
rz(-1.891252) q[3];
sx q[3];
rz(-2.0600187) q[3];
sx q[3];
rz(-1.7699522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.94613218) q[2];
sx q[2];
rz(-1.5798502) q[2];
sx q[2];
rz(-1.3004318) q[2];
rz(-2.2293633) q[3];
sx q[3];
rz(-1.2765086) q[3];
sx q[3];
rz(0.31360489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4881956) q[0];
sx q[0];
rz(-2.6618239) q[0];
sx q[0];
rz(-2.524014) q[0];
rz(3.0600582) q[1];
sx q[1];
rz(-0.50059861) q[1];
sx q[1];
rz(2.4542782) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0162553) q[0];
sx q[0];
rz(-2.9000686) q[0];
sx q[0];
rz(1.3708273) q[0];
rz(-0.26726802) q[2];
sx q[2];
rz(-1.6785926) q[2];
sx q[2];
rz(1.2085309) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0159434) q[1];
sx q[1];
rz(-1.1563035) q[1];
sx q[1];
rz(-1.260956) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0275484) q[3];
sx q[3];
rz(-1.6121658) q[3];
sx q[3];
rz(-2.4376805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0240137) q[2];
sx q[2];
rz(-1.964183) q[2];
sx q[2];
rz(-3.0702316) q[2];
rz(1.9060382) q[3];
sx q[3];
rz(-0.33696285) q[3];
sx q[3];
rz(2.5753042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.5736893) q[0];
sx q[0];
rz(-2.8901143) q[0];
sx q[0];
rz(-3.0036744) q[0];
rz(2.5714286) q[1];
sx q[1];
rz(-1.0823931) q[1];
sx q[1];
rz(-0.015448419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041067657) q[0];
sx q[0];
rz(-2.3599632) q[0];
sx q[0];
rz(-1.1055205) q[0];
rz(-pi) q[1];
x q[1];
rz(2.040898) q[2];
sx q[2];
rz(-0.86720556) q[2];
sx q[2];
rz(-2.4936465) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9271157) q[1];
sx q[1];
rz(-1.0815911) q[1];
sx q[1];
rz(1.7577043) q[1];
rz(-pi) q[2];
rz(-0.055850765) q[3];
sx q[3];
rz(-1.7540364) q[3];
sx q[3];
rz(-0.88411546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40921679) q[2];
sx q[2];
rz(-0.78745431) q[2];
sx q[2];
rz(-2.709205) q[2];
rz(-2.2435097) q[3];
sx q[3];
rz(-0.87528527) q[3];
sx q[3];
rz(1.557365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4846004) q[0];
sx q[0];
rz(-1.1285755) q[0];
sx q[0];
rz(-2.3791671) q[0];
rz(0.55105974) q[1];
sx q[1];
rz(-1.8012128) q[1];
sx q[1];
rz(2.6691379) q[1];
rz(-2.77825) q[2];
sx q[2];
rz(-1.5283593) q[2];
sx q[2];
rz(-1.9178393) q[2];
rz(-2.2021709) q[3];
sx q[3];
rz(-2.4236854) q[3];
sx q[3];
rz(-3.051306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

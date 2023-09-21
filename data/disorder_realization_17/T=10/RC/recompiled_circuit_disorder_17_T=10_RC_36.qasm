OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9085812) q[0];
sx q[0];
rz(-1.9549978) q[0];
sx q[0];
rz(-1.1458122) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(-0.29247984) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16380331) q[0];
sx q[0];
rz(-1.5153432) q[0];
sx q[0];
rz(-1.9652912) q[0];
rz(-pi) q[1];
rz(-2.3621759) q[2];
sx q[2];
rz(-1.6072818) q[2];
sx q[2];
rz(-1.2365637) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8279743) q[1];
sx q[1];
rz(-0.59578124) q[1];
sx q[1];
rz(1.4890563) q[1];
x q[2];
rz(1.5745893) q[3];
sx q[3];
rz(-0.47382254) q[3];
sx q[3];
rz(-2.4592318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.73074377) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(-1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(-2.399562) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(1.3630294) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5308967) q[0];
sx q[0];
rz(-2.1801821) q[0];
sx q[0];
rz(-1.0884398) q[0];
x q[1];
rz(-2.734191) q[2];
sx q[2];
rz(-2.6371837) q[2];
sx q[2];
rz(-0.87770578) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9957673) q[1];
sx q[1];
rz(-1.2312504) q[1];
sx q[1];
rz(-1.5911566) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3046673) q[3];
sx q[3];
rz(-1.7566655) q[3];
sx q[3];
rz(2.4974125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30119511) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.2724686) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(-3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(-0.2628251) q[0];
rz(-2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(-1.9062818) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6567898) q[0];
sx q[0];
rz(-0.018253837) q[0];
sx q[0];
rz(0.59469964) q[0];
rz(-3.1149408) q[2];
sx q[2];
rz(-1.3304552) q[2];
sx q[2];
rz(-1.9453366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23535145) q[1];
sx q[1];
rz(-1.6760577) q[1];
sx q[1];
rz(-0.32354849) q[1];
rz(0.34706195) q[3];
sx q[3];
rz(-2.2357781) q[3];
sx q[3];
rz(-2.5393328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-0.18033218) q[2];
rz(-2.5056433) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(-0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3746049) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(2.4225127) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(1.0346574) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661515) q[0];
sx q[0];
rz(-3.0794641) q[0];
sx q[0];
rz(-0.32390578) q[0];
rz(-pi) q[1];
rz(1.2345418) q[2];
sx q[2];
rz(-0.14306919) q[2];
sx q[2];
rz(2.2640995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0532916) q[1];
sx q[1];
rz(-2.3846855) q[1];
sx q[1];
rz(0.41815586) q[1];
rz(-pi) q[2];
rz(0.47311584) q[3];
sx q[3];
rz(-1.4149727) q[3];
sx q[3];
rz(-1.6254049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2085312) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(1.9832206) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(-2.2896144) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1119969) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(2.239256) q[0];
rz(-0.92102712) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-0.62228084) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8181994) q[0];
sx q[0];
rz(-0.33938956) q[0];
sx q[0];
rz(-2.1470977) q[0];
x q[1];
rz(3.0201966) q[2];
sx q[2];
rz(-2.4897794) q[2];
sx q[2];
rz(-2.7898942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8497148) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(-0.090976322) q[1];
rz(-pi) q[2];
rz(-1.9641818) q[3];
sx q[3];
rz(-2.0188361) q[3];
sx q[3];
rz(-3.0569227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92419147) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(1.6993258) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570046) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(3.0969627) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(3.1071641) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444367) q[0];
sx q[0];
rz(-1.3855055) q[0];
sx q[0];
rz(-0.013884355) q[0];
x q[1];
rz(-1.351864) q[2];
sx q[2];
rz(-2.6282675) q[2];
sx q[2];
rz(-0.22253144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62528246) q[1];
sx q[1];
rz(-1.0370266) q[1];
sx q[1];
rz(-2.7620402) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7354292) q[3];
sx q[3];
rz(-0.87222404) q[3];
sx q[3];
rz(-2.5005831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3520711) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(-2.1748523) q[2];
rz(-3.1294075) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(-0.99408856) q[0];
rz(3.0864691) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(-2.4023043) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2544884) q[0];
sx q[0];
rz(-1.6657889) q[0];
sx q[0];
rz(1.7083005) q[0];
rz(-pi) q[1];
rz(-1.0871068) q[2];
sx q[2];
rz(-1.0488044) q[2];
sx q[2];
rz(-0.99624485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90034396) q[1];
sx q[1];
rz(-0.56300357) q[1];
sx q[1];
rz(2.7808933) q[1];
rz(1.8592632) q[3];
sx q[3];
rz(-0.97086421) q[3];
sx q[3];
rz(0.44655061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5560975) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(0.56376702) q[2];
rz(0.24027763) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(-1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(-2.5575496) q[0];
rz(2.7208327) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(-2.8631794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42785545) q[0];
sx q[0];
rz(-0.58361485) q[0];
sx q[0];
rz(-1.7996656) q[0];
rz(0.39324795) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(0.72088036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6555357) q[1];
sx q[1];
rz(-1.9209407) q[1];
sx q[1];
rz(-2.9085367) q[1];
rz(-1.0233381) q[3];
sx q[3];
rz(-1.7806782) q[3];
sx q[3];
rz(1.7367712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(-2.753624) q[2];
rz(1.7913473) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(-0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(-2.4097089) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(3.0678715) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6442935) q[0];
sx q[0];
rz(-1.6719581) q[0];
sx q[0];
rz(2.0072719) q[0];
rz(-pi) q[1];
rz(-1.6905626) q[2];
sx q[2];
rz(-1.4904067) q[2];
sx q[2];
rz(0.62270852) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8158055) q[1];
sx q[1];
rz(-2.1539638) q[1];
sx q[1];
rz(1.5259685) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57314408) q[3];
sx q[3];
rz(-1.4775652) q[3];
sx q[3];
rz(-1.7236818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5471197) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(0.42868844) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4093032) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(-1.0428585) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(-1.3483378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6207328) q[0];
sx q[0];
rz(-3.0254786) q[0];
sx q[0];
rz(1.7945047) q[0];
rz(-1.6595608) q[2];
sx q[2];
rz(-2.4030493) q[2];
sx q[2];
rz(2.2752787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9546982) q[1];
sx q[1];
rz(-2.8098626) q[1];
sx q[1];
rz(1.3518672) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4753307) q[3];
sx q[3];
rz(-0.30295918) q[3];
sx q[3];
rz(3.016824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(3.0467765) q[2];
rz(1.9684277) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6288347) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(0.082967233) q[2];
sx q[2];
rz(-1.7832179) q[2];
sx q[2];
rz(-0.28945343) q[2];
rz(0.81090609) q[3];
sx q[3];
rz(-1.1699642) q[3];
sx q[3];
rz(-1.5311833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
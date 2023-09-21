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
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(-1.9519238) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0931041) q[0];
sx q[0];
rz(-0.63360533) q[0];
sx q[0];
rz(2.5392169) q[0];
rz(-pi) q[1];
rz(2.4742545) q[2];
sx q[2];
rz(-0.22775209) q[2];
sx q[2];
rz(1.3078794) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78566879) q[1];
sx q[1];
rz(-0.34468109) q[1];
sx q[1];
rz(2.0211401) q[1];
rz(2.7351904) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(-1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(-2.544196) q[2];
rz(1.3655837) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7146724) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(-0.33357099) q[0];
rz(2.0479653) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(-3.0283668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0229867) q[0];
sx q[0];
rz(-1.9959873) q[0];
sx q[0];
rz(3.0990764) q[0];
x q[1];
rz(3.0683238) q[2];
sx q[2];
rz(-1.5511302) q[2];
sx q[2];
rz(2.427223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51616878) q[1];
sx q[1];
rz(-2.5516769) q[1];
sx q[1];
rz(2.3708458) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96666386) q[3];
sx q[3];
rz(-1.3738487) q[3];
sx q[3];
rz(-2.7316824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-2.9193027) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(2.0942028) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(3.0139794) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022284431) q[0];
sx q[0];
rz(-2.2420609) q[0];
sx q[0];
rz(0.22247252) q[0];
x q[1];
rz(0.67655501) q[2];
sx q[2];
rz(-2.1495719) q[2];
sx q[2];
rz(2.3253331) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81880169) q[1];
sx q[1];
rz(-1.4885508) q[1];
sx q[1];
rz(1.8104042) q[1];
x q[2];
rz(2.8887799) q[3];
sx q[3];
rz(-0.95889927) q[3];
sx q[3];
rz(1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3601274) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(2.7919853) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.67868245) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(-2.7754916) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(2.7691832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.303064) q[0];
sx q[0];
rz(-1.4493363) q[0];
sx q[0];
rz(-0.86304201) q[0];
x q[1];
rz(-2.8065368) q[2];
sx q[2];
rz(-1.0154361) q[2];
sx q[2];
rz(2.7044538) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9875033) q[1];
sx q[1];
rz(-1.9470864) q[1];
sx q[1];
rz(-2.980568) q[1];
rz(-pi) q[2];
rz(-0.72752556) q[3];
sx q[3];
rz(-1.1480867) q[3];
sx q[3];
rz(-0.62733516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7238414) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(0.27238971) q[2];
rz(0.90879905) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(0.4549543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(1.1313261) q[0];
rz(-2.5436026) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(0.67684832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8880496) q[0];
sx q[0];
rz(-1.6497668) q[0];
sx q[0];
rz(0.44102863) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7860239) q[2];
sx q[2];
rz(-0.85126801) q[2];
sx q[2];
rz(2.1446251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4557867) q[1];
sx q[1];
rz(-1.5965441) q[1];
sx q[1];
rz(2.8867433) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5464209) q[3];
sx q[3];
rz(-1.1773603) q[3];
sx q[3];
rz(-2.3778621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1260219) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.4591249) q[0];
rz(2.3964264) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(0.93313342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81290302) q[0];
sx q[0];
rz(-1.110154) q[0];
sx q[0];
rz(2.8115389) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60284166) q[2];
sx q[2];
rz(-1.2614998) q[2];
sx q[2];
rz(1.7939292) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7416523) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(-1.9386577) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91306367) q[3];
sx q[3];
rz(-1.7115895) q[3];
sx q[3];
rz(-1.4482244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3559945) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(2.8549426) q[2];
rz(-2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(1.02007) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(2.1405623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9065549) q[0];
sx q[0];
rz(-0.28309238) q[0];
sx q[0];
rz(-1.7049768) q[0];
rz(-pi) q[1];
rz(-2.311003) q[2];
sx q[2];
rz(-0.48362728) q[2];
sx q[2];
rz(2.4110576) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.475358) q[1];
sx q[1];
rz(-1.1821786) q[1];
sx q[1];
rz(2.9734441) q[1];
rz(-0.27490297) q[3];
sx q[3];
rz(-0.66920815) q[3];
sx q[3];
rz(-1.4698524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(-0.10350791) q[2];
rz(-0.59182709) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(-1.7991964) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(0.38988316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99446699) q[0];
sx q[0];
rz(-1.9305221) q[0];
sx q[0];
rz(2.5229182) q[0];
rz(-pi) q[1];
rz(-0.50750081) q[2];
sx q[2];
rz(-1.2050873) q[2];
sx q[2];
rz(-2.4909004) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4854359) q[1];
sx q[1];
rz(-1.1404783) q[1];
sx q[1];
rz(2.7780611) q[1];
rz(-pi) q[2];
rz(2.249986) q[3];
sx q[3];
rz(-1.7779751) q[3];
sx q[3];
rz(-0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.3020246) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(0.71425444) q[2];
rz(-0.8977302) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(0.77734787) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(-2.8651967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24647507) q[0];
sx q[0];
rz(-0.44024375) q[0];
sx q[0];
rz(-0.38245364) q[0];
x q[1];
rz(2.4403964) q[2];
sx q[2];
rz(-2.8465135) q[2];
sx q[2];
rz(1.3311177) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.249835) q[1];
sx q[1];
rz(-1.3502305) q[1];
sx q[1];
rz(2.9794429) q[1];
x q[2];
rz(-1.8180088) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2892264) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(-3.1048807) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-2.7888443) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(-1.030285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093826483) q[0];
sx q[0];
rz(-0.86137912) q[0];
sx q[0];
rz(-1.8295656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68600168) q[2];
sx q[2];
rz(-1.3091607) q[2];
sx q[2];
rz(-2.5930282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5489588) q[1];
sx q[1];
rz(-1.1332129) q[1];
sx q[1];
rz(-1.2195107) q[1];
rz(-pi) q[2];
rz(-1.9029721) q[3];
sx q[3];
rz(-0.98364753) q[3];
sx q[3];
rz(1.4004933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
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
rz(0.43351602) q[2];
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

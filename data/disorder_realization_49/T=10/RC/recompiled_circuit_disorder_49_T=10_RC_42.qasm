OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(3.1363857) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(4.2978573) q[1];
sx q[1];
rz(8.2351091) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7548646) q[0];
sx q[0];
rz(-1.0611738) q[0];
sx q[0];
rz(1.9652363) q[0];
rz(1.7132684) q[2];
sx q[2];
rz(-1.7490897) q[2];
sx q[2];
rz(-1.9879736) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3559239) q[1];
sx q[1];
rz(-0.34468109) q[1];
sx q[1];
rz(-1.1204526) q[1];
rz(-pi) q[2];
rz(-0.94250836) q[3];
sx q[3];
rz(-1.2357467) q[3];
sx q[3];
rz(0.059700746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(-0.5973967) q[2];
rz(-1.776009) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4269203) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(-0.11322583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43464381) q[0];
sx q[0];
rz(-1.5320677) q[0];
sx q[0];
rz(-1.9963272) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5905154) q[2];
sx q[2];
rz(-1.644051) q[2];
sx q[2];
rz(-0.8578701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6254239) q[1];
sx q[1];
rz(-2.5516769) q[1];
sx q[1];
rz(-0.77074681) q[1];
x q[2];
rz(-1.232997) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(0.88463569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.018628) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(0.22228995) q[2];
rz(-0.22953454) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(0.88071841) q[0];
rz(-2.0942028) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-0.12761322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532967) q[0];
sx q[0];
rz(-1.3971546) q[0];
sx q[0];
rz(2.2542473) q[0];
rz(-pi) q[1];
rz(0.80673809) q[2];
sx q[2];
rz(-2.2819937) q[2];
sx q[2];
rz(-2.9850609) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0763921) q[1];
sx q[1];
rz(-0.25307357) q[1];
sx q[1];
rz(1.2364926) q[1];
rz(-pi) q[2];
rz(0.94362887) q[3];
sx q[3];
rz(-1.7769995) q[3];
sx q[3];
rz(-2.9507153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(2.8313417) q[2];
rz(2.7919853) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4629102) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(2.983685) q[0];
rz(-2.7754916) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(2.7691832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87293738) q[0];
sx q[0];
rz(-2.4252709) q[0];
sx q[0];
rz(-1.7563845) q[0];
rz(-0.33505586) q[2];
sx q[2];
rz(-1.0154361) q[2];
sx q[2];
rz(0.43713883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1540893) q[1];
sx q[1];
rz(-1.1945063) q[1];
sx q[1];
rz(-0.16102468) q[1];
rz(-pi) q[2];
rz(1.0286721) q[3];
sx q[3];
rz(-0.91915932) q[3];
sx q[3];
rz(1.8478912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(-0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(-2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(2.0102665) q[0];
rz(0.5979901) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-0.67684832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35447435) q[0];
sx q[0];
rz(-2.0103543) q[0];
sx q[0];
rz(-1.6580824) q[0];
rz(-0.73111515) q[2];
sx q[2];
rz(-1.7321246) q[2];
sx q[2];
rz(-2.4246755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4557867) q[1];
sx q[1];
rz(-1.5965441) q[1];
sx q[1];
rz(2.8867433) q[1];
rz(0.058651794) q[3];
sx q[3];
rz(-0.39415112) q[3];
sx q[3];
rz(-0.7002206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(1.1452902) q[2];
rz(-2.9925313) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(-0.74516621) q[1];
sx q[1];
rz(-2.7909653) q[1];
sx q[1];
rz(2.2084592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695078) q[0];
sx q[0];
rz(-2.5818995) q[0];
sx q[0];
rz(0.99225386) q[0];
x q[1];
rz(-1.9408579) q[2];
sx q[2];
rz(-1.0002631) q[2];
sx q[2];
rz(2.7119315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16490368) q[1];
sx q[1];
rz(-0.36786825) q[1];
sx q[1];
rz(1.5644172) q[1];
rz(-pi) q[2];
rz(-2.228529) q[3];
sx q[3];
rz(-1.4300031) q[3];
sx q[3];
rz(1.6933683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3559945) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(0.24946985) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(-1.4666784) q[0];
rz(-1.02007) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(1.0010304) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46465835) q[0];
sx q[0];
rz(-1.6081728) q[0];
sx q[0];
rz(-1.8514762) q[0];
rz(-0.34044388) q[2];
sx q[2];
rz(-1.921244) q[2];
sx q[2];
rz(0.0705138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96879362) q[1];
sx q[1];
rz(-1.4152923) q[1];
sx q[1];
rz(1.1771727) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4909199) q[3];
sx q[3];
rz(-1.7400029) q[3];
sx q[3];
rz(-2.8229439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.137407) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(3.0380847) q[2];
rz(-0.59182709) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-2.3039968) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(-1.3423963) q[1];
sx q[1];
rz(-1.9806769) q[1];
sx q[1];
rz(2.7517095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33078157) q[0];
sx q[0];
rz(-2.1446052) q[0];
sx q[0];
rz(-2.0033037) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6673442) q[2];
sx q[2];
rz(-0.61605011) q[2];
sx q[2];
rz(-1.6499856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.08323616) q[1];
sx q[1];
rz(-2.5857158) q[1];
sx q[1];
rz(-0.91169375) q[1];
x q[2];
rz(0.26384683) q[3];
sx q[3];
rz(-2.2328394) q[3];
sx q[3];
rz(-0.86316934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8395681) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(2.4273382) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(2.3642448) q[0];
rz(-2.2946987) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(-2.8651967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8951176) q[0];
sx q[0];
rz(-0.44024375) q[0];
sx q[0];
rz(-2.759139) q[0];
rz(-pi) q[1];
rz(2.9133965) q[2];
sx q[2];
rz(-1.7595292) q[2];
sx q[2];
rz(2.2224094) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.249835) q[1];
sx q[1];
rz(-1.7913622) q[1];
sx q[1];
rz(-0.16214976) q[1];
rz(-pi) q[2];
rz(1.3235839) q[3];
sx q[3];
rz(-0.97655481) q[3];
sx q[3];
rz(2.0127399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(0.12750553) q[2];
rz(3.1048807) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(-1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(2.7888443) q[0];
rz(-0.57669512) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83537) q[0];
sx q[0];
rz(-1.3754002) q[0];
sx q[0];
rz(2.4154001) q[0];
rz(-pi) q[1];
rz(-2.7416517) q[2];
sx q[2];
rz(-0.72657864) q[2];
sx q[2];
rz(-0.71619294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5489588) q[1];
sx q[1];
rz(-2.0083798) q[1];
sx q[1];
rz(1.9220819) q[1];
rz(-2.528245) q[3];
sx q[3];
rz(-1.8457335) q[3];
sx q[3];
rz(3.1230694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0673922) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(-2.4009005) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-3.0497131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108903) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(-2.4717992) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(0.40217051) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(-2.4767247) q[3];
sx q[3];
rz(-2.3330199) q[3];
sx q[3];
rz(-1.4521269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

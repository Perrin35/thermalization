OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18093827) q[0];
sx q[0];
rz(-3.0783983) q[0];
sx q[0];
rz(0.59046459) q[0];
rz(-0.2904627) q[1];
sx q[1];
rz(-1.3757179) q[1];
sx q[1];
rz(-3.1380993) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1549415) q[0];
sx q[0];
rz(-1.9738324) q[0];
sx q[0];
rz(2.0986845) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8344049) q[2];
sx q[2];
rz(-2.2279539) q[2];
sx q[2];
rz(-0.023935723) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0776135) q[1];
sx q[1];
rz(-0.031647041) q[1];
sx q[1];
rz(-0.74546234) q[1];
x q[2];
rz(0.9736075) q[3];
sx q[3];
rz(-1.9499717) q[3];
sx q[3];
rz(1.3446863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2575839) q[2];
sx q[2];
rz(-0.85870063) q[2];
sx q[2];
rz(-2.1250471) q[2];
rz(-0.86333418) q[3];
sx q[3];
rz(-0.038067929) q[3];
sx q[3];
rz(1.6421002) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6739552) q[0];
sx q[0];
rz(-1.6683945) q[0];
sx q[0];
rz(0.13482811) q[0];
rz(-0.22678953) q[1];
sx q[1];
rz(-0.062954523) q[1];
sx q[1];
rz(-0.24980587) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9591208) q[0];
sx q[0];
rz(-1.9879436) q[0];
sx q[0];
rz(-0.66138256) q[0];
rz(-pi) q[1];
rz(-1.4748361) q[2];
sx q[2];
rz(-2.6393106) q[2];
sx q[2];
rz(1.407215) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61714482) q[1];
sx q[1];
rz(-3.063319) q[1];
sx q[1];
rz(-3.0039178) q[1];
rz(-pi) q[2];
rz(1.9365015) q[3];
sx q[3];
rz(-1.6964751) q[3];
sx q[3];
rz(2.2179536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9164385) q[2];
sx q[2];
rz(-0.068000451) q[2];
sx q[2];
rz(-2.1615243) q[2];
rz(-2.2465536) q[3];
sx q[3];
rz(-0.74256247) q[3];
sx q[3];
rz(-1.9612954) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5535683) q[0];
sx q[0];
rz(-0.50694412) q[0];
sx q[0];
rz(1.6059599) q[0];
rz(-1.5060679) q[1];
sx q[1];
rz(-2.3160544) q[1];
sx q[1];
rz(2.1855386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6076304) q[0];
sx q[0];
rz(-1.7844943) q[0];
sx q[0];
rz(-1.6713118) q[0];
rz(2.4105775) q[2];
sx q[2];
rz(-1.6201978) q[2];
sx q[2];
rz(1.4412944) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3967665) q[1];
sx q[1];
rz(-1.464572) q[1];
sx q[1];
rz(0.83858072) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6272229) q[3];
sx q[3];
rz(-1.8764495) q[3];
sx q[3];
rz(-0.54600588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8354336) q[2];
sx q[2];
rz(-0.71194887) q[2];
sx q[2];
rz(0.63208675) q[2];
rz(0.88405526) q[3];
sx q[3];
rz(-3.1243117) q[3];
sx q[3];
rz(-0.89068762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4759091) q[0];
sx q[0];
rz(-1.2983687) q[0];
sx q[0];
rz(2.9744398) q[0];
rz(1.2627603) q[1];
sx q[1];
rz(-2.1985168) q[1];
sx q[1];
rz(1.7335588) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3638813) q[0];
sx q[0];
rz(-1.9929307) q[0];
sx q[0];
rz(1.2516096) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0037065119) q[2];
sx q[2];
rz(-1.7825244) q[2];
sx q[2];
rz(2.0293411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2307407) q[1];
sx q[1];
rz(-0.28225275) q[1];
sx q[1];
rz(-0.7464472) q[1];
rz(1.8620328) q[3];
sx q[3];
rz(-2.3197745) q[3];
sx q[3];
rz(0.90437168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35272804) q[2];
sx q[2];
rz(-0.015268607) q[2];
sx q[2];
rz(-0.35066476) q[2];
rz(-2.8777425) q[3];
sx q[3];
rz(-0.00084547384) q[3];
sx q[3];
rz(1.9381757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8769787) q[0];
sx q[0];
rz(-0.74540859) q[0];
sx q[0];
rz(1.718148) q[0];
rz(-0.28165948) q[1];
sx q[1];
rz(-1.5080844) q[1];
sx q[1];
rz(1.3196094) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5277313) q[0];
sx q[0];
rz(-2.0855122) q[0];
sx q[0];
rz(2.6524833) q[0];
x q[1];
rz(2.2036668) q[2];
sx q[2];
rz(-0.034215052) q[2];
sx q[2];
rz(-2.2020415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7235683) q[1];
sx q[1];
rz(-0.71360525) q[1];
sx q[1];
rz(-2.6221971) q[1];
rz(-pi) q[2];
rz(-1.9044456) q[3];
sx q[3];
rz(-1.9556442) q[3];
sx q[3];
rz(-2.9001482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.264297) q[2];
sx q[2];
rz(-0.046961203) q[2];
sx q[2];
rz(1.4287255) q[2];
rz(-1.9127539) q[3];
sx q[3];
rz(-3.0890586) q[3];
sx q[3];
rz(-3.0534237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0838715) q[0];
sx q[0];
rz(-3.0241522) q[0];
sx q[0];
rz(2.591326) q[0];
rz(2.8394207) q[1];
sx q[1];
rz(-1.7487532) q[1];
sx q[1];
rz(-2.7064533) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36212039) q[0];
sx q[0];
rz(-1.605662) q[0];
sx q[0];
rz(0.21135862) q[0];
rz(-pi) q[1];
rz(-1.5698293) q[2];
sx q[2];
rz(-1.5782326) q[2];
sx q[2];
rz(-0.32583562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8543431) q[1];
sx q[1];
rz(-1.6222427) q[1];
sx q[1];
rz(1.1456942) q[1];
rz(0.76558785) q[3];
sx q[3];
rz(-0.88589719) q[3];
sx q[3];
rz(1.1455889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2415194) q[2];
sx q[2];
rz(-3.1406431) q[2];
sx q[2];
rz(2.1692236) q[2];
rz(-2.1420245) q[3];
sx q[3];
rz(-3.1344423) q[3];
sx q[3];
rz(-2.0991367) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480176) q[0];
sx q[0];
rz(-1.8997718) q[0];
sx q[0];
rz(1.0663363) q[0];
rz(-1.713133) q[1];
sx q[1];
rz(-0.26930535) q[1];
sx q[1];
rz(1.2880633) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6664827) q[0];
sx q[0];
rz(-1.2550071) q[0];
sx q[0];
rz(-0.72986365) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7884862) q[2];
sx q[2];
rz(-2.2688008) q[2];
sx q[2];
rz(1.526213) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24780547) q[1];
sx q[1];
rz(-1.5310186) q[1];
sx q[1];
rz(-3.1282148) q[1];
rz(1.9729105) q[3];
sx q[3];
rz(-2.5743604) q[3];
sx q[3];
rz(0.81517787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0224033) q[2];
sx q[2];
rz(-3.0722805) q[2];
sx q[2];
rz(-0.27528396) q[2];
rz(-0.22661181) q[3];
sx q[3];
rz(-0.35609069) q[3];
sx q[3];
rz(2.7039458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76448035) q[0];
sx q[0];
rz(-0.28700101) q[0];
sx q[0];
rz(-0.58718938) q[0];
rz(1.6181234) q[1];
sx q[1];
rz(-2.0240929) q[1];
sx q[1];
rz(1.5709741) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.068741) q[0];
sx q[0];
rz(-2.223461) q[0];
sx q[0];
rz(-2.1101584) q[0];
rz(-pi) q[1];
x q[1];
rz(0.067085548) q[2];
sx q[2];
rz(-2.0338981) q[2];
sx q[2];
rz(2.9912134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5636922) q[1];
sx q[1];
rz(-1.6017388) q[1];
sx q[1];
rz(1.6146496) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32125485) q[3];
sx q[3];
rz(-1.0215838) q[3];
sx q[3];
rz(-2.3556701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1945343) q[2];
sx q[2];
rz(-0.18320601) q[2];
sx q[2];
rz(1.3018695) q[2];
rz(2.7070847) q[3];
sx q[3];
rz(-3.1012111) q[3];
sx q[3];
rz(-2.6935008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7942363) q[0];
sx q[0];
rz(-0.11581049) q[0];
sx q[0];
rz(1.7606803) q[0];
rz(-1.5538838) q[1];
sx q[1];
rz(-2.1788308) q[1];
sx q[1];
rz(-0.10996058) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23672432) q[0];
sx q[0];
rz(-2.3873608) q[0];
sx q[0];
rz(1.2931025) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7919586) q[2];
sx q[2];
rz(-2.261637) q[2];
sx q[2];
rz(-0.12689374) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0011855652) q[1];
sx q[1];
rz(-2.3275536) q[1];
sx q[1];
rz(-0.034239727) q[1];
rz(2.8300555) q[3];
sx q[3];
rz(-1.4185801) q[3];
sx q[3];
rz(1.0675586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3515557) q[2];
sx q[2];
rz(-2.07708) q[2];
sx q[2];
rz(-2.7891187) q[2];
rz(-0.6414837) q[3];
sx q[3];
rz(-0.04647579) q[3];
sx q[3];
rz(2.7789153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8620257) q[0];
sx q[0];
rz(-2.8643705) q[0];
sx q[0];
rz(-2.5773881) q[0];
rz(-1.5203681) q[1];
sx q[1];
rz(-1.0313326) q[1];
sx q[1];
rz(-0.079781562) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.736981) q[0];
sx q[0];
rz(-0.68697646) q[0];
sx q[0];
rz(-1.8115329) q[0];
rz(1.5230082) q[2];
sx q[2];
rz(-1.523368) q[2];
sx q[2];
rz(-3.0509867) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3120739) q[1];
sx q[1];
rz(-0.79527277) q[1];
sx q[1];
rz(-1.9327232) q[1];
rz(-1.6542808) q[3];
sx q[3];
rz(-1.8125839) q[3];
sx q[3];
rz(2.7906362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0148049) q[2];
sx q[2];
rz(-3.1316275) q[2];
sx q[2];
rz(-2.0662181) q[2];
rz(-0.77438313) q[3];
sx q[3];
rz(-0.024024809) q[3];
sx q[3];
rz(2.8703441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423691) q[0];
sx q[0];
rz(-1.5689701) q[0];
sx q[0];
rz(-1.5690621) q[0];
rz(2.6161999) q[1];
sx q[1];
rz(-0.071594302) q[1];
sx q[1];
rz(-0.20803861) q[1];
rz(2.138562) q[2];
sx q[2];
rz(-1.3896639) q[2];
sx q[2];
rz(1.9682503) q[2];
rz(-1.0925068) q[3];
sx q[3];
rz(-2.1989965) q[3];
sx q[3];
rz(0.60529706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

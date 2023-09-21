OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2274269) q[0];
sx q[0];
rz(-1.5072855) q[0];
sx q[0];
rz(-2.8660197) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17272858) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(-1.2781065) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62568134) q[1];
sx q[1];
rz(-1.4206919) q[1];
sx q[1];
rz(-0.73966571) q[1];
rz(-3.1148881) q[3];
sx q[3];
rz(-1.6741802) q[3];
sx q[3];
rz(1.0552693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(2.9075918) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4734128) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(1.0043253) q[0];
rz(1.6218119) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53045814) q[0];
sx q[0];
rz(-2.3728752) q[0];
sx q[0];
rz(0.64972933) q[0];
rz(-0.1728671) q[2];
sx q[2];
rz(-1.7555408) q[2];
sx q[2];
rz(-0.52344054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4304632) q[1];
sx q[1];
rz(-1.7257479) q[1];
sx q[1];
rz(-0.98053812) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6456804) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(-1.0752614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8721547) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(0.088236563) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(-2.2580106) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040366216) q[0];
sx q[0];
rz(-1.4063615) q[0];
sx q[0];
rz(0.40584392) q[0];
rz(-1.9687037) q[2];
sx q[2];
rz(-1.8430317) q[2];
sx q[2];
rz(-0.18323252) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59144943) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(-1.8557465) q[1];
rz(-pi) q[2];
rz(-1.4530573) q[3];
sx q[3];
rz(-2.2265834) q[3];
sx q[3];
rz(-1.1586787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3016004) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.600425) q[0];
sx q[0];
rz(-2.4006002) q[0];
sx q[0];
rz(-0.11359544) q[0];
rz(-pi) q[1];
rz(-1.6646531) q[2];
sx q[2];
rz(-1.8641194) q[2];
sx q[2];
rz(-2.9039127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95851129) q[1];
sx q[1];
rz(-1.6861048) q[1];
sx q[1];
rz(0.99460852) q[1];
rz(-0.22927852) q[3];
sx q[3];
rz(-1.5203272) q[3];
sx q[3];
rz(-0.63427395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52175534) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(0.27553976) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(2.2648947) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(-2.2263288) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1268069) q[0];
sx q[0];
rz(-2.1984587) q[0];
sx q[0];
rz(-3.0833551) q[0];
x q[1];
rz(1.5259597) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(2.5068138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60587864) q[1];
sx q[1];
rz(-1.6026346) q[1];
sx q[1];
rz(-1.4153626) q[1];
x q[2];
rz(-2.4004585) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(-0.65628624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9203732) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-0.39658305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.840344) q[0];
sx q[0];
rz(-2.6123971) q[0];
sx q[0];
rz(1.0521786) q[0];
rz(-1.8814538) q[2];
sx q[2];
rz(-2.0271745) q[2];
sx q[2];
rz(-3.0820456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.086416883) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(-0.1187101) q[1];
rz(2.7820884) q[3];
sx q[3];
rz(-0.61509575) q[3];
sx q[3];
rz(2.2744327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(-1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(-2.0152337) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267245) q[0];
sx q[0];
rz(-1.65979) q[0];
sx q[0];
rz(0.29863775) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0790714) q[2];
sx q[2];
rz(-1.4462496) q[2];
sx q[2];
rz(1.5073656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4451051) q[1];
sx q[1];
rz(-1.3828053) q[1];
sx q[1];
rz(0.45447116) q[1];
rz(-pi) q[2];
rz(1.3133615) q[3];
sx q[3];
rz(-0.35850393) q[3];
sx q[3];
rz(-0.076971019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(3.0287108) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8772802) q[0];
sx q[0];
rz(-2.7240629) q[0];
sx q[0];
rz(2.2420922) q[0];
x q[1];
rz(-0.44185408) q[2];
sx q[2];
rz(-2.0896857) q[2];
sx q[2];
rz(-1.2043124) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7777562) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(-1.9883518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7912346) q[3];
sx q[3];
rz(-2.7535097) q[3];
sx q[3];
rz(-0.9006075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(2.5194871) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(-0.48535767) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982614) q[0];
sx q[0];
rz(-1.4364463) q[0];
sx q[0];
rz(1.5166548) q[0];
x q[1];
rz(0.59818563) q[2];
sx q[2];
rz(-1.7028114) q[2];
sx q[2];
rz(0.26192947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7945054) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(1.1361213) q[1];
x q[2];
rz(-0.63391179) q[3];
sx q[3];
rz(-2.4571107) q[3];
sx q[3];
rz(-2.4661635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-3.1260417) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(-1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(0.5982582) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3529417) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
x q[1];
rz(2.5355849) q[2];
sx q[2];
rz(-0.25642828) q[2];
sx q[2];
rz(-2.3026349) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88565247) q[1];
sx q[1];
rz(-1.2068032) q[1];
sx q[1];
rz(0.10140681) q[1];
rz(-pi) q[2];
rz(-0.80501276) q[3];
sx q[3];
rz(-2.9152438) q[3];
sx q[3];
rz(-1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.839529) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(2.5219953) q[2];
sx q[2];
rz(-1.4705114) q[2];
sx q[2];
rz(0.32497139) q[2];
rz(2.0486352) q[3];
sx q[3];
rz(-0.7328877) q[3];
sx q[3];
rz(2.4536798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
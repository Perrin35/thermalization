OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(2.7603005) q[1];
sx q[1];
rz(-2.5420904) q[1];
sx q[1];
rz(-1.376027) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24297548) q[0];
sx q[0];
rz(-0.78865047) q[0];
sx q[0];
rz(-2.2126161) q[0];
rz(-2.0618093) q[2];
sx q[2];
rz(-0.91677374) q[2];
sx q[2];
rz(0.71066463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0378117) q[1];
sx q[1];
rz(-1.722464) q[1];
sx q[1];
rz(-1.2396461) q[1];
rz(-0.052470603) q[3];
sx q[3];
rz(-0.82004181) q[3];
sx q[3];
rz(0.47580556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(2.1526745) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20724021) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(-0.99769366) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(2.4172799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99740072) q[0];
sx q[0];
rz(-1.2574982) q[0];
sx q[0];
rz(-0.86667592) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19947796) q[2];
sx q[2];
rz(-1.4885159) q[2];
sx q[2];
rz(-2.2897838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7608632) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(-2.1058583) q[1];
rz(-pi) q[2];
rz(2.4507387) q[3];
sx q[3];
rz(-0.30403954) q[3];
sx q[3];
rz(0.090304852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(-0.13606717) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(-0.46229258) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.2190855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8273979) q[0];
sx q[0];
rz(-0.44566804) q[0];
sx q[0];
rz(0.92339869) q[0];
rz(-pi) q[1];
rz(2.2723324) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(-1.9321835) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29875007) q[1];
sx q[1];
rz(-1.0423653) q[1];
sx q[1];
rz(-0.28038402) q[1];
x q[2];
rz(1.3027906) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(-2.5915495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7200155) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(-3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(2.6233327) q[0];
rz(2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(2.3148361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087238833) q[0];
sx q[0];
rz(-2.0728489) q[0];
sx q[0];
rz(-1.0994083) q[0];
rz(-pi) q[1];
rz(-2.3739359) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(-1.3354288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0747448) q[1];
sx q[1];
rz(-1.2859935) q[1];
sx q[1];
rz(0.39798255) q[1];
rz(-pi) q[2];
rz(-1.8825674) q[3];
sx q[3];
rz(-2.2707553) q[3];
sx q[3];
rz(-1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-2.0157053) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1681686) q[0];
sx q[0];
rz(-1.0963206) q[0];
sx q[0];
rz(-3.0939177) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4777849) q[2];
sx q[2];
rz(-1.2592578) q[2];
sx q[2];
rz(-2.2500492) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83708977) q[1];
sx q[1];
rz(-1.883852) q[1];
sx q[1];
rz(-0.00043991107) q[1];
rz(-pi) q[2];
rz(-0.40163715) q[3];
sx q[3];
rz(-1.4816227) q[3];
sx q[3];
rz(3.0107486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29331648) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(0.69333386) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(1.7745811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54129823) q[0];
sx q[0];
rz(-1.5026662) q[0];
sx q[0];
rz(-3.1085204) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96732803) q[2];
sx q[2];
rz(-0.8562932) q[2];
sx q[2];
rz(-0.29010233) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94783084) q[1];
sx q[1];
rz(-1.281922) q[1];
sx q[1];
rz(3.0572901) q[1];
x q[2];
rz(1.4671765) q[3];
sx q[3];
rz(-1.2195671) q[3];
sx q[3];
rz(-0.87408376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-2.7381251) q[2];
rz(-0.48163313) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(2.6223555) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.9708995) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3160352) q[0];
sx q[0];
rz(-0.023840126) q[0];
sx q[0];
rz(2.8799879) q[0];
rz(3.0253719) q[2];
sx q[2];
rz(-1.6214317) q[2];
sx q[2];
rz(-1.2223787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2611744) q[1];
sx q[1];
rz(-0.92217991) q[1];
sx q[1];
rz(-2.3748114) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6550754) q[3];
sx q[3];
rz(-2.5518637) q[3];
sx q[3];
rz(3.1004578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(-2.5308385) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(0.92774123) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(1.7469453) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0028249) q[0];
sx q[0];
rz(-1.7768304) q[0];
sx q[0];
rz(-0.10319184) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6693194) q[2];
sx q[2];
rz(-1.7211282) q[2];
sx q[2];
rz(2.3538102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7568126) q[1];
sx q[1];
rz(-1.2287178) q[1];
sx q[1];
rz(0.53201075) q[1];
x q[2];
rz(2.2046702) q[3];
sx q[3];
rz(-0.86808944) q[3];
sx q[3];
rz(-3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(2.6867552) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4421473) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(3.033175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3603044) q[0];
sx q[0];
rz(-2.4753248) q[0];
sx q[0];
rz(-1.6643307) q[0];
rz(1.7494781) q[2];
sx q[2];
rz(-1.5249426) q[2];
sx q[2];
rz(2.1669441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11549599) q[1];
sx q[1];
rz(-2.8671088) q[1];
sx q[1];
rz(0.91067578) q[1];
rz(-pi) q[2];
rz(-0.18687825) q[3];
sx q[3];
rz(-0.87516057) q[3];
sx q[3];
rz(0.38284341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(2.4627731) q[0];
rz(2.7774096) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-0.055158786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77627742) q[0];
sx q[0];
rz(-0.3826097) q[0];
sx q[0];
rz(-2.9676653) q[0];
rz(-pi) q[1];
rz(0.55969413) q[2];
sx q[2];
rz(-1.7494546) q[2];
sx q[2];
rz(2.6924804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2436279) q[1];
sx q[1];
rz(-2.1568255) q[1];
sx q[1];
rz(0.90378739) q[1];
rz(-1.1147898) q[3];
sx q[3];
rz(-1.5378693) q[3];
sx q[3];
rz(1.1820716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(2.2035051) q[3];
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
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162083) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(1.6677042) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(1.5260074) q[3];
sx q[3];
rz(-1.3309892) q[3];
sx q[3];
rz(0.25467024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
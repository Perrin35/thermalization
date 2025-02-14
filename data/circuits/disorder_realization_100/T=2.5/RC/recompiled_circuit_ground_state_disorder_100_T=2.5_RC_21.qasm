OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7772726) q[0];
sx q[0];
rz(-0.84982818) q[0];
sx q[0];
rz(-1.538869) q[0];
rz(-0.76378769) q[1];
sx q[1];
rz(-0.053215947) q[1];
sx q[1];
rz(0.8344714) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6042701) q[0];
sx q[0];
rz(-1.6217082) q[0];
sx q[0];
rz(-1.0983653) q[0];
x q[1];
rz(2.5080483) q[2];
sx q[2];
rz(-2.0764645) q[2];
sx q[2];
rz(-1.0552931) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38307962) q[1];
sx q[1];
rz(-2.5316409) q[1];
sx q[1];
rz(1.8083231) q[1];
x q[2];
rz(2.2356307) q[3];
sx q[3];
rz(-0.49091456) q[3];
sx q[3];
rz(2.5708446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6795464) q[2];
sx q[2];
rz(-1.5507689) q[2];
sx q[2];
rz(0.66590434) q[2];
rz(2.7146961) q[3];
sx q[3];
rz(-0.96564966) q[3];
sx q[3];
rz(0.20195937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17215984) q[0];
sx q[0];
rz(-0.96123022) q[0];
sx q[0];
rz(2.6184671) q[0];
rz(-0.41921774) q[1];
sx q[1];
rz(-2.9328177) q[1];
sx q[1];
rz(2.8255393) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2284477) q[0];
sx q[0];
rz(-0.70437594) q[0];
sx q[0];
rz(0.26972847) q[0];
rz(-pi) q[1];
rz(1.5957556) q[2];
sx q[2];
rz(-1.7169368) q[2];
sx q[2];
rz(-1.6294831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16686116) q[1];
sx q[1];
rz(-0.63878794) q[1];
sx q[1];
rz(2.0898934) q[1];
rz(0.42556588) q[3];
sx q[3];
rz(-1.4168973) q[3];
sx q[3];
rz(-2.9776412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6009964) q[2];
sx q[2];
rz(-0.54800802) q[2];
sx q[2];
rz(1.126464) q[2];
rz(0.90538853) q[3];
sx q[3];
rz(-2.0784046) q[3];
sx q[3];
rz(-1.5409013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.044540502) q[0];
sx q[0];
rz(-1.4690761) q[0];
sx q[0];
rz(1.6497521) q[0];
rz(2.7482694) q[1];
sx q[1];
rz(-1.2231188) q[1];
sx q[1];
rz(0.16600674) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1309083) q[0];
sx q[0];
rz(-1.8839399) q[0];
sx q[0];
rz(-0.74096276) q[0];
rz(-pi) q[1];
rz(-2.2691378) q[2];
sx q[2];
rz(-1.4373682) q[2];
sx q[2];
rz(-0.39233559) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9505009) q[1];
sx q[1];
rz(-2.3984004) q[1];
sx q[1];
rz(-2.8673001) q[1];
rz(0.29434926) q[3];
sx q[3];
rz(-2.6913319) q[3];
sx q[3];
rz(3.1401304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7321221) q[2];
sx q[2];
rz(-2.8889443) q[2];
sx q[2];
rz(2.9518993) q[2];
rz(2.7267552) q[3];
sx q[3];
rz(-1.8388351) q[3];
sx q[3];
rz(0.10492575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.6787427) q[0];
sx q[0];
rz(-0.78063709) q[0];
sx q[0];
rz(-1.8829874) q[0];
rz(1.3858093) q[1];
sx q[1];
rz(-0.36711991) q[1];
sx q[1];
rz(2.3585336) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5393128) q[0];
sx q[0];
rz(-0.040518213) q[0];
sx q[0];
rz(-2.4862073) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5747397) q[2];
sx q[2];
rz(-2.3794305) q[2];
sx q[2];
rz(2.5944425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.051797051) q[1];
sx q[1];
rz(-2.0690114) q[1];
sx q[1];
rz(-0.55023043) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7854758) q[3];
sx q[3];
rz(-1.5443729) q[3];
sx q[3];
rz(1.6230447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.9328402) q[2];
sx q[2];
rz(-1.1609266) q[2];
sx q[2];
rz(-1.0932385) q[2];
rz(2.7161963) q[3];
sx q[3];
rz(-1.930611) q[3];
sx q[3];
rz(-1.0824599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071674034) q[0];
sx q[0];
rz(-2.2921083) q[0];
sx q[0];
rz(-1.7639532) q[0];
rz(-2.2397974) q[1];
sx q[1];
rz(-1.8865562) q[1];
sx q[1];
rz(-0.76595438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1655584) q[0];
sx q[0];
rz(-1.9355825) q[0];
sx q[0];
rz(2.4976612) q[0];
x q[1];
rz(1.8661921) q[2];
sx q[2];
rz(-0.89933853) q[2];
sx q[2];
rz(1.7341933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6038725) q[1];
sx q[1];
rz(-1.0728379) q[1];
sx q[1];
rz(-1.7437976) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28217989) q[3];
sx q[3];
rz(-1.9934142) q[3];
sx q[3];
rz(0.19920345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9506266) q[2];
sx q[2];
rz(-2.2910255) q[2];
sx q[2];
rz(1.4858049) q[2];
rz(-1.1367669) q[3];
sx q[3];
rz(-2.6906689) q[3];
sx q[3];
rz(-0.76751149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9039827) q[0];
sx q[0];
rz(-0.5466277) q[0];
sx q[0];
rz(0.23608635) q[0];
rz(1.6646594) q[1];
sx q[1];
rz(-1.6496941) q[1];
sx q[1];
rz(-1.1621071) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88173819) q[0];
sx q[0];
rz(-0.82603541) q[0];
sx q[0];
rz(-1.8709454) q[0];
rz(-pi) q[1];
rz(-2.7323551) q[2];
sx q[2];
rz(-1.6820568) q[2];
sx q[2];
rz(1.8429327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2047822) q[1];
sx q[1];
rz(-2.6814662) q[1];
sx q[1];
rz(2.2304753) q[1];
x q[2];
rz(0.44329537) q[3];
sx q[3];
rz(-0.05847769) q[3];
sx q[3];
rz(-0.33340461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7621925) q[2];
sx q[2];
rz(-2.6683932) q[2];
sx q[2];
rz(-3.1326478) q[2];
rz(-1.3116607) q[3];
sx q[3];
rz(-1.0511755) q[3];
sx q[3];
rz(1.8478307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0084429) q[0];
sx q[0];
rz(-0.75859469) q[0];
sx q[0];
rz(1.5593418) q[0];
rz(2.3186231) q[1];
sx q[1];
rz(-0.4987078) q[1];
sx q[1];
rz(1.71436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8860075) q[0];
sx q[0];
rz(-2.581121) q[0];
sx q[0];
rz(2.1917402) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4983369) q[2];
sx q[2];
rz(-1.7417041) q[2];
sx q[2];
rz(2.3115445) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8135951) q[1];
sx q[1];
rz(-0.57472052) q[1];
sx q[1];
rz(-1.5411472) q[1];
x q[2];
rz(-2.8893647) q[3];
sx q[3];
rz(-0.42740932) q[3];
sx q[3];
rz(-0.073425882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83671776) q[2];
sx q[2];
rz(-2.0676401) q[2];
sx q[2];
rz(0.3982591) q[2];
rz(0.38175976) q[3];
sx q[3];
rz(-1.441381) q[3];
sx q[3];
rz(1.0086584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037982382) q[0];
sx q[0];
rz(-1.4382265) q[0];
sx q[0];
rz(-1.2305228) q[0];
rz(-1.1331753) q[1];
sx q[1];
rz(-0.76825348) q[1];
sx q[1];
rz(-0.59655985) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8997079) q[0];
sx q[0];
rz(-1.5615789) q[0];
sx q[0];
rz(2.8901744) q[0];
rz(1.9279154) q[2];
sx q[2];
rz(-0.51110635) q[2];
sx q[2];
rz(2.7765283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79227263) q[1];
sx q[1];
rz(-1.3111924) q[1];
sx q[1];
rz(2.1672897) q[1];
x q[2];
rz(-1.9539709) q[3];
sx q[3];
rz(-0.98742332) q[3];
sx q[3];
rz(-0.72900984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8489428) q[2];
sx q[2];
rz(-2.5299447) q[2];
sx q[2];
rz(-1.4673648) q[2];
rz(-2.8749706) q[3];
sx q[3];
rz(-2.1610503) q[3];
sx q[3];
rz(-2.9149741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.6999231) q[0];
sx q[0];
rz(-0.00049463153) q[0];
sx q[0];
rz(-2.0200404) q[0];
rz(0.17555217) q[1];
sx q[1];
rz(-0.67070621) q[1];
sx q[1];
rz(0.49351722) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.314007) q[0];
sx q[0];
rz(-0.092484154) q[0];
sx q[0];
rz(2.2580763) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7337903) q[2];
sx q[2];
rz(-0.85821292) q[2];
sx q[2];
rz(-1.4607371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9914876) q[1];
sx q[1];
rz(-2.5511523) q[1];
sx q[1];
rz(-1.1439249) q[1];
x q[2];
rz(1.4021691) q[3];
sx q[3];
rz(-2.6624749) q[3];
sx q[3];
rz(1.428133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4404122) q[2];
sx q[2];
rz(-0.36786011) q[2];
sx q[2];
rz(-2.5555447) q[2];
rz(1.690257) q[3];
sx q[3];
rz(-1.3380932) q[3];
sx q[3];
rz(0.51226789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.119809) q[0];
sx q[0];
rz(-2.1202987) q[0];
sx q[0];
rz(2.9721066) q[0];
rz(-0.64300621) q[1];
sx q[1];
rz(-2.2702859) q[1];
sx q[1];
rz(-0.51173425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62637037) q[0];
sx q[0];
rz(-3.0128509) q[0];
sx q[0];
rz(-0.66855742) q[0];
rz(-pi) q[1];
rz(-0.24004748) q[2];
sx q[2];
rz(-2.8721951) q[2];
sx q[2];
rz(0.9123957) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3822126) q[1];
sx q[1];
rz(-2.4649715) q[1];
sx q[1];
rz(-1.8130019) q[1];
x q[2];
rz(-2.0848537) q[3];
sx q[3];
rz(-1.3680046) q[3];
sx q[3];
rz(-2.6904484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52342755) q[2];
sx q[2];
rz(-0.7343505) q[2];
sx q[2];
rz(-3.1132474) q[2];
rz(0.47447765) q[3];
sx q[3];
rz(-0.35895434) q[3];
sx q[3];
rz(2.1019782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2059373) q[0];
sx q[0];
rz(-1.3858929) q[0];
sx q[0];
rz(-0.53336597) q[0];
rz(0.43279303) q[1];
sx q[1];
rz(-1.5846033) q[1];
sx q[1];
rz(-0.63176647) q[1];
rz(-2.9262056) q[2];
sx q[2];
rz(-1.4949847) q[2];
sx q[2];
rz(1.1204213) q[2];
rz(1.778936) q[3];
sx q[3];
rz(-1.647462) q[3];
sx q[3];
rz(0.36858222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

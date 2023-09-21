OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(6.8350514) q[0];
sx q[0];
rz(9.4466136) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(2.9266761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41479933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(0.8154072) q[0];
x q[1];
rz(0.46918842) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(-0.56433041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51337459) q[1];
sx q[1];
rz(-2.3623423) q[1];
sx q[1];
rz(0.90374225) q[1];
rz(-pi) q[2];
rz(-2.5932556) q[3];
sx q[3];
rz(-2.3468446) q[3];
sx q[3];
rz(-2.7693975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(-2.2136097) q[0];
rz(-1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-2.2448418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0052764) q[0];
sx q[0];
rz(-1.7668084) q[0];
sx q[0];
rz(-2.2851903) q[0];
x q[1];
rz(-2.9727544) q[2];
sx q[2];
rz(-1.4988006) q[2];
sx q[2];
rz(0.26693401) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0838544) q[1];
sx q[1];
rz(-1.6946304) q[1];
sx q[1];
rz(0.79353516) q[1];
x q[2];
rz(-1.1284626) q[3];
sx q[3];
rz(-2.3395174) q[3];
sx q[3];
rz(-2.341552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26560489) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(2.7868328) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-2.1133912) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(0.43513402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1895771) q[0];
sx q[0];
rz(-0.36882419) q[0];
sx q[0];
rz(0.30216218) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8912796) q[2];
sx q[2];
rz(-1.7134943) q[2];
sx q[2];
rz(-0.38562361) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55089009) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(-2.9571556) q[1];
x q[2];
rz(1.2262672) q[3];
sx q[3];
rz(-2.3039654) q[3];
sx q[3];
rz(-1.6813577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(-3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(0.91745013) q[0];
rz(0.67287412) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53752758) q[0];
sx q[0];
rz(-0.91188216) q[0];
sx q[0];
rz(1.5327246) q[0];
rz(-pi) q[1];
rz(-0.8226383) q[2];
sx q[2];
rz(-0.3974786) q[2];
sx q[2];
rz(2.9646404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59135624) q[1];
sx q[1];
rz(-1.1657506) q[1];
sx q[1];
rz(2.3637799) q[1];
x q[2];
rz(0.78062765) q[3];
sx q[3];
rz(-2.160191) q[3];
sx q[3];
rz(-1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(-2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-1.1013793) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6081776) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(1.7081225) q[0];
rz(-pi) q[1];
rz(1.9031992) q[2];
sx q[2];
rz(-1.9661511) q[2];
sx q[2];
rz(-0.3619286) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2561803) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(2.1229565) q[1];
rz(-0.080658241) q[3];
sx q[3];
rz(-2.6622052) q[3];
sx q[3];
rz(-0.56848923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7197363) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(2.6110113) q[2];
rz(-1.4060219) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(-1.6249599) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96777746) q[0];
sx q[0];
rz(-0.41668188) q[0];
sx q[0];
rz(1.5370876) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2595915) q[2];
sx q[2];
rz(-1.1278369) q[2];
sx q[2];
rz(-1.358658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1291618) q[1];
sx q[1];
rz(-0.42814246) q[1];
sx q[1];
rz(1.3151602) q[1];
rz(-pi) q[2];
rz(0.63038007) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(-2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(-2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(0.75659928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3698759) q[0];
sx q[0];
rz(-1.6635832) q[0];
sx q[0];
rz(-1.7232399) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43811626) q[2];
sx q[2];
rz(-2.6926059) q[2];
sx q[2];
rz(0.24030906) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8826897) q[1];
sx q[1];
rz(-1.8001302) q[1];
sx q[1];
rz(-0.16112666) q[1];
x q[2];
rz(-0.19161253) q[3];
sx q[3];
rz(-1.5125456) q[3];
sx q[3];
rz(0.91764698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0044272) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(0.19443092) q[2];
rz(-2.2284609) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(2.1527122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522404) q[0];
sx q[0];
rz(-1.1712495) q[0];
sx q[0];
rz(0.35895343) q[0];
rz(-pi) q[1];
rz(-2.9724389) q[2];
sx q[2];
rz(-3.0468468) q[2];
sx q[2];
rz(2.6731238) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5576396) q[1];
sx q[1];
rz(-1.156731) q[1];
sx q[1];
rz(3.0686892) q[1];
rz(2.057468) q[3];
sx q[3];
rz(-1.1063965) q[3];
sx q[3];
rz(-1.407479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6797592) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(0.87337714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4163923) q[0];
sx q[0];
rz(-1.7500449) q[0];
sx q[0];
rz(-3.0730625) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4344425) q[2];
sx q[2];
rz(-1.7169723) q[2];
sx q[2];
rz(1.2863976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5492591) q[1];
sx q[1];
rz(-2.1570286) q[1];
sx q[1];
rz(1.1997644) q[1];
x q[2];
rz(-2.14823) q[3];
sx q[3];
rz(-1.1953029) q[3];
sx q[3];
rz(1.7253699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(2.7673289) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69797126) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(2.1886254) q[0];
rz(-0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.3964765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8499334) q[0];
sx q[0];
rz(-2.8072661) q[0];
sx q[0];
rz(2.8773984) q[0];
rz(0.65812494) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(2.4326774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6203306) q[1];
sx q[1];
rz(-1.8552823) q[1];
sx q[1];
rz(-1.2465338) q[1];
rz(-0.20874899) q[3];
sx q[3];
rz(-0.26780805) q[3];
sx q[3];
rz(0.97313125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(2.9818025) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(-2.8543499) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0970584) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-1.3143905) q[2];
sx q[2];
rz(-2.495043) q[2];
sx q[2];
rz(-1.3174353) q[2];
rz(1.521048) q[3];
sx q[3];
rz(-2.1038901) q[3];
sx q[3];
rz(-0.62705561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

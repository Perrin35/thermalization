OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5069198) q[0];
sx q[0];
rz(-1.4022175) q[0];
sx q[0];
rz(0.75048598) q[0];
rz(-3.0783202) q[1];
sx q[1];
rz(-0.81875357) q[1];
sx q[1];
rz(-1.0190581) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69438231) q[0];
sx q[0];
rz(-1.3827494) q[0];
sx q[0];
rz(-2.5716554) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2804777) q[2];
sx q[2];
rz(-1.3657331) q[2];
sx q[2];
rz(0.67180639) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77324919) q[1];
sx q[1];
rz(-1.7071525) q[1];
sx q[1];
rz(0.58500599) q[1];
rz(1.9492703) q[3];
sx q[3];
rz(-1.0590388) q[3];
sx q[3];
rz(3.0853693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84746209) q[2];
sx q[2];
rz(-2.1475466) q[2];
sx q[2];
rz(-1.4685941) q[2];
rz(-2.9351506) q[3];
sx q[3];
rz(-1.3279746) q[3];
sx q[3];
rz(0.65202057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.9701397) q[0];
sx q[0];
rz(-1.7296706) q[0];
sx q[0];
rz(-3.0385802) q[0];
rz(-2.7439694) q[1];
sx q[1];
rz(-2.7151974) q[1];
sx q[1];
rz(-2.2893589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641076) q[0];
sx q[0];
rz(-1.0380569) q[0];
sx q[0];
rz(-3.0716776) q[0];
rz(-1.909341) q[2];
sx q[2];
rz(-1.3024835) q[2];
sx q[2];
rz(1.2635096) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0956079) q[1];
sx q[1];
rz(-1.5569481) q[1];
sx q[1];
rz(1.0471837) q[1];
rz(-2.9217124) q[3];
sx q[3];
rz(-0.19769719) q[3];
sx q[3];
rz(1.8085331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6323866) q[2];
sx q[2];
rz(-1.0078603) q[2];
sx q[2];
rz(2.3243135) q[2];
rz(-0.20538524) q[3];
sx q[3];
rz(-0.21586625) q[3];
sx q[3];
rz(0.5425905) q[3];
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
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3096454) q[0];
sx q[0];
rz(-2.0537856) q[0];
sx q[0];
rz(0.52016869) q[0];
rz(2.4179516) q[1];
sx q[1];
rz(-1.2836722) q[1];
sx q[1];
rz(2.9626194) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13143927) q[0];
sx q[0];
rz(-2.4418533) q[0];
sx q[0];
rz(0.72314163) q[0];
x q[1];
rz(-0.34988259) q[2];
sx q[2];
rz(-1.3777108) q[2];
sx q[2];
rz(-1.3909222) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3368028) q[1];
sx q[1];
rz(-1.8095142) q[1];
sx q[1];
rz(1.6126339) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8427027) q[3];
sx q[3];
rz(-0.93122175) q[3];
sx q[3];
rz(-3.0323882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.2554864) q[2];
sx q[2];
rz(-0.28033689) q[2];
sx q[2];
rz(-2.172016) q[2];
rz(2.7995836) q[3];
sx q[3];
rz(-1.0229144) q[3];
sx q[3];
rz(-1.8146993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.904357) q[0];
sx q[0];
rz(-0.80191737) q[0];
sx q[0];
rz(3.0661769) q[0];
rz(-2.8341809) q[1];
sx q[1];
rz(-1.1630029) q[1];
sx q[1];
rz(2.7154162) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89875561) q[0];
sx q[0];
rz(-0.97427827) q[0];
sx q[0];
rz(-0.92208029) q[0];
x q[1];
rz(0.46838872) q[2];
sx q[2];
rz(-1.3184183) q[2];
sx q[2];
rz(0.017367432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4170917) q[1];
sx q[1];
rz(-1.0001273) q[1];
sx q[1];
rz(0.23391926) q[1];
x q[2];
rz(0.05699031) q[3];
sx q[3];
rz(-0.40613031) q[3];
sx q[3];
rz(-0.95088357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9093466) q[2];
sx q[2];
rz(-1.987395) q[2];
sx q[2];
rz(-2.2171891) q[2];
rz(2.038548) q[3];
sx q[3];
rz(-2.0274935) q[3];
sx q[3];
rz(-1.2801142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1396609) q[0];
sx q[0];
rz(-0.49134555) q[0];
sx q[0];
rz(-2.370148) q[0];
rz(1.8292684) q[1];
sx q[1];
rz(-1.3714561) q[1];
sx q[1];
rz(-1.2244474) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0010013) q[0];
sx q[0];
rz(-0.13990211) q[0];
sx q[0];
rz(1.3683659) q[0];
rz(0.61842711) q[2];
sx q[2];
rz(-1.5723516) q[2];
sx q[2];
rz(-2.8032982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6836732) q[1];
sx q[1];
rz(-0.39622849) q[1];
sx q[1];
rz(1.8093682) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1125074) q[3];
sx q[3];
rz(-0.84064129) q[3];
sx q[3];
rz(2.2793913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9839342) q[2];
sx q[2];
rz(-0.57571405) q[2];
sx q[2];
rz(-2.0470587) q[2];
rz(-3.0648699) q[3];
sx q[3];
rz(-0.50944296) q[3];
sx q[3];
rz(-2.4491687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43129608) q[0];
sx q[0];
rz(-0.65839473) q[0];
sx q[0];
rz(-0.21507138) q[0];
rz(-0.98945016) q[1];
sx q[1];
rz(-0.96885252) q[1];
sx q[1];
rz(1.4873803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3475018) q[0];
sx q[0];
rz(-0.54301822) q[0];
sx q[0];
rz(0.31859441) q[0];
rz(-pi) q[1];
rz(0.62724386) q[2];
sx q[2];
rz(-2.6656642) q[2];
sx q[2];
rz(1.5679596) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.038092117) q[1];
sx q[1];
rz(-1.9327144) q[1];
sx q[1];
rz(2.017971) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0607674) q[3];
sx q[3];
rz(-0.96749159) q[3];
sx q[3];
rz(1.5321466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4789077) q[2];
sx q[2];
rz(-1.4143133) q[2];
sx q[2];
rz(1.5773845) q[2];
rz(2.1936737) q[3];
sx q[3];
rz(-1.0736059) q[3];
sx q[3];
rz(-1.7387559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.3941512) q[0];
sx q[0];
rz(-1.489137) q[0];
sx q[0];
rz(-2.5519651) q[0];
rz(1.6472752) q[1];
sx q[1];
rz(-1.3811771) q[1];
sx q[1];
rz(0.078701198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704961) q[0];
sx q[0];
rz(-1.9645623) q[0];
sx q[0];
rz(-1.257819) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3316657) q[2];
sx q[2];
rz(-2.0239365) q[2];
sx q[2];
rz(-1.5117642) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5906395) q[1];
sx q[1];
rz(-1.3627684) q[1];
sx q[1];
rz(2.8880638) q[1];
rz(2.3669993) q[3];
sx q[3];
rz(-2.4586185) q[3];
sx q[3];
rz(-2.3247064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2176723) q[2];
sx q[2];
rz(-2.2675026) q[2];
sx q[2];
rz(3.0295642) q[2];
rz(2.7364386) q[3];
sx q[3];
rz(-1.3978037) q[3];
sx q[3];
rz(-2.9356094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9236295) q[0];
sx q[0];
rz(-2.1681652) q[0];
sx q[0];
rz(-1.4189036) q[0];
rz(2.0619242) q[1];
sx q[1];
rz(-0.38366145) q[1];
sx q[1];
rz(1.4071646) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16033835) q[0];
sx q[0];
rz(-1.8869074) q[0];
sx q[0];
rz(-0.42686064) q[0];
rz(-pi) q[1];
rz(0.17282829) q[2];
sx q[2];
rz(-0.43718526) q[2];
sx q[2];
rz(0.83807105) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8720334) q[1];
sx q[1];
rz(-2.5692794) q[1];
sx q[1];
rz(0.78929269) q[1];
rz(-pi) q[2];
rz(-1.2958762) q[3];
sx q[3];
rz(-1.3310199) q[3];
sx q[3];
rz(1.4942684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.63274038) q[2];
sx q[2];
rz(-0.036157046) q[2];
sx q[2];
rz(0.73842326) q[2];
rz(-2.9449055) q[3];
sx q[3];
rz(-0.6179215) q[3];
sx q[3];
rz(-0.64016199) q[3];
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
rz(-pi/2) q[3];
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
rz(2.7554373) q[0];
sx q[0];
rz(-0.25179395) q[0];
sx q[0];
rz(-3.1016896) q[0];
rz(2.9207322) q[1];
sx q[1];
rz(-1.618229) q[1];
sx q[1];
rz(1.9853282) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45922088) q[0];
sx q[0];
rz(-0.30293884) q[0];
sx q[0];
rz(-2.0762864) q[0];
x q[1];
rz(-0.21577452) q[2];
sx q[2];
rz(-1.3808983) q[2];
sx q[2];
rz(-0.18958651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58312341) q[1];
sx q[1];
rz(-1.236318) q[1];
sx q[1];
rz(-0.98812798) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74083225) q[3];
sx q[3];
rz(-0.15578905) q[3];
sx q[3];
rz(0.039905101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7463189) q[2];
sx q[2];
rz(-0.36397448) q[2];
sx q[2];
rz(3.086997) q[2];
rz(2.3553081) q[3];
sx q[3];
rz(-1.1966642) q[3];
sx q[3];
rz(-1.1597077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7603124) q[0];
sx q[0];
rz(-0.78802839) q[0];
sx q[0];
rz(0.81083167) q[0];
rz(0.52931085) q[1];
sx q[1];
rz(-1.7050754) q[1];
sx q[1];
rz(-1.1108105) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62880077) q[0];
sx q[0];
rz(-1.1847825) q[0];
sx q[0];
rz(2.5284472) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6535437) q[2];
sx q[2];
rz(-2.077861) q[2];
sx q[2];
rz(1.6973073) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6880092) q[1];
sx q[1];
rz(-2.3613259) q[1];
sx q[1];
rz(0.23201402) q[1];
rz(-pi) q[2];
rz(-1.1485898) q[3];
sx q[3];
rz(-0.65311382) q[3];
sx q[3];
rz(-0.70142581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3436053) q[2];
sx q[2];
rz(-0.56861773) q[2];
sx q[2];
rz(-2.9000751) q[2];
rz(-2.8999515) q[3];
sx q[3];
rz(-1.4915978) q[3];
sx q[3];
rz(-0.1514761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80736137) q[0];
sx q[0];
rz(-1.7932899) q[0];
sx q[0];
rz(0.67247969) q[0];
rz(-0.45202759) q[1];
sx q[1];
rz(-0.94257911) q[1];
sx q[1];
rz(-0.51530757) q[1];
rz(-1.4448901) q[2];
sx q[2];
rz(-0.44582146) q[2];
sx q[2];
rz(1.3498342) q[2];
rz(0.47361278) q[3];
sx q[3];
rz(-2.0423741) q[3];
sx q[3];
rz(-3.1033677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

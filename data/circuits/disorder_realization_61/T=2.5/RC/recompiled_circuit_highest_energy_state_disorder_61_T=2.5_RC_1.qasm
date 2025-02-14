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
rz(-2.6391368) q[0];
sx q[0];
rz(-1.1075736) q[0];
sx q[0];
rz(-1.3753608) q[0];
rz(-0.31887588) q[1];
sx q[1];
rz(-1.640929) q[1];
sx q[1];
rz(0.73135102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7773444) q[0];
sx q[0];
rz(-0.70667446) q[0];
sx q[0];
rz(-1.0663435) q[0];
rz(0.66197239) q[2];
sx q[2];
rz(-1.5675002) q[2];
sx q[2];
rz(-0.73288435) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.13818422) q[1];
sx q[1];
rz(-0.78552442) q[1];
sx q[1];
rz(-1.6351321) q[1];
x q[2];
rz(-2.6177789) q[3];
sx q[3];
rz(-0.18693811) q[3];
sx q[3];
rz(-0.65526774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88087624) q[2];
sx q[2];
rz(-1.2141576) q[2];
sx q[2];
rz(2.1991849) q[2];
rz(2.2830394) q[3];
sx q[3];
rz(-2.5311354) q[3];
sx q[3];
rz(0.71949351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(3.0829561) q[0];
rz(3.0851641) q[1];
sx q[1];
rz(-1.3899048) q[1];
sx q[1];
rz(-1.1172392) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53848828) q[0];
sx q[0];
rz(-1.9637917) q[0];
sx q[0];
rz(-0.80727838) q[0];
rz(-1.2710167) q[2];
sx q[2];
rz(-1.4477056) q[2];
sx q[2];
rz(1.8649776) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5176425) q[1];
sx q[1];
rz(-2.8887667) q[1];
sx q[1];
rz(-1.9060017) q[1];
x q[2];
rz(-1.4499217) q[3];
sx q[3];
rz(-1.0981005) q[3];
sx q[3];
rz(-0.25388381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3650018) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(-0.4064202) q[2];
rz(-2.8255919) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(-0.78280226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2773748) q[0];
sx q[0];
rz(-0.45806956) q[0];
sx q[0];
rz(1.3360485) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.3986992) q[1];
sx q[1];
rz(-2.1727402) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1135546) q[0];
sx q[0];
rz(-2.4062706) q[0];
sx q[0];
rz(0.57082973) q[0];
rz(-2.5843262) q[2];
sx q[2];
rz(-0.87199482) q[2];
sx q[2];
rz(3.0521309) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3826461) q[1];
sx q[1];
rz(-2.8422656) q[1];
sx q[1];
rz(-1.4274197) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72520599) q[3];
sx q[3];
rz(-1.7659597) q[3];
sx q[3];
rz(2.8233087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1442673) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(-2.9664795) q[2];
rz(1.3052321) q[3];
sx q[3];
rz(-2.1988018) q[3];
sx q[3];
rz(-1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99821943) q[0];
sx q[0];
rz(-2.0046736) q[0];
sx q[0];
rz(1.0293707) q[0];
rz(-2.3221305) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(2.3447461) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734755) q[0];
sx q[0];
rz(-1.1321018) q[0];
sx q[0];
rz(1.7832157) q[0];
rz(-0.0060002319) q[2];
sx q[2];
rz(-1.0428535) q[2];
sx q[2];
rz(2.1354298) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5863492) q[1];
sx q[1];
rz(-1.2995016) q[1];
sx q[1];
rz(2.9284412) q[1];
rz(-2.0011033) q[3];
sx q[3];
rz(-2.1874551) q[3];
sx q[3];
rz(-2.5466621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13021079) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(1.2288564) q[2];
rz(-2.7025488) q[3];
sx q[3];
rz(-1.565758) q[3];
sx q[3];
rz(-1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(2.3332692) q[0];
rz(1.3574379) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(0.70560169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1418512) q[0];
sx q[0];
rz(-1.5361642) q[0];
sx q[0];
rz(-0.3006199) q[0];
x q[1];
rz(-1.2998492) q[2];
sx q[2];
rz(-2.4204905) q[2];
sx q[2];
rz(1.1972028) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8568154) q[1];
sx q[1];
rz(-1.5681303) q[1];
sx q[1];
rz(1.3176012) q[1];
rz(-pi) q[2];
rz(-0.6198778) q[3];
sx q[3];
rz(-2.2408463) q[3];
sx q[3];
rz(-2.1170962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48050532) q[2];
sx q[2];
rz(-1.2837807) q[2];
sx q[2];
rz(-2.535848) q[2];
rz(3.0211966) q[3];
sx q[3];
rz(-1.2984637) q[3];
sx q[3];
rz(2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8814988) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(-0.0023728097) q[0];
rz(2.782605) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(-2.733309) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4160894) q[0];
sx q[0];
rz(-1.444284) q[0];
sx q[0];
rz(0.036065412) q[0];
x q[1];
rz(0.26770182) q[2];
sx q[2];
rz(-0.98065871) q[2];
sx q[2];
rz(-0.4877643) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6039415) q[1];
sx q[1];
rz(-1.3408325) q[1];
sx q[1];
rz(-1.5020788) q[1];
x q[2];
rz(2.1903135) q[3];
sx q[3];
rz(-0.86182098) q[3];
sx q[3];
rz(-3.0326642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6458873) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(-3.0118946) q[2];
rz(-0.4194704) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(0.81370846) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7523338) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(0.68354052) q[0];
rz(0.93841249) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(1.9155115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83663229) q[0];
sx q[0];
rz(-1.2042521) q[0];
sx q[0];
rz(-2.5185597) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8713538) q[2];
sx q[2];
rz(-2.5242189) q[2];
sx q[2];
rz(0.63948217) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6437373) q[1];
sx q[1];
rz(-1.2853936) q[1];
sx q[1];
rz(2.0394112) q[1];
x q[2];
rz(0.88507401) q[3];
sx q[3];
rz(-2.367691) q[3];
sx q[3];
rz(-2.2507071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.73416) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(0.70719353) q[2];
rz(1.6678984) q[3];
sx q[3];
rz(-0.67631045) q[3];
sx q[3];
rz(-1.428712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4555175) q[0];
sx q[0];
rz(-0.34611836) q[0];
sx q[0];
rz(-0.92380512) q[0];
rz(-1.0651945) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(2.0069897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2475928) q[0];
sx q[0];
rz(-1.1957279) q[0];
sx q[0];
rz(-3.0984237) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85757014) q[2];
sx q[2];
rz(-2.3152707) q[2];
sx q[2];
rz(-0.37304953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75443134) q[1];
sx q[1];
rz(-2.608077) q[1];
sx q[1];
rz(1.6759765) q[1];
rz(-pi) q[2];
rz(1.9530961) q[3];
sx q[3];
rz(-2.144882) q[3];
sx q[3];
rz(-2.7758383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.137546) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(-2.498632) q[2];
rz(-0.17063394) q[3];
sx q[3];
rz(-0.65244397) q[3];
sx q[3];
rz(-2.046106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3951025) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(-0.80628959) q[0];
rz(3.131648) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(-2.8795805) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4777139) q[0];
sx q[0];
rz(-1.0482402) q[0];
sx q[0];
rz(0.19139238) q[0];
x q[1];
rz(2.4960619) q[2];
sx q[2];
rz(-1.4420266) q[2];
sx q[2];
rz(-1.9632531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4993742) q[1];
sx q[1];
rz(-2.0009577) q[1];
sx q[1];
rz(-0.22689928) q[1];
rz(-pi) q[2];
rz(-2.1266009) q[3];
sx q[3];
rz(-0.57325809) q[3];
sx q[3];
rz(-2.3044293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(-0.14287512) q[2];
rz(-2.2376132) q[3];
sx q[3];
rz(-1.6429792) q[3];
sx q[3];
rz(-0.8814632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7664465) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(-2.4850856) q[0];
rz(-0.87908602) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(-0.62969977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64163369) q[0];
sx q[0];
rz(-2.2370506) q[0];
sx q[0];
rz(-2.0883661) q[0];
rz(-1.0378605) q[2];
sx q[2];
rz(-1.3383972) q[2];
sx q[2];
rz(-2.971422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9650642) q[1];
sx q[1];
rz(-1.5429436) q[1];
sx q[1];
rz(-1.5256397) q[1];
x q[2];
rz(1.042243) q[3];
sx q[3];
rz(-1.3011329) q[3];
sx q[3];
rz(0.93664133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62275824) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(0.68010124) q[2];
rz(0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(-2.0228001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685365) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(-0.14280351) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(-1.8129195) q[2];
sx q[2];
rz(-1.6264781) q[2];
sx q[2];
rz(-1.3774058) q[2];
rz(-0.96137267) q[3];
sx q[3];
rz(-1.6037887) q[3];
sx q[3];
rz(-2.9870839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

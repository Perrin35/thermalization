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
rz(-2.4102416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36424822) q[0];
sx q[0];
rz(-2.4349182) q[0];
sx q[0];
rz(1.0663435) q[0];
rz(-3.1362304) q[2];
sx q[2];
rz(-0.66197936) q[2];
sx q[2];
rz(0.84214166) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22905502) q[1];
sx q[1];
rz(-0.78733912) q[1];
sx q[1];
rz(0.064219193) q[1];
rz(-pi) q[2];
rz(-0.16234397) q[3];
sx q[3];
rz(-1.663891) q[3];
sx q[3];
rz(-2.7422991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.88087624) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-0.94240776) q[2];
rz(2.2830394) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(-0.71949351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(-0.058636531) q[0];
rz(-3.0851641) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(-1.1172392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6031044) q[0];
sx q[0];
rz(-1.9637917) q[0];
sx q[0];
rz(-0.80727838) q[0];
rz(-pi) q[1];
rz(1.2710167) q[2];
sx q[2];
rz(-1.4477056) q[2];
sx q[2];
rz(-1.8649776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2721418) q[1];
sx q[1];
rz(-1.488416) q[1];
sx q[1];
rz(1.3314962) q[1];
rz(-0.23157236) q[3];
sx q[3];
rz(-2.6548214) q[3];
sx q[3];
rz(0.51460217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3650018) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(-0.4064202) q[2];
rz(2.8255919) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(-2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2773748) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(1.8055441) q[0];
rz(-2.863073) q[1];
sx q[1];
rz(-1.3986992) q[1];
sx q[1];
rz(2.1727402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8271885) q[0];
sx q[0];
rz(-2.1705856) q[0];
sx q[0];
rz(1.1161854) q[0];
rz(-pi) q[1];
rz(-2.1325705) q[2];
sx q[2];
rz(-2.2780905) q[2];
sx q[2];
rz(0.85894859) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67476701) q[1];
sx q[1];
rz(-1.5286501) q[1];
sx q[1];
rz(1.274363) q[1];
x q[2];
rz(-2.4163867) q[3];
sx q[3];
rz(-1.7659597) q[3];
sx q[3];
rz(2.8233087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99732533) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(2.9664795) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99821943) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(1.0293707) q[0];
rz(2.3221305) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(0.79684657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5433079) q[0];
sx q[0];
rz(-2.6572022) q[0];
sx q[0];
rz(-0.422307) q[0];
x q[1];
rz(3.1355924) q[2];
sx q[2];
rz(-2.0987392) q[2];
sx q[2];
rz(1.0061629) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5863492) q[1];
sx q[1];
rz(-1.2995016) q[1];
sx q[1];
rz(0.21315141) q[1];
rz(-1.1404893) q[3];
sx q[3];
rz(-2.1874551) q[3];
sx q[3];
rz(-0.59493055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13021079) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(-1.9127362) q[2];
rz(2.7025488) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(-1.9108093) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6586886) q[0];
sx q[0];
rz(-1.5331601) q[0];
sx q[0];
rz(-2.3332692) q[0];
rz(-1.7841548) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(0.70560169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1418512) q[0];
sx q[0];
rz(-1.6054285) q[0];
sx q[0];
rz(-0.3006199) q[0];
x q[1];
rz(0.23106261) q[2];
sx q[2];
rz(-0.88132826) q[2];
sx q[2];
rz(-1.5901515) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8452705) q[1];
sx q[1];
rz(-0.25320881) q[1];
sx q[1];
rz(-1.5814387) q[1];
rz(-2.5217149) q[3];
sx q[3];
rz(-0.90074632) q[3];
sx q[3];
rz(1.0244964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6610873) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(2.535848) q[2];
rz(-0.12039603) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(-2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8814988) q[0];
sx q[0];
rz(-0.81387481) q[0];
sx q[0];
rz(-3.1392198) q[0];
rz(2.782605) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(0.40828362) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15015442) q[0];
sx q[0];
rz(-1.6065734) q[0];
sx q[0];
rz(-1.6973901) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8738908) q[2];
sx q[2];
rz(-0.98065871) q[2];
sx q[2];
rz(0.4877643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5376512) q[1];
sx q[1];
rz(-1.8007601) q[1];
sx q[1];
rz(-1.5020788) q[1];
rz(-pi) q[2];
rz(-0.8114641) q[3];
sx q[3];
rz(-2.0272019) q[3];
sx q[3];
rz(1.8965669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4957054) q[2];
sx q[2];
rz(-2.3453823) q[2];
sx q[2];
rz(-3.0118946) q[2];
rz(-2.7221223) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(2.3278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3892589) q[0];
sx q[0];
rz(-2.3663754) q[0];
sx q[0];
rz(2.4580521) q[0];
rz(2.2031802) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(-1.9155115) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968564) q[0];
sx q[0];
rz(-2.4312907) q[0];
sx q[0];
rz(-0.5819178) q[0];
rz(-pi) q[1];
rz(2.8713538) q[2];
sx q[2];
rz(-2.5242189) q[2];
sx q[2];
rz(-0.63948217) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9270806) q[1];
sx q[1];
rz(-2.0190372) q[1];
sx q[1];
rz(2.8238676) q[1];
rz(-pi) q[2];
rz(-0.92323233) q[3];
sx q[3];
rz(-1.1123163) q[3];
sx q[3];
rz(1.9322559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40743264) q[2];
sx q[2];
rz(-1.4340883) q[2];
sx q[2];
rz(-0.70719353) q[2];
rz(1.6678984) q[3];
sx q[3];
rz(-2.4652822) q[3];
sx q[3];
rz(1.428712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68607512) q[0];
sx q[0];
rz(-0.34611836) q[0];
sx q[0];
rz(2.2177875) q[0];
rz(-1.0651945) q[1];
sx q[1];
rz(-1.9520452) q[1];
sx q[1];
rz(-2.0069897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4489733) q[0];
sx q[0];
rz(-1.53063) q[0];
sx q[0];
rz(1.9461826) q[0];
x q[1];
rz(-0.85757014) q[2];
sx q[2];
rz(-2.3152707) q[2];
sx q[2];
rz(-2.7685431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3871613) q[1];
sx q[1];
rz(-0.5335156) q[1];
sx q[1];
rz(-1.6759765) q[1];
rz(-1.1884965) q[3];
sx q[3];
rz(-2.144882) q[3];
sx q[3];
rz(-2.7758383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0040466641) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(-2.498632) q[2];
rz(2.9709587) q[3];
sx q[3];
rz(-0.65244397) q[3];
sx q[3];
rz(-2.046106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3951025) q[0];
sx q[0];
rz(-0.13245067) q[0];
sx q[0];
rz(-0.80628959) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-0.6178304) q[1];
sx q[1];
rz(-0.26201216) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1073955) q[0];
sx q[0];
rz(-2.5881564) q[0];
sx q[0];
rz(1.2518) q[0];
x q[1];
rz(0.64553078) q[2];
sx q[2];
rz(-1.699566) q[2];
sx q[2];
rz(-1.9632531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6422185) q[1];
sx q[1];
rz(-2.0009577) q[1];
sx q[1];
rz(0.22689928) q[1];
x q[2];
rz(-1.069182) q[3];
sx q[3];
rz(-1.2805689) q[3];
sx q[3];
rz(-2.888916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3917196) q[2];
sx q[2];
rz(-0.55730692) q[2];
sx q[2];
rz(0.14287512) q[2];
rz(2.2376132) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(-0.8814632) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7664465) q[0];
sx q[0];
rz(-1.9117993) q[0];
sx q[0];
rz(0.65650702) q[0];
rz(-2.2625066) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(2.5118929) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7555576) q[0];
sx q[0];
rz(-2.3229001) q[0];
sx q[0];
rz(-2.5798803) q[0];
x q[1];
rz(2.8734287) q[2];
sx q[2];
rz(-1.0536453) q[2];
sx q[2];
rz(1.6059396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39552626) q[1];
sx q[1];
rz(-1.5256572) q[1];
sx q[1];
rz(-3.1137115) q[1];
rz(-1.0694169) q[3];
sx q[3];
rz(-2.5540941) q[3];
sx q[3];
rz(2.0794433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62275824) q[2];
sx q[2];
rz(-2.3164985) q[2];
sx q[2];
rz(-2.4614914) q[2];
rz(-2.425219) q[3];
sx q[3];
rz(-3.0383737) q[3];
sx q[3];
rz(-1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685365) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(-0.14280351) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(-0.057351107) q[2];
sx q[2];
rz(-1.3290559) q[2];
sx q[2];
rz(0.17964687) q[2];
rz(-1.6283926) q[3];
sx q[3];
rz(-0.61020281) q[3];
sx q[3];
rz(-1.4635066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2375862) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(1.0884652) q[0];
rz(-pi) q[1];
rz(-1.4346052) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.9905123) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99416713) q[1];
sx q[1];
rz(-2.6174183) q[1];
sx q[1];
rz(-1.0679354) q[1];
x q[2];
rz(2.0099785) q[3];
sx q[3];
rz(-1.7712799) q[3];
sx q[3];
rz(-2.2905614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.819954) q[2];
sx q[2];
rz(-1.7724089) q[2];
sx q[2];
rz(-2.3036172) q[2];
rz(-0.49301246) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74719602) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(-1.863742) q[0];
rz(-2.9648119) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(2.7094254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5837001) q[0];
sx q[0];
rz(-1.0520792) q[0];
sx q[0];
rz(1.8899263) q[0];
rz(-pi) q[1];
rz(-2.3199109) q[2];
sx q[2];
rz(-2.5857918) q[2];
sx q[2];
rz(1.0210212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.00694) q[1];
sx q[1];
rz(-1.4538527) q[1];
sx q[1];
rz(2.0140531) q[1];
rz(0.024859419) q[3];
sx q[3];
rz(-1.0059352) q[3];
sx q[3];
rz(-2.3853175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5923578) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(1.3782079) q[3];
sx q[3];
rz(-1.8816201) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(0.41734636) q[0];
rz(-1.4886645) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(-2.6352077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15784141) q[0];
sx q[0];
rz(-0.39641532) q[0];
sx q[0];
rz(-2.0435964) q[0];
rz(-pi) q[1];
rz(0.26489139) q[2];
sx q[2];
rz(-0.75220097) q[2];
sx q[2];
rz(-2.0843992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80458927) q[1];
sx q[1];
rz(-0.5491921) q[1];
sx q[1];
rz(2.5227491) q[1];
x q[2];
rz(1.9403463) q[3];
sx q[3];
rz(-1.4003716) q[3];
sx q[3];
rz(2.2263262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4425519) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(2.55012) q[2];
rz(0.58602035) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(-1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(2.3024094) q[0];
rz(0.025578586) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2617333) q[0];
sx q[0];
rz(-1.121959) q[0];
sx q[0];
rz(-3.0405322) q[0];
x q[1];
rz(2.2976539) q[2];
sx q[2];
rz(-2.2288449) q[2];
sx q[2];
rz(-2.7968) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38494426) q[1];
sx q[1];
rz(-1.1886016) q[1];
sx q[1];
rz(2.3734943) q[1];
rz(2.7110093) q[3];
sx q[3];
rz(-1.3123371) q[3];
sx q[3];
rz(2.5653198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-3.0419066) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(0.25594041) q[0];
rz(-2.6804965) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(-0.76006132) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5417267) q[0];
sx q[0];
rz(-1.4540298) q[0];
sx q[0];
rz(1.4670502) q[0];
rz(-pi) q[1];
rz(-0.063644479) q[2];
sx q[2];
rz(-1.1568937) q[2];
sx q[2];
rz(2.6457583) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6872245) q[1];
sx q[1];
rz(-1.6974653) q[1];
sx q[1];
rz(-1.6997937) q[1];
x q[2];
rz(-2.7759336) q[3];
sx q[3];
rz(-2.1846111) q[3];
sx q[3];
rz(-2.1706276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5710859) q[2];
sx q[2];
rz(-1.3723624) q[2];
sx q[2];
rz(-2.467353) q[2];
rz(0.21480602) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53133416) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(2.0746322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2444006) q[0];
sx q[0];
rz(-0.025248993) q[0];
sx q[0];
rz(0.61141725) q[0];
rz(-2.3382171) q[2];
sx q[2];
rz(-2.5640045) q[2];
sx q[2];
rz(-0.20882777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3031591) q[1];
sx q[1];
rz(-1.5232956) q[1];
sx q[1];
rz(-0.82364239) q[1];
rz(1.4076036) q[3];
sx q[3];
rz(-0.7810775) q[3];
sx q[3];
rz(1.2490602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.431488) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(2.5816494) q[2];
rz(-0.7263178) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(2.9395318) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(2.2834159) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35408033) q[0];
sx q[0];
rz(-1.6991827) q[0];
sx q[0];
rz(1.4228729) q[0];
rz(-pi) q[1];
rz(0.045072149) q[2];
sx q[2];
rz(-1.9809234) q[2];
sx q[2];
rz(-2.8571667) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.757829) q[1];
sx q[1];
rz(-1.5997636) q[1];
sx q[1];
rz(-1.5860735) q[1];
rz(2.5704727) q[3];
sx q[3];
rz(-1.2763599) q[3];
sx q[3];
rz(-1.4049698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28785607) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(1.3254335) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(0.37471399) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(-1.7920866) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23096202) q[0];
sx q[0];
rz(-0.9696784) q[0];
sx q[0];
rz(1.0966976) q[0];
rz(-1.5937514) q[2];
sx q[2];
rz(-0.85283961) q[2];
sx q[2];
rz(0.1644451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.046557758) q[1];
sx q[1];
rz(-1.7525502) q[1];
sx q[1];
rz(-0.05038105) q[1];
rz(-1.8324864) q[3];
sx q[3];
rz(-1.040254) q[3];
sx q[3];
rz(1.4025276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(-0.17503861) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(-1.5375686) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984734) q[0];
sx q[0];
rz(-2.7612918) q[0];
sx q[0];
rz(0.93912504) q[0];
rz(-pi) q[1];
rz(-0.23065718) q[2];
sx q[2];
rz(-0.42867491) q[2];
sx q[2];
rz(0.32282695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4496085) q[1];
sx q[1];
rz(-2.2712628) q[1];
sx q[1];
rz(1.7193754) q[1];
rz(-pi) q[2];
rz(-0.36112862) q[3];
sx q[3];
rz(-1.1361406) q[3];
sx q[3];
rz(1.9643195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(0.76134479) q[2];
rz(-0.90041655) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(-3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(0.21690579) q[0];
rz(-0.63198173) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023708658) q[0];
sx q[0];
rz(-0.81239359) q[0];
sx q[0];
rz(-2.6931767) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5920792) q[2];
sx q[2];
rz(-1.9351442) q[2];
sx q[2];
rz(-3.0669616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54143822) q[1];
sx q[1];
rz(-2.6791875) q[1];
sx q[1];
rz(0.17141436) q[1];
rz(-pi) q[2];
rz(0.44317742) q[3];
sx q[3];
rz(-0.9842397) q[3];
sx q[3];
rz(0.12936684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0163991) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(2.1968502) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(-0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5007297) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-0.8846994) q[1];
sx q[1];
rz(-2.2365166) q[1];
sx q[1];
rz(2.8832163) q[1];
rz(2.5169218) q[2];
sx q[2];
rz(-0.93745898) q[2];
sx q[2];
rz(-2.6425101) q[2];
rz(2.9363019) q[3];
sx q[3];
rz(-1.1355573) q[3];
sx q[3];
rz(-0.60187403) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
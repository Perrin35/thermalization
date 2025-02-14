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
rz(-0.1470546) q[0];
sx q[0];
rz(-1.4015863) q[0];
sx q[0];
rz(0.92010486) q[0];
rz(3.0609581) q[1];
sx q[1];
rz(-0.57890761) q[1];
sx q[1];
rz(-0.97715598) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1950705) q[0];
sx q[0];
rz(-1.3597288) q[0];
sx q[0];
rz(2.7618183) q[0];
rz(-pi) q[1];
rz(0.13222569) q[2];
sx q[2];
rz(-1.5572845) q[2];
sx q[2];
rz(2.9861762) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5861533) q[1];
sx q[1];
rz(-1.4358913) q[1];
sx q[1];
rz(-0.6026661) q[1];
rz(-pi) q[2];
rz(-1.8771873) q[3];
sx q[3];
rz(-0.9054817) q[3];
sx q[3];
rz(0.09395919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-2.5207632) q[2];
rz(0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(-0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(-2.5625693) q[0];
rz(-0.82194263) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(0.45375219) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4849783) q[0];
sx q[0];
rz(-1.7957644) q[0];
sx q[0];
rz(-0.63240914) q[0];
rz(2.6751509) q[2];
sx q[2];
rz(-1.6714408) q[2];
sx q[2];
rz(-0.60324861) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0099185506) q[1];
sx q[1];
rz(-1.6955396) q[1];
sx q[1];
rz(-2.591955) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9048432) q[3];
sx q[3];
rz(-0.65757759) q[3];
sx q[3];
rz(2.6386054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4072676) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(2.5767051) q[2];
rz(-2.7867553) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1693717) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(-2.5900904) q[0];
rz(0.36477271) q[1];
sx q[1];
rz(-2.8021937) q[1];
sx q[1];
rz(0.50484467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14061117) q[0];
sx q[0];
rz(-1.7303082) q[0];
sx q[0];
rz(-2.7310452) q[0];
x q[1];
rz(-1.4771492) q[2];
sx q[2];
rz(-1.159707) q[2];
sx q[2];
rz(1.1716441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80124679) q[1];
sx q[1];
rz(-0.97724229) q[1];
sx q[1];
rz(2.1562063) q[1];
rz(-1.9471719) q[3];
sx q[3];
rz(-2.0224114) q[3];
sx q[3];
rz(0.87290369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-0.047253963) q[2];
rz(2.3047678) q[3];
sx q[3];
rz(-2.4182726) q[3];
sx q[3];
rz(-2.1164472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4243917) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(1.3264054) q[0];
rz(1.6560417) q[1];
sx q[1];
rz(-2.0573719) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0795006) q[0];
sx q[0];
rz(-1.5695111) q[0];
sx q[0];
rz(1.5527524) q[0];
x q[1];
rz(-0.97135404) q[2];
sx q[2];
rz(-0.58341372) q[2];
sx q[2];
rz(-2.3177878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30764515) q[1];
sx q[1];
rz(-1.2356504) q[1];
sx q[1];
rz(-0.076471523) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67948273) q[3];
sx q[3];
rz(-0.75089291) q[3];
sx q[3];
rz(-0.39284947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9012458) q[2];
sx q[2];
rz(-2.2525658) q[2];
sx q[2];
rz(-1.7690313) q[2];
rz(1.9339804) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(2.8670368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8127301) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(0.10511705) q[0];
rz(0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-2.0786659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8369097) q[0];
sx q[0];
rz(-1.5567008) q[0];
sx q[0];
rz(2.430116) q[0];
rz(-pi) q[1];
rz(2.8811997) q[2];
sx q[2];
rz(-0.77885926) q[2];
sx q[2];
rz(0.078140251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2984933) q[1];
sx q[1];
rz(-1.5122422) q[1];
sx q[1];
rz(1.5850409) q[1];
rz(-pi) q[2];
rz(-0.0023554597) q[3];
sx q[3];
rz(-1.071953) q[3];
sx q[3];
rz(-0.82060196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.062332705) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(1.3366535) q[2];
rz(0.054232728) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(-0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17276758) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(-1.927595) q[0];
rz(2.0955775) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(-0.85375839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19148286) q[0];
sx q[0];
rz(-1.5111418) q[0];
sx q[0];
rz(-1.4364987) q[0];
rz(-pi) q[1];
rz(-0.75550236) q[2];
sx q[2];
rz(-1.1248338) q[2];
sx q[2];
rz(-1.8309636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5750372) q[1];
sx q[1];
rz(-2.2872529) q[1];
sx q[1];
rz(-0.8030007) q[1];
rz(-pi) q[2];
rz(0.3993897) q[3];
sx q[3];
rz(-2.7654184) q[3];
sx q[3];
rz(-1.1272421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0115016) q[2];
sx q[2];
rz(-1.4902196) q[2];
sx q[2];
rz(-2.3940274) q[2];
rz(-0.93303624) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(-2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0316684) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(2.5001496) q[0];
rz(-2.172442) q[1];
sx q[1];
rz(-2.0717924) q[1];
sx q[1];
rz(-1.1136805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9819558) q[0];
sx q[0];
rz(-2.7244719) q[0];
sx q[0];
rz(1.3288767) q[0];
x q[1];
rz(-0.53345726) q[2];
sx q[2];
rz(-1.0048303) q[2];
sx q[2];
rz(0.12803687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2526452) q[1];
sx q[1];
rz(-1.749649) q[1];
sx q[1];
rz(-2.0986544) q[1];
rz(0.05984743) q[3];
sx q[3];
rz(-2.6204797) q[3];
sx q[3];
rz(1.1731565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.574719) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(2.3835772) q[2];
rz(-2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(-2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4161943) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(-2.388227) q[0];
rz(0.92195177) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(-3.0070378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3938569) q[0];
sx q[0];
rz(-1.6574291) q[0];
sx q[0];
rz(-1.2133693) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3165264) q[2];
sx q[2];
rz(-0.94506028) q[2];
sx q[2];
rz(2.3601687) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53271097) q[1];
sx q[1];
rz(-1.5517762) q[1];
sx q[1];
rz(3.0518107) q[1];
rz(-pi) q[2];
rz(-2.4681925) q[3];
sx q[3];
rz(-1.9385425) q[3];
sx q[3];
rz(1.0750225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3677463) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(-1.9473677) q[2];
rz(1.1456683) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(1.0660508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0399748) q[0];
sx q[0];
rz(-2.6823253) q[0];
sx q[0];
rz(-0.38247821) q[0];
rz(0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(1.8006178) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1857698) q[0];
sx q[0];
rz(-0.18322028) q[0];
sx q[0];
rz(-2.630156) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4828311) q[2];
sx q[2];
rz(-2.3260657) q[2];
sx q[2];
rz(1.2620827) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.080440259) q[1];
sx q[1];
rz(-2.3576479) q[1];
sx q[1];
rz(-2.817201) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3519456) q[3];
sx q[3];
rz(-2.7776981) q[3];
sx q[3];
rz(-0.6536676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0535023) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(3.0255393) q[2];
rz(1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5905404) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(2.9569448) q[0];
rz(-0.91839904) q[1];
sx q[1];
rz(-1.5946439) q[1];
sx q[1];
rz(-1.7041697) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0032063403) q[0];
sx q[0];
rz(-2.0420364) q[0];
sx q[0];
rz(2.9127761) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.499648) q[2];
sx q[2];
rz(-2.2916823) q[2];
sx q[2];
rz(-2.3325553) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.559222) q[1];
sx q[1];
rz(-2.5032024) q[1];
sx q[1];
rz(1.8965782) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40664169) q[3];
sx q[3];
rz(-2.4235117) q[3];
sx q[3];
rz(-2.6585916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.664428) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(2.3264558) q[3];
sx q[3];
rz(-0.4709979) q[3];
sx q[3];
rz(2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(0.78085113) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(1.9736171) q[2];
sx q[2];
rz(-2.8388173) q[2];
sx q[2];
rz(-1.7025316) q[2];
rz(1.5315957) q[3];
sx q[3];
rz(-1.8146252) q[3];
sx q[3];
rz(-1.6011325) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

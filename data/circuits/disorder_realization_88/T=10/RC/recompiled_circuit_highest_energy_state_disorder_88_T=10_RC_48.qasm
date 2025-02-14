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
rz(4.881599) q[0];
sx q[0];
rz(10.344883) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(-2.562685) q[1];
sx q[1];
rz(0.97715598) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8590737) q[0];
sx q[0];
rz(-0.4319829) q[0];
sx q[0];
rz(0.52406128) q[0];
rz(-pi) q[1];
rz(0.10213587) q[2];
sx q[2];
rz(-3.0086824) q[2];
sx q[2];
rz(-1.5166211) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5861533) q[1];
sx q[1];
rz(-1.7057014) q[1];
sx q[1];
rz(-2.5389266) q[1];
rz(-pi) q[2];
rz(-2.7746088) q[3];
sx q[3];
rz(-2.4189848) q[3];
sx q[3];
rz(2.5740576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-2.5207632) q[2];
rz(-0.51554716) q[3];
sx q[3];
rz(-1.7190869) q[3];
sx q[3];
rz(-0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8949378) q[0];
sx q[0];
rz(-1.6064914) q[0];
sx q[0];
rz(-2.5625693) q[0];
rz(-0.82194263) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(-2.6878405) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61868405) q[0];
sx q[0];
rz(-0.6660307) q[0];
sx q[0];
rz(-2.7722107) q[0];
rz(-0.46644174) q[2];
sx q[2];
rz(-1.6714408) q[2];
sx q[2];
rz(2.538344) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0099185506) q[1];
sx q[1];
rz(-1.6955396) q[1];
sx q[1];
rz(0.54963767) q[1];
rz(-pi) q[2];
rz(0.9048432) q[3];
sx q[3];
rz(-0.65757759) q[3];
sx q[3];
rz(0.50298728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.73432505) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(-0.56488758) q[2];
rz(2.7867553) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.7179276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722209) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(-2.5900904) q[0];
rz(2.7768199) q[1];
sx q[1];
rz(-2.8021937) q[1];
sx q[1];
rz(2.636748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6423787) q[0];
sx q[0];
rz(-1.9758245) q[0];
sx q[0];
rz(1.7444872) q[0];
rz(-pi) q[1];
rz(1.4771492) q[2];
sx q[2];
rz(-1.159707) q[2];
sx q[2];
rz(1.9699485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4703731) q[1];
sx q[1];
rz(-0.80802901) q[1];
sx q[1];
rz(-2.4554182) q[1];
rz(-pi) q[2];
rz(1.9471719) q[3];
sx q[3];
rz(-2.0224114) q[3];
sx q[3];
rz(2.268689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-0.047253963) q[2];
rz(-2.3047678) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(1.0251454) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4243917) q[0];
sx q[0];
rz(-2.1121139) q[0];
sx q[0];
rz(-1.3264054) q[0];
rz(1.4855509) q[1];
sx q[1];
rz(-2.0573719) q[1];
sx q[1];
rz(-2.5066689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50872749) q[0];
sx q[0];
rz(-1.5527525) q[0];
sx q[0];
rz(0.0012854171) q[0];
x q[1];
rz(0.35648326) q[2];
sx q[2];
rz(-2.0429789) q[2];
sx q[2];
rz(-3.0038578) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0628218) q[1];
sx q[1];
rz(-0.34343849) q[1];
sx q[1];
rz(1.354864) q[1];
rz(-pi) q[2];
rz(2.4621099) q[3];
sx q[3];
rz(-2.3906997) q[3];
sx q[3];
rz(0.39284947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9012458) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(1.7690313) q[2];
rz(-1.9339804) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8127301) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(-3.0364756) q[0];
rz(-0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-1.0629268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369097) q[0];
sx q[0];
rz(-1.5567008) q[0];
sx q[0];
rz(0.71147664) q[0];
rz(-pi) q[1];
rz(-2.8811997) q[2];
sx q[2];
rz(-2.3627334) q[2];
sx q[2];
rz(0.078140251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2731367) q[1];
sx q[1];
rz(-1.5565762) q[1];
sx q[1];
rz(0.058560024) q[1];
rz(-1.0719518) q[3];
sx q[3];
rz(-1.5728647) q[3];
sx q[3];
rz(2.3902714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.062332705) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(-1.3366535) q[2];
rz(3.0873599) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.17276758) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(-1.927595) q[0];
rz(-2.0955775) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(-2.2878343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3873685) q[0];
sx q[0];
rz(-1.7048536) q[0];
sx q[0];
rz(-0.06019528) q[0];
rz(-pi) q[1];
rz(-2.1519203) q[2];
sx q[2];
rz(-0.90384353) q[2];
sx q[2];
rz(3.0158531) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5665555) q[1];
sx q[1];
rz(-0.85433975) q[1];
sx q[1];
rz(-2.3385919) q[1];
x q[2];
rz(0.34900174) q[3];
sx q[3];
rz(-1.7141388) q[3];
sx q[3];
rz(0.81763148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0115016) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(2.3940274) q[2];
rz(-0.93303624) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-0.30195495) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1099243) q[0];
sx q[0];
rz(-2.1878991) q[0];
sx q[0];
rz(0.64144301) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(-1.1136805) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1596368) q[0];
sx q[0];
rz(-0.41712077) q[0];
sx q[0];
rz(-1.812716) q[0];
rz(-pi) q[1];
x q[1];
rz(0.895787) q[2];
sx q[2];
rz(-2.3844344) q[2];
sx q[2];
rz(0.96162187) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1195928) q[1];
sx q[1];
rz(-0.55460677) q[1];
sx q[1];
rz(-1.2261934) q[1];
rz(-pi) q[2];
rz(-3.0817452) q[3];
sx q[3];
rz(-0.52111292) q[3];
sx q[3];
rz(-1.1731565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5668737) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(0.75801545) q[2];
rz(2.8356683) q[3];
sx q[3];
rz(-2.7829792) q[3];
sx q[3];
rz(-2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7253983) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(0.7533657) q[0];
rz(-2.2196409) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(0.13455483) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7477357) q[0];
sx q[0];
rz(-1.6574291) q[0];
sx q[0];
rz(1.9282233) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8250663) q[2];
sx q[2];
rz(-0.94506028) q[2];
sx q[2];
rz(0.78142399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1052195) q[1];
sx q[1];
rz(-1.660562) q[1];
sx q[1];
rz(-1.5516993) q[1];
rz(-pi) q[2];
rz(-0.67340019) q[3];
sx q[3];
rz(-1.2030501) q[3];
sx q[3];
rz(-2.0665702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77384633) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(-1.194225) q[2];
rz(1.1456683) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(-2.0755419) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0399748) q[0];
sx q[0];
rz(-2.6823253) q[0];
sx q[0];
rz(0.38247821) q[0];
rz(0.76599145) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.3409748) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6670973) q[0];
sx q[0];
rz(-1.7303559) q[0];
sx q[0];
rz(1.480353) q[0];
x q[1];
rz(0.75720301) q[2];
sx q[2];
rz(-1.5067889) q[2];
sx q[2];
rz(-2.7724977) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0611524) q[1];
sx q[1];
rz(-0.78394475) q[1];
sx q[1];
rz(-2.817201) q[1];
rz(-pi) q[2];
rz(-0.082499113) q[3];
sx q[3];
rz(-1.9256221) q[3];
sx q[3];
rz(0.88731836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0880903) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(-3.0255393) q[2];
rz(1.8807489) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.457823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(-0.1846479) q[0];
rz(2.2231936) q[1];
sx q[1];
rz(-1.5946439) q[1];
sx q[1];
rz(1.437423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0032063403) q[0];
sx q[0];
rz(-1.0995563) q[0];
sx q[0];
rz(-2.9127761) q[0];
x q[1];
rz(-1.499648) q[2];
sx q[2];
rz(-2.2916823) q[2];
sx q[2];
rz(-2.3325553) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5823707) q[1];
sx q[1];
rz(-2.5032024) q[1];
sx q[1];
rz(1.8965782) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67622306) q[3];
sx q[3];
rz(-1.3075324) q[3];
sx q[3];
rz(0.77419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.664428) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(2.5028382) q[2];
rz(0.81513682) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.12185085) q[2];
sx q[2];
rz(-1.2929299) q[2];
sx q[2];
rz(-1.2826512) q[2];
rz(-0.24400901) q[3];
sx q[3];
rz(-1.5327558) q[3];
sx q[3];
rz(-0.039804606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

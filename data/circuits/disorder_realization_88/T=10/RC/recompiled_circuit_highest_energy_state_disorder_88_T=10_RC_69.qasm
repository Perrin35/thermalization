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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(2.2214878) q[0];
rz(3.0609581) q[1];
sx q[1];
rz(-0.57890761) q[1];
sx q[1];
rz(2.1644367) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1950705) q[0];
sx q[0];
rz(-1.3597288) q[0];
sx q[0];
rz(0.37977438) q[0];
x q[1];
rz(1.5571655) q[2];
sx q[2];
rz(-1.7030099) q[2];
sx q[2];
rz(-1.4135828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5861533) q[1];
sx q[1];
rz(-1.7057014) q[1];
sx q[1];
rz(-2.5389266) q[1];
rz(-pi) q[2];
rz(2.7746088) q[3];
sx q[3];
rz(-2.4189848) q[3];
sx q[3];
rz(0.56753502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5620293) q[2];
sx q[2];
rz(-1.2144438) q[2];
sx q[2];
rz(0.62082949) q[2];
rz(-2.6260455) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(-0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(0.57902336) q[0];
rz(0.82194263) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(-0.45375219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61868405) q[0];
sx q[0];
rz(-2.4755619) q[0];
sx q[0];
rz(0.36938195) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2208925) q[2];
sx q[2];
rz(-2.6651987) q[2];
sx q[2];
rz(1.1645137) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4848043) q[1];
sx q[1];
rz(-2.1156807) q[1];
sx q[1];
rz(1.4247895) q[1];
x q[2];
rz(-0.9048432) q[3];
sx q[3];
rz(-2.4840151) q[3];
sx q[3];
rz(-2.6386054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4072676) q[2];
sx q[2];
rz(-3.0027323) q[2];
sx q[2];
rz(2.5767051) q[2];
rz(2.7867553) q[3];
sx q[3];
rz(-2.1653039) q[3];
sx q[3];
rz(1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1693717) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(2.5900904) q[0];
rz(-2.7768199) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(-0.50484467) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4992139) q[0];
sx q[0];
rz(-1.1657682) q[0];
sx q[0];
rz(1.7444872) q[0];
x q[1];
rz(1.6644434) q[2];
sx q[2];
rz(-1.9818857) q[2];
sx q[2];
rz(1.9699485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3403459) q[1];
sx q[1];
rz(-0.97724229) q[1];
sx q[1];
rz(2.1562063) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9471719) q[3];
sx q[3];
rz(-1.1191812) q[3];
sx q[3];
rz(2.268689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(-3.0943387) q[2];
rz(-2.3047678) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.717201) q[0];
sx q[0];
rz(-1.0294788) q[0];
sx q[0];
rz(-1.3264054) q[0];
rz(1.4855509) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.062092) q[0];
sx q[0];
rz(-1.5695111) q[0];
sx q[0];
rz(-1.5527524) q[0];
rz(-pi) q[1];
rz(-0.35648326) q[2];
sx q[2];
rz(-1.0986137) q[2];
sx q[2];
rz(-3.0038578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30764515) q[1];
sx q[1];
rz(-1.9059423) q[1];
sx q[1];
rz(3.0651211) q[1];
rz(2.1011971) q[3];
sx q[3];
rz(-2.1302855) q[3];
sx q[3];
rz(-0.44236576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2403468) q[2];
sx q[2];
rz(-0.88902688) q[2];
sx q[2];
rz(-1.3725613) q[2];
rz(1.2076123) q[3];
sx q[3];
rz(-0.6438846) q[3];
sx q[3];
rz(-0.27455583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32886252) q[0];
sx q[0];
rz(-2.6565318) q[0];
sx q[0];
rz(0.10511705) q[0];
rz(0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-2.0786659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30468291) q[0];
sx q[0];
rz(-1.5567008) q[0];
sx q[0];
rz(2.430116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.379871) q[2];
sx q[2];
rz(-1.3889379) q[2];
sx q[2];
rz(1.3052502) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2731367) q[1];
sx q[1];
rz(-1.5850164) q[1];
sx q[1];
rz(0.058560024) q[1];
rz(1.5751198) q[3];
sx q[3];
rz(-0.49884847) q[3];
sx q[3];
rz(2.3259142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.062332705) q[2];
sx q[2];
rz(-1.7534813) q[2];
sx q[2];
rz(-1.8049392) q[2];
rz(-3.0873599) q[3];
sx q[3];
rz(-1.1905866) q[3];
sx q[3];
rz(0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9688251) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(1.927595) q[0];
rz(2.0955775) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(2.2878343) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1778422) q[0];
sx q[0];
rz(-0.14687777) q[0];
sx q[0];
rz(-1.9903723) q[0];
x q[1];
rz(2.5327024) q[2];
sx q[2];
rz(-0.85431803) q[2];
sx q[2];
rz(2.4520055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5665555) q[1];
sx q[1];
rz(-0.85433975) q[1];
sx q[1];
rz(-0.8030007) q[1];
rz(1.7231971) q[3];
sx q[3];
rz(-1.2255242) q[3];
sx q[3];
rz(-0.70122805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0115016) q[2];
sx q[2];
rz(-1.4902196) q[2];
sx q[2];
rz(-2.3940274) q[2];
rz(2.2085564) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(2.8396377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1099243) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(-2.5001496) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(-1.1136805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1596368) q[0];
sx q[0];
rz(-0.41712077) q[0];
sx q[0];
rz(1.3288767) q[0];
rz(-0.93514438) q[2];
sx q[2];
rz(-2.0144267) q[2];
sx q[2];
rz(-2.005524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.021999849) q[1];
sx q[1];
rz(-0.55460677) q[1];
sx q[1];
rz(-1.2261934) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0817452) q[3];
sx q[3];
rz(-2.6204797) q[3];
sx q[3];
rz(1.9684362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5668737) q[2];
sx q[2];
rz(-0.69872624) q[2];
sx q[2];
rz(2.3835772) q[2];
rz(2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.4161943) q[0];
sx q[0];
rz(-2.5328126) q[0];
sx q[0];
rz(2.388227) q[0];
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
rz(-3.0908703) q[0];
sx q[0];
rz(-0.36733741) q[0];
sx q[0];
rz(-1.3274756) q[0];
rz(2.5002648) q[2];
sx q[2];
rz(-1.3654764) q[2];
sx q[2];
rz(0.63831282) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0363732) q[1];
sx q[1];
rz(-1.660562) q[1];
sx q[1];
rz(1.5898934) q[1];
rz(2.4681925) q[3];
sx q[3];
rz(-1.9385425) q[3];
sx q[3];
rz(2.0665702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3677463) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(1.194225) q[2];
rz(-1.1456683) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(-2.0755419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1016178) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(2.7591144) q[0];
rz(0.76599145) q[1];
sx q[1];
rz(-1.6440369) q[1];
sx q[1];
rz(1.3409748) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4744953) q[0];
sx q[0];
rz(-1.7303559) q[0];
sx q[0];
rz(1.480353) q[0];
rz(-0.093042298) q[2];
sx q[2];
rz(-0.75936717) q[2];
sx q[2];
rz(-2.0074646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7786918) q[1];
sx q[1];
rz(-2.3039989) q[1];
sx q[1];
rz(1.2630839) q[1];
x q[2];
rz(1.2148592) q[3];
sx q[3];
rz(-1.4934469) q[3];
sx q[3];
rz(2.429395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0880903) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(-3.0255393) q[2];
rz(-1.2608438) q[3];
sx q[3];
rz(-2.0033658) q[3];
sx q[3];
rz(1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(0.1846479) q[0];
rz(-2.2231936) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(1.437423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4686615) q[0];
sx q[0];
rz(-1.7743006) q[0];
sx q[0];
rz(-2.0528021) q[0];
x q[1];
rz(-1.6419446) q[2];
sx q[2];
rz(-0.84991036) q[2];
sx q[2];
rz(0.80903731) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1841764) q[1];
sx q[1];
rz(-0.97089689) q[1];
sx q[1];
rz(-2.9084212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2380794) q[3];
sx q[3];
rz(-2.2196688) q[3];
sx q[3];
rz(1.0024662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47716466) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(2.3264558) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(0.42069978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(2.3607415) q[1];
sx q[1];
rz(-1.4518705) q[1];
sx q[1];
rz(-0.89444583) q[1];
rz(0.12185085) q[2];
sx q[2];
rz(-1.2929299) q[2];
sx q[2];
rz(-1.2826512) q[2];
rz(-2.9853447) q[3];
sx q[3];
rz(-2.8946946) q[3];
sx q[3];
rz(1.3794086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

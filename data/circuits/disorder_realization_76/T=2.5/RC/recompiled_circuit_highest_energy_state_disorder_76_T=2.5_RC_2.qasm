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
rz(0.87585706) q[0];
sx q[0];
rz(2.1729204) q[0];
sx q[0];
rz(7.2090413) q[0];
rz(-2.5246188) q[1];
sx q[1];
rz(-2.4763835) q[1];
sx q[1];
rz(1.8170504) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0449555) q[0];
sx q[0];
rz(-2.5384266) q[0];
sx q[0];
rz(0.11267333) q[0];
x q[1];
rz(-1.6683031) q[2];
sx q[2];
rz(-0.98819369) q[2];
sx q[2];
rz(-2.1511457) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2046537) q[1];
sx q[1];
rz(-2.0273682) q[1];
sx q[1];
rz(0.66557933) q[1];
rz(-pi) q[2];
rz(-0.91574131) q[3];
sx q[3];
rz(-0.91980442) q[3];
sx q[3];
rz(0.30177339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2609743) q[2];
sx q[2];
rz(-1.7558492) q[2];
sx q[2];
rz(-2.3495038) q[2];
rz(-2.2657307) q[3];
sx q[3];
rz(-0.10871092) q[3];
sx q[3];
rz(2.47827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51139128) q[0];
sx q[0];
rz(-1.8237317) q[0];
sx q[0];
rz(-2.200101) q[0];
rz(-2.1425653) q[1];
sx q[1];
rz(-0.92728725) q[1];
sx q[1];
rz(3.0587382) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6951286) q[0];
sx q[0];
rz(-2.3020491) q[0];
sx q[0];
rz(1.7056998) q[0];
x q[1];
rz(0.24225927) q[2];
sx q[2];
rz(-1.2238811) q[2];
sx q[2];
rz(1.0936979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1513843) q[1];
sx q[1];
rz(-0.71077222) q[1];
sx q[1];
rz(0.73020331) q[1];
rz(2.1762455) q[3];
sx q[3];
rz(-2.1248528) q[3];
sx q[3];
rz(0.29747552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89722172) q[2];
sx q[2];
rz(-1.7873849) q[2];
sx q[2];
rz(-1.2574035) q[2];
rz(0.8463549) q[3];
sx q[3];
rz(-0.1592764) q[3];
sx q[3];
rz(2.4634821) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9427247) q[0];
sx q[0];
rz(-1.4582448) q[0];
sx q[0];
rz(0.22920907) q[0];
rz(2.9275059) q[1];
sx q[1];
rz(-0.87132088) q[1];
sx q[1];
rz(-2.7600938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3840528) q[0];
sx q[0];
rz(-2.4017576) q[0];
sx q[0];
rz(-0.53348855) q[0];
rz(-pi) q[1];
rz(1.3573285) q[2];
sx q[2];
rz(-0.74356438) q[2];
sx q[2];
rz(1.5539757) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.72813624) q[1];
sx q[1];
rz(-1.6237139) q[1];
sx q[1];
rz(-0.028832988) q[1];
rz(-pi) q[2];
rz(1.5656785) q[3];
sx q[3];
rz(-1.9178784) q[3];
sx q[3];
rz(-1.0597313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4182338) q[2];
sx q[2];
rz(-2.8853719) q[2];
sx q[2];
rz(2.5733433) q[2];
rz(1.4189643) q[3];
sx q[3];
rz(-1.7122372) q[3];
sx q[3];
rz(2.0645781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028246183) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(-2.8908253) q[0];
rz(0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(1.2006522) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1538382) q[0];
sx q[0];
rz(-2.2306577) q[0];
sx q[0];
rz(0.46100088) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1293639) q[2];
sx q[2];
rz(-1.3462634) q[2];
sx q[2];
rz(-2.6972636) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1657287) q[1];
sx q[1];
rz(-1.4714676) q[1];
sx q[1];
rz(-2.4286859) q[1];
x q[2];
rz(2.0968998) q[3];
sx q[3];
rz(-1.3402437) q[3];
sx q[3];
rz(1.6342722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4667929) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(1.2910845) q[2];
rz(-0.42168266) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(1.5765367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52294937) q[0];
sx q[0];
rz(-0.72294253) q[0];
sx q[0];
rz(1.500754) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5539955) q[1];
sx q[1];
rz(1.1281475) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90445176) q[0];
sx q[0];
rz(-1.6078976) q[0];
sx q[0];
rz(-2.9536329) q[0];
rz(2.7893547) q[2];
sx q[2];
rz(-2.0697429) q[2];
sx q[2];
rz(1.3832472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50666821) q[1];
sx q[1];
rz(-1.4574058) q[1];
sx q[1];
rz(-0.054111295) q[1];
rz(0.87781436) q[3];
sx q[3];
rz(-2.7846309) q[3];
sx q[3];
rz(1.2277384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.03380123) q[2];
sx q[2];
rz(-2.0267603) q[2];
sx q[2];
rz(-2.4647253) q[2];
rz(-2.233861) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(2.7270253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9229729) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(-2.3722017) q[0];
rz(-1.9006624) q[1];
sx q[1];
rz(-0.85378328) q[1];
sx q[1];
rz(0.21533899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32672918) q[0];
sx q[0];
rz(-1.0007326) q[0];
sx q[0];
rz(0.22320052) q[0];
rz(-pi) q[1];
rz(1.5211283) q[2];
sx q[2];
rz(-2.81041) q[2];
sx q[2];
rz(-1.0612203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29500719) q[1];
sx q[1];
rz(-0.93448105) q[1];
sx q[1];
rz(-1.8759439) q[1];
x q[2];
rz(0.39677119) q[3];
sx q[3];
rz(-2.4037529) q[3];
sx q[3];
rz(-0.48860197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26663366) q[2];
sx q[2];
rz(-1.5077488) q[2];
sx q[2];
rz(0.61666644) q[2];
rz(-0.2002317) q[3];
sx q[3];
rz(-0.73035208) q[3];
sx q[3];
rz(-1.1327789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1277593) q[0];
sx q[0];
rz(-2.5795689) q[0];
sx q[0];
rz(-3.1264937) q[0];
rz(-0.12241441) q[1];
sx q[1];
rz(-1.8465123) q[1];
sx q[1];
rz(-1.0282358) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434179) q[0];
sx q[0];
rz(-0.91599303) q[0];
sx q[0];
rz(1.6408987) q[0];
rz(-1.2790658) q[2];
sx q[2];
rz(-2.5030067) q[2];
sx q[2];
rz(1.2995468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0575917) q[1];
sx q[1];
rz(-2.0023953) q[1];
sx q[1];
rz(1.6017385) q[1];
rz(1.3002197) q[3];
sx q[3];
rz(-1.2660053) q[3];
sx q[3];
rz(0.51606015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5369109) q[2];
sx q[2];
rz(-1.9313507) q[2];
sx q[2];
rz(-0.49986419) q[2];
rz(-0.31575051) q[3];
sx q[3];
rz(-2.4707268) q[3];
sx q[3];
rz(2.1226814) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71378088) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(0.94386238) q[0];
rz(-1.4048514) q[1];
sx q[1];
rz(-1.2662788) q[1];
sx q[1];
rz(-2.2199383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93225828) q[0];
sx q[0];
rz(-2.8186322) q[0];
sx q[0];
rz(-2.9684116) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87041847) q[2];
sx q[2];
rz(-2.7590115) q[2];
sx q[2];
rz(2.062881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5068622) q[1];
sx q[1];
rz(-2.1980739) q[1];
sx q[1];
rz(-2.2192713) q[1];
rz(-1.5783903) q[3];
sx q[3];
rz(-0.51285997) q[3];
sx q[3];
rz(-2.9160485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1450119) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(1.5865631) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.8127706) q[3];
sx q[3];
rz(-1.740295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15040511) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(0.64252585) q[0];
rz(1.4379028) q[1];
sx q[1];
rz(-2.3846886) q[1];
sx q[1];
rz(0.14370758) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3966345) q[0];
sx q[0];
rz(-1.2681343) q[0];
sx q[0];
rz(2.55654) q[0];
x q[1];
rz(0.917786) q[2];
sx q[2];
rz(-0.76030234) q[2];
sx q[2];
rz(2.0298634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66389304) q[1];
sx q[1];
rz(-1.8346922) q[1];
sx q[1];
rz(-2.4724835) q[1];
x q[2];
rz(-0.43817839) q[3];
sx q[3];
rz(-1.9219805) q[3];
sx q[3];
rz(2.8170057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5370499) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(0.26270467) q[2];
rz(-0.25775868) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(-2.2885382) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2951374) q[0];
sx q[0];
rz(-1.6520123) q[0];
sx q[0];
rz(1.9211796) q[0];
rz(-2.659761) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-2.2442472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37384826) q[0];
sx q[0];
rz(-1.4486827) q[0];
sx q[0];
rz(-0.51885817) q[0];
rz(2.8209575) q[2];
sx q[2];
rz(-0.33679397) q[2];
sx q[2];
rz(-1.8048665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.022696115) q[1];
sx q[1];
rz(-1.5210946) q[1];
sx q[1];
rz(-1.8291874) q[1];
rz(-pi) q[2];
rz(1.4225716) q[3];
sx q[3];
rz(-0.61398849) q[3];
sx q[3];
rz(-2.2999291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64409488) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(-3.0268055) q[2];
rz(2.3980906) q[3];
sx q[3];
rz(-1.2657575) q[3];
sx q[3];
rz(2.6928597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153618) q[0];
sx q[0];
rz(-0.76114934) q[0];
sx q[0];
rz(2.1834955) q[0];
rz(1.505898) q[1];
sx q[1];
rz(-2.0151357) q[1];
sx q[1];
rz(-1.9912079) q[1];
rz(-0.68810473) q[2];
sx q[2];
rz(-2.2672014) q[2];
sx q[2];
rz(-1.3939569) q[2];
rz(-0.43396797) q[3];
sx q[3];
rz(-1.4186191) q[3];
sx q[3];
rz(-0.71604244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

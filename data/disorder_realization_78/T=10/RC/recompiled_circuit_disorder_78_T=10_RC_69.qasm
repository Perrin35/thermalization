OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(-2.5506033) q[0];
sx q[0];
rz(-0.58340573) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(4.1252131) q[1];
sx q[1];
rz(10.317378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4757724) q[0];
sx q[0];
rz(-1.2736397) q[0];
sx q[0];
rz(-2.1509403) q[0];
rz(-2.4721488) q[2];
sx q[2];
rz(-1.1602243) q[2];
sx q[2];
rz(-2.4759811) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6751911) q[1];
sx q[1];
rz(-2.2474504) q[1];
sx q[1];
rz(2.7387709) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0953232) q[3];
sx q[3];
rz(-0.75152961) q[3];
sx q[3];
rz(-0.7788333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083667) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-2.5657186) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.974568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28536797) q[0];
sx q[0];
rz(-1.8282187) q[0];
sx q[0];
rz(-0.060398922) q[0];
rz(-2.1115233) q[2];
sx q[2];
rz(-1.3582555) q[2];
sx q[2];
rz(0.18268798) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45707073) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(-0.80656273) q[1];
rz(-0.59252177) q[3];
sx q[3];
rz(-0.7624818) q[3];
sx q[3];
rz(-3.1069063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.033826753) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(-3.068148) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(-1.0148369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0600216) q[0];
sx q[0];
rz(-0.9284174) q[0];
sx q[0];
rz(1.5006256) q[0];
x q[1];
rz(-2.0929167) q[2];
sx q[2];
rz(-0.45959696) q[2];
sx q[2];
rz(0.93333474) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9008357) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(0.97297538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88481836) q[3];
sx q[3];
rz(-1.3664477) q[3];
sx q[3];
rz(-0.94061461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(0.69747654) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(2.5090384) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(-0.63582173) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83946562) q[0];
sx q[0];
rz(-1.3913904) q[0];
sx q[0];
rz(2.0749712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84391528) q[2];
sx q[2];
rz(-1.4471495) q[2];
sx q[2];
rz(-2.4713949) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0337452) q[1];
sx q[1];
rz(-2.3305232) q[1];
sx q[1];
rz(-1.8268405) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28663978) q[3];
sx q[3];
rz(-1.3661824) q[3];
sx q[3];
rz(2.3972547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9277966) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-0.68871838) q[2];
rz(-2.8074746) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(2.9978602) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(-0.74514666) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(2.863046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700096) q[0];
sx q[0];
rz(-2.2402813) q[0];
sx q[0];
rz(-2.5894126) q[0];
x q[1];
rz(-1.6804382) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(1.7313752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.46036938) q[1];
sx q[1];
rz(-1.9378098) q[1];
sx q[1];
rz(-2.2296434) q[1];
x q[2];
rz(-0.81551084) q[3];
sx q[3];
rz(-1.1482571) q[3];
sx q[3];
rz(0.58327196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(0.29423514) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(-0.63502216) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(3.0335398) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5393339) q[0];
sx q[0];
rz(-0.56319153) q[0];
sx q[0];
rz(-0.62298933) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4322386) q[2];
sx q[2];
rz(-1.1922622) q[2];
sx q[2];
rz(-1.1059424) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15247004) q[1];
sx q[1];
rz(-1.9580541) q[1];
sx q[1];
rz(2.6168004) q[1];
rz(-0.26515682) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1487427) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(-2.0097849) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-2.2479642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.87003) q[0];
sx q[0];
rz(-1.8989519) q[0];
sx q[0];
rz(0.66332711) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8473179) q[2];
sx q[2];
rz(-2.0308959) q[2];
sx q[2];
rz(2.0725046) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5205198) q[1];
sx q[1];
rz(-2.2398661) q[1];
sx q[1];
rz(1.3135765) q[1];
rz(-pi) q[2];
rz(3.0393533) q[3];
sx q[3];
rz(-0.7623626) q[3];
sx q[3];
rz(-1.2650714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(0.52948362) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(2.7990394) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787518) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(-3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.1539248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.740828) q[0];
sx q[0];
rz(-0.65069288) q[0];
sx q[0];
rz(3.0853737) q[0];
x q[1];
rz(-2.5904028) q[2];
sx q[2];
rz(-1.6973719) q[2];
sx q[2];
rz(1.3431431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0358419) q[1];
sx q[1];
rz(-2.3301947) q[1];
sx q[1];
rz(-1.6559385) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2418601) q[3];
sx q[3];
rz(-1.8814439) q[3];
sx q[3];
rz(1.8048546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5899137) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(0.38044688) q[2];
rz(2.0137265) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-2.5019116) q[0];
rz(1.2387964) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78478783) q[0];
sx q[0];
rz(-2.6492282) q[0];
sx q[0];
rz(-2.3193293) q[0];
rz(-pi) q[1];
rz(1.8413999) q[2];
sx q[2];
rz(-0.81297183) q[2];
sx q[2];
rz(1.1635309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6744086) q[1];
sx q[1];
rz(-1.3896175) q[1];
sx q[1];
rz(2.6513537) q[1];
rz(1.7438668) q[3];
sx q[3];
rz(-0.67865463) q[3];
sx q[3];
rz(2.9242587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-0.38273746) q[2];
rz(-0.9283723) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(0.25892192) q[0];
rz(0.71031538) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(2.6616667) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5878764) q[0];
sx q[0];
rz(-1.237861) q[0];
sx q[0];
rz(-2.2399708) q[0];
rz(-pi) q[1];
rz(2.2357113) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(0.71042176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6595315) q[1];
sx q[1];
rz(-1.2188984) q[1];
sx q[1];
rz(0.54606502) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2009833) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(-1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2991128) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-0.97933979) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(-2.813415) q[2];
sx q[2];
rz(-0.50445088) q[2];
sx q[2];
rz(2.9403461) q[2];
rz(-1.8388207) q[3];
sx q[3];
rz(-1.839073) q[3];
sx q[3];
rz(3.0243235) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

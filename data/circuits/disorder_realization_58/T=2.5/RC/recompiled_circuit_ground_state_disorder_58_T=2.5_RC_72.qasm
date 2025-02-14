OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6405606) q[0];
sx q[0];
rz(-0.80067331) q[0];
sx q[0];
rz(-0.14525695) q[0];
rz(2.2486806) q[1];
sx q[1];
rz(-0.4141663) q[1];
sx q[1];
rz(-1.1069586) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9710179) q[0];
sx q[0];
rz(-0.98890328) q[0];
sx q[0];
rz(-2.2547743) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2824481) q[2];
sx q[2];
rz(-2.7224053) q[2];
sx q[2];
rz(-2.223658) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8095137) q[1];
sx q[1];
rz(-0.65843907) q[1];
sx q[1];
rz(2.0488942) q[1];
rz(-pi) q[2];
rz(1.065973) q[3];
sx q[3];
rz(-0.42718605) q[3];
sx q[3];
rz(1.8452725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.60078159) q[2];
sx q[2];
rz(-1.7822219) q[2];
sx q[2];
rz(0.092279807) q[2];
rz(-1.7021092) q[3];
sx q[3];
rz(-1.9974134) q[3];
sx q[3];
rz(2.2030742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.1321024) q[0];
sx q[0];
rz(-1.9785896) q[0];
sx q[0];
rz(0.60390419) q[0];
rz(1.6101135) q[1];
sx q[1];
rz(-1.8096626) q[1];
sx q[1];
rz(1.5276705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22900362) q[0];
sx q[0];
rz(-2.4618076) q[0];
sx q[0];
rz(-2.0849821) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3983922) q[2];
sx q[2];
rz(-2.024352) q[2];
sx q[2];
rz(-1.989606) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81948796) q[1];
sx q[1];
rz(-0.27399644) q[1];
sx q[1];
rz(0.042983965) q[1];
rz(-0.0068596938) q[3];
sx q[3];
rz(-3.0433972) q[3];
sx q[3];
rz(-2.6650037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70832002) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(-1.0478919) q[2];
rz(-0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(-1.9765114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79353756) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(-1.0070356) q[0];
rz(0.29193613) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(2.4709591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1517253) q[0];
sx q[0];
rz(-1.4833996) q[0];
sx q[0];
rz(-0.94957994) q[0];
x q[1];
rz(0.88460716) q[2];
sx q[2];
rz(-0.90145196) q[2];
sx q[2];
rz(1.455292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5340928) q[1];
sx q[1];
rz(-2.1779446) q[1];
sx q[1];
rz(-1.9837888) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.241983) q[3];
sx q[3];
rz(-1.4610664) q[3];
sx q[3];
rz(1.3789267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2900419) q[2];
sx q[2];
rz(-2.2883577) q[2];
sx q[2];
rz(-0.89149371) q[2];
rz(-1.1024891) q[3];
sx q[3];
rz(-0.99415556) q[3];
sx q[3];
rz(1.7360784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.3201228) q[0];
sx q[0];
rz(-1.7550884) q[0];
sx q[0];
rz(1.0572877) q[0];
rz(0.0013466324) q[1];
sx q[1];
rz(-2.3206382) q[1];
sx q[1];
rz(-0.79426208) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0214572) q[0];
sx q[0];
rz(-1.5524143) q[0];
sx q[0];
rz(-0.6619428) q[0];
x q[1];
rz(2.7873331) q[2];
sx q[2];
rz(-1.2984167) q[2];
sx q[2];
rz(-2.543022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9201607) q[1];
sx q[1];
rz(-2.2756612) q[1];
sx q[1];
rz(-2.1842498) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4554855) q[3];
sx q[3];
rz(-1.274935) q[3];
sx q[3];
rz(1.4747696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0346251) q[2];
sx q[2];
rz(-0.62109533) q[2];
sx q[2];
rz(-2.1194439) q[2];
rz(-1.8978097) q[3];
sx q[3];
rz(-0.94977489) q[3];
sx q[3];
rz(1.3589842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46955243) q[0];
sx q[0];
rz(-1.38009) q[0];
sx q[0];
rz(0.45355466) q[0];
rz(-1.0376616) q[1];
sx q[1];
rz(-2.0062168) q[1];
sx q[1];
rz(2.3496148) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2642518) q[0];
sx q[0];
rz(-1.8487572) q[0];
sx q[0];
rz(-3.0959913) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7466082) q[2];
sx q[2];
rz(-2.1126267) q[2];
sx q[2];
rz(-1.3196368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5663858) q[1];
sx q[1];
rz(-2.7820911) q[1];
sx q[1];
rz(-2.8043037) q[1];
x q[2];
rz(1.7089473) q[3];
sx q[3];
rz(-1.5524716) q[3];
sx q[3];
rz(1.5973171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3114634) q[2];
sx q[2];
rz(-2.4757803) q[2];
sx q[2];
rz(-2.8361481) q[2];
rz(-0.18181248) q[3];
sx q[3];
rz(-1.612855) q[3];
sx q[3];
rz(-0.73604933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5859454) q[0];
sx q[0];
rz(-2.3190627) q[0];
sx q[0];
rz(0.12538759) q[0];
rz(-1.5749982) q[1];
sx q[1];
rz(-1.6927203) q[1];
sx q[1];
rz(-0.032141846) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9956995) q[0];
sx q[0];
rz(-0.75135485) q[0];
sx q[0];
rz(2.8180647) q[0];
rz(1.5015542) q[2];
sx q[2];
rz(-0.73261315) q[2];
sx q[2];
rz(0.68147269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7508538) q[1];
sx q[1];
rz(-1.097569) q[1];
sx q[1];
rz(2.3473032) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12539668) q[3];
sx q[3];
rz(-2.3920357) q[3];
sx q[3];
rz(-2.9182485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0142168) q[2];
sx q[2];
rz(-1.2140423) q[2];
sx q[2];
rz(1.1227013) q[2];
rz(-3.0681916) q[3];
sx q[3];
rz(-2.1843036) q[3];
sx q[3];
rz(-0.61298031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70925322) q[0];
sx q[0];
rz(-1.657635) q[0];
sx q[0];
rz(2.6313229) q[0];
rz(-0.36422745) q[1];
sx q[1];
rz(-0.42114708) q[1];
sx q[1];
rz(1.4998923) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071935805) q[0];
sx q[0];
rz(-0.47635117) q[0];
sx q[0];
rz(-1.8029965) q[0];
rz(-3.0232863) q[2];
sx q[2];
rz(-1.4139172) q[2];
sx q[2];
rz(-1.1821234) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81699521) q[1];
sx q[1];
rz(-2.5909068) q[1];
sx q[1];
rz(-2.703642) q[1];
rz(-0.73697258) q[3];
sx q[3];
rz(-1.7993357) q[3];
sx q[3];
rz(1.4481973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69786543) q[2];
sx q[2];
rz(-1.5550193) q[2];
sx q[2];
rz(-0.86722428) q[2];
rz(-1.5602268) q[3];
sx q[3];
rz(-2.9428704) q[3];
sx q[3];
rz(1.3377415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7241868) q[0];
sx q[0];
rz(-3.0751808) q[0];
sx q[0];
rz(2.7253286) q[0];
rz(1.9789713) q[1];
sx q[1];
rz(-1.1888209) q[1];
sx q[1];
rz(-2.3847041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0666733) q[0];
sx q[0];
rz(-0.52353379) q[0];
sx q[0];
rz(-2.174211) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1882638) q[2];
sx q[2];
rz(-2.4673415) q[2];
sx q[2];
rz(-1.9566388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50584201) q[1];
sx q[1];
rz(-1.381784) q[1];
sx q[1];
rz(-2.0057949) q[1];
x q[2];
rz(2.5060095) q[3];
sx q[3];
rz(-1.1867689) q[3];
sx q[3];
rz(1.7049198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4550712) q[2];
sx q[2];
rz(-1.4342118) q[2];
sx q[2];
rz(-1.7835468) q[2];
rz(2.3250735) q[3];
sx q[3];
rz(-1.9101382) q[3];
sx q[3];
rz(2.9787298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002976) q[0];
sx q[0];
rz(-1.6930027) q[0];
sx q[0];
rz(-1.7342389) q[0];
rz(0.74527144) q[1];
sx q[1];
rz(-2.407357) q[1];
sx q[1];
rz(-1.5596681) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9055515) q[0];
sx q[0];
rz(-0.92285448) q[0];
sx q[0];
rz(2.7929162) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4938131) q[2];
sx q[2];
rz(-1.5967007) q[2];
sx q[2];
rz(-1.4454973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7551387) q[1];
sx q[1];
rz(-1.6999131) q[1];
sx q[1];
rz(-2.0172202) q[1];
rz(1.4780864) q[3];
sx q[3];
rz(-0.57766908) q[3];
sx q[3];
rz(-0.81837624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4173296) q[2];
sx q[2];
rz(-1.686692) q[2];
sx q[2];
rz(2.0261185) q[2];
rz(-0.25350246) q[3];
sx q[3];
rz(-1.8799672) q[3];
sx q[3];
rz(-3.1326262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215006) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(2.3790835) q[0];
rz(2.1992042) q[1];
sx q[1];
rz(-2.0768879) q[1];
sx q[1];
rz(1.2368088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.002703) q[0];
sx q[0];
rz(-1.5093439) q[0];
sx q[0];
rz(-2.6451254) q[0];
x q[1];
rz(-0.89434172) q[2];
sx q[2];
rz(-2.333775) q[2];
sx q[2];
rz(1.7730912) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0227388) q[1];
sx q[1];
rz(-0.035422649) q[1];
sx q[1];
rz(0.78411786) q[1];
rz(-1.9754161) q[3];
sx q[3];
rz(-1.5050833) q[3];
sx q[3];
rz(-0.031377553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2463871) q[2];
sx q[2];
rz(-0.095194101) q[2];
sx q[2];
rz(-3.134356) q[2];
rz(0.17035189) q[3];
sx q[3];
rz(-1.6536313) q[3];
sx q[3];
rz(-0.65792221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687179) q[0];
sx q[0];
rz(-1.7418516) q[0];
sx q[0];
rz(-2.0280784) q[0];
rz(-0.089182236) q[1];
sx q[1];
rz(-1.6249648) q[1];
sx q[1];
rz(-1.9393495) q[1];
rz(-1.4349277) q[2];
sx q[2];
rz(-1.7262222) q[2];
sx q[2];
rz(-1.8651967) q[2];
rz(0.54696541) q[3];
sx q[3];
rz(-2.7934358) q[3];
sx q[3];
rz(-1.1478333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(0.068109186) q[0];
sx q[0];
rz(13.343233) q[0];
rz(-1.3070973) q[1];
sx q[1];
rz(-0.96944648) q[1];
sx q[1];
rz(-1.8196655) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.043524) q[0];
sx q[0];
rz(-2.4130445) q[0];
sx q[0];
rz(1.3619955) q[0];
x q[1];
rz(-1.2250617) q[2];
sx q[2];
rz(-1.099473) q[2];
sx q[2];
rz(2.2668348) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.894134) q[1];
sx q[1];
rz(-1.1620635) q[1];
sx q[1];
rz(-0.031886727) q[1];
rz(0.98422756) q[3];
sx q[3];
rz(-2.7968458) q[3];
sx q[3];
rz(-1.00351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(2.784101) q[2];
rz(-0.56719559) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(-0.43958694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95328632) q[0];
sx q[0];
rz(-0.28306857) q[0];
sx q[0];
rz(1.8081007) q[0];
rz(-0.4370583) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(-0.83998799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6463722) q[0];
sx q[0];
rz(-0.53421181) q[0];
sx q[0];
rz(-1.2665126) q[0];
x q[1];
rz(-2.5469668) q[2];
sx q[2];
rz(-1.42808) q[2];
sx q[2];
rz(0.12446257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4639791) q[1];
sx q[1];
rz(-2.2871454) q[1];
sx q[1];
rz(2.3542627) q[1];
rz(-pi) q[2];
rz(-1.5072823) q[3];
sx q[3];
rz(-0.5507142) q[3];
sx q[3];
rz(-0.39001071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4240894) q[2];
sx q[2];
rz(-1.6987897) q[2];
sx q[2];
rz(-1.4060414) q[2];
rz(-0.16215912) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9185987) q[0];
sx q[0];
rz(-3.0610237) q[0];
sx q[0];
rz(0.42174569) q[0];
rz(-0.84284198) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(-2.8657894) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3746412) q[0];
sx q[0];
rz(-1.8948484) q[0];
sx q[0];
rz(1.3270017) q[0];
rz(-pi) q[1];
rz(2.6709072) q[2];
sx q[2];
rz(-1.3725201) q[2];
sx q[2];
rz(0.40775611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0581857) q[1];
sx q[1];
rz(-0.32489932) q[1];
sx q[1];
rz(0.16375457) q[1];
rz(-0.86902501) q[3];
sx q[3];
rz(-1.6459673) q[3];
sx q[3];
rz(-1.4787256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80321035) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(-0.69022834) q[2];
rz(-0.35987443) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(0.38351044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279219) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(2.2699455) q[0];
rz(-0.40920416) q[1];
sx q[1];
rz(-0.49574655) q[1];
sx q[1];
rz(-0.31235487) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86996468) q[0];
sx q[0];
rz(-1.6883435) q[0];
sx q[0];
rz(0.44674504) q[0];
rz(2.6957507) q[2];
sx q[2];
rz(-0.22805691) q[2];
sx q[2];
rz(-2.7331309) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1725268) q[1];
sx q[1];
rz(-1.2102264) q[1];
sx q[1];
rz(1.7883302) q[1];
x q[2];
rz(1.7658224) q[3];
sx q[3];
rz(-0.93149501) q[3];
sx q[3];
rz(-2.2488058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13545869) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(0.82497605) q[2];
rz(1.4247591) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(0.43535522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41998) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(-0.87442526) q[0];
rz(2.8127316) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(-1.0985451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.335137) q[0];
sx q[0];
rz(-2.7493434) q[0];
sx q[0];
rz(-0.29229887) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1055036) q[2];
sx q[2];
rz(-2.3054625) q[2];
sx q[2];
rz(-0.058908894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85531536) q[1];
sx q[1];
rz(-0.5868656) q[1];
sx q[1];
rz(-1.0081968) q[1];
x q[2];
rz(-2.1608814) q[3];
sx q[3];
rz(-1.3385337) q[3];
sx q[3];
rz(1.6764318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.141779) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(-1.4874124) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(-0.28356799) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11033002) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(1.3856101) q[0];
rz(-2.022187) q[1];
sx q[1];
rz(-1.4675354) q[1];
sx q[1];
rz(1.3853692) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0273756) q[0];
sx q[0];
rz(-1.8226591) q[0];
sx q[0];
rz(-2.2172539) q[0];
rz(-pi) q[1];
rz(-0.15689705) q[2];
sx q[2];
rz(-2.6378535) q[2];
sx q[2];
rz(-2.4258326) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1839528) q[1];
sx q[1];
rz(-1.7639177) q[1];
sx q[1];
rz(-2.814358) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9870944) q[3];
sx q[3];
rz(-1.6067926) q[3];
sx q[3];
rz(-2.846623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0783656) q[2];
sx q[2];
rz(-0.93457064) q[2];
sx q[2];
rz(2.7654977) q[2];
rz(-2.0070576) q[3];
sx q[3];
rz(-0.36342707) q[3];
sx q[3];
rz(-2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1217839) q[0];
sx q[0];
rz(-2.5406295) q[0];
sx q[0];
rz(0.10251775) q[0];
rz(0.58492297) q[1];
sx q[1];
rz(-0.94766098) q[1];
sx q[1];
rz(-1.1835416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2823806) q[0];
sx q[0];
rz(-1.4420274) q[0];
sx q[0];
rz(2.9984401) q[0];
rz(-pi) q[1];
rz(-2.4790484) q[2];
sx q[2];
rz(-1.8383611) q[2];
sx q[2];
rz(-0.91987687) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2940604) q[1];
sx q[1];
rz(-2.4449663) q[1];
sx q[1];
rz(1.1678334) q[1];
rz(0.4966708) q[3];
sx q[3];
rz(-2.6210945) q[3];
sx q[3];
rz(-0.48540533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3465603) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(0.19980508) q[2];
rz(1.0605158) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(0.04960355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878433) q[0];
sx q[0];
rz(-1.6596154) q[0];
sx q[0];
rz(-0.31295452) q[0];
rz(2.1921659) q[1];
sx q[1];
rz(-1.3525617) q[1];
sx q[1];
rz(-1.4535646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67927406) q[0];
sx q[0];
rz(-1.7574213) q[0];
sx q[0];
rz(-1.6629205) q[0];
rz(-pi) q[1];
rz(1.6498927) q[2];
sx q[2];
rz(-2.1888615) q[2];
sx q[2];
rz(1.6343509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.955757) q[1];
sx q[1];
rz(-1.794853) q[1];
sx q[1];
rz(-2.4332389) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9096776) q[3];
sx q[3];
rz(-2.0032855) q[3];
sx q[3];
rz(-0.56604715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7877385) q[2];
sx q[2];
rz(-0.9757897) q[2];
sx q[2];
rz(-1.3346416) q[2];
rz(-1.3245964) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(1.2142396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036309328) q[0];
sx q[0];
rz(-0.61573017) q[0];
sx q[0];
rz(-0.69806725) q[0];
rz(-0.87567323) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(-0.36044136) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9746285) q[0];
sx q[0];
rz(-0.95050838) q[0];
sx q[0];
rz(-2.1168769) q[0];
x q[1];
rz(-1.5064042) q[2];
sx q[2];
rz(-1.9528104) q[2];
sx q[2];
rz(0.40859336) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0086509) q[1];
sx q[1];
rz(-1.1906149) q[1];
sx q[1];
rz(0.93312414) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4029986) q[3];
sx q[3];
rz(-2.1133964) q[3];
sx q[3];
rz(0.29796539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91161072) q[2];
sx q[2];
rz(-1.4260099) q[2];
sx q[2];
rz(-0.81725517) q[2];
rz(0.83003712) q[3];
sx q[3];
rz(-2.6006112) q[3];
sx q[3];
rz(0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5068186) q[0];
sx q[0];
rz(-0.84832484) q[0];
sx q[0];
rz(-3.1095374) q[0];
rz(-1.3196779) q[1];
sx q[1];
rz(-1.7664884) q[1];
sx q[1];
rz(-0.65418902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.74943) q[0];
sx q[0];
rz(-1.4512968) q[0];
sx q[0];
rz(0.197535) q[0];
rz(2.3532969) q[2];
sx q[2];
rz(-1.239778) q[2];
sx q[2];
rz(-0.85110474) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54683751) q[1];
sx q[1];
rz(-1.0887301) q[1];
sx q[1];
rz(1.0422816) q[1];
x q[2];
rz(-0.93864949) q[3];
sx q[3];
rz(-1.3397927) q[3];
sx q[3];
rz(-1.3579901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0111982) q[2];
sx q[2];
rz(-2.0989213) q[2];
sx q[2];
rz(-1.0065669) q[2];
rz(-0.69342962) q[3];
sx q[3];
rz(-1.8208241) q[3];
sx q[3];
rz(1.2815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4257767) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(0.88059942) q[1];
sx q[1];
rz(-2.0279299) q[1];
sx q[1];
rz(0.50150064) q[1];
rz(2.2702552) q[2];
sx q[2];
rz(-1.9819145) q[2];
sx q[2];
rz(-2.869538) q[2];
rz(1.952259) q[3];
sx q[3];
rz(-2.2305924) q[3];
sx q[3];
rz(2.3371405) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

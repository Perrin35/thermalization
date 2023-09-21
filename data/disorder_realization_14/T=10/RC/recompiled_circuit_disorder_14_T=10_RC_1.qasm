OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.68552652) q[0];
sx q[0];
rz(-2.752562) q[0];
sx q[0];
rz(-2.2580137) q[0];
rz(-0.0097302516) q[1];
sx q[1];
rz(-1.4571804) q[1];
sx q[1];
rz(1.943346) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6057974) q[0];
sx q[0];
rz(-1.3942766) q[0];
sx q[0];
rz(-1.6685041) q[0];
rz(2.3660907) q[2];
sx q[2];
rz(-1.8047793) q[2];
sx q[2];
rz(0.50228679) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4858422) q[1];
sx q[1];
rz(-1.2409004) q[1];
sx q[1];
rz(-0.54539036) q[1];
rz(1.4000113) q[3];
sx q[3];
rz(-1.6163974) q[3];
sx q[3];
rz(3.1053229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1951695) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(-0.18134376) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.2657335) q[3];
sx q[3];
rz(0.75631022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(-2.7547577) q[0];
rz(2.6392) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5997255) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72774784) q[0];
sx q[0];
rz(-1.9808597) q[0];
sx q[0];
rz(-2.4983665) q[0];
rz(0.16427152) q[2];
sx q[2];
rz(-1.1195682) q[2];
sx q[2];
rz(2.3300366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.220037) q[1];
sx q[1];
rz(-1.5350635) q[1];
sx q[1];
rz(1.5608556) q[1];
rz(-pi) q[2];
rz(-2.6582791) q[3];
sx q[3];
rz(-1.9938853) q[3];
sx q[3];
rz(-1.9139293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.60275045) q[2];
sx q[2];
rz(-0.87783146) q[2];
sx q[2];
rz(1.7269469) q[2];
rz(0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8787815) q[0];
sx q[0];
rz(-1.6101863) q[0];
sx q[0];
rz(0.61022726) q[0];
rz(1.2894851) q[1];
sx q[1];
rz(-2.162343) q[1];
sx q[1];
rz(-2.1496225) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6863166) q[0];
sx q[0];
rz(-2.1222097) q[0];
sx q[0];
rz(0.21689143) q[0];
rz(0.66833468) q[2];
sx q[2];
rz(-1.4112345) q[2];
sx q[2];
rz(0.72232407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5132644) q[1];
sx q[1];
rz(-0.4497512) q[1];
sx q[1];
rz(-0.20723923) q[1];
rz(-0.45735995) q[3];
sx q[3];
rz(-1.7914346) q[3];
sx q[3];
rz(1.8822576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38301864) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(-0.26829159) q[2];
rz(2.7475083) q[3];
sx q[3];
rz(-1.9106617) q[3];
sx q[3];
rz(-0.18850732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2878993) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(-2.7365141) q[0];
rz(-2.4515117) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(-2.4437723) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47397428) q[0];
sx q[0];
rz(-1.6606497) q[0];
sx q[0];
rz(3.1180179) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9525098) q[2];
sx q[2];
rz(-1.0587947) q[2];
sx q[2];
rz(2.2005759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78742541) q[1];
sx q[1];
rz(-2.0560871) q[1];
sx q[1];
rz(2.9334958) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29699765) q[3];
sx q[3];
rz(-1.9747435) q[3];
sx q[3];
rz(1.7904074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71965376) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.6323803) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.6633165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(-1.0255381) q[0];
rz(-2.569596) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(2.5122723) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2239265) q[0];
sx q[0];
rz(-1.1610982) q[0];
sx q[0];
rz(1.8976589) q[0];
x q[1];
rz(2.8867678) q[2];
sx q[2];
rz(-0.91081589) q[2];
sx q[2];
rz(-0.13435907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1197966) q[1];
sx q[1];
rz(-2.6152059) q[1];
sx q[1];
rz(0.22967931) q[1];
x q[2];
rz(-1.945433) q[3];
sx q[3];
rz(-2.5064838) q[3];
sx q[3];
rz(2.0600968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5489674) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(-0.34234753) q[2];
rz(-1.7957548) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(-2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65258566) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(0.18519369) q[0];
rz(1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(-1.3669744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2469251) q[0];
sx q[0];
rz(-1.3498382) q[0];
sx q[0];
rz(2.5470683) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.30079) q[2];
sx q[2];
rz(-1.7992939) q[2];
sx q[2];
rz(-0.92323869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5609315) q[1];
sx q[1];
rz(-1.7226189) q[1];
sx q[1];
rz(3.0488357) q[1];
rz(-pi) q[2];
rz(1.1432511) q[3];
sx q[3];
rz(-2.098009) q[3];
sx q[3];
rz(-2.9851819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24484816) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(0.883376) q[2];
rz(-1.4128489) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25154034) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(1.863377) q[0];
rz(-0.037840769) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(-1.7657123) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8517075) q[0];
sx q[0];
rz(-0.87806784) q[0];
sx q[0];
rz(0.38498199) q[0];
rz(-pi) q[1];
rz(0.059869754) q[2];
sx q[2];
rz(-2.0436358) q[2];
sx q[2];
rz(1.5470031) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.059767698) q[1];
sx q[1];
rz(-1.4653112) q[1];
sx q[1];
rz(-0.11591537) q[1];
rz(-pi) q[2];
rz(-1.901058) q[3];
sx q[3];
rz(-1.5152647) q[3];
sx q[3];
rz(-2.8020669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4574796) q[2];
sx q[2];
rz(-2.423968) q[2];
sx q[2];
rz(-0.73105556) q[2];
rz(-3.030792) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(0.69537648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59657997) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(-2.4080283) q[0];
rz(2.5336174) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(0.2342934) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83345862) q[0];
sx q[0];
rz(-1.3123543) q[0];
sx q[0];
rz(1.0779557) q[0];
x q[1];
rz(-0.73306002) q[2];
sx q[2];
rz(-2.2736533) q[2];
sx q[2];
rz(-1.6349863) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56470796) q[1];
sx q[1];
rz(-1.3815834) q[1];
sx q[1];
rz(2.3343759) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33741823) q[3];
sx q[3];
rz(-1.0113582) q[3];
sx q[3];
rz(-2.9448201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0155448) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(1.7162494) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5995246) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(1.9482127) q[0];
rz(1.2606196) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(0.9448005) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1856857) q[0];
sx q[0];
rz(-2.3392896) q[0];
sx q[0];
rz(0.54922744) q[0];
rz(-0.94061942) q[2];
sx q[2];
rz(-2.665328) q[2];
sx q[2];
rz(-2.1911052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8560801) q[1];
sx q[1];
rz(-1.7417007) q[1];
sx q[1];
rz(0.28275615) q[1];
rz(-2.9444669) q[3];
sx q[3];
rz(-1.2211495) q[3];
sx q[3];
rz(2.8417348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9251359) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(0.28820583) q[2];
rz(-2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.5079927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(2.2286041) q[0];
rz(-0.37462014) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(-0.8909117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83715314) q[0];
sx q[0];
rz(-1.3225978) q[0];
sx q[0];
rz(2.7960143) q[0];
rz(2.8577523) q[2];
sx q[2];
rz(-2.1643057) q[2];
sx q[2];
rz(0.80988353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76162321) q[1];
sx q[1];
rz(-2.7459811) q[1];
sx q[1];
rz(-0.53669866) q[1];
rz(-pi) q[2];
rz(2.9534146) q[3];
sx q[3];
rz(-1.6450226) q[3];
sx q[3];
rz(1.0159514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23641071) q[2];
sx q[2];
rz(-2.2832401) q[2];
sx q[2];
rz(2.6399844) q[2];
rz(-1.2891399) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(2.6954209) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5065153) q[0];
sx q[0];
rz(-1.7000533) q[0];
sx q[0];
rz(0.62361367) q[0];
rz(1.1322017) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(0.26666178) q[2];
sx q[2];
rz(-1.7186135) q[2];
sx q[2];
rz(2.7415822) q[2];
rz(3.0922079) q[3];
sx q[3];
rz(-2.1168843) q[3];
sx q[3];
rz(-2.4612853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
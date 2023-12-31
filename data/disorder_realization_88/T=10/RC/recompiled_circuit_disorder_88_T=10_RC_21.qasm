OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(2.8136301) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(-2.9405138) q[1];
sx q[1];
rz(-0.091436401) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51973242) q[0];
sx q[0];
rz(-0.73017263) q[0];
sx q[0];
rz(-0.61868389) q[0];
rz(-pi) q[1];
rz(-2.6684127) q[2];
sx q[2];
rz(-0.28943974) q[2];
sx q[2];
rz(2.6602886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58127296) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(-1.9074349) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8880635) q[3];
sx q[3];
rz(-2.5228365) q[3];
sx q[3];
rz(-1.7537102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(-1.1260024) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(2.2385105) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098410957) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(-0.19031659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8961401) q[0];
sx q[0];
rz(-1.3724694) q[0];
sx q[0];
rz(-1.1597) q[0];
x q[1];
rz(-2.5311567) q[2];
sx q[2];
rz(-0.62610859) q[2];
sx q[2];
rz(0.49952835) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39365921) q[1];
sx q[1];
rz(-2.2927082) q[1];
sx q[1];
rz(-3.1118803) q[1];
x q[2];
rz(-2.042949) q[3];
sx q[3];
rz(-0.90001366) q[3];
sx q[3];
rz(2.9428279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(-2.3213342) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4688063) q[0];
sx q[0];
rz(-2.5350223) q[0];
sx q[0];
rz(2.9787105) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6349995) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.9931672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51860147) q[1];
sx q[1];
rz(-0.90497436) q[1];
sx q[1];
rz(1.4856505) q[1];
x q[2];
rz(0.30001254) q[3];
sx q[3];
rz(-2.1999151) q[3];
sx q[3];
rz(2.574207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(1.2134264) q[2];
rz(-0.16472566) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(-0.18606342) q[0];
rz(-2.9371254) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064917795) q[0];
sx q[0];
rz(-3.0041822) q[0];
sx q[0];
rz(1.2491559) q[0];
rz(0.16764955) q[2];
sx q[2];
rz(-1.285458) q[2];
sx q[2];
rz(-0.97380762) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0399196) q[1];
sx q[1];
rz(-2.1664201) q[1];
sx q[1];
rz(-2.4458829) q[1];
rz(3.1277539) q[3];
sx q[3];
rz(-1.7503563) q[3];
sx q[3];
rz(-1.9765215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(0.91919351) q[2];
rz(0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(1.2478158) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677251) q[0];
sx q[0];
rz(-2.6350832) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.426288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844937) q[0];
sx q[0];
rz(-1.6528659) q[0];
sx q[0];
rz(-0.9135855) q[0];
rz(-1.8833141) q[2];
sx q[2];
rz(-1.7973571) q[2];
sx q[2];
rz(-2.8001919) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9702455) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(-1.4096178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3705195) q[3];
sx q[3];
rz(-1.8043451) q[3];
sx q[3];
rz(-2.0296823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8695801) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(-0.4367691) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(1.0304931) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(2.5700263) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58493462) q[0];
sx q[0];
rz(-0.99781636) q[0];
sx q[0];
rz(-0.95227382) q[0];
rz(-pi) q[1];
rz(1.4095441) q[2];
sx q[2];
rz(-2.469392) q[2];
sx q[2];
rz(1.6598998) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4042582) q[1];
sx q[1];
rz(-1.3306276) q[1];
sx q[1];
rz(2.0391697) q[1];
x q[2];
rz(-0.84047517) q[3];
sx q[3];
rz(-2.2395036) q[3];
sx q[3];
rz(-0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17343865) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(-0.56419939) q[2];
rz(3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-2.8469767) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0903704) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-0.66551048) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3411361) q[0];
sx q[0];
rz(-1.8130214) q[0];
sx q[0];
rz(-3.1405297) q[0];
rz(-pi) q[1];
rz(2.268157) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(0.21542491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29384089) q[1];
sx q[1];
rz(-1.7898702) q[1];
sx q[1];
rz(-2.0391383) q[1];
rz(2.1036759) q[3];
sx q[3];
rz(-2.8623192) q[3];
sx q[3];
rz(-0.52475196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15726382) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.7187913) q[2];
rz(3.026399) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-0.19083047) q[0];
rz(0.62675369) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-2.802882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8341634) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(2.9817581) q[0];
rz(1.3391101) q[2];
sx q[2];
rz(-2.0121687) q[2];
sx q[2];
rz(1.0362792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6845477) q[1];
sx q[1];
rz(-1.0984048) q[1];
sx q[1];
rz(0.12462516) q[1];
x q[2];
rz(-0.87266463) q[3];
sx q[3];
rz(-0.41655311) q[3];
sx q[3];
rz(-1.0636898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-2.4411566) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(-2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.609628) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(-2.530653) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(0.13959612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1698648) q[0];
sx q[0];
rz(-1.5507878) q[0];
sx q[0];
rz(1.6318897) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81462461) q[2];
sx q[2];
rz(-1.5280208) q[2];
sx q[2];
rz(-0.2864366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2792564) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(0.52211296) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.052501909) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(-0.98464636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(2.810478) q[2];
rz(-0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-0.087879114) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.4846444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907222) q[0];
sx q[0];
rz(-1.8403247) q[0];
sx q[0];
rz(1.1140633) q[0];
rz(-pi) q[1];
rz(2.9160203) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(3.0849506) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(-1.9468716) q[1];
rz(-1.5649892) q[3];
sx q[3];
rz(-2.4366637) q[3];
sx q[3];
rz(-2.4095636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2075656) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375297) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(2.9539625) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(0.27992237) q[2];
sx q[2];
rz(-1.4301849) q[2];
sx q[2];
rz(2.4591597) q[2];
rz(-0.88541661) q[3];
sx q[3];
rz(-1.0049184) q[3];
sx q[3];
rz(0.64941209) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

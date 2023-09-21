OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(4.3972754) q[0];
sx q[0];
rz(9.7527405) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5782195) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(-0.63011516) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7056866) q[2];
sx q[2];
rz(-1.3139259) q[2];
sx q[2];
rz(-0.97193064) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4788471) q[1];
sx q[1];
rz(-1.4884243) q[1];
sx q[1];
rz(-1.3326416) q[1];
rz(0.21858469) q[3];
sx q[3];
rz(-0.9871452) q[3];
sx q[3];
rz(-1.0047131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(-0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098410957) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(-0.69357187) q[0];
rz(-1.0961078) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(-0.19031659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2454525) q[0];
sx q[0];
rz(-1.3724694) q[0];
sx q[0];
rz(-1.1597) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61043592) q[2];
sx q[2];
rz(-2.5154841) q[2];
sx q[2];
rz(2.6420643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39365921) q[1];
sx q[1];
rz(-2.2927082) q[1];
sx q[1];
rz(3.1118803) q[1];
rz(-2.4137647) q[3];
sx q[3];
rz(-1.2065294) q[3];
sx q[3];
rz(-1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(-0.1427342) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(-0.90863168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.2480199) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4753715) q[0];
sx q[0];
rz(-2.1682122) q[0];
sx q[0];
rz(-1.682838) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6349995) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(1.9931672) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51860147) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(1.4856505) q[1];
rz(-pi) q[2];
rz(0.30001254) q[3];
sx q[3];
rz(-0.94167751) q[3];
sx q[3];
rz(0.56738561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32039207) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(0.16472566) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(-2.9555292) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545249) q[0];
sx q[0];
rz(-1.6141119) q[0];
sx q[0];
rz(-1.4403507) q[0];
rz(-pi) q[1];
rz(2.0879891) q[2];
sx q[2];
rz(-2.811811) q[2];
sx q[2];
rz(1.6264548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.101673) q[1];
sx q[1];
rz(-0.97517255) q[1];
sx q[1];
rz(-0.69570978) q[1];
x q[2];
rz(1.4947055) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(-2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(2.2223991) q[2];
rz(-2.8202608) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(-2.2633973) q[0];
rz(1.325266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.7153046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8910687) q[0];
sx q[0];
rz(-2.2254125) q[0];
sx q[0];
rz(-3.0380681) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8833141) q[2];
sx q[2];
rz(-1.7973571) q[2];
sx q[2];
rz(-2.8001919) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1713472) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(1.7319748) q[1];
rz(-pi) q[2];
rz(-0.32894965) q[3];
sx q[3];
rz(-0.79862404) q[3];
sx q[3];
rz(-2.9165099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(-2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3549266) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(1.0304931) q[0];
rz(-2.4018535) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(0.57156634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7874998) q[0];
sx q[0];
rz(-2.0797074) q[0];
sx q[0];
rz(2.4718168) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90494855) q[2];
sx q[2];
rz(-1.4706503) q[2];
sx q[2];
rz(0.03749321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7476269) q[1];
sx q[1];
rz(-2.6193301) q[1];
sx q[1];
rz(2.0678492) q[1];
rz(-pi) q[2];
rz(2.4404293) q[3];
sx q[3];
rz(-0.94651604) q[3];
sx q[3];
rz(3.0091156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(-3.0155904) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(-2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(-2.4760822) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8004566) q[0];
sx q[0];
rz(-1.3285713) q[0];
sx q[0];
rz(-0.0010629396) q[0];
x q[1];
rz(-0.87343563) q[2];
sx q[2];
rz(-1.6849815) q[2];
sx q[2];
rz(2.9261677) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29384089) q[1];
sx q[1];
rz(-1.3517225) q[1];
sx q[1];
rz(-2.0391383) q[1];
rz(-pi) q[2];
rz(2.9969278) q[3];
sx q[3];
rz(-1.8105227) q[3];
sx q[3];
rz(-3.1160115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15726382) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(1.4228014) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(0.19083047) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-2.802882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8341634) q[0];
sx q[0];
rz(-1.9142262) q[0];
sx q[0];
rz(-0.15983454) q[0];
x q[1];
rz(-1.8024826) q[2];
sx q[2];
rz(-2.0121687) q[2];
sx q[2];
rz(1.0362792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.056811) q[1];
sx q[1];
rz(-1.6817131) q[1];
sx q[1];
rz(-2.0463498) q[1];
rz(-pi) q[2];
rz(-2.8645105) q[3];
sx q[3];
rz(-1.8859366) q[3];
sx q[3];
rz(-1.806123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(-0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(-2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53196466) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-3.0019965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7394373) q[0];
sx q[0];
rz(-1.6318775) q[0];
sx q[0];
rz(-3.1215467) q[0];
rz(-3.082824) q[2];
sx q[2];
rz(-2.3261056) q[2];
sx q[2];
rz(1.2440484) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4715251) q[1];
sx q[1];
rz(-1.0499665) q[1];
sx q[1];
rz(-1.4937558) q[1];
rz(-1.5009874) q[3];
sx q[3];
rz(-2.4959292) q[3];
sx q[3];
rz(-1.071969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6167986) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(2.810478) q[2];
rz(2.3838499) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(-3.0537135) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(-0.18558003) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.4846444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907222) q[0];
sx q[0];
rz(-1.301268) q[0];
sx q[0];
rz(2.0275293) q[0];
x q[1];
rz(2.9160203) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(-0.05664209) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1440925) q[1];
sx q[1];
rz(-1.9449241) q[1];
sx q[1];
rz(3.0348026) q[1];
rz(-pi) q[2];
rz(-3.1366523) q[3];
sx q[3];
rz(-0.86588174) q[3];
sx q[3];
rz(-2.4019394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93402702) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-2.8616703) q[2];
sx q[2];
rz(-1.4301849) q[2];
sx q[2];
rz(2.4591597) q[2];
rz(-2.4545112) q[3];
sx q[3];
rz(-2.1344746) q[3];
sx q[3];
rz(1.8070756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
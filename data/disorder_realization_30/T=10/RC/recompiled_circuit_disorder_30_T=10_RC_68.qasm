OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(2.8770652) q[0];
sx q[0];
rz(9.8192083) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(-1.2000097) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.278468) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(1.6186884) q[0];
rz(-pi) q[1];
rz(-2.1985487) q[2];
sx q[2];
rz(-2.606751) q[2];
sx q[2];
rz(-0.56378555) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29772273) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(-2.2754199) q[1];
x q[2];
rz(1.6148189) q[3];
sx q[3];
rz(-1.1929973) q[3];
sx q[3];
rz(-1.4461185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(2.8519894) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(-3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263379) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(-3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.3551691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2750759) q[0];
sx q[0];
rz(-1.5059024) q[0];
sx q[0];
rz(-2.1823723) q[0];
x q[1];
rz(1.4488892) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(0.09300692) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.720571) q[1];
sx q[1];
rz(-2.6877626) q[1];
sx q[1];
rz(1.4245016) q[1];
x q[2];
rz(-1.7543206) q[3];
sx q[3];
rz(-1.6471383) q[3];
sx q[3];
rz(-0.21218382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(0.74869853) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(-2.0764988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0899723) q[0];
sx q[0];
rz(-1.211477) q[0];
sx q[0];
rz(-3.1348455) q[0];
rz(-pi) q[1];
rz(-2.2421354) q[2];
sx q[2];
rz(-0.81842917) q[2];
sx q[2];
rz(0.98086548) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9948431) q[1];
sx q[1];
rz(-0.61770505) q[1];
sx q[1];
rz(-1.6330994) q[1];
x q[2];
rz(-0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(0.041785985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.76434) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-0.011118523) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5549094) q[0];
sx q[0];
rz(-1.1197829) q[0];
sx q[0];
rz(-0.45278544) q[0];
rz(0.63979062) q[2];
sx q[2];
rz(-2.1961229) q[2];
sx q[2];
rz(-1.7196136) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74777714) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(-2.087489) q[1];
rz(-1.3783781) q[3];
sx q[3];
rz(-1.8835861) q[3];
sx q[3];
rz(2.6421412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6461688) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(2.8584976) q[2];
rz(-0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(-2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.8146851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6700232) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(2.3877386) q[0];
rz(-pi) q[1];
rz(-2.5577776) q[2];
sx q[2];
rz(-2.2540255) q[2];
sx q[2];
rz(-0.77606397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2111645) q[1];
sx q[1];
rz(-1.8976364) q[1];
sx q[1];
rz(1.8153166) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96552403) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.999324) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(-1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(2.4601049) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-2.025827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9003446) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(-2.3099398) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7115466) q[2];
sx q[2];
rz(-2.1950245) q[2];
sx q[2];
rz(-0.66832322) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0259077) q[1];
sx q[1];
rz(-0.58870643) q[1];
sx q[1];
rz(0.031850978) q[1];
rz(-pi) q[2];
rz(1.3917771) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(2.9359407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(0.3195233) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091244) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(-1.1122423) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(-2.5792714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494941) q[0];
sx q[0];
rz(-1.5856165) q[0];
sx q[0];
rz(2.2390319) q[0];
x q[1];
rz(-0.33195915) q[2];
sx q[2];
rz(-1.6662285) q[2];
sx q[2];
rz(-1.7711668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.454969) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(-1.7186233) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.087248487) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(2.0283386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(-3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.6202392) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6103219) q[0];
sx q[0];
rz(-1.6390419) q[0];
sx q[0];
rz(0.1971498) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7395758) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(1.6059665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8926881) q[1];
sx q[1];
rz(-2.3612594) q[1];
sx q[1];
rz(-0.08463879) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89150724) q[3];
sx q[3];
rz(-1.4818947) q[3];
sx q[3];
rz(2.8241983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(2.8709581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5198869) q[0];
sx q[0];
rz(-1.6019551) q[0];
sx q[0];
rz(-3.1057538) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80910271) q[2];
sx q[2];
rz(-1.2247694) q[2];
sx q[2];
rz(1.5580633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2808025) q[1];
sx q[1];
rz(-1.6165015) q[1];
sx q[1];
rz(1.526282) q[1];
rz(-pi) q[2];
rz(-0.6587894) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(-0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(-2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(2.5922095) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(-2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-0.49490067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1586944) q[0];
sx q[0];
rz(-1.532908) q[0];
sx q[0];
rz(-1.0844896) q[0];
x q[1];
rz(-3.1206563) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(-2.7431969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7086664) q[1];
sx q[1];
rz(-2.4927944) q[1];
sx q[1];
rz(-1.7556346) q[1];
rz(-pi) q[2];
rz(1.6845735) q[3];
sx q[3];
rz(-2.4416231) q[3];
sx q[3];
rz(-1.9663119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(-0.32785329) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(-1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-1.0275502) q[2];
sx q[2];
rz(-1.7474568) q[2];
sx q[2];
rz(-2.2251868) q[2];
rz(-0.76392382) q[3];
sx q[3];
rz(-0.36119701) q[3];
sx q[3];
rz(1.2154538) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(4.2098213) q[0];
sx q[0];
rz(9.8888483) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-2.6455542) q[0];
sx q[0];
rz(0.92143671) q[0];
x q[1];
rz(-1.60381) q[2];
sx q[2];
rz(-1.0812949) q[2];
sx q[2];
rz(2.2061493) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1844993) q[1];
sx q[1];
rz(-0.59809369) q[1];
sx q[1];
rz(-2.4615272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9908882) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(2.5022262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(0.67392504) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4085061) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(0.79663509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63882534) q[0];
sx q[0];
rz(-0.728038) q[0];
sx q[0];
rz(-2.0730221) q[0];
rz(-2.8393306) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(0.1421393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14256515) q[1];
sx q[1];
rz(-1.5579281) q[1];
sx q[1];
rz(0.60728118) q[1];
rz(-pi) q[2];
rz(2.2315352) q[3];
sx q[3];
rz(-2.4168192) q[3];
sx q[3];
rz(-1.7006601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.18156302) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(2.6484047) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(1.440381) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-0.70297855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36613208) q[0];
sx q[0];
rz(-1.6876939) q[0];
sx q[0];
rz(-1.2890105) q[0];
x q[1];
rz(-1.9269283) q[2];
sx q[2];
rz(-1.0991569) q[2];
sx q[2];
rz(0.65442649) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86537251) q[1];
sx q[1];
rz(-0.11905383) q[1];
sx q[1];
rz(-0.40465506) q[1];
rz(-pi) q[2];
rz(2.1941575) q[3];
sx q[3];
rz(-1.4294335) q[3];
sx q[3];
rz(2.3879104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.253809) q[0];
sx q[0];
rz(-2.1953708) q[0];
sx q[0];
rz(-0.33872351) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2345384) q[2];
sx q[2];
rz(-0.90869892) q[2];
sx q[2];
rz(1.7788356) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.024151) q[1];
sx q[1];
rz(-1.5643331) q[1];
sx q[1];
rz(0.27177377) q[1];
rz(-1.8654278) q[3];
sx q[3];
rz(-2.2664824) q[3];
sx q[3];
rz(2.5254315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(2.9768067) q[2];
rz(0.22848836) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(-1.9794827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1253818) q[0];
sx q[0];
rz(-1.0291161) q[0];
sx q[0];
rz(-2.0340232) q[0];
rz(-0.39761333) q[2];
sx q[2];
rz(-1.0828472) q[2];
sx q[2];
rz(-2.5401126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4967242) q[1];
sx q[1];
rz(-0.86360303) q[1];
sx q[1];
rz(1.0134539) q[1];
rz(-2.4200053) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(1.23502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(-3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-2.390958) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(1.3060588) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40697843) q[0];
sx q[0];
rz(-0.66440551) q[0];
sx q[0];
rz(-2.2218496) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3681709) q[2];
sx q[2];
rz(-1.1154004) q[2];
sx q[2];
rz(0.95603285) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8514511) q[1];
sx q[1];
rz(-0.72039225) q[1];
sx q[1];
rz(0.011522567) q[1];
x q[2];
rz(-1.7197051) q[3];
sx q[3];
rz(-0.9871261) q[3];
sx q[3];
rz(0.35051171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(0.60069096) q[2];
rz(-1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(-2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-1.0466928) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(2.7244862) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313275) q[0];
sx q[0];
rz(-1.350704) q[0];
sx q[0];
rz(-0.025028153) q[0];
rz(-1.276937) q[2];
sx q[2];
rz(-1.2760578) q[2];
sx q[2];
rz(-2.7851832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48737803) q[1];
sx q[1];
rz(-1.1114792) q[1];
sx q[1];
rz(-0.76141255) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7583371) q[3];
sx q[3];
rz(-0.68985046) q[3];
sx q[3];
rz(-0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(2.588429) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(-0.52545351) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32507867) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(-1.4272383) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.534879) q[0];
sx q[0];
rz(-2.0438072) q[0];
sx q[0];
rz(0.19434778) q[0];
x q[1];
rz(0.33371146) q[2];
sx q[2];
rz(-1.6534272) q[2];
sx q[2];
rz(0.73252788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6205794) q[1];
sx q[1];
rz(-2.3304906) q[1];
sx q[1];
rz(2.051342) q[1];
rz(1.1376082) q[3];
sx q[3];
rz(-0.28372753) q[3];
sx q[3];
rz(1.197049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444721) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(-0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-0.91167489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104886) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(-1.2937806) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5981204) q[2];
sx q[2];
rz(-1.6795571) q[2];
sx q[2];
rz(0.47765884) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.72494353) q[1];
sx q[1];
rz(-1.1133725) q[1];
sx q[1];
rz(1.7979421) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0802286) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(2.3232943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-2.8923477) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(-0.58139873) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6045195) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(-1.3634691) q[0];
rz(2.084311) q[2];
sx q[2];
rz(-2.1272749) q[2];
sx q[2];
rz(1.8884115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.33503867) q[1];
sx q[1];
rz(-1.6400669) q[1];
sx q[1];
rz(-0.075242234) q[1];
rz(-2.1595702) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(2.2002937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-0.20467219) q[2];
rz(-1.7278016) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(-1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(0.12116184) q[2];
sx q[2];
rz(-2.0222752) q[2];
sx q[2];
rz(-3.0565699) q[2];
rz(-0.99863573) q[3];
sx q[3];
rz(-1.5012267) q[3];
sx q[3];
rz(-2.5517626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

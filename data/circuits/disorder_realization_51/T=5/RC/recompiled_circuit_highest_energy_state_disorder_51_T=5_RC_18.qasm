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
rz(0.60740745) q[0];
sx q[0];
rz(2.6976801) q[0];
sx q[0];
rz(9.1023268) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(-2.9236887) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8810824) q[0];
sx q[0];
rz(-2.5985607) q[0];
sx q[0];
rz(0.61123993) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71157168) q[2];
sx q[2];
rz(-0.24453881) q[2];
sx q[2];
rz(-2.0681579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3609546) q[1];
sx q[1];
rz(-1.3129957) q[1];
sx q[1];
rz(1.4514489) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8508223) q[3];
sx q[3];
rz(-1.2655846) q[3];
sx q[3];
rz(1.206516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0107062) q[2];
sx q[2];
rz(-2.71038) q[2];
sx q[2];
rz(-2.9475589) q[2];
rz(2.7583097) q[3];
sx q[3];
rz(-0.88519874) q[3];
sx q[3];
rz(-1.6398199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5633157) q[0];
sx q[0];
rz(-1.0551772) q[0];
sx q[0];
rz(-1.1042327) q[0];
rz(-2.4597994) q[1];
sx q[1];
rz(-1.7742523) q[1];
sx q[1];
rz(-1.404748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8786312) q[0];
sx q[0];
rz(-1.6978953) q[0];
sx q[0];
rz(-2.8841579) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66281976) q[2];
sx q[2];
rz(-0.79403764) q[2];
sx q[2];
rz(2.8341849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4841744) q[1];
sx q[1];
rz(-0.60125105) q[1];
sx q[1];
rz(-2.5674393) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6854755) q[3];
sx q[3];
rz(-1.8240098) q[3];
sx q[3];
rz(-2.6270783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76689395) q[2];
sx q[2];
rz(-2.1465492) q[2];
sx q[2];
rz(1.2332756) q[2];
rz(-0.79902664) q[3];
sx q[3];
rz(-0.34501758) q[3];
sx q[3];
rz(-0.12921216) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98053539) q[0];
sx q[0];
rz(-2.0687456) q[0];
sx q[0];
rz(-0.55915731) q[0];
rz(-2.0237538) q[1];
sx q[1];
rz(-1.5562156) q[1];
sx q[1];
rz(-3.0303755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2905544) q[0];
sx q[0];
rz(-1.9469598) q[0];
sx q[0];
rz(-1.4567503) q[0];
rz(2.5170277) q[2];
sx q[2];
rz(-1.3497769) q[2];
sx q[2];
rz(-2.2159037) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.049455482) q[1];
sx q[1];
rz(-0.73284819) q[1];
sx q[1];
rz(-1.8548099) q[1];
rz(-2.7484557) q[3];
sx q[3];
rz(-1.7066147) q[3];
sx q[3];
rz(-0.47752646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3126276) q[2];
sx q[2];
rz(-1.4451005) q[2];
sx q[2];
rz(1.9435389) q[2];
rz(1.7184006) q[3];
sx q[3];
rz(-1.6953902) q[3];
sx q[3];
rz(-0.47553441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6296185) q[0];
sx q[0];
rz(-2.8422575) q[0];
sx q[0];
rz(-0.83628118) q[0];
rz(-1.5760999) q[1];
sx q[1];
rz(-2.0471768) q[1];
sx q[1];
rz(-2.2907168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353514) q[0];
sx q[0];
rz(-1.5374105) q[0];
sx q[0];
rz(0.99793156) q[0];
rz(-pi) q[1];
rz(-0.81340547) q[2];
sx q[2];
rz(-2.4508173) q[2];
sx q[2];
rz(-1.1479614) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5961413) q[1];
sx q[1];
rz(-0.79789466) q[1];
sx q[1];
rz(-1.376838) q[1];
rz(0.20765813) q[3];
sx q[3];
rz(-1.9586143) q[3];
sx q[3];
rz(-2.4352899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1026844) q[2];
sx q[2];
rz(-2.2816198) q[2];
sx q[2];
rz(2.6343708) q[2];
rz(0.66633362) q[3];
sx q[3];
rz(-2.3907876) q[3];
sx q[3];
rz(-0.91608086) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136593) q[0];
sx q[0];
rz(-1.4480696) q[0];
sx q[0];
rz(-1.6130945) q[0];
rz(-1.5799892) q[1];
sx q[1];
rz(-2.0785619) q[1];
sx q[1];
rz(-2.9772421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30475475) q[0];
sx q[0];
rz(-1.4753818) q[0];
sx q[0];
rz(-1.4897904) q[0];
rz(-2.3037203) q[2];
sx q[2];
rz(-0.83014402) q[2];
sx q[2];
rz(1.4655419) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0009672) q[1];
sx q[1];
rz(-0.78636175) q[1];
sx q[1];
rz(-2.0669851) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23779121) q[3];
sx q[3];
rz(-1.2850518) q[3];
sx q[3];
rz(-1.3956436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5932172) q[2];
sx q[2];
rz(-0.84725738) q[2];
sx q[2];
rz(-2.4753921) q[2];
rz(0.14144746) q[3];
sx q[3];
rz(-1.5902767) q[3];
sx q[3];
rz(2.1470054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2675466) q[0];
sx q[0];
rz(-2.5135437) q[0];
sx q[0];
rz(0.48126599) q[0];
rz(-0.73356837) q[1];
sx q[1];
rz(-1.2679408) q[1];
sx q[1];
rz(3.0487294) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5294612) q[0];
sx q[0];
rz(-2.4680637) q[0];
sx q[0];
rz(-2.3271534) q[0];
rz(-3.0888482) q[2];
sx q[2];
rz(-0.03869066) q[2];
sx q[2];
rz(0.82088137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9290849) q[1];
sx q[1];
rz(-1.9783535) q[1];
sx q[1];
rz(-0.69504884) q[1];
rz(-pi) q[2];
rz(0.95488432) q[3];
sx q[3];
rz(-0.43514565) q[3];
sx q[3];
rz(-2.8762115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4197454) q[2];
sx q[2];
rz(-1.5200619) q[2];
sx q[2];
rz(1.9160371) q[2];
rz(-0.078977481) q[3];
sx q[3];
rz(-1.1342528) q[3];
sx q[3];
rz(1.2823766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6135898) q[0];
sx q[0];
rz(-2.4230175) q[0];
sx q[0];
rz(0.80793107) q[0];
rz(-1.6174053) q[1];
sx q[1];
rz(-0.15847358) q[1];
sx q[1];
rz(-0.70404109) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28939279) q[0];
sx q[0];
rz(-1.5796229) q[0];
sx q[0];
rz(1.5509964) q[0];
rz(-pi) q[1];
rz(-1.5624244) q[2];
sx q[2];
rz(-1.5666448) q[2];
sx q[2];
rz(-0.81134398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0320183) q[1];
sx q[1];
rz(-2.5201511) q[1];
sx q[1];
rz(-1.6098609) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1857587) q[3];
sx q[3];
rz(-1.4356239) q[3];
sx q[3];
rz(1.2379902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.070270553) q[2];
sx q[2];
rz(-2.4532048) q[2];
sx q[2];
rz(-2.2570611) q[2];
rz(-2.6279348) q[3];
sx q[3];
rz(-1.1762041) q[3];
sx q[3];
rz(2.6562712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1270444) q[0];
sx q[0];
rz(-2.3764648) q[0];
sx q[0];
rz(2.2807518) q[0];
rz(3.1061106) q[1];
sx q[1];
rz(-1.2018964) q[1];
sx q[1];
rz(1.6532345) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4299773) q[0];
sx q[0];
rz(-1.9693144) q[0];
sx q[0];
rz(-2.7048955) q[0];
x q[1];
rz(-2.9126337) q[2];
sx q[2];
rz(-1.8081852) q[2];
sx q[2];
rz(1.6473532) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79566222) q[1];
sx q[1];
rz(-1.3952858) q[1];
sx q[1];
rz(1.1815168) q[1];
rz(-pi) q[2];
rz(2.866167) q[3];
sx q[3];
rz(-1.4793494) q[3];
sx q[3];
rz(-0.25642727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.615768) q[2];
sx q[2];
rz(-1.1349698) q[2];
sx q[2];
rz(-1.6305249) q[2];
rz(-2.3903971) q[3];
sx q[3];
rz(-0.54371756) q[3];
sx q[3];
rz(0.80120075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57665956) q[0];
sx q[0];
rz(-1.3929921) q[0];
sx q[0];
rz(2.9826953) q[0];
rz(-1.502602) q[1];
sx q[1];
rz(-1.331012) q[1];
sx q[1];
rz(-1.9749036) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5432271) q[0];
sx q[0];
rz(-1.0002213) q[0];
sx q[0];
rz(1.2735418) q[0];
rz(-pi) q[1];
rz(-1.8966271) q[2];
sx q[2];
rz(-1.3052757) q[2];
sx q[2];
rz(2.7139783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37196463) q[1];
sx q[1];
rz(-1.5329016) q[1];
sx q[1];
rz(0.016788646) q[1];
x q[2];
rz(1.6004815) q[3];
sx q[3];
rz(-1.5497713) q[3];
sx q[3];
rz(2.0898745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73929536) q[2];
sx q[2];
rz(-2.4810956) q[2];
sx q[2];
rz(-2.8759549) q[2];
rz(1.9499251) q[3];
sx q[3];
rz(-1.3375514) q[3];
sx q[3];
rz(2.5875097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3091076) q[0];
sx q[0];
rz(-0.6159679) q[0];
sx q[0];
rz(0.68514222) q[0];
rz(-2.5849672) q[1];
sx q[1];
rz(-1.4645422) q[1];
sx q[1];
rz(-2.305078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292824) q[0];
sx q[0];
rz(-2.5959281) q[0];
sx q[0];
rz(-0.87581234) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30843251) q[2];
sx q[2];
rz(-2.0698202) q[2];
sx q[2];
rz(-0.67828808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3565607) q[1];
sx q[1];
rz(-1.1379269) q[1];
sx q[1];
rz(2.4487316) q[1];
rz(-pi) q[2];
rz(-2.0835288) q[3];
sx q[3];
rz(-0.31121507) q[3];
sx q[3];
rz(-1.0070359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4461925) q[2];
sx q[2];
rz(-1.4245028) q[2];
sx q[2];
rz(2.2851473) q[2];
rz(0.19927464) q[3];
sx q[3];
rz(-0.98237413) q[3];
sx q[3];
rz(-0.91163951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0740501) q[0];
sx q[0];
rz(-2.1265125) q[0];
sx q[0];
rz(1.4763005) q[0];
rz(0.53786565) q[1];
sx q[1];
rz(-1.270351) q[1];
sx q[1];
rz(-1.7958633) q[1];
rz(0.79277586) q[2];
sx q[2];
rz(-2.1332827) q[2];
sx q[2];
rz(2.573043) q[2];
rz(2.2945401) q[3];
sx q[3];
rz(-1.77917) q[3];
sx q[3];
rz(2.8875881) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

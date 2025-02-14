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
rz(-2.8430804) q[0];
sx q[0];
rz(-1.2895583) q[0];
sx q[0];
rz(3.1182752) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(-2.4040931) q[1];
sx q[1];
rz(-1.4758543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1064426) q[0];
sx q[0];
rz(-1.5272836) q[0];
sx q[0];
rz(-1.50108) q[0];
rz(-1.2405922) q[2];
sx q[2];
rz(-1.6088952) q[2];
sx q[2];
rz(0.35534258) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45881328) q[1];
sx q[1];
rz(-2.2662357) q[1];
sx q[1];
rz(2.0942445) q[1];
x q[2];
rz(-0.51835097) q[3];
sx q[3];
rz(-2.1795594) q[3];
sx q[3];
rz(-1.6776372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4618571) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(-0.30065817) q[2];
rz(-2.4833208) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7268426) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(-2.9657189) q[0];
rz(0.56100065) q[1];
sx q[1];
rz(-0.70204061) q[1];
sx q[1];
rz(0.19634136) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25717673) q[0];
sx q[0];
rz(-2.0129343) q[0];
sx q[0];
rz(-0.37143645) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11949338) q[2];
sx q[2];
rz(-1.6867078) q[2];
sx q[2];
rz(-1.3136065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2113116) q[1];
sx q[1];
rz(-2.086211) q[1];
sx q[1];
rz(0.0023018259) q[1];
rz(2.7815468) q[3];
sx q[3];
rz(-0.36491766) q[3];
sx q[3];
rz(1.5555738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0520797) q[2];
sx q[2];
rz(-2.5255346) q[2];
sx q[2];
rz(1.6717795) q[2];
rz(-0.52516627) q[3];
sx q[3];
rz(-1.7142121) q[3];
sx q[3];
rz(2.4770881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1216275) q[0];
sx q[0];
rz(-2.2522734) q[0];
sx q[0];
rz(0.63051939) q[0];
rz(-1.1446674) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(-3.0847881) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6358179) q[0];
sx q[0];
rz(-0.56418252) q[0];
sx q[0];
rz(-0.83585541) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88448435) q[2];
sx q[2];
rz(-1.1911885) q[2];
sx q[2];
rz(1.1432709) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5377429) q[1];
sx q[1];
rz(-1.7122388) q[1];
sx q[1];
rz(0.46700041) q[1];
rz(-2.013754) q[3];
sx q[3];
rz(-2.6747799) q[3];
sx q[3];
rz(1.6396322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8593665) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(0.43914208) q[2];
rz(-1.3965083) q[3];
sx q[3];
rz(-1.2625932) q[3];
sx q[3];
rz(2.5631574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0948931) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(-1.0149581) q[0];
rz(-3.0789442) q[1];
sx q[1];
rz(-2.1868314) q[1];
sx q[1];
rz(-2.8573724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64062041) q[0];
sx q[0];
rz(-2.5729502) q[0];
sx q[0];
rz(-0.68822648) q[0];
rz(-pi) q[1];
rz(1.7106871) q[2];
sx q[2];
rz(-1.0770742) q[2];
sx q[2];
rz(-1.3700652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1831467) q[1];
sx q[1];
rz(-1.5827521) q[1];
sx q[1];
rz(-0.35838106) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.013945097) q[3];
sx q[3];
rz(-1.317465) q[3];
sx q[3];
rz(2.0512745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8404428) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(-2.3811049) q[2];
rz(0.31989756) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29192057) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(-2.736295) q[0];
rz(0.48794508) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(-2.778756) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10754171) q[0];
sx q[0];
rz(-2.0011609) q[0];
sx q[0];
rz(1.910188) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9486729) q[2];
sx q[2];
rz(-2.123017) q[2];
sx q[2];
rz(-3.0966126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1048442) q[1];
sx q[1];
rz(-1.3499773) q[1];
sx q[1];
rz(-1.5209496) q[1];
rz(1.2651578) q[3];
sx q[3];
rz(-0.89857093) q[3];
sx q[3];
rz(-1.6068939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1399347) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(-0.72251594) q[2];
rz(-0.51269382) q[3];
sx q[3];
rz(-0.86944681) q[3];
sx q[3];
rz(1.9374013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2672511) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(-2.784506) q[1];
sx q[1];
rz(-1.6386702) q[1];
sx q[1];
rz(1.1579317) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6053033) q[0];
sx q[0];
rz(-2.2380296) q[0];
sx q[0];
rz(2.8298004) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.544417) q[2];
sx q[2];
rz(-0.9998601) q[2];
sx q[2];
rz(1.939029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56933438) q[1];
sx q[1];
rz(-0.63779325) q[1];
sx q[1];
rz(-2.7946154) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80145349) q[3];
sx q[3];
rz(-2.0870001) q[3];
sx q[3];
rz(3.0954697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5623515) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(2.5518899) q[2];
rz(-1.8288745) q[3];
sx q[3];
rz(-1.4785942) q[3];
sx q[3];
rz(-2.5780799) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24484672) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(-0.27031159) q[0];
rz(-2.7145794) q[1];
sx q[1];
rz(-1.3753563) q[1];
sx q[1];
rz(-2.7332773) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9664551) q[0];
sx q[0];
rz(-1.2027367) q[0];
sx q[0];
rz(1.4428663) q[0];
x q[1];
rz(-1.4502268) q[2];
sx q[2];
rz(-0.70115108) q[2];
sx q[2];
rz(2.1599959) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85874365) q[1];
sx q[1];
rz(-0.74902422) q[1];
sx q[1];
rz(-0.94094679) q[1];
rz(2.4189878) q[3];
sx q[3];
rz(-1.6480903) q[3];
sx q[3];
rz(-2.3376436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2030877) q[2];
sx q[2];
rz(-2.0577343) q[2];
sx q[2];
rz(-2.2173524) q[2];
rz(2.3067394) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(2.1448962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860109) q[0];
sx q[0];
rz(-0.38014933) q[0];
sx q[0];
rz(-2.7741449) q[0];
rz(0.5683178) q[1];
sx q[1];
rz(-1.3325997) q[1];
sx q[1];
rz(-0.41499358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85431722) q[0];
sx q[0];
rz(-1.1181896) q[0];
sx q[0];
rz(-0.29083473) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0239564) q[2];
sx q[2];
rz(-1.1073565) q[2];
sx q[2];
rz(-1.5228423) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.69232443) q[1];
sx q[1];
rz(-2.0776761) q[1];
sx q[1];
rz(-2.436422) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1363963) q[3];
sx q[3];
rz(-1.2250161) q[3];
sx q[3];
rz(2.7246662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7830398) q[2];
sx q[2];
rz(-0.8873322) q[2];
sx q[2];
rz(-2.3084194) q[2];
rz(0.35946515) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(1.6370714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9678765) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(-3.1267082) q[0];
rz(-1.124565) q[1];
sx q[1];
rz(-2.896307) q[1];
sx q[1];
rz(0.10717779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5456148) q[0];
sx q[0];
rz(-2.2920485) q[0];
sx q[0];
rz(-0.19907339) q[0];
rz(-pi) q[1];
rz(2.31865) q[2];
sx q[2];
rz(-1.9230033) q[2];
sx q[2];
rz(2.4885881) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1746782) q[1];
sx q[1];
rz(-1.0941097) q[1];
sx q[1];
rz(-1.4033474) q[1];
rz(-pi) q[2];
rz(0.26132432) q[3];
sx q[3];
rz(-1.3991941) q[3];
sx q[3];
rz(1.2444519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(2.2819819) q[2];
rz(0.25935069) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49941007) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(-1.1580178) q[0];
rz(-2.6462789) q[1];
sx q[1];
rz(-1.9751578) q[1];
sx q[1];
rz(-0.29702979) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756145) q[0];
sx q[0];
rz(-2.8225464) q[0];
sx q[0];
rz(-1.5837529) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0115602) q[2];
sx q[2];
rz(-0.30270019) q[2];
sx q[2];
rz(1.5720194) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.90290711) q[1];
sx q[1];
rz(-0.33329646) q[1];
sx q[1];
rz(0.500228) q[1];
rz(-pi) q[2];
rz(-2.2333651) q[3];
sx q[3];
rz(-1.55791) q[3];
sx q[3];
rz(2.9290309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.872252) q[2];
sx q[2];
rz(-2.7597235) q[2];
sx q[2];
rz(-2.2807109) q[2];
rz(-2.6779029) q[3];
sx q[3];
rz(-0.90353614) q[3];
sx q[3];
rz(-1.1369107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1562445) q[0];
sx q[0];
rz(-1.5736268) q[0];
sx q[0];
rz(-1.9851984) q[0];
rz(1.8078177) q[1];
sx q[1];
rz(-1.540624) q[1];
sx q[1];
rz(-2.435138) q[1];
rz(-2.1352519) q[2];
sx q[2];
rz(-2.153572) q[2];
sx q[2];
rz(-1.5700271) q[2];
rz(-1.9577338) q[3];
sx q[3];
rz(-1.1087742) q[3];
sx q[3];
rz(-2.5840989) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

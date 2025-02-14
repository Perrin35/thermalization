OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0526643) q[0];
sx q[0];
rz(-2.0117691) q[0];
sx q[0];
rz(-3.1273754) q[0];
rz(-3.0832503) q[1];
sx q[1];
rz(-0.74906936) q[1];
sx q[1];
rz(-0.60325375) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7503643) q[0];
sx q[0];
rz(-2.3846013) q[0];
sx q[0];
rz(2.6177177) q[0];
rz(3.0809359) q[2];
sx q[2];
rz(-0.91691518) q[2];
sx q[2];
rz(2.4210409) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64390427) q[1];
sx q[1];
rz(-0.7086646) q[1];
sx q[1];
rz(-0.24073118) q[1];
rz(-pi) q[2];
rz(-1.353065) q[3];
sx q[3];
rz(-2.4231632) q[3];
sx q[3];
rz(-1.4361385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3732036) q[2];
sx q[2];
rz(-1.6037805) q[2];
sx q[2];
rz(2.0453889) q[2];
rz(-2.0632035) q[3];
sx q[3];
rz(-0.53467852) q[3];
sx q[3];
rz(2.6064742) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815149) q[0];
sx q[0];
rz(-2.09477) q[0];
sx q[0];
rz(-0.4775508) q[0];
rz(-1.011147) q[1];
sx q[1];
rz(-2.5944581) q[1];
sx q[1];
rz(2.0571041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3318429) q[0];
sx q[0];
rz(-1.4536152) q[0];
sx q[0];
rz(1.940889) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1126419) q[2];
sx q[2];
rz(-1.1626557) q[2];
sx q[2];
rz(-0.52157079) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4762759) q[1];
sx q[1];
rz(-1.3937104) q[1];
sx q[1];
rz(-0.55426532) q[1];
x q[2];
rz(-2.4020477) q[3];
sx q[3];
rz(-1.4149743) q[3];
sx q[3];
rz(2.9062944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61887211) q[2];
sx q[2];
rz(-2.7174157) q[2];
sx q[2];
rz(-1.0901394) q[2];
rz(-2.5777396) q[3];
sx q[3];
rz(-0.68958759) q[3];
sx q[3];
rz(-0.032698154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.76145935) q[0];
sx q[0];
rz(-0.021012336) q[0];
sx q[0];
rz(-1.6047961) q[0];
rz(1.8236632) q[1];
sx q[1];
rz(-2.047796) q[1];
sx q[1];
rz(-0.39753786) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1434162) q[0];
sx q[0];
rz(-0.93188028) q[0];
sx q[0];
rz(0.27235106) q[0];
rz(-pi) q[1];
rz(0.2961646) q[2];
sx q[2];
rz(-2.1284747) q[2];
sx q[2];
rz(0.22685834) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1405903) q[1];
sx q[1];
rz(-2.4574124) q[1];
sx q[1];
rz(-1.160847) q[1];
rz(2.0769172) q[3];
sx q[3];
rz(-1.1017513) q[3];
sx q[3];
rz(-0.367638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2354108) q[2];
sx q[2];
rz(-1.1778888) q[2];
sx q[2];
rz(-2.5664491) q[2];
rz(2.4681674) q[3];
sx q[3];
rz(-1.5410825) q[3];
sx q[3];
rz(1.5448236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6509318) q[0];
sx q[0];
rz(-1.1216102) q[0];
sx q[0];
rz(-0.90890539) q[0];
rz(0.017223651) q[1];
sx q[1];
rz(-0.52809087) q[1];
sx q[1];
rz(3.1294894) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7229268) q[0];
sx q[0];
rz(-1.8932492) q[0];
sx q[0];
rz(0.61005436) q[0];
x q[1];
rz(1.034017) q[2];
sx q[2];
rz(-2.1488948) q[2];
sx q[2];
rz(0.54852099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3084532) q[1];
sx q[1];
rz(-2.1670682) q[1];
sx q[1];
rz(-0.27383974) q[1];
rz(-0.34194591) q[3];
sx q[3];
rz(-0.13445915) q[3];
sx q[3];
rz(0.97245261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0688613) q[2];
sx q[2];
rz(-0.27421811) q[2];
sx q[2];
rz(2.7552628) q[2];
rz(-2.7522411) q[3];
sx q[3];
rz(-1.6081622) q[3];
sx q[3];
rz(-1.0079591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35578457) q[0];
sx q[0];
rz(-0.72440994) q[0];
sx q[0];
rz(0.32387787) q[0];
rz(0.38662275) q[1];
sx q[1];
rz(-1.7522248) q[1];
sx q[1];
rz(1.5547543) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4938439) q[0];
sx q[0];
rz(-0.84058981) q[0];
sx q[0];
rz(-1.6573919) q[0];
rz(2.1065169) q[2];
sx q[2];
rz(-1.8958099) q[2];
sx q[2];
rz(-2.6485788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8530121) q[1];
sx q[1];
rz(-1.1094571) q[1];
sx q[1];
rz(2.4235759) q[1];
x q[2];
rz(2.0915786) q[3];
sx q[3];
rz(-2.2987597) q[3];
sx q[3];
rz(-1.7075552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3891478) q[2];
sx q[2];
rz(-2.2989595) q[2];
sx q[2];
rz(3.0774934) q[2];
rz(-2.5373503) q[3];
sx q[3];
rz(-1.3925545) q[3];
sx q[3];
rz(1.909168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391249) q[0];
sx q[0];
rz(-2.5305643) q[0];
sx q[0];
rz(0.87919277) q[0];
rz(-2.7912256) q[1];
sx q[1];
rz(-1.8354974) q[1];
sx q[1];
rz(-0.15215692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50987303) q[0];
sx q[0];
rz(-1.841396) q[0];
sx q[0];
rz(-2.4218049) q[0];
rz(-pi) q[1];
rz(2.0538267) q[2];
sx q[2];
rz(-1.3994875) q[2];
sx q[2];
rz(-0.47836253) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85206807) q[1];
sx q[1];
rz(-1.7729323) q[1];
sx q[1];
rz(0.44092559) q[1];
x q[2];
rz(-2.426126) q[3];
sx q[3];
rz(-0.56583929) q[3];
sx q[3];
rz(0.51763207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0396314) q[2];
sx q[2];
rz(-2.8122718) q[2];
sx q[2];
rz(-2.9015818) q[2];
rz(0.65252423) q[3];
sx q[3];
rz(-2.0034761) q[3];
sx q[3];
rz(1.8664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0535102) q[0];
sx q[0];
rz(-1.2367915) q[0];
sx q[0];
rz(2.2834593) q[0];
rz(1.3872604) q[1];
sx q[1];
rz(-0.45448449) q[1];
sx q[1];
rz(2.8087356) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80477205) q[0];
sx q[0];
rz(-2.4737278) q[0];
sx q[0];
rz(1.0475545) q[0];
rz(-pi) q[1];
rz(1.5120686) q[2];
sx q[2];
rz(-1.3212122) q[2];
sx q[2];
rz(1.6317612) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12182799) q[1];
sx q[1];
rz(-2.6693444) q[1];
sx q[1];
rz(-0.46544238) q[1];
x q[2];
rz(0.44705963) q[3];
sx q[3];
rz(-0.70526988) q[3];
sx q[3];
rz(1.7763607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7994999) q[2];
sx q[2];
rz(-2.9436593) q[2];
sx q[2];
rz(-1.0510772) q[2];
rz(1.6445232) q[3];
sx q[3];
rz(-1.1727419) q[3];
sx q[3];
rz(2.3433949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368211) q[0];
sx q[0];
rz(-1.0268651) q[0];
sx q[0];
rz(0.89212242) q[0];
rz(0.46961531) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(0.77286744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5839856) q[0];
sx q[0];
rz(-1.0482885) q[0];
sx q[0];
rz(0.64675348) q[0];
x q[1];
rz(0.30420423) q[2];
sx q[2];
rz(-2.4149564) q[2];
sx q[2];
rz(1.3192504) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9443431) q[1];
sx q[1];
rz(-0.57749236) q[1];
sx q[1];
rz(-2.5275699) q[1];
rz(-pi) q[2];
rz(-0.77247932) q[3];
sx q[3];
rz(-0.33145839) q[3];
sx q[3];
rz(-2.9951688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5272687) q[2];
sx q[2];
rz(-1.2246776) q[2];
sx q[2];
rz(-0.69250715) q[2];
rz(2.6912189) q[3];
sx q[3];
rz(-1.9362484) q[3];
sx q[3];
rz(-1.6374121) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5699128) q[0];
sx q[0];
rz(-0.79140651) q[0];
sx q[0];
rz(-2.1929542) q[0];
rz(2.0750849) q[1];
sx q[1];
rz(-0.98943168) q[1];
sx q[1];
rz(-1.8705286) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61425754) q[0];
sx q[0];
rz(-1.7760881) q[0];
sx q[0];
rz(3.0468449) q[0];
x q[1];
rz(2.1461851) q[2];
sx q[2];
rz(-0.040608309) q[2];
sx q[2];
rz(0.81088582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5258057) q[1];
sx q[1];
rz(-0.64103014) q[1];
sx q[1];
rz(1.3435783) q[1];
rz(0.44235559) q[3];
sx q[3];
rz(-1.6887553) q[3];
sx q[3];
rz(-0.30260146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5991685) q[2];
sx q[2];
rz(-1.8343265) q[2];
sx q[2];
rz(-0.14031169) q[2];
rz(3.1091651) q[3];
sx q[3];
rz(-2.2166538) q[3];
sx q[3];
rz(-1.9353346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5133544) q[0];
sx q[0];
rz(-1.4435377) q[0];
sx q[0];
rz(-2.3640609) q[0];
rz(-2.2041722) q[1];
sx q[1];
rz(-1.2317069) q[1];
sx q[1];
rz(0.44313988) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6339119) q[0];
sx q[0];
rz(-0.11735317) q[0];
sx q[0];
rz(-0.76274921) q[0];
rz(2.037077) q[2];
sx q[2];
rz(-2.4906922) q[2];
sx q[2];
rz(-0.69233957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4037343) q[1];
sx q[1];
rz(-0.95721204) q[1];
sx q[1];
rz(-1.1627083) q[1];
rz(-pi) q[2];
rz(2.9786213) q[3];
sx q[3];
rz(-1.2633557) q[3];
sx q[3];
rz(-0.93443894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0810658) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(3.013986) q[2];
rz(0.29868948) q[3];
sx q[3];
rz(-1.9689711) q[3];
sx q[3];
rz(-2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5304607) q[0];
sx q[0];
rz(-1.8831384) q[0];
sx q[0];
rz(0.25181121) q[0];
rz(0.61027377) q[1];
sx q[1];
rz(-2.2563969) q[1];
sx q[1];
rz(1.0856249) q[1];
rz(-1.6444141) q[2];
sx q[2];
rz(-2.4215019) q[2];
sx q[2];
rz(-0.99046594) q[2];
rz(1.1518703) q[3];
sx q[3];
rz(-0.48696951) q[3];
sx q[3];
rz(1.3305668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

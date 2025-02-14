OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17570198) q[0];
sx q[0];
rz(5.0134563) q[0];
sx q[0];
rz(11.391517) q[0];
rz(-1.084561) q[1];
sx q[1];
rz(-1.4718055) q[1];
sx q[1];
rz(-2.5785799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7783053) q[0];
sx q[0];
rz(-1.7666923) q[0];
sx q[0];
rz(1.6250074) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47663759) q[2];
sx q[2];
rz(-0.99299201) q[2];
sx q[2];
rz(3.1057242) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7864816) q[1];
sx q[1];
rz(-1.4110089) q[1];
sx q[1];
rz(0.33336877) q[1];
x q[2];
rz(-2.7607542) q[3];
sx q[3];
rz(-1.9054753) q[3];
sx q[3];
rz(-3.1370441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53434831) q[2];
sx q[2];
rz(-1.4899985) q[2];
sx q[2];
rz(-1.6687757) q[2];
rz(1.2502753) q[3];
sx q[3];
rz(-1.5273124) q[3];
sx q[3];
rz(-0.97243029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8984574) q[0];
sx q[0];
rz(-0.055483015) q[0];
sx q[0];
rz(0.45795101) q[0];
rz(3.0398439) q[1];
sx q[1];
rz(-1.9513756) q[1];
sx q[1];
rz(-2.8153458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1622489) q[0];
sx q[0];
rz(-1.152713) q[0];
sx q[0];
rz(0.3808977) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3497222) q[2];
sx q[2];
rz(-1.862163) q[2];
sx q[2];
rz(-1.86369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40405289) q[1];
sx q[1];
rz(-2.0758325) q[1];
sx q[1];
rz(-2.2801599) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2539034) q[3];
sx q[3];
rz(-2.2884048) q[3];
sx q[3];
rz(-2.9760391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.132167) q[2];
sx q[2];
rz(-2.1475809) q[2];
sx q[2];
rz(1.8780635) q[2];
rz(0.020091232) q[3];
sx q[3];
rz(-0.76954904) q[3];
sx q[3];
rz(-1.2300389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(3.1175213) q[0];
sx q[0];
rz(-3.0776403) q[0];
sx q[0];
rz(-0.57620311) q[0];
rz(-1.4441215) q[1];
sx q[1];
rz(-1.2497808) q[1];
sx q[1];
rz(-1.8796657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44718868) q[0];
sx q[0];
rz(-1.879829) q[0];
sx q[0];
rz(0.033559994) q[0];
x q[1];
rz(2.0327225) q[2];
sx q[2];
rz(-2.0512159) q[2];
sx q[2];
rz(-1.2269208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.539725) q[1];
sx q[1];
rz(-0.66942184) q[1];
sx q[1];
rz(0.4008534) q[1];
rz(-0.71407302) q[3];
sx q[3];
rz(-2.0104109) q[3];
sx q[3];
rz(0.022865828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5835517) q[2];
sx q[2];
rz(-2.5463107) q[2];
sx q[2];
rz(-2.9884647) q[2];
rz(2.2554452) q[3];
sx q[3];
rz(-1.359442) q[3];
sx q[3];
rz(1.9396293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8837638) q[0];
sx q[0];
rz(-1.9572636) q[0];
sx q[0];
rz(2.9828239) q[0];
rz(-1.0348882) q[1];
sx q[1];
rz(-2.1885469) q[1];
sx q[1];
rz(0.75611702) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0949815) q[0];
sx q[0];
rz(-1.5205407) q[0];
sx q[0];
rz(3.0836225) q[0];
x q[1];
rz(1.8240575) q[2];
sx q[2];
rz(-1.54414) q[2];
sx q[2];
rz(2.4566112) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19263527) q[1];
sx q[1];
rz(-1.3657943) q[1];
sx q[1];
rz(-2.8947049) q[1];
x q[2];
rz(1.7773553) q[3];
sx q[3];
rz(-1.3404003) q[3];
sx q[3];
rz(-1.2588866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54255828) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(-0.8503882) q[2];
rz(1.6709857) q[3];
sx q[3];
rz(-1.8149523) q[3];
sx q[3];
rz(0.82408041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4849512) q[0];
sx q[0];
rz(-1.739772) q[0];
sx q[0];
rz(2.9630419) q[0];
rz(-1.8932331) q[1];
sx q[1];
rz(-1.5310042) q[1];
sx q[1];
rz(-2.607883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8781136) q[0];
sx q[0];
rz(-1.7537259) q[0];
sx q[0];
rz(0.35533743) q[0];
rz(-pi) q[1];
rz(2.6421946) q[2];
sx q[2];
rz(-2.207107) q[2];
sx q[2];
rz(0.15313521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0977749) q[1];
sx q[1];
rz(-2.2378451) q[1];
sx q[1];
rz(-0.066867877) q[1];
x q[2];
rz(0.097919252) q[3];
sx q[3];
rz(-1.0684425) q[3];
sx q[3];
rz(0.97508783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8278213) q[2];
sx q[2];
rz(-1.8905819) q[2];
sx q[2];
rz(-1.0531462) q[2];
rz(0.18038067) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(-0.36261305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5798222) q[0];
sx q[0];
rz(-1.0012015) q[0];
sx q[0];
rz(-1.2286105) q[0];
rz(-2.6749532) q[1];
sx q[1];
rz(-1.9231611) q[1];
sx q[1];
rz(3.012588) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6138959) q[0];
sx q[0];
rz(-1.6718699) q[0];
sx q[0];
rz(1.4368426) q[0];
rz(-pi) q[1];
rz(-2.9797709) q[2];
sx q[2];
rz(-2.2941049) q[2];
sx q[2];
rz(-0.79297334) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5552703) q[1];
sx q[1];
rz(-0.96809371) q[1];
sx q[1];
rz(-1.9307677) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0971076) q[3];
sx q[3];
rz(-2.8420181) q[3];
sx q[3];
rz(0.64944525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0064319) q[2];
sx q[2];
rz(-2.5401523) q[2];
sx q[2];
rz(2.637291) q[2];
rz(-2.048061) q[3];
sx q[3];
rz(-1.3239599) q[3];
sx q[3];
rz(2.1849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4920376) q[0];
sx q[0];
rz(-1.858359) q[0];
sx q[0];
rz(1.5851703) q[0];
rz(-2.5632312) q[1];
sx q[1];
rz(-1.2588986) q[1];
sx q[1];
rz(-3.0392821) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226782) q[0];
sx q[0];
rz(-0.75531976) q[0];
sx q[0];
rz(0.97277555) q[0];
rz(-pi) q[1];
x q[1];
rz(1.341171) q[2];
sx q[2];
rz(-1.9582677) q[2];
sx q[2];
rz(1.6796527) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7150517) q[1];
sx q[1];
rz(-0.3449966) q[1];
sx q[1];
rz(-2.662602) q[1];
rz(-pi) q[2];
rz(-2.2874114) q[3];
sx q[3];
rz(-0.9583756) q[3];
sx q[3];
rz(-1.8068594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1954701) q[2];
sx q[2];
rz(-1.2637694) q[2];
sx q[2];
rz(2.5275687) q[2];
rz(-0.92769235) q[3];
sx q[3];
rz(-1.5427019) q[3];
sx q[3];
rz(-0.6774925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55352655) q[0];
sx q[0];
rz(-2.8480242) q[0];
sx q[0];
rz(1.199383) q[0];
rz(-1.9623494) q[1];
sx q[1];
rz(-1.0284547) q[1];
sx q[1];
rz(2.1882122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48483585) q[0];
sx q[0];
rz(-2.0868882) q[0];
sx q[0];
rz(-0.94658312) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0327263) q[2];
sx q[2];
rz(-2.1061387) q[2];
sx q[2];
rz(1.1606729) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39722193) q[1];
sx q[1];
rz(-0.83852856) q[1];
sx q[1];
rz(-0.1485172) q[1];
x q[2];
rz(1.4918421) q[3];
sx q[3];
rz(-2.051226) q[3];
sx q[3];
rz(-1.7460268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.74932468) q[2];
sx q[2];
rz(-0.12568036) q[2];
sx q[2];
rz(0.18088642) q[2];
rz(-2.0285897) q[3];
sx q[3];
rz(-2.1235316) q[3];
sx q[3];
rz(0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0522633) q[0];
sx q[0];
rz(-0.90173975) q[0];
sx q[0];
rz(2.3774636) q[0];
rz(-2.1185421) q[1];
sx q[1];
rz(-1.5307531) q[1];
sx q[1];
rz(-1.4952362) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4039072) q[0];
sx q[0];
rz(-1.7244743) q[0];
sx q[0];
rz(0.30074461) q[0];
rz(-pi) q[1];
rz(0.97025292) q[2];
sx q[2];
rz(-1.8489328) q[2];
sx q[2];
rz(0.32626611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0937522) q[1];
sx q[1];
rz(-1.6421659) q[1];
sx q[1];
rz(-0.38629766) q[1];
rz(1.1782618) q[3];
sx q[3];
rz(-1.9826673) q[3];
sx q[3];
rz(0.28366006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19227795) q[2];
sx q[2];
rz(-1.01769) q[2];
sx q[2];
rz(-2.6940572) q[2];
rz(0.13628515) q[3];
sx q[3];
rz(-2.648573) q[3];
sx q[3];
rz(2.5435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78551453) q[0];
sx q[0];
rz(-2.1165753) q[0];
sx q[0];
rz(-0.4726952) q[0];
rz(-0.39974943) q[1];
sx q[1];
rz(-0.70416299) q[1];
sx q[1];
rz(-0.8185111) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21441701) q[0];
sx q[0];
rz(-2.7821988) q[0];
sx q[0];
rz(2.6136616) q[0];
rz(-pi) q[1];
rz(-2.4651756) q[2];
sx q[2];
rz(-1.2349326) q[2];
sx q[2];
rz(-1.5458653) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0014035066) q[1];
sx q[1];
rz(-1.2480987) q[1];
sx q[1];
rz(3.138209) q[1];
rz(-1.3977526) q[3];
sx q[3];
rz(-1.9132083) q[3];
sx q[3];
rz(0.69998031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5588536) q[2];
sx q[2];
rz(-1.0820729) q[2];
sx q[2];
rz(-1.2782798) q[2];
rz(-1.1996972) q[3];
sx q[3];
rz(-1.4408828) q[3];
sx q[3];
rz(0.27967683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77547913) q[0];
sx q[0];
rz(-1.7178602) q[0];
sx q[0];
rz(-1.3445509) q[0];
rz(1.4115502) q[1];
sx q[1];
rz(-2.330214) q[1];
sx q[1];
rz(-0.65705962) q[1];
rz(3.1362202) q[2];
sx q[2];
rz(-3.0006144) q[2];
sx q[2];
rz(-0.49868546) q[2];
rz(0.33956066) q[3];
sx q[3];
rz(-0.93501285) q[3];
sx q[3];
rz(1.2502678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

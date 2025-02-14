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
rz(-1.7186681) q[0];
sx q[0];
rz(-1.0942425) q[0];
sx q[0];
rz(-2.8835468) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(0.25564495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5244183) q[0];
sx q[0];
rz(-0.99992311) q[0];
sx q[0];
rz(-0.85321315) q[0];
rz(-pi) q[1];
rz(-2.7529703) q[2];
sx q[2];
rz(-1.6291233) q[2];
sx q[2];
rz(-1.3921392) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2110677) q[1];
sx q[1];
rz(-1.3881359) q[1];
sx q[1];
rz(-0.069932368) q[1];
rz(-pi) q[2];
rz(-0.49778865) q[3];
sx q[3];
rz(-1.3827795) q[3];
sx q[3];
rz(-1.9888442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(1.1154491) q[2];
rz(-0.014178064) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(-0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9344591) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(1.2472664) q[0];
rz(2.6994052) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(0.8173379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984262) q[0];
sx q[0];
rz(-1.048537) q[0];
sx q[0];
rz(2.9199615) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1139718) q[2];
sx q[2];
rz(-2.6386542) q[2];
sx q[2];
rz(2.3523112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25774792) q[1];
sx q[1];
rz(-0.93755975) q[1];
sx q[1];
rz(-1.7764938) q[1];
rz(2.3807008) q[3];
sx q[3];
rz(-2.0756654) q[3];
sx q[3];
rz(0.99772108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3680129) q[2];
sx q[2];
rz(-0.37675884) q[2];
sx q[2];
rz(0.22029857) q[2];
rz(-1.1822654) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(0.063145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910864) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(-0.41172045) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(2.128111) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0385376) q[0];
sx q[0];
rz(-0.25347175) q[0];
sx q[0];
rz(-0.059132476) q[0];
x q[1];
rz(2.011843) q[2];
sx q[2];
rz(-1.6905606) q[2];
sx q[2];
rz(2.3074367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0358082) q[1];
sx q[1];
rz(-0.86584751) q[1];
sx q[1];
rz(-2.1501599) q[1];
x q[2];
rz(2.2974422) q[3];
sx q[3];
rz(-1.6682079) q[3];
sx q[3];
rz(-1.0507492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1608405) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(1.8864924) q[2];
rz(2.8694425) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.18144064) q[0];
sx q[0];
rz(-2.20521) q[0];
sx q[0];
rz(-0.61035672) q[0];
rz(2.4329674) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(1.5707387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2529396) q[0];
sx q[0];
rz(-1.4331237) q[0];
sx q[0];
rz(2.987079) q[0];
rz(-2.7566986) q[2];
sx q[2];
rz(-2.8334624) q[2];
sx q[2];
rz(-2.5727814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29250568) q[1];
sx q[1];
rz(-1.6189112) q[1];
sx q[1];
rz(-0.27166697) q[1];
rz(1.9746154) q[3];
sx q[3];
rz(-2.32956) q[3];
sx q[3];
rz(-0.62893516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2107971) q[2];
sx q[2];
rz(-2.1288629) q[2];
sx q[2];
rz(-0.55541682) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.9925502) q[3];
sx q[3];
rz(-0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50437462) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(0.36886886) q[0];
rz(2.8952307) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(1.4341644) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61489366) q[0];
sx q[0];
rz(-1.5656359) q[0];
sx q[0];
rz(2.3847488) q[0];
x q[1];
rz(-2.7857615) q[2];
sx q[2];
rz(-2.6705378) q[2];
sx q[2];
rz(-1.5233153) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1885438) q[1];
sx q[1];
rz(-1.5750139) q[1];
sx q[1];
rz(-0.40178816) q[1];
x q[2];
rz(1.1617817) q[3];
sx q[3];
rz(-2.8944765) q[3];
sx q[3];
rz(-2.1638526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.65115923) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(1.0820214) q[2];
rz(2.3467482) q[3];
sx q[3];
rz(-1.669603) q[3];
sx q[3];
rz(-1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.865888) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(0.24359447) q[0];
rz(0.98681915) q[1];
sx q[1];
rz(-0.74816626) q[1];
sx q[1];
rz(0.93596828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075386062) q[0];
sx q[0];
rz(-2.6200326) q[0];
sx q[0];
rz(-2.2564476) q[0];
rz(-pi) q[1];
rz(-1.9452106) q[2];
sx q[2];
rz(-2.3585359) q[2];
sx q[2];
rz(0.20314344) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82560357) q[1];
sx q[1];
rz(-2.2554419) q[1];
sx q[1];
rz(1.0427834) q[1];
rz(-pi) q[2];
rz(-0.97096393) q[3];
sx q[3];
rz(-1.5884627) q[3];
sx q[3];
rz(2.3223557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.662107) q[2];
sx q[2];
rz(-2.0027436) q[2];
sx q[2];
rz(-0.34995079) q[2];
rz(0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(-0.85132712) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418942) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(-0.30174524) q[0];
rz(-0.31983495) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(2.005827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58590305) q[0];
sx q[0];
rz(-1.6246038) q[0];
sx q[0];
rz(0.64482989) q[0];
x q[1];
rz(0.93923969) q[2];
sx q[2];
rz(-1.9983091) q[2];
sx q[2];
rz(-1.7981426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37138501) q[1];
sx q[1];
rz(-1.6504297) q[1];
sx q[1];
rz(2.1712028) q[1];
x q[2];
rz(-3.0301827) q[3];
sx q[3];
rz(-2.3476331) q[3];
sx q[3];
rz(-0.68796989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7319506) q[2];
sx q[2];
rz(-2.4739517) q[2];
sx q[2];
rz(-2.9988334) q[2];
rz(1.3879294) q[3];
sx q[3];
rz(-1.1889435) q[3];
sx q[3];
rz(-0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9614354) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(0.55602443) q[0];
rz(-0.22932886) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(-1.75846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62263238) q[0];
sx q[0];
rz(-0.83370249) q[0];
sx q[0];
rz(-1.6181158) q[0];
x q[1];
rz(1.6298619) q[2];
sx q[2];
rz(-2.0286273) q[2];
sx q[2];
rz(0.77564592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30547678) q[1];
sx q[1];
rz(-0.51199847) q[1];
sx q[1];
rz(2.192537) q[1];
rz(0.45342584) q[3];
sx q[3];
rz(-1.2229706) q[3];
sx q[3];
rz(-2.3330101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6250299) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(-1.9005091) q[2];
rz(0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(-2.3144498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0144219) q[0];
sx q[0];
rz(-1.4143455) q[0];
sx q[0];
rz(2.993809) q[0];
rz(-2.6328909) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(0.37857372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5145549) q[0];
sx q[0];
rz(-2.6588024) q[0];
sx q[0];
rz(-2.618769) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97855391) q[2];
sx q[2];
rz(-1.8752974) q[2];
sx q[2];
rz(-2.1155807) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7853773) q[1];
sx q[1];
rz(-1.2357724) q[1];
sx q[1];
rz(2.8122693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76102961) q[3];
sx q[3];
rz(-1.4806403) q[3];
sx q[3];
rz(-0.48616274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1754237) q[2];
sx q[2];
rz(-2.1899026) q[2];
sx q[2];
rz(-1.8625205) q[2];
rz(-1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(0.14555791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(1.3928423) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(3.0774935) q[0];
rz(-0.52866689) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(-0.97506964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3885766) q[0];
sx q[0];
rz(-2.574769) q[0];
sx q[0];
rz(-2.6334769) q[0];
rz(-2.1553173) q[2];
sx q[2];
rz(-2.0615675) q[2];
sx q[2];
rz(-1.402439) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8243557) q[1];
sx q[1];
rz(-1.1350766) q[1];
sx q[1];
rz(-1.3292666) q[1];
rz(-pi) q[2];
rz(-0.24419489) q[3];
sx q[3];
rz(-1.7932442) q[3];
sx q[3];
rz(0.73058587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14572445) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(-0.60689849) q[2];
rz(2.8912344) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(-2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470227) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(2.6429214) q[1];
sx q[1];
rz(-1.5892727) q[1];
sx q[1];
rz(-1.7472063) q[1];
rz(-2.278419) q[2];
sx q[2];
rz(-1.0959846) q[2];
sx q[2];
rz(-1.6352996) q[2];
rz(1.7767033) q[3];
sx q[3];
rz(-1.6905224) q[3];
sx q[3];
rz(-0.57991309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

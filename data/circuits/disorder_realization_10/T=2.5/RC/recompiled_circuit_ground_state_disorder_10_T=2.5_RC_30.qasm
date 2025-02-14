OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(-0.40894142) q[0];
sx q[0];
rz(2.1483108) q[0];
rz(-0.14708695) q[1];
sx q[1];
rz(-0.99647254) q[1];
sx q[1];
rz(-1.7239404) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.111747) q[0];
sx q[0];
rz(-2.2897215) q[0];
sx q[0];
rz(-1.0875888) q[0];
rz(2.1314012) q[2];
sx q[2];
rz(-0.94807887) q[2];
sx q[2];
rz(1.6499008) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4864145) q[1];
sx q[1];
rz(-1.298212) q[1];
sx q[1];
rz(0.64719871) q[1];
x q[2];
rz(3.1088022) q[3];
sx q[3];
rz(-1.7049978) q[3];
sx q[3];
rz(2.5658432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19109569) q[2];
sx q[2];
rz(-1.82093) q[2];
sx q[2];
rz(-1.0547868) q[2];
rz(2.6705006) q[3];
sx q[3];
rz(-1.6867009) q[3];
sx q[3];
rz(0.86827046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3926113) q[0];
sx q[0];
rz(-1.6845717) q[0];
sx q[0];
rz(-2.9066322) q[0];
rz(-1.8670392) q[1];
sx q[1];
rz(-1.8320558) q[1];
sx q[1];
rz(0.68769208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4859568) q[0];
sx q[0];
rz(-1.5036426) q[0];
sx q[0];
rz(-0.32819076) q[0];
rz(-0.48414064) q[2];
sx q[2];
rz(-2.2233367) q[2];
sx q[2];
rz(0.86838858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5357757) q[1];
sx q[1];
rz(-1.5652085) q[1];
sx q[1];
rz(1.5582915) q[1];
rz(-pi) q[2];
rz(-3.0491203) q[3];
sx q[3];
rz(-2.1824129) q[3];
sx q[3];
rz(3.0870952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4253652) q[2];
sx q[2];
rz(-1.3024412) q[2];
sx q[2];
rz(0.20935527) q[2];
rz(-0.9730722) q[3];
sx q[3];
rz(-2.2416185) q[3];
sx q[3];
rz(2.5488241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3000325) q[0];
sx q[0];
rz(-2.7356) q[0];
sx q[0];
rz(-5/(11*pi)) q[0];
rz(-1.2380098) q[1];
sx q[1];
rz(-2.369945) q[1];
sx q[1];
rz(-3.0498116) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14427139) q[0];
sx q[0];
rz(-2.1073807) q[0];
sx q[0];
rz(-0.40383213) q[0];
rz(0.068938418) q[2];
sx q[2];
rz(-1.4747915) q[2];
sx q[2];
rz(1.7092741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8615177) q[1];
sx q[1];
rz(-0.28805486) q[1];
sx q[1];
rz(-0.83863284) q[1];
rz(-pi) q[2];
rz(-0.067236891) q[3];
sx q[3];
rz(-1.1394258) q[3];
sx q[3];
rz(-0.19156415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38516513) q[2];
sx q[2];
rz(-1.5084718) q[2];
sx q[2];
rz(0.37919322) q[2];
rz(0.82342974) q[3];
sx q[3];
rz(-2.7354666) q[3];
sx q[3];
rz(2.8520975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.32651153) q[0];
sx q[0];
rz(-2.264475) q[0];
sx q[0];
rz(-0.99863482) q[0];
rz(2.6594992) q[1];
sx q[1];
rz(-2.1908052) q[1];
sx q[1];
rz(-1.3179717) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7319152) q[0];
sx q[0];
rz(-1.240726) q[0];
sx q[0];
rz(0.56185742) q[0];
rz(-pi) q[1];
rz(0.36366574) q[2];
sx q[2];
rz(-0.61700706) q[2];
sx q[2];
rz(-1.7638159) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6819344) q[1];
sx q[1];
rz(-2.1005957) q[1];
sx q[1];
rz(-1.8145241) q[1];
x q[2];
rz(-0.79341268) q[3];
sx q[3];
rz(-1.1856836) q[3];
sx q[3];
rz(-2.957893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92129293) q[2];
sx q[2];
rz(-2.2396294) q[2];
sx q[2];
rz(-0.48545066) q[2];
rz(2.7576533) q[3];
sx q[3];
rz(-1.5555614) q[3];
sx q[3];
rz(0.78401047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081472814) q[0];
sx q[0];
rz(-3.0373242) q[0];
sx q[0];
rz(3.1299348) q[0];
rz(-3.0043789) q[1];
sx q[1];
rz(-1.5882746) q[1];
sx q[1];
rz(-1.8395909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0948167) q[0];
sx q[0];
rz(-1.4621549) q[0];
sx q[0];
rz(-2.957445) q[0];
rz(-2.6259093) q[2];
sx q[2];
rz(-2.1123675) q[2];
sx q[2];
rz(-2.6560419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80655386) q[1];
sx q[1];
rz(-1.8825899) q[1];
sx q[1];
rz(-0.11374082) q[1];
x q[2];
rz(-1.9309031) q[3];
sx q[3];
rz(-1.9269639) q[3];
sx q[3];
rz(2.3115013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7234708) q[2];
sx q[2];
rz(-1.301441) q[2];
sx q[2];
rz(2.8483025) q[2];
rz(-0.063118525) q[3];
sx q[3];
rz(-1.2226353) q[3];
sx q[3];
rz(2.2466834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0224595) q[0];
sx q[0];
rz(-2.8526511) q[0];
sx q[0];
rz(-3.1030848) q[0];
rz(-1.0844082) q[1];
sx q[1];
rz(-2.7190828) q[1];
sx q[1];
rz(2.1748621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9733676) q[0];
sx q[0];
rz(-1.7666139) q[0];
sx q[0];
rz(1.7652579) q[0];
rz(0.013029702) q[2];
sx q[2];
rz(-1.8378864) q[2];
sx q[2];
rz(-2.992127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7679784) q[1];
sx q[1];
rz(-1.2315589) q[1];
sx q[1];
rz(-1.9935196) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8916675) q[3];
sx q[3];
rz(-2.2062613) q[3];
sx q[3];
rz(2.9912586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66095573) q[2];
sx q[2];
rz(-2.9065242) q[2];
sx q[2];
rz(-2.1601775) q[2];
rz(-1.6437982) q[3];
sx q[3];
rz(-1.8794329) q[3];
sx q[3];
rz(-1.5463382) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82866955) q[0];
sx q[0];
rz(-2.4607615) q[0];
sx q[0];
rz(-2.8269826) q[0];
rz(-2.9529052) q[1];
sx q[1];
rz(-0.43470946) q[1];
sx q[1];
rz(0.46868086) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4026827) q[0];
sx q[0];
rz(-2.8651868) q[0];
sx q[0];
rz(2.6457647) q[0];
rz(3.1397083) q[2];
sx q[2];
rz(-1.4030289) q[2];
sx q[2];
rz(-2.014007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0464161) q[1];
sx q[1];
rz(-0.79907387) q[1];
sx q[1];
rz(2.0002736) q[1];
x q[2];
rz(-0.46814274) q[3];
sx q[3];
rz(-1.8624412) q[3];
sx q[3];
rz(-2.1848752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20588747) q[2];
sx q[2];
rz(-2.2998655) q[2];
sx q[2];
rz(-2.1708798) q[2];
rz(-2.6054221) q[3];
sx q[3];
rz(-1.9325117) q[3];
sx q[3];
rz(0.71603388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69990528) q[0];
sx q[0];
rz(-0.55460414) q[0];
sx q[0];
rz(1.4458789) q[0];
rz(1.0384167) q[1];
sx q[1];
rz(-1.1035792) q[1];
sx q[1];
rz(-3.0605002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4202284) q[0];
sx q[0];
rz(-2.2330771) q[0];
sx q[0];
rz(-1.2032979) q[0];
rz(-pi) q[1];
rz(2.7620188) q[2];
sx q[2];
rz(-2.4502769) q[2];
sx q[2];
rz(-0.33502455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75540076) q[1];
sx q[1];
rz(-0.31495783) q[1];
sx q[1];
rz(-1.1909199) q[1];
x q[2];
rz(0.66488336) q[3];
sx q[3];
rz(-0.68782998) q[3];
sx q[3];
rz(2.6874264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9424092) q[2];
sx q[2];
rz(-2.3641391) q[2];
sx q[2];
rz(-0.24277631) q[2];
rz(1.5822004) q[3];
sx q[3];
rz(-2.715761) q[3];
sx q[3];
rz(-1.2529713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9026069) q[0];
sx q[0];
rz(-2.1098397) q[0];
sx q[0];
rz(1.8865939) q[0];
rz(1.7651419) q[1];
sx q[1];
rz(-2.1170728) q[1];
sx q[1];
rz(-2.4741516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7783035) q[0];
sx q[0];
rz(-0.61909715) q[0];
sx q[0];
rz(0.97514345) q[0];
rz(1.6133283) q[2];
sx q[2];
rz(-1.0746403) q[2];
sx q[2];
rz(-2.4537078) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.046958663) q[1];
sx q[1];
rz(-2.486645) q[1];
sx q[1];
rz(1.0818693) q[1];
rz(0.38102229) q[3];
sx q[3];
rz(-1.8748611) q[3];
sx q[3];
rz(-1.1883433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5218375) q[2];
sx q[2];
rz(-2.4027368) q[2];
sx q[2];
rz(2.6832306) q[2];
rz(1.6058263) q[3];
sx q[3];
rz(-1.0777487) q[3];
sx q[3];
rz(-0.88466907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89852029) q[0];
sx q[0];
rz(-0.5721108) q[0];
sx q[0];
rz(2.7761053) q[0];
rz(-0.99994031) q[1];
sx q[1];
rz(-1.702407) q[1];
sx q[1];
rz(-1.4568636) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92577584) q[0];
sx q[0];
rz(-2.121637) q[0];
sx q[0];
rz(-2.9047545) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3028872) q[2];
sx q[2];
rz(-1.678907) q[2];
sx q[2];
rz(-1.2488169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0914336) q[1];
sx q[1];
rz(-0.83178565) q[1];
sx q[1];
rz(1.3350525) q[1];
rz(-pi) q[2];
rz(-2.6014464) q[3];
sx q[3];
rz(-0.77278256) q[3];
sx q[3];
rz(-0.55071044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0008056) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(-1.0055536) q[2];
rz(-0.65070659) q[3];
sx q[3];
rz(-1.4395827) q[3];
sx q[3];
rz(-0.85056359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477796) q[0];
sx q[0];
rz(-0.78421264) q[0];
sx q[0];
rz(-1.0017851) q[0];
rz(1.7306937) q[1];
sx q[1];
rz(-1.7041364) q[1];
sx q[1];
rz(-2.0028353) q[1];
rz(-3.0244556) q[2];
sx q[2];
rz(-1.8979372) q[2];
sx q[2];
rz(-2.7330782) q[2];
rz(-1.6639573) q[3];
sx q[3];
rz(-1.6747083) q[3];
sx q[3];
rz(-1.2887521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

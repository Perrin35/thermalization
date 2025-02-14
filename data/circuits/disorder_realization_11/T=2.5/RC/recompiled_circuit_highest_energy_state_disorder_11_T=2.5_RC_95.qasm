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
rz(-3.0885347) q[0];
sx q[0];
rz(-0.36202708) q[0];
sx q[0];
rz(-1.9757353) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(-1.7133763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4025164) q[0];
sx q[0];
rz(-1.5105643) q[0];
sx q[0];
rz(2.6291558) q[0];
rz(-pi) q[1];
rz(2.8642902) q[2];
sx q[2];
rz(-0.26099527) q[2];
sx q[2];
rz(0.1467315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5184091) q[1];
sx q[1];
rz(-1.3819547) q[1];
sx q[1];
rz(0.13407003) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0936419) q[3];
sx q[3];
rz(-2.6029498) q[3];
sx q[3];
rz(-0.25687309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4521788) q[2];
sx q[2];
rz(-3.1248326) q[2];
sx q[2];
rz(-2.9407799) q[2];
rz(0.14874841) q[3];
sx q[3];
rz(-0.0047618682) q[3];
sx q[3];
rz(-2.8301921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1397322) q[0];
sx q[0];
rz(-2.5480324) q[0];
sx q[0];
rz(1.0396022) q[0];
rz(-0.014558583) q[1];
sx q[1];
rz(-1.9083551) q[1];
sx q[1];
rz(1.5537517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717465) q[0];
sx q[0];
rz(-0.87939191) q[0];
sx q[0];
rz(2.2798674) q[0];
x q[1];
rz(-1.3755685) q[2];
sx q[2];
rz(-0.075610925) q[2];
sx q[2];
rz(1.6603254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.55744104) q[1];
sx q[1];
rz(-1.5366035) q[1];
sx q[1];
rz(-1.8299915) q[1];
rz(-pi) q[2];
rz(-2.8087141) q[3];
sx q[3];
rz(-1.7529704) q[3];
sx q[3];
rz(1.5150439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62850922) q[2];
sx q[2];
rz(-1.6078948) q[2];
sx q[2];
rz(1.3879363) q[2];
rz(-1.3758818) q[3];
sx q[3];
rz(-1.0466156) q[3];
sx q[3];
rz(-2.8468813) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3917711) q[0];
sx q[0];
rz(-2.9048558) q[0];
sx q[0];
rz(2.5340875) q[0];
rz(-1.5982184) q[1];
sx q[1];
rz(-2.9608455) q[1];
sx q[1];
rz(-2.1764596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.799149) q[0];
sx q[0];
rz(-1.5913054) q[0];
sx q[0];
rz(-0.0040702013) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34334646) q[2];
sx q[2];
rz(-2.008524) q[2];
sx q[2];
rz(-2.154532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17441347) q[1];
sx q[1];
rz(-1.7153746) q[1];
sx q[1];
rz(0.036199613) q[1];
rz(-pi) q[2];
x q[2];
rz(0.083250982) q[3];
sx q[3];
rz(-1.4752582) q[3];
sx q[3];
rz(-2.9527309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33660108) q[2];
sx q[2];
rz(-2.470863) q[2];
sx q[2];
rz(2.2750308) q[2];
rz(-2.0349515) q[3];
sx q[3];
rz(-1.5525147) q[3];
sx q[3];
rz(-1.470587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57257819) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(-1.6160075) q[0];
rz(-3.1309879) q[1];
sx q[1];
rz(-0.0037825982) q[1];
sx q[1];
rz(-2.3912281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.01975) q[0];
sx q[0];
rz(-1.6656309) q[0];
sx q[0];
rz(2.4767843) q[0];
rz(-1.343614) q[2];
sx q[2];
rz(-2.1641693) q[2];
sx q[2];
rz(-2.5636755) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.039783104) q[1];
sx q[1];
rz(-1.0652115) q[1];
sx q[1];
rz(-1.702448) q[1];
x q[2];
rz(0.31171215) q[3];
sx q[3];
rz(-1.4501713) q[3];
sx q[3];
rz(-1.8256622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41287199) q[2];
sx q[2];
rz(-1.0888638) q[2];
sx q[2];
rz(-1.2415761) q[2];
rz(-0.0056754644) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(2.2976105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64549696) q[0];
sx q[0];
rz(-3.0912919) q[0];
sx q[0];
rz(2.2182933) q[0];
rz(-0.80054379) q[1];
sx q[1];
rz(-3.1381331) q[1];
sx q[1];
rz(2.961535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411856) q[0];
sx q[0];
rz(-1.6123471) q[0];
sx q[0];
rz(1.8096258) q[0];
rz(-2.4550986) q[2];
sx q[2];
rz(-3.0149934) q[2];
sx q[2];
rz(2.0172269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99111588) q[1];
sx q[1];
rz(-1.2186945) q[1];
sx q[1];
rz(-1.3153418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8534142) q[3];
sx q[3];
rz(-0.37943951) q[3];
sx q[3];
rz(1.5490378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8397612) q[2];
sx q[2];
rz(-1.8415035) q[2];
sx q[2];
rz(-1.7054455) q[2];
rz(-1.2561717) q[3];
sx q[3];
rz(-1.6585645) q[3];
sx q[3];
rz(-3.0459611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65211463) q[0];
sx q[0];
rz(-2.5697932) q[0];
sx q[0];
rz(-2.6982464) q[0];
rz(-1.1706785) q[1];
sx q[1];
rz(-0.0010633855) q[1];
sx q[1];
rz(0.43177691) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5988749) q[0];
sx q[0];
rz(-0.86158991) q[0];
sx q[0];
rz(1.1559691) q[0];
rz(-pi) q[1];
x q[1];
rz(3.058039) q[2];
sx q[2];
rz(-1.7423034) q[2];
sx q[2];
rz(0.87860859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1237632) q[1];
sx q[1];
rz(-0.69849724) q[1];
sx q[1];
rz(0.61459728) q[1];
x q[2];
rz(-1.8006936) q[3];
sx q[3];
rz(-1.2480924) q[3];
sx q[3];
rz(1.475212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7410437) q[2];
sx q[2];
rz(-2.2811175) q[2];
sx q[2];
rz(-1.7575556) q[2];
rz(2.4436229) q[3];
sx q[3];
rz(-2.3845086) q[3];
sx q[3];
rz(2.82011) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8292002) q[0];
sx q[0];
rz(-1.7990524) q[0];
sx q[0];
rz(-2.7685352) q[0];
rz(-0.27698764) q[1];
sx q[1];
rz(-3.1412558) q[1];
sx q[1];
rz(0.75609797) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5146268) q[0];
sx q[0];
rz(-1.1921765) q[0];
sx q[0];
rz(-0.65831229) q[0];
x q[1];
rz(1.2038846) q[2];
sx q[2];
rz(-1.1808625) q[2];
sx q[2];
rz(-0.45118794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71891975) q[1];
sx q[1];
rz(-2.4399099) q[1];
sx q[1];
rz(0.91928457) q[1];
rz(-pi) q[2];
rz(2.8384865) q[3];
sx q[3];
rz(-1.167093) q[3];
sx q[3];
rz(-1.2991326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9296391) q[2];
sx q[2];
rz(-0.55277199) q[2];
sx q[2];
rz(2.2201404) q[2];
rz(-0.067342162) q[3];
sx q[3];
rz(-1.208297) q[3];
sx q[3];
rz(-1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9901554) q[0];
sx q[0];
rz(-2.86148) q[0];
sx q[0];
rz(-0.13023278) q[0];
rz(0.79365927) q[1];
sx q[1];
rz(-3.1399813) q[1];
sx q[1];
rz(0.28009716) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1156688) q[0];
sx q[0];
rz(-1.6243132) q[0];
sx q[0];
rz(-1.491886) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.75861) q[2];
sx q[2];
rz(-0.23442253) q[2];
sx q[2];
rz(1.2146666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7982095) q[1];
sx q[1];
rz(-1.9771132) q[1];
sx q[1];
rz(1.0465606) q[1];
rz(1.0956826) q[3];
sx q[3];
rz(-0.91714782) q[3];
sx q[3];
rz(-0.16330367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.085999504) q[2];
sx q[2];
rz(-1.4693825) q[2];
sx q[2];
rz(2.3386193) q[2];
rz(-1.6064074) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(-2.3101961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11237535) q[0];
sx q[0];
rz(-3.1387098) q[0];
sx q[0];
rz(0.10920864) q[0];
rz(0.373492) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(-2.5901897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5686865) q[0];
sx q[0];
rz(-2.150642) q[0];
sx q[0];
rz(2.2115179) q[0];
rz(-pi) q[1];
rz(0.18420561) q[2];
sx q[2];
rz(-1.4478683) q[2];
sx q[2];
rz(-0.21609989) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.051259816) q[1];
sx q[1];
rz(-0.96104927) q[1];
sx q[1];
rz(2.4734797) q[1];
x q[2];
rz(-1.8165473) q[3];
sx q[3];
rz(-1.9634523) q[3];
sx q[3];
rz(-2.806854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6138844) q[2];
sx q[2];
rz(-1.8324499) q[2];
sx q[2];
rz(-1.8180397) q[2];
rz(-1.2644794) q[3];
sx q[3];
rz(-1.8529961) q[3];
sx q[3];
rz(3.1356139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5962113) q[0];
sx q[0];
rz(-2.5061506) q[0];
sx q[0];
rz(0.7363466) q[0];
rz(2.9296854) q[1];
sx q[1];
rz(-1.0286464) q[1];
sx q[1];
rz(-1.597499) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422393) q[0];
sx q[0];
rz(-2.8618715) q[0];
sx q[0];
rz(-0.065252462) q[0];
rz(0.030363516) q[2];
sx q[2];
rz(-1.5734451) q[2];
sx q[2];
rz(-2.9752258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5527344) q[1];
sx q[1];
rz(-0.97179669) q[1];
sx q[1];
rz(-1.0042648) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5251901) q[3];
sx q[3];
rz(-2.6129301) q[3];
sx q[3];
rz(-1.566949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.836901) q[2];
sx q[2];
rz(-0.83536124) q[2];
sx q[2];
rz(-1.2738073) q[2];
rz(-1.6972208) q[3];
sx q[3];
rz(-0.074051753) q[3];
sx q[3];
rz(-1.4950289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16784167) q[0];
sx q[0];
rz(-1.583562) q[0];
sx q[0];
rz(-1.2927443) q[0];
rz(1.6043067) q[1];
sx q[1];
rz(-2.2289386) q[1];
sx q[1];
rz(-2.9569721) q[1];
rz(-0.040097728) q[2];
sx q[2];
rz(-1.5617149) q[2];
sx q[2];
rz(2.8446773) q[2];
rz(-0.96327412) q[3];
sx q[3];
rz(-1.4701075) q[3];
sx q[3];
rz(-0.95095271) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

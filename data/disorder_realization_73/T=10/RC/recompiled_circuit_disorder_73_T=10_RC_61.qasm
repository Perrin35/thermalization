OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(-0.72416645) q[0];
rz(0.63996285) q[1];
sx q[1];
rz(-0.53007403) q[1];
sx q[1];
rz(-0.78483265) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1278348) q[0];
sx q[0];
rz(-2.7779967) q[0];
sx q[0];
rz(-2.512393) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41854026) q[2];
sx q[2];
rz(-1.6633908) q[2];
sx q[2];
rz(1.6469524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9468294) q[1];
sx q[1];
rz(-1.0461055) q[1];
sx q[1];
rz(2.8835433) q[1];
rz(-1.1562528) q[3];
sx q[3];
rz(-1.539955) q[3];
sx q[3];
rz(2.8950092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(-3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.2215866) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(-2.8979229) q[0];
rz(2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.3557281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49039098) q[0];
sx q[0];
rz(-2.8848007) q[0];
sx q[0];
rz(-1.6780361) q[0];
rz(-pi) q[1];
rz(-2.7849814) q[2];
sx q[2];
rz(-1.2043673) q[2];
sx q[2];
rz(0.29417843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6527378) q[1];
sx q[1];
rz(-1.4421717) q[1];
sx q[1];
rz(-1.1210404) q[1];
x q[2];
rz(1.6173108) q[3];
sx q[3];
rz(-0.97913137) q[3];
sx q[3];
rz(-3.068963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-2.8919354) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(-2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24519414) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(0.89843345) q[0];
rz(1.8067182) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(-1.2737087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35667426) q[0];
sx q[0];
rz(-1.9264364) q[0];
sx q[0];
rz(0.26892923) q[0];
x q[1];
rz(2.8286335) q[2];
sx q[2];
rz(-0.3096146) q[2];
sx q[2];
rz(-1.6329873) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18759218) q[1];
sx q[1];
rz(-2.2224269) q[1];
sx q[1];
rz(-1.1263532) q[1];
x q[2];
rz(-0.7234296) q[3];
sx q[3];
rz(-1.9577141) q[3];
sx q[3];
rz(-0.5667516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0818103) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(-0.95345062) q[2];
rz(3.1070784) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2801441) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-0.074137069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5384597) q[0];
sx q[0];
rz(-2.5939301) q[0];
sx q[0];
rz(-2.0752226) q[0];
rz(-pi) q[1];
rz(2.2494227) q[2];
sx q[2];
rz(-1.8461508) q[2];
sx q[2];
rz(-0.25472578) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.292865) q[1];
sx q[1];
rz(-1.25602) q[1];
sx q[1];
rz(-3.0175812) q[1];
rz(-1.7449964) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(2.1679227) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(1.779153) q[0];
rz(2.3249987) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.978925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8487932) q[0];
sx q[0];
rz(-1.6822364) q[0];
sx q[0];
rz(3.1116027) q[0];
x q[1];
rz(-2.0432203) q[2];
sx q[2];
rz(-2.782151) q[2];
sx q[2];
rz(0.098509468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77046493) q[1];
sx q[1];
rz(-0.41509291) q[1];
sx q[1];
rz(-2.0626555) q[1];
rz(2.707162) q[3];
sx q[3];
rz(-1.0940922) q[3];
sx q[3];
rz(1.7720951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6267307) q[2];
sx q[2];
rz(-1.0802439) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(0.90406117) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(-2.8335559) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338617) q[0];
sx q[0];
rz(-0.15796414) q[0];
sx q[0];
rz(-0.64361848) q[0];
x q[1];
rz(0.61200895) q[2];
sx q[2];
rz(-2.0888121) q[2];
sx q[2];
rz(-2.3305364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6412515) q[1];
sx q[1];
rz(-0.89741035) q[1];
sx q[1];
rz(0.25232368) q[1];
rz(-pi) q[2];
rz(-2.8940053) q[3];
sx q[3];
rz(-2.7831804) q[3];
sx q[3];
rz(-0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77928153) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(1.9936838) q[2];
rz(-2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(2.4961297) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(-0.1299783) q[0];
rz(3.1107483) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(-2.470509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4422465) q[0];
sx q[0];
rz(-1.6244349) q[0];
sx q[0];
rz(-1.6058812) q[0];
rz(-pi) q[1];
rz(-3.1110711) q[2];
sx q[2];
rz(-1.3866716) q[2];
sx q[2];
rz(-2.7437999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56747251) q[1];
sx q[1];
rz(-1.6735055) q[1];
sx q[1];
rz(-1.9436388) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6299134) q[3];
sx q[3];
rz(-2.4654508) q[3];
sx q[3];
rz(-2.6638871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.122763) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(-3.1130062) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.9776483) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(2.7291765) q[0];
rz(-1.4498129) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(1.9746045) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1056571) q[0];
sx q[0];
rz(-1.783839) q[0];
sx q[0];
rz(-0.94353326) q[0];
rz(-pi) q[1];
rz(-1.251986) q[2];
sx q[2];
rz(-0.73482162) q[2];
sx q[2];
rz(-0.97359818) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7653212) q[1];
sx q[1];
rz(-1.5404535) q[1];
sx q[1];
rz(2.9529851) q[1];
rz(2.9917631) q[3];
sx q[3];
rz(-1.0458046) q[3];
sx q[3];
rz(2.6168952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2010487) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(-2.9947301) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(-1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-2.2145859) q[0];
rz(1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(-1.4896726) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.888436) q[0];
sx q[0];
rz(-0.90421593) q[0];
sx q[0];
rz(2.1783834) q[0];
rz(1.6696879) q[2];
sx q[2];
rz(-0.47404587) q[2];
sx q[2];
rz(-2.3941819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5628964) q[1];
sx q[1];
rz(-1.7637196) q[1];
sx q[1];
rz(-1.3955411) q[1];
x q[2];
rz(0.82724039) q[3];
sx q[3];
rz(-2.4628371) q[3];
sx q[3];
rz(2.7620706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5902517) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(-0.57724214) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.5883314) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(2.8531895) q[0];
rz(-2.6092031) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(0.14702252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9592181) q[0];
sx q[0];
rz(-1.9518513) q[0];
sx q[0];
rz(0.94840886) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23994259) q[2];
sx q[2];
rz(-1.6389756) q[2];
sx q[2];
rz(-2.6477674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8204931) q[1];
sx q[1];
rz(-1.5309146) q[1];
sx q[1];
rz(0.95221968) q[1];
rz(-pi) q[2];
rz(1.4149425) q[3];
sx q[3];
rz(-0.999756) q[3];
sx q[3];
rz(1.5172289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0816575) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(-0.75731164) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(2.6894765) q[2];
sx q[2];
rz(-1.5099667) q[2];
sx q[2];
rz(-2.920334) q[2];
rz(-0.13027262) q[3];
sx q[3];
rz(-1.0126922) q[3];
sx q[3];
rz(-1.0425413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
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
rz(-3.0191874) q[0];
sx q[0];
rz(-2.2221017) q[0];
sx q[0];
rz(-1.9424196) q[0];
rz(0.1872669) q[1];
sx q[1];
rz(6.8254677) q[1];
sx q[1];
rz(10.944278) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30795112) q[0];
sx q[0];
rz(-2.184377) q[0];
sx q[0];
rz(-1.475901) q[0];
rz(-pi) q[1];
rz(2.7678732) q[2];
sx q[2];
rz(-1.1454095) q[2];
sx q[2];
rz(-0.56217867) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28827661) q[1];
sx q[1];
rz(-1.9717934) q[1];
sx q[1];
rz(-0.53944352) q[1];
rz(-pi) q[2];
rz(1.6193797) q[3];
sx q[3];
rz(-1.3268927) q[3];
sx q[3];
rz(2.3590306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2970994) q[2];
sx q[2];
rz(-2.5145734) q[2];
sx q[2];
rz(1.7681047) q[2];
rz(0.80396906) q[3];
sx q[3];
rz(-1.6425902) q[3];
sx q[3];
rz(1.1871626) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74541575) q[0];
sx q[0];
rz(-2.4276623) q[0];
sx q[0];
rz(-2.7667238) q[0];
rz(2.6238341) q[1];
sx q[1];
rz(-1.9381783) q[1];
sx q[1];
rz(-1.3688603) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8027975) q[0];
sx q[0];
rz(-1.5707004) q[0];
sx q[0];
rz(3.1406458) q[0];
x q[1];
rz(-0.19646074) q[2];
sx q[2];
rz(-1.5972553) q[2];
sx q[2];
rz(2.485254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9437286) q[1];
sx q[1];
rz(-1.9790589) q[1];
sx q[1];
rz(0.54710435) q[1];
rz(-pi) q[2];
rz(0.13624713) q[3];
sx q[3];
rz(-1.7833424) q[3];
sx q[3];
rz(-0.023879026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.123473) q[2];
sx q[2];
rz(-0.70088434) q[2];
sx q[2];
rz(-0.71462053) q[2];
rz(-0.2054275) q[3];
sx q[3];
rz(-1.3222062) q[3];
sx q[3];
rz(1.2986758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648723) q[0];
sx q[0];
rz(-1.0370075) q[0];
sx q[0];
rz(0.49764693) q[0];
rz(2.5045555) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(3.0348437) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3489368) q[0];
sx q[0];
rz(-1.7829527) q[0];
sx q[0];
rz(1.8769708) q[0];
rz(3.1349772) q[2];
sx q[2];
rz(-0.90108904) q[2];
sx q[2];
rz(-1.0158599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52912583) q[1];
sx q[1];
rz(-2.4262316) q[1];
sx q[1];
rz(-1.6523408) q[1];
x q[2];
rz(-1.2744997) q[3];
sx q[3];
rz(-1.955964) q[3];
sx q[3];
rz(0.22814685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5367624) q[2];
sx q[2];
rz(-0.75211516) q[2];
sx q[2];
rz(0.19482782) q[2];
rz(0.005006494) q[3];
sx q[3];
rz(-0.61401335) q[3];
sx q[3];
rz(2.3495638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1057338) q[0];
sx q[0];
rz(-0.53426131) q[0];
sx q[0];
rz(-2.1015097) q[0];
rz(-0.19373521) q[1];
sx q[1];
rz(-1.5015142) q[1];
sx q[1];
rz(-0.74660444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6116007) q[0];
sx q[0];
rz(-2.009729) q[0];
sx q[0];
rz(-3.0143987) q[0];
rz(-0.87073054) q[2];
sx q[2];
rz(-1.680003) q[2];
sx q[2];
rz(0.93217119) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6317217) q[1];
sx q[1];
rz(-1.6225909) q[1];
sx q[1];
rz(-2.447261) q[1];
x q[2];
rz(-2.901731) q[3];
sx q[3];
rz(-1.4647533) q[3];
sx q[3];
rz(0.21940809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99183434) q[2];
sx q[2];
rz(-1.6712302) q[2];
sx q[2];
rz(2.9370918) q[2];
rz(1.1931984) q[3];
sx q[3];
rz(-0.85719332) q[3];
sx q[3];
rz(0.96113718) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1292773) q[0];
sx q[0];
rz(-0.17450541) q[0];
sx q[0];
rz(-1.0804863) q[0];
rz(0.37733817) q[1];
sx q[1];
rz(-1.5107379) q[1];
sx q[1];
rz(2.2468755) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.13131) q[0];
sx q[0];
rz(-1.2169516) q[0];
sx q[0];
rz(1.3283587) q[0];
rz(-2.4885043) q[2];
sx q[2];
rz(-2.2342367) q[2];
sx q[2];
rz(-1.1764248) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94863588) q[1];
sx q[1];
rz(-2.0564449) q[1];
sx q[1];
rz(-0.57106496) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92833251) q[3];
sx q[3];
rz(-1.8683507) q[3];
sx q[3];
rz(-1.8326108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6696024) q[2];
sx q[2];
rz(-0.44595465) q[2];
sx q[2];
rz(0.10776821) q[2];
rz(0.17555155) q[3];
sx q[3];
rz(-1.2551509) q[3];
sx q[3];
rz(-0.56041437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948697) q[0];
sx q[0];
rz(-2.6304498) q[0];
sx q[0];
rz(-1.0119337) q[0];
rz(-1.6366421) q[1];
sx q[1];
rz(-0.71957809) q[1];
sx q[1];
rz(1.0171657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0295231) q[0];
sx q[0];
rz(-1.3991303) q[0];
sx q[0];
rz(1.8708234) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1332729) q[2];
sx q[2];
rz(-2.2485844) q[2];
sx q[2];
rz(0.00046367292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6222657) q[1];
sx q[1];
rz(-0.19058386) q[1];
sx q[1];
rz(-3.0960073) q[1];
rz(-pi) q[2];
rz(1.8886376) q[3];
sx q[3];
rz(-2.6719465) q[3];
sx q[3];
rz(-1.3263248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73198685) q[2];
sx q[2];
rz(-2.4797532) q[2];
sx q[2];
rz(-2.1448263) q[2];
rz(-0.064920001) q[3];
sx q[3];
rz(-1.4109979) q[3];
sx q[3];
rz(1.9916649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053059269) q[0];
sx q[0];
rz(-0.94045883) q[0];
sx q[0];
rz(0.9084107) q[0];
rz(2.9216595) q[1];
sx q[1];
rz(-1.5426153) q[1];
sx q[1];
rz(1.251108) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.661099) q[0];
sx q[0];
rz(-1.1014043) q[0];
sx q[0];
rz(2.8493657) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32749356) q[2];
sx q[2];
rz(-2.8411189) q[2];
sx q[2];
rz(-0.97285482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0323769) q[1];
sx q[1];
rz(-2.2046979) q[1];
sx q[1];
rz(-0.79297592) q[1];
rz(-pi) q[2];
rz(0.12507579) q[3];
sx q[3];
rz(-0.8440869) q[3];
sx q[3];
rz(-2.7419326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0327586) q[2];
sx q[2];
rz(-2.3075576) q[2];
sx q[2];
rz(-1.6081107) q[2];
rz(1.9882625) q[3];
sx q[3];
rz(-2.0861552) q[3];
sx q[3];
rz(1.4920894) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.15025) q[0];
sx q[0];
rz(-0.38337502) q[0];
sx q[0];
rz(0.94069329) q[0];
rz(3.0062145) q[1];
sx q[1];
rz(-1.4043413) q[1];
sx q[1];
rz(-1.0248331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0198284) q[0];
sx q[0];
rz(-1.262895) q[0];
sx q[0];
rz(2.204049) q[0];
rz(-1.631279) q[2];
sx q[2];
rz(-2.1776458) q[2];
sx q[2];
rz(-1.902193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35541818) q[1];
sx q[1];
rz(-0.23705951) q[1];
sx q[1];
rz(-0.60225822) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5190184) q[3];
sx q[3];
rz(-1.9815784) q[3];
sx q[3];
rz(0.67922986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3465053) q[2];
sx q[2];
rz(-0.89007178) q[2];
sx q[2];
rz(0.65004641) q[2];
rz(2.5551689) q[3];
sx q[3];
rz(-1.4805877) q[3];
sx q[3];
rz(2.39095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745678) q[0];
sx q[0];
rz(-0.50849193) q[0];
sx q[0];
rz(-1.1453999) q[0];
rz(-1.8264495) q[1];
sx q[1];
rz(-1.9086842) q[1];
sx q[1];
rz(-2.4536536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75632771) q[0];
sx q[0];
rz(-1.1996128) q[0];
sx q[0];
rz(0.20275499) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1150949) q[2];
sx q[2];
rz(-2.470068) q[2];
sx q[2];
rz(-1.4948927) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2738139) q[1];
sx q[1];
rz(-1.7325337) q[1];
sx q[1];
rz(-2.383197) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2086772) q[3];
sx q[3];
rz(-2.5543) q[3];
sx q[3];
rz(-1.6577394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.73100662) q[2];
sx q[2];
rz(-2.035391) q[2];
sx q[2];
rz(-2.2753687) q[2];
rz(-2.2335562) q[3];
sx q[3];
rz(-2.3767545) q[3];
sx q[3];
rz(-1.0959371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662125) q[0];
sx q[0];
rz(-1.1840273) q[0];
sx q[0];
rz(-2.7259735) q[0];
rz(-0.79187727) q[1];
sx q[1];
rz(-1.2245347) q[1];
sx q[1];
rz(-0.22722879) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.613406) q[0];
sx q[0];
rz(-3.1216756) q[0];
sx q[0];
rz(2.5100757) q[0];
rz(-pi) q[1];
rz(1.4899026) q[2];
sx q[2];
rz(-2.3757907) q[2];
sx q[2];
rz(-1.3882989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1243678) q[1];
sx q[1];
rz(-1.827888) q[1];
sx q[1];
rz(0.59558792) q[1];
rz(-pi) q[2];
rz(1.0978183) q[3];
sx q[3];
rz(-1.1036345) q[3];
sx q[3];
rz(0.43097365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8772584) q[2];
sx q[2];
rz(-0.72089583) q[2];
sx q[2];
rz(-2.7016675) q[2];
rz(-0.59709966) q[3];
sx q[3];
rz(-0.66949451) q[3];
sx q[3];
rz(-1.3574903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.6482342) q[0];
sx q[0];
rz(-1.5536722) q[0];
sx q[0];
rz(-1.8552725) q[0];
rz(-1.413912) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(-0.4072047) q[2];
sx q[2];
rz(-1.5199678) q[2];
sx q[2];
rz(0.54553568) q[2];
rz(2.3233933) q[3];
sx q[3];
rz(-1.7530734) q[3];
sx q[3];
rz(-1.659041) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

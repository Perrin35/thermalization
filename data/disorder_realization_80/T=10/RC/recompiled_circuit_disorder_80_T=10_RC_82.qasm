OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(2.9449129) q[0];
sx q[0];
rz(11.37698) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(2.300188) q[1];
sx q[1];
rz(9.0471164) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6807251) q[0];
sx q[0];
rz(-1.6922249) q[0];
sx q[0];
rz(-0.073343883) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6701277) q[2];
sx q[2];
rz(-0.28495312) q[2];
sx q[2];
rz(3.1379267) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.764896) q[1];
sx q[1];
rz(-0.88250676) q[1];
sx q[1];
rz(2.1535758) q[1];
rz(-pi) q[2];
rz(-0.30794413) q[3];
sx q[3];
rz(-0.43711284) q[3];
sx q[3];
rz(-1.8969769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.010014023) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(2.4689891) q[2];
rz(0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(-2.9130274) q[0];
rz(0.16054343) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(2.8536318) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8156793) q[0];
sx q[0];
rz(-2.8796112) q[0];
sx q[0];
rz(-2.322305) q[0];
rz(-pi) q[1];
rz(-2.7520913) q[2];
sx q[2];
rz(-2.1283538) q[2];
sx q[2];
rz(2.8745289) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49111734) q[1];
sx q[1];
rz(-2.3820291) q[1];
sx q[1];
rz(2.526545) q[1];
x q[2];
rz(-3.037022) q[3];
sx q[3];
rz(-0.35211709) q[3];
sx q[3];
rz(-0.62192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(1.4734369) q[2];
rz(2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(2.766585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32524747) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-2.8741799) q[0];
rz(-1.422241) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(-2.1898988) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62943447) q[0];
sx q[0];
rz(-0.21244563) q[0];
sx q[0];
rz(-2.0149219) q[0];
rz(-1.1859602) q[2];
sx q[2];
rz(-2.5927612) q[2];
sx q[2];
rz(1.8897111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3141331) q[1];
sx q[1];
rz(-2.3312807) q[1];
sx q[1];
rz(2.00287) q[1];
rz(-0.88055196) q[3];
sx q[3];
rz(-0.75062245) q[3];
sx q[3];
rz(2.7504138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13016985) q[2];
sx q[2];
rz(-1.6458076) q[2];
sx q[2];
rz(-2.7974131) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(-1.244506) q[0];
rz(0.1291153) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(0.37277645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0817954) q[0];
sx q[0];
rz(-2.3764388) q[0];
sx q[0];
rz(-2.9761936) q[0];
rz(-1.9984841) q[2];
sx q[2];
rz(-0.83979411) q[2];
sx q[2];
rz(-1.0243624) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4200538) q[1];
sx q[1];
rz(-2.5905847) q[1];
sx q[1];
rz(-1.0934456) q[1];
x q[2];
rz(2.0471441) q[3];
sx q[3];
rz(-1.8653449) q[3];
sx q[3];
rz(1.095872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(2.2367031) q[2];
rz(2.7010226) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47914094) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.8537846) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(2.6073661) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68191093) q[0];
sx q[0];
rz(-0.43897438) q[0];
sx q[0];
rz(0.79553332) q[0];
rz(2.062837) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(-0.505503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.6492669) q[1];
sx q[1];
rz(-1.6345134) q[1];
sx q[1];
rz(-1.4953514) q[1];
rz(-pi) q[2];
rz(-0.71877919) q[3];
sx q[3];
rz(-1.6703509) q[3];
sx q[3];
rz(-1.0574785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(-1.7701373) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(-0.21970704) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4390398) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-2.904536) q[0];
rz(-1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(1.7664849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9000589) q[0];
sx q[0];
rz(-2.7068479) q[0];
sx q[0];
rz(-0.66381201) q[0];
x q[1];
rz(-1.6429971) q[2];
sx q[2];
rz(-0.82319665) q[2];
sx q[2];
rz(-2.2526134) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68723893) q[1];
sx q[1];
rz(-2.7308186) q[1];
sx q[1];
rz(3.0803068) q[1];
x q[2];
rz(-0.70011219) q[3];
sx q[3];
rz(-0.60986076) q[3];
sx q[3];
rz(1.8240579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53283006) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(2.9928845) q[2];
rz(-0.016629774) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6454813) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(0.87316978) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-1.1154122) q[1];
sx q[1];
rz(-0.08392863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5556363) q[0];
sx q[0];
rz(-1.3480098) q[0];
sx q[0];
rz(2.1369834) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7883045) q[2];
sx q[2];
rz(-1.366426) q[2];
sx q[2];
rz(-0.7962966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88372701) q[1];
sx q[1];
rz(-2.1818433) q[1];
sx q[1];
rz(-2.3962767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.300802) q[3];
sx q[3];
rz(-1.6222686) q[3];
sx q[3];
rz(1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7579047) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(2.596358) q[2];
rz(2.7111354) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(1.2114725) q[0];
rz(-2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0154008) q[0];
sx q[0];
rz(-1.3500824) q[0];
sx q[0];
rz(-1.6507694) q[0];
x q[1];
rz(-2.321645) q[2];
sx q[2];
rz(-0.87265271) q[2];
sx q[2];
rz(1.1993711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1434053) q[1];
sx q[1];
rz(-0.37969509) q[1];
sx q[1];
rz(3.0642964) q[1];
rz(-pi) q[2];
rz(0.66611992) q[3];
sx q[3];
rz(-1.564581) q[3];
sx q[3];
rz(2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(-2.5366606) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873001) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(2.5262685) q[0];
rz(0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(-2.5659134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3424123) q[0];
sx q[0];
rz(-0.77501446) q[0];
sx q[0];
rz(2.4393625) q[0];
x q[1];
rz(-1.766878) q[2];
sx q[2];
rz(-1.4596645) q[2];
sx q[2];
rz(2.9642504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5513735) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(-1.3361841) q[1];
rz(-pi) q[2];
rz(-1.0300893) q[3];
sx q[3];
rz(-1.4095777) q[3];
sx q[3];
rz(-2.2459523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(-0.021961948) q[2];
rz(2.9711376) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4158674) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(-3.1126853) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(-0.96910563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983738) q[0];
sx q[0];
rz(-1.6341097) q[0];
sx q[0];
rz(-1.6427342) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74659851) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(1.5310841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9358878) q[1];
sx q[1];
rz(-1.3921326) q[1];
sx q[1];
rz(-0.67358576) q[1];
x q[2];
rz(2.7006847) q[3];
sx q[3];
rz(-1.8280067) q[3];
sx q[3];
rz(2.5901026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9344899) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-3.0974401) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(-1.953223) q[2];
sx q[2];
rz(-1.8845176) q[2];
sx q[2];
rz(0.11558576) q[2];
rz(0.58892693) q[3];
sx q[3];
rz(-0.52994655) q[3];
sx q[3];
rz(0.51701057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

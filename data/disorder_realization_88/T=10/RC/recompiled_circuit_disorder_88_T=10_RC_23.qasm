OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(2.8136301) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(0.091436401) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51973242) q[0];
sx q[0];
rz(-2.41142) q[0];
sx q[0];
rz(2.5229088) q[0];
rz(-pi) q[1];
rz(0.47317998) q[2];
sx q[2];
rz(-2.8521529) q[2];
sx q[2];
rz(0.48130408) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4788471) q[1];
sx q[1];
rz(-1.4884243) q[1];
sx q[1];
rz(1.3326416) q[1];
rz(0.9760194) q[3];
sx q[3];
rz(-1.3888437) q[3];
sx q[3];
rz(-2.6973157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66449195) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-1.1260024) q[2];
rz(-0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(-2.0454848) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(0.19031659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7305671) q[0];
sx q[0];
rz(-1.1682296) q[0];
sx q[0];
rz(-2.9257724) q[0];
rz(-0.61043592) q[2];
sx q[2];
rz(-2.5154841) q[2];
sx q[2];
rz(-2.6420643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9840934) q[1];
sx q[1];
rz(-1.5484973) q[1];
sx q[1];
rz(-0.84866546) q[1];
rz(-pi) q[2];
rz(2.042949) q[3];
sx q[3];
rz(-2.241579) q[3];
sx q[3];
rz(-0.19876476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(-0.1427342) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74293566) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(2.8495158) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4688063) q[0];
sx q[0];
rz(-2.5350223) q[0];
sx q[0];
rz(-0.16288217) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5065932) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.1484255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6229912) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(1.6559421) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1850584) q[3];
sx q[3];
rz(-0.68813656) q[3];
sx q[3];
rz(-3.0582173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32039207) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.9281663) q[2];
rz(-2.976867) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(-0.18606342) q[0];
rz(-2.9371254) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870677) q[0];
sx q[0];
rz(-1.5274807) q[0];
sx q[0];
rz(1.7012419) q[0];
x q[1];
rz(0.16764955) q[2];
sx q[2];
rz(-1.8561346) q[2];
sx q[2];
rz(-2.167785) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0399196) q[1];
sx q[1];
rz(-2.1664201) q[1];
sx q[1];
rz(0.69570978) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1277539) q[3];
sx q[3];
rz(-1.3912364) q[3];
sx q[3];
rz(-1.9765215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(0.32133189) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(-1.2478158) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(-2.2633973) q[0];
rz(1.325266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(-1.7153046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8910687) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(-3.0380681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9274474) q[2];
sx q[2];
rz(-0.38380917) q[2];
sx q[2];
rz(-2.5197033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6225699) q[1];
sx q[1];
rz(-1.6791324) q[1];
sx q[1];
rz(-2.3035754) q[1];
rz(-0.32894965) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(2.9165099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3549266) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(-2.1110995) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(2.5700263) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8067236) q[0];
sx q[0];
rz(-2.3248701) q[0];
sx q[0];
rz(-2.4094894) q[0];
rz(-pi) q[1];
rz(-1.7320485) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(1.4816928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.188365) q[1];
sx q[1];
rz(-2.0247012) q[1];
sx q[1];
rz(2.8737349) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70116331) q[3];
sx q[3];
rz(-2.1950766) q[3];
sx q[3];
rz(0.13247709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(-0.56419939) q[2];
rz(3.0155904) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(-0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(-2.4760822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8004566) q[0];
sx q[0];
rz(-1.3285713) q[0];
sx q[0];
rz(3.1405297) q[0];
rz(1.3940784) q[2];
sx q[2];
rz(-0.70509796) q[2];
sx q[2];
rz(-1.9213898) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29384089) q[1];
sx q[1];
rz(-1.7898702) q[1];
sx q[1];
rz(-2.0391383) q[1];
rz(0.14466488) q[3];
sx q[3];
rz(-1.33107) q[3];
sx q[3];
rz(-3.1160115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(1.4228014) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(-0.62675369) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(2.802882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20913798) q[0];
sx q[0];
rz(-1.4203686) q[0];
sx q[0];
rz(-1.223279) q[0];
x q[1];
rz(-2.68967) q[2];
sx q[2];
rz(-1.7799313) q[2];
sx q[2];
rz(-0.6349596) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.056811) q[1];
sx q[1];
rz(-1.6817131) q[1];
sx q[1];
rz(2.0463498) q[1];
rz(-1.2440153) q[3];
sx q[3];
rz(-1.8339001) q[3];
sx q[3];
rz(-2.9941878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64615858) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-2.4411566) q[2];
rz(-2.2436079) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(0.65834808) q[0];
rz(-2.530653) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(0.13959612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9717279) q[0];
sx q[0];
rz(-1.5507878) q[0];
sx q[0];
rz(1.5097029) q[0];
x q[1];
rz(0.058768674) q[2];
sx q[2];
rz(-0.81548703) q[2];
sx q[2];
rz(1.8975443) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2792564) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(-0.52211296) q[1];
rz(-pi) q[2];
rz(2.2152882) q[3];
sx q[3];
rz(-1.6127805) q[3];
sx q[3];
rz(0.55461649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6167986) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(-0.75774276) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(-3.0537135) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(-0.18558003) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.4846444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907222) q[0];
sx q[0];
rz(-1.301268) q[0];
sx q[0];
rz(2.0275293) q[0];
rz(-2.9160203) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(0.05664209) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6757322) q[1];
sx q[1];
rz(-1.4714186) q[1];
sx q[1];
rz(-1.9468716) q[1];
rz(-1.5649892) q[3];
sx q[3];
rz(-2.4366637) q[3];
sx q[3];
rz(-2.4095636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-0.55220848) q[2];
rz(0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26375297) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(2.6681343) q[2];
sx q[2];
rz(-2.8291694) q[2];
sx q[2];
rz(1.3419801) q[2];
rz(-2.4545112) q[3];
sx q[3];
rz(-2.1344746) q[3];
sx q[3];
rz(1.8070756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

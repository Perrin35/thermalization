OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(-0.52596337) q[0];
sx q[0];
rz(-2.9232803) q[0];
rz(1.4766308) q[1];
sx q[1];
rz(-2.6638439) q[1];
sx q[1];
rz(-0.55396095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6787348) q[0];
sx q[0];
rz(-0.7174055) q[0];
sx q[0];
rz(2.3166902) q[0];
x q[1];
rz(-0.5289107) q[2];
sx q[2];
rz(-1.5195091) q[2];
sx q[2];
rz(1.2085268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3240149) q[1];
sx q[1];
rz(-0.87000532) q[1];
sx q[1];
rz(0.27591095) q[1];
rz(2.6036156) q[3];
sx q[3];
rz(-1.2012641) q[3];
sx q[3];
rz(1.6224976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37609425) q[2];
sx q[2];
rz(-0.29355294) q[2];
sx q[2];
rz(1.2676839) q[2];
rz(2.3119161) q[3];
sx q[3];
rz(-1.5209578) q[3];
sx q[3];
rz(-0.42384306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.053452881) q[0];
sx q[0];
rz(-1.3351048) q[0];
sx q[0];
rz(-1.394519) q[0];
rz(2.246619) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(1.5825533) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532063) q[0];
sx q[0];
rz(-1.0307285) q[0];
sx q[0];
rz(0.77629838) q[0];
rz(-2.7855273) q[2];
sx q[2];
rz(-1.1741956) q[2];
sx q[2];
rz(-0.44677904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8445721) q[1];
sx q[1];
rz(-2.3851042) q[1];
sx q[1];
rz(-0.46539657) q[1];
rz(-0.96554324) q[3];
sx q[3];
rz(-0.50600921) q[3];
sx q[3];
rz(3.0288896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6643657) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(0.030755432) q[2];
rz(2.6886046) q[3];
sx q[3];
rz(-2.9121297) q[3];
sx q[3];
rz(-0.17791137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396616) q[0];
sx q[0];
rz(-3.0406096) q[0];
sx q[0];
rz(-0.82823753) q[0];
rz(-3.0939057) q[1];
sx q[1];
rz(-2.2779155) q[1];
sx q[1];
rz(-1.9140859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9253474) q[0];
sx q[0];
rz(-1.1309393) q[0];
sx q[0];
rz(-0.60103215) q[0];
x q[1];
rz(1.1647878) q[2];
sx q[2];
rz(-1.7696524) q[2];
sx q[2];
rz(2.5420497) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3674254) q[1];
sx q[1];
rz(-1.8313421) q[1];
sx q[1];
rz(-1.114964) q[1];
rz(-pi) q[2];
rz(1.5121721) q[3];
sx q[3];
rz(-0.72949648) q[3];
sx q[3];
rz(2.5945455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0487655) q[2];
sx q[2];
rz(-2.2266677) q[2];
sx q[2];
rz(-2.0406593) q[2];
rz(-0.01384211) q[3];
sx q[3];
rz(-1.3661386) q[3];
sx q[3];
rz(2.3069416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.58532995) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(-2.8072667) q[0];
rz(0.73257929) q[1];
sx q[1];
rz(-0.94894797) q[1];
sx q[1];
rz(3.0308731) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77483141) q[0];
sx q[0];
rz(-0.80935055) q[0];
sx q[0];
rz(-1.5917718) q[0];
x q[1];
rz(-2.6312073) q[2];
sx q[2];
rz(-0.18922999) q[2];
sx q[2];
rz(0.58923474) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1724989) q[1];
sx q[1];
rz(-0.90414541) q[1];
sx q[1];
rz(-2.4039925) q[1];
rz(2.3663051) q[3];
sx q[3];
rz(-1.2167769) q[3];
sx q[3];
rz(-0.86590761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.4577786) q[2];
sx q[2];
rz(-0.28202287) q[2];
sx q[2];
rz(1.4078183) q[2];
rz(1.1553361) q[3];
sx q[3];
rz(-1.1563533) q[3];
sx q[3];
rz(-0.19481625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69283501) q[0];
sx q[0];
rz(-2.223707) q[0];
sx q[0];
rz(-0.99739972) q[0];
rz(1.6150486) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(1.7220727) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8527232) q[0];
sx q[0];
rz(-2.8670681) q[0];
sx q[0];
rz(-2.2122266) q[0];
rz(2.3676374) q[2];
sx q[2];
rz(-1.4593235) q[2];
sx q[2];
rz(2.39448) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0816272) q[1];
sx q[1];
rz(-0.94725376) q[1];
sx q[1];
rz(2.7452031) q[1];
rz(3.0180198) q[3];
sx q[3];
rz(-1.5383178) q[3];
sx q[3];
rz(-2.6952254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.851696) q[2];
sx q[2];
rz(-3.0471314) q[2];
sx q[2];
rz(2.8186901) q[2];
rz(1.1139392) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(0.41306257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776529) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(0.80663484) q[0];
rz(2.5576162) q[1];
sx q[1];
rz(-2.0239425) q[1];
sx q[1];
rz(-1.9516099) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89878824) q[0];
sx q[0];
rz(-1.4477647) q[0];
sx q[0];
rz(1.2205628) q[0];
x q[1];
rz(2.6072845) q[2];
sx q[2];
rz(-1.3883841) q[2];
sx q[2];
rz(1.4774496) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8998225) q[1];
sx q[1];
rz(-2.9651151) q[1];
sx q[1];
rz(1.2736257) q[1];
x q[2];
rz(-1.7403931) q[3];
sx q[3];
rz(-0.84753321) q[3];
sx q[3];
rz(-2.4270428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40818647) q[2];
sx q[2];
rz(-1.5032282) q[2];
sx q[2];
rz(-2.9564986) q[2];
rz(-1.6144276) q[3];
sx q[3];
rz(-1.4027169) q[3];
sx q[3];
rz(0.68814284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1899034) q[0];
sx q[0];
rz(-0.95174319) q[0];
sx q[0];
rz(0.21251799) q[0];
rz(-2.3566133) q[1];
sx q[1];
rz(-1.6467983) q[1];
sx q[1];
rz(0.77883887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3991645) q[0];
sx q[0];
rz(-2.140283) q[0];
sx q[0];
rz(0.67230255) q[0];
rz(-1.1981702) q[2];
sx q[2];
rz(-2.1261132) q[2];
sx q[2];
rz(0.057387847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6952371) q[1];
sx q[1];
rz(-1.8255705) q[1];
sx q[1];
rz(-0.22344113) q[1];
x q[2];
rz(-2.1206585) q[3];
sx q[3];
rz(-2.3519197) q[3];
sx q[3];
rz(2.0146973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51398858) q[2];
sx q[2];
rz(-2.6476761) q[2];
sx q[2];
rz(-0.40840515) q[2];
rz(-2.4397395) q[3];
sx q[3];
rz(-0.93188325) q[3];
sx q[3];
rz(-1.2592038) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7931165) q[0];
sx q[0];
rz(-1.2154673) q[0];
sx q[0];
rz(3.075573) q[0];
rz(-1.5090212) q[1];
sx q[1];
rz(-2.0965818) q[1];
sx q[1];
rz(0.95796934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45487472) q[0];
sx q[0];
rz(-1.8544079) q[0];
sx q[0];
rz(0.4768178) q[0];
rz(1.0220549) q[2];
sx q[2];
rz(-0.83344747) q[2];
sx q[2];
rz(-1.1481783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1537452) q[1];
sx q[1];
rz(-1.8829131) q[1];
sx q[1];
rz(-0.69674833) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10029467) q[3];
sx q[3];
rz(-1.7236606) q[3];
sx q[3];
rz(-2.4748442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8336739) q[2];
sx q[2];
rz(-1.6180399) q[2];
sx q[2];
rz(-1.7306805) q[2];
rz(-2.446567) q[3];
sx q[3];
rz(-1.6618988) q[3];
sx q[3];
rz(0.01785774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0628292) q[0];
sx q[0];
rz(-1.6595027) q[0];
sx q[0];
rz(2.1516946) q[0];
rz(2.6784189) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(1.90082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1027066) q[0];
sx q[0];
rz(-0.30065824) q[0];
sx q[0];
rz(-1.3749529) q[0];
rz(-pi) q[1];
rz(1.979901) q[2];
sx q[2];
rz(-0.4767524) q[2];
sx q[2];
rz(-2.7331405) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6077784) q[1];
sx q[1];
rz(-1.3097714) q[1];
sx q[1];
rz(-0.29284524) q[1];
x q[2];
rz(-0.97384805) q[3];
sx q[3];
rz(-2.429109) q[3];
sx q[3];
rz(-1.2013916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3339633) q[2];
sx q[2];
rz(-1.3796076) q[2];
sx q[2];
rz(-0.71869746) q[2];
rz(1.1527609) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(-1.0144455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5929247) q[0];
sx q[0];
rz(-0.34807006) q[0];
sx q[0];
rz(-0.80192178) q[0];
rz(-2.0536664) q[1];
sx q[1];
rz(-1.8049003) q[1];
sx q[1];
rz(-2.4235639) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.021026) q[0];
sx q[0];
rz(-0.88827288) q[0];
sx q[0];
rz(-1.1521856) q[0];
rz(-0.74943351) q[2];
sx q[2];
rz(-0.66808703) q[2];
sx q[2];
rz(-1.4632478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5637569) q[1];
sx q[1];
rz(-0.33591649) q[1];
sx q[1];
rz(-1.4116686) q[1];
rz(-pi) q[2];
rz(2.2949785) q[3];
sx q[3];
rz(-2.1606956) q[3];
sx q[3];
rz(-2.5835832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7903018) q[2];
sx q[2];
rz(-0.21778926) q[2];
sx q[2];
rz(-0.51266074) q[2];
rz(2.3971108) q[3];
sx q[3];
rz(-1.0267886) q[3];
sx q[3];
rz(0.7640394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9164593) q[0];
sx q[0];
rz(-1.2703348) q[0];
sx q[0];
rz(-1.2402007) q[0];
rz(-0.71612877) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(1.345558) q[2];
sx q[2];
rz(-1.6788531) q[2];
sx q[2];
rz(-2.0137871) q[2];
rz(1.2342831) q[3];
sx q[3];
rz(-1.4740305) q[3];
sx q[3];
rz(-2.8421845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

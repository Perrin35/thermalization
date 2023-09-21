OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479552) q[0];
sx q[0];
rz(-1.4076828) q[0];
sx q[0];
rz(-1.9440584) q[0];
rz(-pi) q[1];
rz(-1.6127365) q[2];
sx q[2];
rz(-1.1239927) q[2];
sx q[2];
rz(0.97181335) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4692987) q[1];
sx q[1];
rz(-2.0672332) q[1];
sx q[1];
rz(1.4908355) q[1];
x q[2];
rz(-1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(1.0306851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(-2.3922065) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-2.170927) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(-2.326139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.095123) q[0];
sx q[0];
rz(-2.2182584) q[0];
sx q[0];
rz(-1.6460653) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72950659) q[2];
sx q[2];
rz(-1.365005) q[2];
sx q[2];
rz(-1.6934998) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26317877) q[1];
sx q[1];
rz(-2.8385332) q[1];
sx q[1];
rz(-0.11942272) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64579441) q[3];
sx q[3];
rz(-1.3919953) q[3];
sx q[3];
rz(1.7088695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6796391) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(-2.5088076) q[2];
rz(-1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3011424) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(-1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777269) q[0];
sx q[0];
rz(-2.8054872) q[0];
sx q[0];
rz(-1.1019812) q[0];
rz(-0.61727662) q[2];
sx q[2];
rz(-1.0059788) q[2];
sx q[2];
rz(-1.1519943) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-1.0832936) q[1];
sx q[1];
rz(1.8600149) q[1];
x q[2];
rz(-1.4286832) q[3];
sx q[3];
rz(-2.6283773) q[3];
sx q[3];
rz(-1.7538479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6040566) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(0.58004722) q[2];
rz(2.3245658) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(1.9918611) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(-0.27483637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9592322) q[0];
sx q[0];
rz(-1.4920456) q[0];
sx q[0];
rz(-1.6596646) q[0];
x q[1];
rz(-2.0312112) q[2];
sx q[2];
rz(-1.5693671) q[2];
sx q[2];
rz(-1.187385) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8400164) q[1];
sx q[1];
rz(-0.77743545) q[1];
sx q[1];
rz(0.90228723) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87423012) q[3];
sx q[3];
rz(-1.6682373) q[3];
sx q[3];
rz(1.8271354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.794902) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(0.99779469) q[0];
rz(-2.9580341) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.6246187) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51380101) q[0];
sx q[0];
rz(-0.9538981) q[0];
sx q[0];
rz(-2.6898726) q[0];
x q[1];
rz(1.9618271) q[2];
sx q[2];
rz(-2.070825) q[2];
sx q[2];
rz(-0.21991877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5114054) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(1.1955111) q[1];
x q[2];
rz(-0.98832163) q[3];
sx q[3];
rz(-1.054793) q[3];
sx q[3];
rz(-2.2367246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(-0.3240164) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76535392) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.8776241) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(-2.8009159) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1893038) q[0];
sx q[0];
rz(-1.5049107) q[0];
sx q[0];
rz(-0.031866372) q[0];
rz(0.75603007) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(1.1083958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50960474) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(1.8297086) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86990279) q[3];
sx q[3];
rz(-2.5330336) q[3];
sx q[3];
rz(-0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(0.56345338) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1691549) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-2.6766434) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60663215) q[0];
sx q[0];
rz(-2.6323942) q[0];
sx q[0];
rz(1.9726994) q[0];
x q[1];
rz(-2.1452227) q[2];
sx q[2];
rz(-1.9949706) q[2];
sx q[2];
rz(0.043957274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.834224) q[1];
sx q[1];
rz(-2.8982179) q[1];
sx q[1];
rz(0.4537531) q[1];
rz(-pi) q[2];
rz(1.0602337) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(0.75170654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6222318) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(-0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(3.0632339) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0067622234) q[0];
sx q[0];
rz(-1.6400596) q[0];
sx q[0];
rz(-0.80471054) q[0];
rz(-2.9044754) q[2];
sx q[2];
rz(-1.3318921) q[2];
sx q[2];
rz(2.4799926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6137177) q[1];
sx q[1];
rz(-0.67786874) q[1];
sx q[1];
rz(-1.2916958) q[1];
x q[2];
rz(1.3571635) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(-0.3233288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-2.890214) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(0.87402469) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4955935) q[0];
sx q[0];
rz(-2.3388303) q[0];
sx q[0];
rz(-1.5828703) q[0];
rz(-pi) q[1];
rz(-1.9782412) q[2];
sx q[2];
rz(-1.904084) q[2];
sx q[2];
rz(-2.8826706) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8928788) q[1];
sx q[1];
rz(-1.6422179) q[1];
sx q[1];
rz(-2.0715269) q[1];
rz(-pi) q[2];
rz(-1.2285352) q[3];
sx q[3];
rz(-0.98981333) q[3];
sx q[3];
rz(0.57623219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(-1.9571346) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(-1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1353564) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(2.4172879) q[0];
rz(0.18877098) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.1788517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282384) q[0];
sx q[0];
rz(-1.8431219) q[0];
sx q[0];
rz(-0.21270919) q[0];
rz(-pi) q[1];
rz(0.63209052) q[2];
sx q[2];
rz(-0.1569911) q[2];
sx q[2];
rz(2.0096411) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8995754) q[1];
sx q[1];
rz(-1.7574818) q[1];
sx q[1];
rz(0.94355299) q[1];
x q[2];
rz(-1.7970656) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(0.98480485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5205004) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-1.1307217) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(1.422613) q[3];
sx q[3];
rz(-1.2978745) q[3];
sx q[3];
rz(0.10980448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
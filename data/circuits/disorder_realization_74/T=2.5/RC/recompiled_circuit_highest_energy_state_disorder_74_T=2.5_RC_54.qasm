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
rz(-2.3251301) q[0];
sx q[0];
rz(-0.10179585) q[0];
sx q[0];
rz(-0.53959227) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(0.53909477) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2236299) q[0];
sx q[0];
rz(-1.5453891) q[0];
sx q[0];
rz(-0.44141234) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4194054) q[2];
sx q[2];
rz(-1.2149723) q[2];
sx q[2];
rz(-2.9265938) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2869563) q[1];
sx q[1];
rz(-0.35321924) q[1];
sx q[1];
rz(1.1562721) q[1];
x q[2];
rz(-0.081955628) q[3];
sx q[3];
rz(-0.92343077) q[3];
sx q[3];
rz(-2.820197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3794136) q[2];
sx q[2];
rz(-2.2654686) q[2];
sx q[2];
rz(-1.1890821) q[2];
rz(1.9573697) q[3];
sx q[3];
rz(-2.2370179) q[3];
sx q[3];
rz(1.5949465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9963843) q[0];
sx q[0];
rz(-1.5897911) q[0];
sx q[0];
rz(-1.8183964) q[0];
rz(2.6595751) q[1];
sx q[1];
rz(-2.2299485) q[1];
sx q[1];
rz(-0.97420305) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4344113) q[0];
sx q[0];
rz(-1.392258) q[0];
sx q[0];
rz(-2.2212127) q[0];
rz(-1.4575794) q[2];
sx q[2];
rz(-1.6694476) q[2];
sx q[2];
rz(1.1455918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5017101) q[1];
sx q[1];
rz(-2.5474265) q[1];
sx q[1];
rz(-2.7133248) q[1];
rz(1.3311576) q[3];
sx q[3];
rz(-2.1240892) q[3];
sx q[3];
rz(-0.026913337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1366068) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(2.3919487) q[2];
rz(-2.7096115) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5169446) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(-0.6063478) q[0];
rz(-1.8079405) q[1];
sx q[1];
rz(-1.9275815) q[1];
sx q[1];
rz(0.25951728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57588803) q[0];
sx q[0];
rz(-1.4153) q[0];
sx q[0];
rz(-1.3208273) q[0];
rz(2.9664842) q[2];
sx q[2];
rz(-1.3710183) q[2];
sx q[2];
rz(1.7410884) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3644581) q[1];
sx q[1];
rz(-1.9455457) q[1];
sx q[1];
rz(-1.115429) q[1];
rz(-0.43478888) q[3];
sx q[3];
rz(-1.0514976) q[3];
sx q[3];
rz(0.10552191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8848662) q[2];
sx q[2];
rz(-1.7776411) q[2];
sx q[2];
rz(-1.5474896) q[2];
rz(-2.3571842) q[3];
sx q[3];
rz(-1.514785) q[3];
sx q[3];
rz(-1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49482685) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(2.8821017) q[0];
rz(1.2207813) q[1];
sx q[1];
rz(-0.44993284) q[1];
sx q[1];
rz(-1.4422013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4239765) q[0];
sx q[0];
rz(-2.0956796) q[0];
sx q[0];
rz(2.5931326) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7624315) q[2];
sx q[2];
rz(-2.537832) q[2];
sx q[2];
rz(-1.1549283) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.759023) q[1];
sx q[1];
rz(-1.8174531) q[1];
sx q[1];
rz(2.7794629) q[1];
rz(-pi) q[2];
rz(-0.88860926) q[3];
sx q[3];
rz(-1.6251792) q[3];
sx q[3];
rz(2.3171484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1059025) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(-2.7316459) q[2];
rz(1.2116872) q[3];
sx q[3];
rz(-0.68710059) q[3];
sx q[3];
rz(1.3338026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164417) q[0];
sx q[0];
rz(-2.4254159) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(0.43824276) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(-0.73572198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11392191) q[0];
sx q[0];
rz(-0.7589853) q[0];
sx q[0];
rz(0.12667292) q[0];
x q[1];
rz(0.76916285) q[2];
sx q[2];
rz(-0.20649466) q[2];
sx q[2];
rz(-1.002305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3408518) q[1];
sx q[1];
rz(-1.8634808) q[1];
sx q[1];
rz(-0.26009788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92145958) q[3];
sx q[3];
rz(-1.9917352) q[3];
sx q[3];
rz(-0.73900797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2030486) q[2];
sx q[2];
rz(-1.4578578) q[2];
sx q[2];
rz(-2.5277444) q[2];
rz(-2.7326873) q[3];
sx q[3];
rz(-0.96858612) q[3];
sx q[3];
rz(0.74234211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
rz(-2.3530389) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(-1.9045389) q[0];
rz(1.6963814) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(-2.9663185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69578457) q[0];
sx q[0];
rz(-2.9668167) q[0];
sx q[0];
rz(1.4265027) q[0];
rz(-1.0277599) q[2];
sx q[2];
rz(-2.1547085) q[2];
sx q[2];
rz(1.8915382) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1239667) q[1];
sx q[1];
rz(-1.062289) q[1];
sx q[1];
rz(-2.070049) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57717429) q[3];
sx q[3];
rz(-2.1041311) q[3];
sx q[3];
rz(1.2791025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1767629) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(0.97243398) q[2];
rz(-1.6644299) q[3];
sx q[3];
rz(-1.1494613) q[3];
sx q[3];
rz(2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85457388) q[0];
sx q[0];
rz(-3.119097) q[0];
sx q[0];
rz(-0.35312411) q[0];
rz(-2.1137721) q[1];
sx q[1];
rz(-1.9408344) q[1];
sx q[1];
rz(-1.0110528) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42732692) q[0];
sx q[0];
rz(-1.2135226) q[0];
sx q[0];
rz(1.5861804) q[0];
x q[1];
rz(-0.86507823) q[2];
sx q[2];
rz(-1.0891682) q[2];
sx q[2];
rz(0.31246802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2680929) q[1];
sx q[1];
rz(-2.756739) q[1];
sx q[1];
rz(1.9995688) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8500438) q[3];
sx q[3];
rz(-1.8729775) q[3];
sx q[3];
rz(-0.19552375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8280243) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(-1.774452) q[2];
rz(1.8732871) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(-2.2936599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1139514) q[0];
sx q[0];
rz(-0.60751644) q[0];
sx q[0];
rz(2.0116346) q[0];
rz(-0.21895151) q[1];
sx q[1];
rz(-1.4066701) q[1];
sx q[1];
rz(0.91046441) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51800358) q[0];
sx q[0];
rz(-1.2388889) q[0];
sx q[0];
rz(-2.5123572) q[0];
rz(-0.16436188) q[2];
sx q[2];
rz(-0.48381643) q[2];
sx q[2];
rz(1.1587032) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3748951) q[1];
sx q[1];
rz(-1.8378705) q[1];
sx q[1];
rz(-2.1007936) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72515709) q[3];
sx q[3];
rz(-1.2598383) q[3];
sx q[3];
rz(2.6127315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8255446) q[2];
sx q[2];
rz(-2.3393708) q[2];
sx q[2];
rz(2.3160589) q[2];
rz(-0.46142203) q[3];
sx q[3];
rz(-1.8749219) q[3];
sx q[3];
rz(0.44574827) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5017186) q[0];
sx q[0];
rz(-1.8032782) q[0];
sx q[0];
rz(2.3097532) q[0];
rz(-0.72744751) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(-0.096253455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1883586) q[0];
sx q[0];
rz(-1.112845) q[0];
sx q[0];
rz(-1.0012549) q[0];
rz(0.23230884) q[2];
sx q[2];
rz(-2.8806928) q[2];
sx q[2];
rz(2.1292854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0301275) q[1];
sx q[1];
rz(-0.22997738) q[1];
sx q[1];
rz(1.9783918) q[1];
x q[2];
rz(-2.2489266) q[3];
sx q[3];
rz(-0.24991194) q[3];
sx q[3];
rz(-1.3877317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40621296) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(1.5926788) q[2];
rz(0.63273543) q[3];
sx q[3];
rz(-1.9097208) q[3];
sx q[3];
rz(0.24937853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525472) q[0];
sx q[0];
rz(-2.1010375) q[0];
sx q[0];
rz(-0.35300514) q[0];
rz(-0.32265916) q[1];
sx q[1];
rz(-0.32538515) q[1];
sx q[1];
rz(-2.7696612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4034525) q[0];
sx q[0];
rz(-0.36309013) q[0];
sx q[0];
rz(-0.9842666) q[0];
rz(-pi) q[1];
rz(-0.32970365) q[2];
sx q[2];
rz(-1.7404544) q[2];
sx q[2];
rz(0.28126954) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10807762) q[1];
sx q[1];
rz(-1.7308427) q[1];
sx q[1];
rz(-0.33744244) q[1];
x q[2];
rz(-1.6149893) q[3];
sx q[3];
rz(-2.4165476) q[3];
sx q[3];
rz(-2.8738662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17452621) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(-1.2971499) q[2];
rz(-2.2938812) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(-2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90947718) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(-1.7535946) q[1];
sx q[1];
rz(-1.9276062) q[1];
sx q[1];
rz(3.066317) q[1];
rz(-2.4405932) q[2];
sx q[2];
rz(-2.1908115) q[2];
sx q[2];
rz(0.68344982) q[2];
rz(-1.9907436) q[3];
sx q[3];
rz(-2.3227878) q[3];
sx q[3];
rz(-0.25521758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

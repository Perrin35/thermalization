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
rz(-0.78200114) q[0];
sx q[0];
rz(-1.2712103) q[0];
rz(-2.864569) q[1];
sx q[1];
rz(-2.6695873) q[1];
sx q[1];
rz(-0.0013874887) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822617) q[0];
sx q[0];
rz(-1.9388655) q[0];
sx q[0];
rz(0.1749392) q[0];
rz(-2.6944461) q[2];
sx q[2];
rz(-1.6086173) q[2];
sx q[2];
rz(-2.56074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5056155) q[1];
sx q[1];
rz(-0.50230366) q[1];
sx q[1];
rz(0.14640267) q[1];
rz(-pi) q[2];
rz(2.1625159) q[3];
sx q[3];
rz(-2.7032529) q[3];
sx q[3];
rz(-1.0867659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-2.7367676) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(-0.81545365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1708508) q[0];
sx q[0];
rz(-2.4903957) q[0];
sx q[0];
rz(-3.0424776) q[0];
rz(-pi) q[1];
x q[1];
rz(1.843812) q[2];
sx q[2];
rz(-2.2815939) q[2];
sx q[2];
rz(3.0836011) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7533469) q[1];
sx q[1];
rz(-1.2699632) q[1];
sx q[1];
rz(-1.6080329) q[1];
x q[2];
rz(-0.29173298) q[3];
sx q[3];
rz(-0.66666224) q[3];
sx q[3];
rz(3.047903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6796391) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(0.30953428) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3011424) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32854983) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(-2.9850328) q[0];
x q[1];
rz(-0.91018422) q[2];
sx q[2];
rz(-1.0599469) q[2];
sx q[2];
rz(0.78188932) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98993436) q[1];
sx q[1];
rz(-2.5807568) q[1];
sx q[1];
rz(-0.49353091) q[1];
rz(-3.0619377) q[3];
sx q[3];
rz(-1.0632535) q[3];
sx q[3];
rz(1.5505276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(0.58004722) q[2];
rz(-0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(0.27483637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7461473) q[0];
sx q[0];
rz(-1.6593885) q[0];
sx q[0];
rz(-0.079061411) q[0];
rz(1.574013) q[2];
sx q[2];
rz(-0.46041691) q[2];
sx q[2];
rz(-0.38052961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3596613) q[1];
sx q[1];
rz(-2.0205824) q[1];
sx q[1];
rz(-2.2284472) q[1];
rz(-2.2673625) q[3];
sx q[3];
rz(-1.6682373) q[3];
sx q[3];
rz(1.8271354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(-1.5083195) q[2];
rz(1.9968962) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(-0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(-2.143798) q[0];
rz(-0.18355852) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(-1.6246187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.330634) q[0];
sx q[0];
rz(-1.2067544) q[0];
sx q[0];
rz(-0.90322687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5337358) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(1.5460207) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6301873) q[1];
sx q[1];
rz(-1.9747707) q[1];
sx q[1];
rz(1.1955111) q[1];
x q[2];
rz(-2.5449104) q[3];
sx q[3];
rz(-2.0697069) q[3];
sx q[3];
rz(-0.98017207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(2.8175763) q[2];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(-1.8776241) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(-2.8009159) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7621988) q[0];
sx q[0];
rz(-1.6025935) q[0];
sx q[0];
rz(-1.5048774) q[0];
x q[1];
rz(1.864205) q[2];
sx q[2];
rz(-0.83653203) q[2];
sx q[2];
rz(0.66213995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.016547116) q[1];
sx q[1];
rz(-0.53495896) q[1];
sx q[1];
rz(2.6782481) q[1];
x q[2];
rz(-0.86990279) q[3];
sx q[3];
rz(-2.5330336) q[3];
sx q[3];
rz(-0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(0.54780444) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1691549) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-2.6766434) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3195254) q[0];
sx q[0];
rz(-1.7626581) q[0];
sx q[0];
rz(2.0454387) q[0];
rz(2.1452227) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(-3.0976354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.834224) q[1];
sx q[1];
rz(-0.24337473) q[1];
sx q[1];
rz(-2.6878396) q[1];
rz(-pi) q[2];
rz(1.5246478) q[3];
sx q[3];
rz(-1.5449636) q[3];
sx q[3];
rz(-1.8125364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(-0.078358738) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1348304) q[0];
sx q[0];
rz(-1.6400596) q[0];
sx q[0];
rz(2.3368821) q[0];
rz(2.3381091) q[2];
sx q[2];
rz(-2.8065971) q[2];
sx q[2];
rz(-1.4575046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8790508) q[1];
sx q[1];
rz(-1.7444376) q[1];
sx q[1];
rz(0.91211984) q[1];
rz(-pi) q[2];
rz(3.1073242) q[3];
sx q[3];
rz(-1.4141603) q[3];
sx q[3];
rz(2.6019415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(-2.890214) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.4096227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79779977) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(-0.87402469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286205) q[0];
sx q[0];
rz(-2.3734833) q[0];
sx q[0];
rz(-0.012499768) q[0];
rz(1.1633515) q[2];
sx q[2];
rz(-1.2375087) q[2];
sx q[2];
rz(-0.25892205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8585426) q[1];
sx q[1];
rz(-2.0701323) q[1];
sx q[1];
rz(-3.0602171) q[1];
rz(-pi) q[2];
rz(-2.5328818) q[3];
sx q[3];
rz(-1.2864283) q[3];
sx q[3];
rz(1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-0.7243048) q[0];
rz(0.18877098) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(-1.1788517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282384) q[0];
sx q[0];
rz(-1.2984707) q[0];
sx q[0];
rz(2.9288835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4775425) q[2];
sx q[2];
rz(-1.6972731) q[2];
sx q[2];
rz(1.3716413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5621592) q[1];
sx q[1];
rz(-2.4907618) q[1];
sx q[1];
rz(-1.8821554) q[1];
rz(-1.3445271) q[3];
sx q[3];
rz(-2.0824277) q[3];
sx q[3];
rz(-2.1567878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-1.4964676) q[2];
sx q[2];
rz(-2.7004514) q[2];
sx q[2];
rz(-2.2218291) q[2];
rz(-1.422613) q[3];
sx q[3];
rz(-1.8437181) q[3];
sx q[3];
rz(-3.0317882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(-1.707466) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2375862) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(2.0531274) q[0];
rz(-1.7069874) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.1510804) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5601555) q[1];
sx q[1];
rz(-1.1168224) q[1];
sx q[1];
rz(-2.8698688) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1316142) q[3];
sx q[3];
rz(-1.3703128) q[3];
sx q[3];
rz(2.2905614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3216386) q[2];
sx q[2];
rz(-1.7724089) q[2];
sx q[2];
rz(0.83797541) q[2];
rz(-2.6485802) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-0.078991927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74719602) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.863742) q[0];
rz(2.9648119) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(2.7094254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85053274) q[0];
sx q[0];
rz(-1.8467554) q[0];
sx q[0];
rz(2.6002797) q[0];
rz(-pi) q[1];
rz(-2.7414397) q[2];
sx q[2];
rz(-1.1740985) q[2];
sx q[2];
rz(-0.19043365) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.00694) q[1];
sx q[1];
rz(-1.4538527) q[1];
sx q[1];
rz(-2.0140531) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6099986) q[3];
sx q[3];
rz(-2.5762442) q[3];
sx q[3];
rz(-0.70985868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.8816201) q[3];
sx q[3];
rz(2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040722672) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(1.6529282) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(2.6352077) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9837512) q[0];
sx q[0];
rz(-2.7451773) q[0];
sx q[0];
rz(-1.0979963) q[0];
rz(-pi) q[1];
rz(0.73451368) q[2];
sx q[2];
rz(-1.3909512) q[2];
sx q[2];
rz(-0.70914662) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.312078) q[1];
sx q[1];
rz(-1.8784338) q[1];
sx q[1];
rz(2.679146) q[1];
x q[2];
rz(2.0154325) q[3];
sx q[3];
rz(-2.7362842) q[3];
sx q[3];
rz(-1.0682378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4425519) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(-1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716924) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(0.83918321) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(-1.5485839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73491053) q[0];
sx q[0];
rz(-1.4797749) q[0];
sx q[0];
rz(2.0216366) q[0];
rz(0.80231248) q[2];
sx q[2];
rz(-2.1244086) q[2];
sx q[2];
rz(-2.4137036) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6099445) q[1];
sx q[1];
rz(-0.87024401) q[1];
sx q[1];
rz(2.080426) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2876835) q[3];
sx q[3];
rz(-1.1554171) q[3];
sx q[3];
rz(2.0302041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-0.099686064) q[2];
rz(-2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(-1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(0.25594041) q[0];
rz(2.6804965) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(0.76006132) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5998659) q[0];
sx q[0];
rz(-1.4540298) q[0];
sx q[0];
rz(1.6745425) q[0];
rz(-pi) q[1];
rz(-0.063644479) q[2];
sx q[2];
rz(-1.1568937) q[2];
sx q[2];
rz(2.6457583) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.252994) q[1];
sx q[1];
rz(-2.9610486) q[1];
sx q[1];
rz(-0.79043364) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.101196) q[3];
sx q[3];
rz(-2.4393775) q[3];
sx q[3];
rz(1.5839674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5710859) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(0.6742397) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(-0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6102585) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.1791139) q[0];
rz(0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(-1.0669605) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2039295) q[0];
sx q[0];
rz(-1.5852889) q[0];
sx q[0];
rz(-0.020676215) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80337556) q[2];
sx q[2];
rz(-2.5640045) q[2];
sx q[2];
rz(-0.20882777) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31159196) q[1];
sx q[1];
rz(-0.82468669) q[1];
sx q[1];
rz(-3.0768865) q[1];
rz(-pi) q[2];
rz(-0.79640572) q[3];
sx q[3];
rz(-1.685433) q[3];
sx q[3];
rz(-0.20533268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71010464) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(0.55994326) q[2];
rz(-0.7263178) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(-0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85754919) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(2.2834159) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1976397) q[0];
sx q[0];
rz(-1.7174935) q[0];
sx q[0];
rz(-3.0118045) q[0];
rz(-1.6740587) q[2];
sx q[2];
rz(-0.41245663) q[2];
sx q[2];
rz(-0.39706424) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3837636) q[1];
sx q[1];
rz(-1.5997636) q[1];
sx q[1];
rz(-1.5555192) q[1];
rz(1.9167561) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(-0.35017761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.28785607) q[2];
sx q[2];
rz(-0.47984543) q[2];
sx q[2];
rz(-1.8161592) q[2];
rz(-0.89007968) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(-2.7668787) q[0];
rz(-2.162714) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(-1.3495061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193072) q[0];
sx q[0];
rz(-1.95682) q[0];
sx q[0];
rz(0.65667721) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4235054) q[2];
sx q[2];
rz(-1.5880843) q[2];
sx q[2];
rz(-1.3912488) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0950349) q[1];
sx q[1];
rz(-1.3890424) q[1];
sx q[1];
rz(0.05038105) q[1];
rz(-pi) q[2];
rz(2.7262444) q[3];
sx q[3];
rz(-0.58598622) q[3];
sx q[3];
rz(-1.8893482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(1.2619031) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.6040241) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81048548) q[0];
sx q[0];
rz(-1.2665505) q[0];
sx q[0];
rz(-0.23181339) q[0];
rz(-pi) q[1];
rz(0.41861694) q[2];
sx q[2];
rz(-1.6659684) q[2];
sx q[2];
rz(-1.0375432) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21737145) q[1];
sx q[1];
rz(-1.4573759) q[1];
sx q[1];
rz(0.70593112) q[1];
rz(-pi) q[2];
rz(-2.780464) q[3];
sx q[3];
rz(-1.1361406) q[3];
sx q[3];
rz(-1.9643195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1016772) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(-2.3802479) q[2];
rz(-2.2411761) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(-0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(0.21690579) q[0];
rz(-2.5096109) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(2.1868618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.117884) q[0];
sx q[0];
rz(-2.3291991) q[0];
sx q[0];
rz(-0.448416) q[0];
rz(-pi) q[1];
rz(0.63072272) q[2];
sx q[2];
rz(-0.64881697) q[2];
sx q[2];
rz(-1.1184675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6001544) q[1];
sx q[1];
rz(-2.6791875) q[1];
sx q[1];
rz(-2.9701783) q[1];
rz(0.99777625) q[3];
sx q[3];
rz(-2.4224835) q[3];
sx q[3];
rz(-2.5620808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0163991) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(2.1968502) q[2];
rz(0.38481209) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5007297) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(0.8846994) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-0.62467081) q[2];
sx q[2];
rz(-0.93745898) q[2];
sx q[2];
rz(-2.6425101) q[2];
rz(-1.1576204) q[3];
sx q[3];
rz(-0.4784085) q[3];
sx q[3];
rz(2.9984409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

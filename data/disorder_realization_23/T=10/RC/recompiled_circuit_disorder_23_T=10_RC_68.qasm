OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8795348) q[0];
sx q[0];
rz(-1.4095925) q[0];
sx q[0];
rz(-1.4341266) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2375862) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(-2.0531274) q[0];
x q[1];
rz(-1.4346052) q[2];
sx q[2];
rz(-1.681466) q[2];
sx q[2];
rz(1.1510804) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.13222209) q[1];
sx q[1];
rz(-1.3271866) q[1];
sx q[1];
rz(-1.1019215) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0099785) q[3];
sx q[3];
rz(-1.7712799) q[3];
sx q[3];
rz(2.2905614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.819954) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-0.83797541) q[2];
rz(0.49301246) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-0.078991927) q[3];
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
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74719602) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.2778506) q[0];
rz(-2.9648119) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(0.4321672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9958187) q[0];
sx q[0];
rz(-0.60129014) q[0];
sx q[0];
rz(-2.6390618) q[0];
rz(-pi) q[1];
x q[1];
rz(0.400153) q[2];
sx q[2];
rz(-1.9674941) q[2];
sx q[2];
rz(0.19043365) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1951616) q[1];
sx q[1];
rz(-2.6841607) q[1];
sx q[1];
rz(-1.8381579) q[1];
rz(-pi) q[2];
rz(1.0057955) q[3];
sx q[3];
rz(-1.5917935) q[3];
sx q[3];
rz(0.82783031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5923578) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(-2.6548927) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(-1.6529282) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(-2.6352077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7931472) q[0];
sx q[0];
rz(-1.9217102) q[0];
sx q[0];
rz(-2.9532414) q[0];
rz(-0.73451368) q[2];
sx q[2];
rz(-1.3909512) q[2];
sx q[2];
rz(0.70914662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.312078) q[1];
sx q[1];
rz(-1.8784338) q[1];
sx q[1];
rz(2.679146) q[1];
x q[2];
rz(-1.9403463) q[3];
sx q[3];
rz(-1.7412211) q[3];
sx q[3];
rz(2.2263262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69904077) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(0.58602035) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1699003) q[0];
sx q[0];
rz(-0.62830347) q[0];
sx q[0];
rz(0.83918321) q[0];
rz(3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(-1.5485839) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87985932) q[0];
sx q[0];
rz(-2.0196336) q[0];
sx q[0];
rz(-3.0405322) q[0];
x q[1];
rz(-2.4315005) q[2];
sx q[2];
rz(-0.93821628) q[2];
sx q[2];
rz(-1.8284423) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7566484) q[1];
sx q[1];
rz(-1.1886016) q[1];
sx q[1];
rz(2.3734943) q[1];
rz(-1.8539092) q[3];
sx q[3];
rz(-1.9861756) q[3];
sx q[3];
rz(-2.0302041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(-0.099686064) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-2.8856522) q[0];
rz(-2.6804965) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(-0.76006132) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041199112) q[0];
sx q[0];
rz(-1.6738335) q[0];
sx q[0];
rz(-0.11739199) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9854457) q[2];
sx q[2];
rz(-1.6290602) q[2];
sx q[2];
rz(2.0410048) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.252994) q[1];
sx q[1];
rz(-0.18054403) q[1];
sx q[1];
rz(-0.79043364) q[1];
rz(0.36565904) q[3];
sx q[3];
rz(-2.1846111) q[3];
sx q[3];
rz(0.97096503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5710859) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(-0.6742397) q[2];
rz(0.21480602) q[3];
sx q[3];
rz(-2.6847697) q[3];
sx q[3];
rz(-3.1242483) q[3];
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
sx q[0];
rz(pi/2) q[0];
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
rz(-2.9367661) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(-2.0746322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5087591) q[0];
sx q[0];
rz(-1.5914704) q[0];
sx q[0];
rz(1.5563006) q[0];
x q[1];
rz(-2.7166769) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(2.0770819) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9251717) q[1];
sx q[1];
rz(-2.3932218) q[1];
sx q[1];
rz(1.5009576) q[1];
rz(-2.3451869) q[3];
sx q[3];
rz(-1.4561597) q[3];
sx q[3];
rz(-0.20533268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.431488) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(-2.2834159) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9265384) q[0];
sx q[0];
rz(-2.9460322) q[0];
sx q[0];
rz(-0.8514479) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.467534) q[2];
sx q[2];
rz(-0.41245663) q[2];
sx q[2];
rz(-2.7445284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8691751) q[1];
sx q[1];
rz(-0.032748001) q[1];
sx q[1];
rz(2.6564024) q[1];
rz(-pi) q[2];
rz(0.51123294) q[3];
sx q[3];
rz(-2.5066262) q[3];
sx q[3];
rz(2.883203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.3254335) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(2.1530698) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(-2.7668787) q[0];
rz(0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.7920866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193072) q[0];
sx q[0];
rz(-1.95682) q[0];
sx q[0];
rz(-2.4849154) q[0];
x q[1];
rz(0.026272341) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(0.19933137) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.3186036) q[1];
sx q[1];
rz(-0.18853304) q[1];
sx q[1];
rz(1.3033426) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5958063) q[3];
sx q[3];
rz(-1.7958399) q[3];
sx q[3];
rz(-3.1080266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-0.46869579) q[2];
rz(-1.9474585) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(1.3635925) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63012183) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.2619031) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(-1.5375686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3311072) q[0];
sx q[0];
rz(-1.8750422) q[0];
sx q[0];
rz(2.9097793) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41861694) q[2];
sx q[2];
rz(-1.6659684) q[2];
sx q[2];
rz(2.1040495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2214339) q[1];
sx q[1];
rz(-0.71343525) q[1];
sx q[1];
rz(2.9677797) q[1];
x q[2];
rz(-0.92026199) q[3];
sx q[3];
rz(-0.55763054) q[3];
sx q[3];
rz(-0.44616163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-0.76134479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(-3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(-2.5096109) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(0.95473081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023708658) q[0];
sx q[0];
rz(-2.3291991) q[0];
sx q[0];
rz(-0.448416) q[0];
rz(-pi) q[1];
rz(-2.5108699) q[2];
sx q[2];
rz(-2.4927757) q[2];
sx q[2];
rz(-2.0231252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6001544) q[1];
sx q[1];
rz(-2.6791875) q[1];
sx q[1];
rz(-2.9701783) q[1];
rz(-pi) q[2];
rz(0.99777625) q[3];
sx q[3];
rz(-0.7191092) q[3];
sx q[3];
rz(-0.57951187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0163991) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(-0.94474244) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64086296) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(0.83512886) q[2];
sx q[2];
rz(-2.0616812) q[2];
sx q[2];
rz(2.4731935) q[2];
rz(0.20529071) q[3];
sx q[3];
rz(-2.0060354) q[3];
sx q[3];
rz(2.5397186) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
